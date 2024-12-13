use std::collections::HashMap;
use ngrams::Ngram;
use seq_io::fasta::{Reader, Record};
use std::io::{self, BufReader};
use std_dev::standard_deviation;

const KMER_SIZE: usize = 13;
const OVLP_WINDOWS_SIZE: usize = 100;
const MAX_ACCEPTED_MISMATCH: usize = 0;

struct SequenceRecord {
    name: String,
    sequence: String,
}

struct KmerScore {
    kmer: String,
    score: f64,
}

fn to_records(src: Vec<u8>) -> io::Result<Vec<SequenceRecord>> {
    let mut reader = Reader::new(BufReader::new(src.as_slice()));
    let mut records = Vec::new();

    while let Some(result) = reader.next() {
        let record = result.unwrap();
        let name = record.id().unwrap().to_string();
        let sequence = String::from_utf8(record.full_seq().to_vec())
            .unwrap()
            .to_uppercase()
            .replace("U", "T");
        records.push(SequenceRecord { name, sequence });
    }
    Ok(records)
}

/**
 * Aligns sequences using MAFFT
 */
fn align_sequences(filepath: String) -> Result<Vec<u8>, io::Error> {
    let output = std::process::Command::new("mafft")
        .args([
            "--auto",
            "--quiet",
            "--thread",
            "-1",
            &filepath.clone(),
        ])
        .output()
        .expect("failed to execute MAFFT");

    Ok(output.stdout)
}

fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'U' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => c
    }).collect()
}

fn get_sequence_kmers(sequence: &str, kmer_size: usize) -> Vec<String> {
    sequence.chars().ngrams(kmer_size)
        .filter(|kmer| kmer.iter().all(|c| !"Nn- ".contains(*c)))
        .map(|kmer| kmer.iter().collect())
        .collect()
}

/**
 * Get the related kmers of a given kmer
 *
 * For example, given a kmer "ATCG_T", the related kmer(key) in the segment are:
 * - "ATCGAT"
 * - "ATCGGT"
 * - "ATCGCT"
 * - "ATCGTT"
 */
fn inc_related_primers_freq(freqs: &mut HashMap<String, u16>, kmer: &String) {
    freqs.entry(kmer.clone()).and_modify(|e| *e += 1);
    if MAX_ACCEPTED_MISMATCH == 0 {
        return;
    }

    let keys: Vec<String> = freqs.keys().cloned().collect();
    'outer: for key in keys {
        let mut mismatch_count = 0;
        for (i, c) in key.chars().enumerate() {
            if c != kmer.chars().nth(i).unwrap() {
                mismatch_count += 1;
            }
            if mismatch_count > MAX_ACCEPTED_MISMATCH {
                continue 'outer;
            }
        }
        if kmer.ne(&key) {
            log::debug!("Found related kmer for {} is {}", kmer, key);
            freqs.entry(key).and_modify(|e| *e += 1);
        }
    }
}

fn get_most_freq_kmers(segment_kmers: &Vec<String>) -> Vec<KmerScore> {
    // pre-build the frequency hashmap
    let mut freqs: HashMap<String, u16> = HashMap::new();
    for kmer in segment_kmers {
        freqs.insert(kmer.clone(), 0);
    }

    for kmer in segment_kmers {
        // freqs.entry(kmer.clone()).and_modify(|e| *e += 1);
        inc_related_primers_freq(&mut freqs, kmer);
    }
    let mut scores = Vec::new();
    for (kmer, count) in freqs {
        scores.push(KmerScore {
            kmer: kmer.clone(),
            score: count as f64,
        });
    }
    scores.sort_by_key(|k| -1 * (k.score as i64));
    scores
}

fn partitioning_sequence(sequence: &str, ovlp_windows_size: usize) -> Vec<String> {
    sequence.chars().collect::<Vec<char>>().windows(ovlp_windows_size * 2)
        .step_by(ovlp_windows_size)
        .map(|window| window.iter().collect())
        .collect()
}

fn get_segments(records: &Vec<SequenceRecord>) -> (HashMap<u16, Vec<String>>, HashMap<u16, Vec<String>>) {
    let mut fwd_segments: HashMap<u16, Vec<String>> = HashMap::new();
    let mut rev_segments: HashMap<u16, Vec<String>> = HashMap::new();

    for rec in records.iter() {
        // 1. Partitioning the sequence into segments
        let seq_segments = partitioning_sequence(&rec.sequence, OVLP_WINDOWS_SIZE);
        // 2. Extracting forward k-mers from each segment
        for (segment_idx, segment) in seq_segments.iter().enumerate() {
            let fwd_k_mers = get_sequence_kmers(segment, KMER_SIZE);
            // 3. Push the forward k-mers into a hashmap
            let entry = fwd_segments.entry(segment_idx as u16).or_insert(Vec::new());
            for k_mer in fwd_k_mers {
                entry.push(k_mer);
            }
            // 4. Push the reverse k-mers into a hashmap
            let rev_k_mers = get_sequence_kmers(&reverse_complement(segment), KMER_SIZE);
            let entry = rev_segments.entry(segment_idx as u16).or_insert(Vec::new());
            for k_mer in rev_k_mers {
                entry.push(k_mer);
            }
        }

    }

    (fwd_segments, rev_segments)
}

fn get_gc_percent(kmer: &String) -> f32 {
    let mut gc_count: f32 = 0.0;
    for c in kmer.chars() {
        if c == 'G' || c == 'C' {
            gc_count += 1.0;
        }
    }
    (gc_count / kmer.len() as f32) * 100.0
}

fn get_tm(kmer: &String) -> f32 {
    let mut gc_count: f32 = 0.0;
    let mut at_count: f32 = 0.0;
    for c in kmer.chars() {
        if c == 'G' || c == 'C' {
            gc_count += 1.0;
        }
        if c == 'A' || c == 'T' {
            at_count += 1.0;
        }

    }
    64.9 + (41.0 * (gc_count - 16.4) / (gc_count + at_count))
}

/**
 * Get the threshold for Tm
 *
 * Calculated by find (2*sd(Tm)) + mean(Tm) of the primers
 */
fn get_tm_threshold(primers: Vec<String>) -> f32 {
    let mut tm_values: Vec<f32> = Vec::new();
    for primer in primers {
        tm_values.push(get_tm(&primer));
    }
    let mean_tm = tm_values.iter().sum::<f32>() / tm_values.len() as f32;
    let sd_tm = standard_deviation(&tm_values);
    (2.0 * sd_tm.standard_deviation) + mean_tm
}

fn in_tm_threshold(kmer: &String, threshold: f32) -> bool {
    get_tm(kmer) < threshold
}

fn is_homopolymer(kmer: &String) -> bool {
    let mut prev_char = ' ';
    let mut count = 0;
    for c in kmer.chars() {
        if prev_char == ' ' || c == prev_char {
            count += 1;
        } else {
            count = 0;
        }
        if count >= 5 {
            return true;
        }
        prev_char = c;
    }
    false
}

fn filter_tm_and_homopolymers(kmers: Vec<KmerScore>, threshold: f32) -> Vec<KmerScore> {
    kmers.iter()
        .filter(|k| in_tm_threshold(&k.kmer, threshold) && !is_homopolymer(&k.kmer))
        .map(|k| KmerScore {
            kmer: k.kmer.clone(),
            score: k.score,
        })
        .collect()
}

fn main() -> io::Result<()> {
    env_logger::init();

    let filename = String::from("samples.fasta");

    // 1. Align sequences
    log::debug!("Aligning sequences...");
    let records = match align_sequences(filename) {
        Ok(records) => to_records(records),
        Err(e) => {
            panic!("Error aligning sequences: {}", e);
        }
    };
    log::debug!("Done aligning sequences");
    if records.iter().len() == 0 {
        panic!("No sequences found in the input file");
    }

    // 2. Extracting n-grams from each sequence segments
    log::debug!("Extracting n-grams from each sequence segments...");
    let (fwd_segments, rev_segments) = get_segments(&records?);
    log::debug!("Total segments: fwd:{}, rev:{}", fwd_segments.len(), rev_segments.len());

    // 3. Calculate frequencies of n-grams for each segment both forward/reverse
    log::debug!("Calculating frequencies of n-grams for each segment...");
    let mut fwd_kmers_freqs: Vec<KmerScore> = Vec::new();
    let mut rev_kmers_freqs: Vec<KmerScore> = Vec::new();
    for (segment_idx, segment) in fwd_segments.iter() {
        let fwd_scores = get_most_freq_kmers(segment);
        let best_fwd_score = fwd_scores.first().unwrap();
        fwd_kmers_freqs.push(KmerScore {
            kmer: best_fwd_score.kmer.clone(),
            score: best_fwd_score.score,
        });

        let rev_scores = get_most_freq_kmers(rev_segments.get(segment_idx).unwrap());
        let best_rev_score = rev_scores.first().unwrap();
        rev_kmers_freqs.push(KmerScore {
            kmer: best_rev_score.kmer.clone(),
            score: best_rev_score.score,
        });
    }
    log::debug!("Done calculating frequencies of n-grams for each segment");
    log::debug!("Total forward kmers: {}", fwd_kmers_freqs.len());
    for kmer in &fwd_kmers_freqs[0..5] {
        log::debug!("Forward kmer: {} with score={}, tm={}, is_homopolymers={}", kmer.kmer, kmer.score, get_tm(&kmer.kmer), is_homopolymer(&kmer.kmer));
    }
    for kmer in &rev_kmers_freqs[0..5] {
        log::debug!("Reverse kmer: {} with score={}, tm={}, is_homopolymers={}", kmer.kmer, kmer.score, get_tm(&kmer.kmer), is_homopolymer(&kmer.kmer));
    }
    log::debug!("Total reverse kmers: {}", rev_kmers_freqs.len());

    // Filtering out unmatch Tm and >5nt repeats
    log::debug!("Filtering out unmatch Tm and >5nt repeats...");
    let fwd_tm_threshold = get_tm_threshold(fwd_kmers_freqs.iter().map(|k| k.kmer.clone()).collect());
    let rev_tm_threshold = get_tm_threshold(rev_kmers_freqs.iter().map(|k| k.kmer.clone()).collect());
    log::debug!("Forward Tm threshold: {}", fwd_tm_threshold);
    log::debug!("Reverse Tm threshold: {}", rev_tm_threshold);
    let mut fwd_kmers_filtered: Vec<KmerScore> = filter_tm_and_homopolymers(fwd_kmers_freqs, fwd_tm_threshold);
    let mut rev_kmers_filtered: Vec<KmerScore> = filter_tm_and_homopolymers(rev_kmers_freqs, rev_tm_threshold);
    log::debug!("Done filtering out unmatch Tm and >5nt repeats");

    // 4. Output the primers
    log::debug!("Outputting the primers...");
    for (i, kmer) in fwd_kmers_filtered.iter().enumerate() {
        println!("Forward primer {}: kmer={}, score={}, tm={}, is_homopolymer={}", i, kmer.kmer, kmer.score, get_tm(&kmer.kmer), is_homopolymer(&kmer.kmer));
    }
    for (i, kmer) in rev_kmers_filtered.iter().enumerate() {
        println!("Reverse primer {}: kmer={}, score={}, tm={}, is_homopolymer={}", i, kmer.kmer, kmer.score, get_tm(&kmer.kmer), is_homopolymer(&kmer.kmer));
    }

    fwd_kmers_filtered.sort_by_key(|k| -1 * (k.score as i64));
    rev_kmers_filtered.sort_by_key(|k| -1 * (k.score as i64));

    for (i, item) in fwd_kmers_filtered.iter().enumerate() {
        let gc_percent = get_gc_percent(&item.kmer);
        if 40.0 > gc_percent || gc_percent > 60.0 {
            continue;
        }
        println!("FWD_{},{},Zika virus,64320", i, item.kmer);
    }
    for (i, item) in rev_kmers_filtered.iter().enumerate() {
        let gc_percent = get_gc_percent(&item.kmer);
        if 40.0 > gc_percent || gc_percent > 60.0 {
            continue;
        }
        println!("REV_{},{},Zika virus,64320", i, item.kmer);
    }

    Ok(())
}