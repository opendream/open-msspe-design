mod delta_g;

use ngrams::Ngram;
use seq_io::fasta::{Reader, Record};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::{self, BufReader};
use itertools::Itertools;
use std_dev::standard_deviation;

const KMER_SIZE: usize = 13;
const OVLP_WINDOWS_SIZE: usize = 250;
const MAX_MISMATCH_SEQUENCES: usize = 0;
const MAX_ITERATIONS: usize = 100;

const MV_CONC: f64 = 50.0; // Monovalent cation concentration (mM)
const DV_CONC: f64 = 0.0; // Divalent cation concentration (mM)
const DNTP_CONC: f64 = 0.8; // dNTP concentration (mM)
const DNA_CONC: f64 = 50.0; // Primer concentration (nM)
const ANNEALING_TEMP: f64 = 45.0; // Annealing temperature (°C)

struct SequenceRecord {
    name: String,
    sequence: String,
}

const SEQ_DIR_FWD: u8 = 0x00;
const SEQ_DIR_REV: u8 = 0x01;

#[derive(Hash, Clone)]
struct KmerRecord {
    word: String,
    direction: u8,
    seq_id: usize,
}

impl PartialEq for KmerRecord {
    fn eq(&self, other: &Self) -> bool {
        self.word == other.word && self.direction == other.direction
    }
}

impl Eq for KmerRecord {}

#[derive(Hash, Clone)]
struct KmerFrequency {
    kmer: KmerRecord,
    frequency: usize,
}

impl PartialEq for KmerFrequency {
    fn eq(&self, other: &Self) -> bool {
        self.kmer.word == other.kmer.word
    }
}

impl Eq for KmerFrequency {}

struct KmerStat {
    word: String,
    direction: u8,
    frequency: usize,
    gc_percent: f32,
    tm: f32,
    tm_ok: bool,
    repeats: bool,
    runs: bool,
    delta_g: f32,
    hairpin: bool,
}

struct Segment {
    index: u16,
    kmers: HashMap<KmerRecord, usize>,
}

type KmerSegmentMapping = HashMap<String, HashMap<usize, u32>>;

struct FindingCandidateKmersContext<'a> {
    segments: &'a HashMap<usize, Segment>,
    skipped_segment_indexes: &'a mut HashSet<usize>,
    kmer_segments_mapping: &'a KmerSegmentMapping,
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
        .args(["--auto", "--quiet", "--thread", "-1", &filepath.clone()])
        .output()
        .expect("failed to execute MAFFT");

    Ok(output.stdout)
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

fn find_kmers(sequence: &str, kmer_size: usize) -> Vec<String> {
    sequence
        .chars()
        .ngrams(kmer_size)
        .filter(|kmer| kmer.iter().all(|c| !"Nn- ".contains(*c)))
        .map(|kmer| kmer.iter().collect())
        .unique()
        .collect()
}

fn partitioning_sequence(sequence: &str, ovlp_windows_size: usize) -> Vec<String> {
    sequence
        .chars()
        .collect::<Vec<char>>()
        .windows(ovlp_windows_size * 2)
        .step_by(ovlp_windows_size)
        .map(|window| window.iter().collect())
        .collect()
}

fn get_segments(records: &Vec<SequenceRecord>) -> HashMap<usize, Segment> {
    let mut segments: HashMap<usize, Segment> = HashMap::new();

    for (seq_id, rec) in records.iter().enumerate() {
        // 1. Partitioning the sequence into segments
        let partitions = partitioning_sequence(&rec.sequence, OVLP_WINDOWS_SIZE);
        log::debug!("seq_id={}, partitions={}, first_part_len={}", seq_id, partitions.len(), partitions[0].len());
        // 2. Extracting forward k-mers from each segment
        for (idx, partition) in partitions.iter().enumerate() {
            // Find the segment at idx, if not exists, create a new one
            let segment = segments.entry(idx).or_insert(Segment {
                index: idx as u16,
                kmers: HashMap::new(),
            });
            // Push the k-mers into the list
            let fwd_kmers: Vec<KmerRecord> = find_kmers(partition, KMER_SIZE)
                .iter()
                .map(|kmer| KmerRecord {
                    word: kmer.to_string(),
                    direction: SEQ_DIR_FWD,
                    seq_id,
                })
                .collect();
            let rev_kmers: Vec<KmerRecord> = find_kmers(&reverse_complement(partition), KMER_SIZE)
                .iter()
                .map(|kmer| KmerRecord {
                    word: kmer.to_string(),
                    direction: SEQ_DIR_REV,
                    seq_id,
                })
                .collect();

            for kmer in fwd_kmers.iter() {
                segment
                    .kmers
                    .entry(kmer.clone())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
            for kmer in rev_kmers.iter() {
                segment
                    .kmers
                    .entry(kmer.clone())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
        }
    }

    segments
}

fn make_kmer_segments_mapping(segments: &HashMap<usize, Segment>) -> KmerSegmentMapping {
    let mut kmer_segments_mapping: KmerSegmentMapping = HashMap::new();
    for (i, segment) in segments.iter() {
        for (kmer, freq) in segment.kmers.iter() {
            let entry = KmerRecord {
                word: kmer.word.clone(),
                direction: kmer.direction,
                seq_id: kmer.seq_id,
            };
            kmer_segments_mapping
                .entry(entry.word)
                .or_insert(HashMap::new())
                .insert(*i, *freq as u32);
        }
    }
    kmer_segments_mapping
}

fn find_most_freq_kmer<'a>(
    segments: &'a HashMap<usize, Segment>,
    kmer_segments_mapping: &'a KmerSegmentMapping,
    skipped_segment_indexes: HashSet<usize>,
) -> (KmerFrequency, HashSet<usize>) {
    let mut kmer_total_counts: HashMap<&KmerRecord, usize> = HashMap::new();
    let mut candidate_kmers: HashMap<usize, KmerFrequency> = HashMap::new();
    for (_, segment) in segments.iter() {
        if skipped_segment_indexes.contains(&(segment.index as usize)) {
            continue;
        }

        let mut kmer_map: HashMap<String, &KmerRecord> = HashMap::new();
        let mut kmer_counts: HashMap<String, usize> = HashMap::new();
        for (kmer_record, counts) in segment.kmers.iter() {
            kmer_map.entry(kmer_record.word.clone()).or_insert(kmer_record);
            kmer_counts.entry(kmer_record.word.clone()).and_modify(|e| *e += counts).or_insert(*counts);
            kmer_total_counts
                .entry(kmer_record)
                .and_modify(|e| *e += counts)
                .or_insert(*counts);
        }
        let (word, freq) = kmer_counts.into_iter().max_by(|(_, a), (_, b)| a.cmp(b)).unwrap();
        let winner_kmer = kmer_map.get(&word).unwrap();
        candidate_kmers.insert(
            segment.index as usize,
            KmerFrequency {
                kmer: KmerRecord {
                    word: word.clone(),
                    direction: winner_kmer.direction,
                    seq_id: winner_kmer.seq_id,
                },
                frequency: freq,
            },
        );
    }
    let (_, winner) = candidate_kmers
        .iter()
        .max_by(|(_, a), (_, b)| {
            let a_freq = kmer_total_counts.get(&a.kmer).unwrap();
            let b_freq = kmer_total_counts.get(&b.kmer).unwrap();
            a_freq.cmp(b_freq)
        })
        .unwrap();

    let winner_kmer_segments = kmer_segments_mapping.get(&winner.kmer.word.clone());
    let matched_segment_indexes: Vec<usize> = winner_kmer_segments.unwrap().keys().map(|k| *k).collect();
    let matched_segment_indexes_rev: Vec<usize> = kmer_segments_mapping
        .get(&reverse_complement(&winner.kmer.word.clone()))
        .unwrap()
        .keys()
        .map(|k| *k)
        .collect();

    let mut skipped_segment_indexes: HashSet<usize> = HashSet::new();
    for idx in matched_segment_indexes {
        skipped_segment_indexes.insert(idx);
    }
    for idx in matched_segment_indexes_rev {
        skipped_segment_indexes.insert(idx);
    }

    (
        KmerFrequency {
            kmer: KmerRecord {
                word: winner.kmer.word.clone(),
                direction: winner.kmer.direction,
                seq_id: winner.kmer.seq_id,
            },
            frequency: winner.frequency,
        },
        skipped_segment_indexes,
    )
}

mod tests {
    use super::*;

    #[test]
    fn test_find_most_freq_kmer() {
        let mut segments: HashMap<usize, Segment> = HashMap::new();
        let mut kmer_segments_mapping: KmerSegmentMapping = HashMap::new();
        let mut skipped_segment_indexes: HashSet<usize> = HashSet::new();

        let segment = Segment {
            index: 0,
            kmers: HashMap::new(),
        };
        segments.insert(0, segment);

        let (kmer, skipped) = find_most_freq_kmer(&segments, &kmer_segments_mapping, skipped_segment_indexes);
        assert_eq!(kmer.frequency, 0);
        assert_eq!(skipped.len(), 0);
    }
}

fn get_candidates_kmers<'a>(
    records: &'a Vec<SequenceRecord>,
    segments: &'a HashMap<usize, Segment>,
) -> Vec<KmerFrequency> {
    let mut skipped_segment_indexes: HashSet<usize> = HashSet::new();
    let mut candidate_kmers: Vec<KmerFrequency> = Vec::new();
    let kmer_segment_mappings = make_kmer_segments_mapping(&segments);

    let total_sequences = records.iter().len();
    let min_remaining_sequences = total_sequences - MAX_MISMATCH_SEQUENCES;

    for iter_no in 0..MAX_ITERATIONS {
        let (iter_winner, new_skipped_indexes) = find_most_freq_kmer(
            &segments,
            &kmer_segment_mappings,
            skipped_segment_indexes.clone()
        );
        log::debug!("Iteration={}, winner: {}, freq={}, skips={}", iter_no, iter_winner.kmer.word, iter_winner.frequency, new_skipped_indexes.len());
        for idx in new_skipped_indexes.iter() {
            skipped_segment_indexes.insert(*idx);
        }
        // if remaining missing segments are less than MAX_MISMATCH_SEQUENCES, break
        let x = kmer_segment_mappings.get(&iter_winner.kmer.word);
        log::debug!("Iteration={}, missing={}, min_segments={}", iter_no, x.iter().len(), min_remaining_sequences);
        if kmer_segment_mappings.get(&iter_winner.kmer.word).iter().len() < min_remaining_sequences {
            // break;
        }
        candidate_kmers.push(iter_winner);
    }

    candidate_kmers.iter().map(|k| KmerFrequency {
        kmer: KmerRecord {
            word: k.kmer.word.clone(),
            direction: k.kmer.direction,
            seq_id: k.kmer.seq_id,
        },
        frequency: k.frequency,
    }).collect()
}

fn get_kmer_stats(kmer_records: Vec<KmerFrequency>) -> Vec<KmerStat> {
    // first, finding the threshold for Tm
    let primers: Vec<String> = kmer_records.iter().map(|k| k.kmer.word.clone()).collect();
    let tm_threshold = get_tm_threshold(primers);

    kmer_records
        .iter()
        .map(|kmer_freq| {
            let freq = kmer_freq.frequency;
            let tm = get_tm(kmer_freq.kmer.word.clone());
            let delta_g = 0.0;
            KmerStat {
                word: kmer_freq.kmer.word.clone(),
                direction: kmer_freq.kmer.direction,
                frequency: freq,
                gc_percent: get_gc_percent(kmer_freq.kmer.word.clone()),
                tm,
                tm_ok: in_tm_threshold(kmer_freq.kmer.word.clone(), tm_threshold),
                repeats: is_repeats(kmer_freq.kmer.word.clone()),
                runs: is_run(kmer_freq.kmer.word.clone()),
                delta_g,
                hairpin: delta_g < -9.0,
            }
        })
        .collect()
}

fn get_gc_percent(kmer: String) -> f32 {
    let mut gc_count: f32 = 0.0;
    for c in kmer.chars() {
        if c == 'G' || c == 'C' {
            gc_count += 1.0;
        }
    }
    (gc_count / kmer.len() as f32) * 100.0
}

fn get_tm(kmer: String) -> f32 {
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
        tm_values.push(get_tm(primer));
    }
    let mean_tm = tm_values.iter().sum::<f32>() / tm_values.len() as f32;
    let sd_tm = standard_deviation(&tm_values);
    (2.0 * sd_tm.standard_deviation) + mean_tm
}

fn in_tm_threshold(kmer: String, threshold: f32) -> bool {
    get_tm(kmer) < threshold
}

/**
 * Check if the kmer  has >=5nt di-nucleotide repeats
 *
 * For example, ATATATATATGG is too many AT repeats, then return `true`.
 */
fn is_repeats(kmer: String) -> bool {
    let kmer_ngrams = kmer.chars().ngrams(2);
    let mut repeats = 0;
    let mut last_chunk: Vec<char> = Vec::new();
    for chunk in kmer_ngrams {
        if chunk == last_chunk {
            repeats += 1;
        } else {
            repeats = 0;
        }
        last_chunk = chunk;
    }
    repeats >= 5
}

fn is_run(kmer: String) -> bool {
    let mut runs = 0;
    let mut last_char = ' ';
    for c in kmer.chars() {
        if c == last_char {
            runs += 1;
        } else {
            runs = 0;
        }
        last_char = c;
    }
    runs >= 5
}

fn main() -> io::Result<()> {
    env_logger::init();

    let filename = String::from("samples.fasta");

    // 1. Align sequences
    log::debug!("Aligning sequences...");
    let records = match align_sequences(filename) {
        Ok(records) => to_records(records)?,
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
    let segments = get_segments(&records);
    log::debug!("Done, Total segments: {}", segments.len());

    // 3. Calculate frequencies of n-grams for each segment both forward/reverse
    log::debug!("Calculating frequencies of k-mer for all segments...");
    let candidate_kmers = get_candidates_kmers(&records, &segments);
    log::debug!(
        "Done calculating, Total k-mers left: {}",
        candidate_kmers.len()
    );

    // 4. Filtering out unmatched criteria
    log::debug!("Filtering out unmatched criteria (Tm and >5nt repeats, runs...)");
    let kmer_stats = get_kmer_stats(candidate_kmers);
    let candidate_primers: Vec<&KmerStat> = kmer_stats
        .iter()
        .filter(|kmer_stat| {
            kmer_stat.tm_ok
                && !kmer_stat.repeats
                && !kmer_stat.runs
                && kmer_stat.gc_percent >= 40.0
                && kmer_stat.gc_percent <= 60.0
        })
        .collect();
    log::debug!("Done filtering out unmatched");

    // 5. Output the primers
    log::debug!("Outputting primers...");
    println!("Total primers: {}", candidate_primers.len());

    Ok(())
}
