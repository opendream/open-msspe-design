mod delta_g;
mod graphdb;
mod primer;

use crate::delta_g::{run_ntthal, NtthalOptions};
use crate::primer::{check_primers, CheckPrimerParams, PrimerInfo, DEFAULT_MAX_TM, DEFAULT_MIN_TM};
use clap::Parser;
use graphdb::Edge;
use itertools::Itertools;
use ngrams::Ngram;
use seq_io::fasta::{Reader, Record};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::{self, BufReader};
use std_dev::standard_deviation;

const KMER_SIZE: usize = 13;
const WINDOW_SIZE: usize = 500;
const OVERLAP_SIZE: usize = 250;
const MAX_MISMATCH_SEGMENTS: usize = 1;
const MAX_ITERATIONS: usize = 1000;
const SEARCH_WINDOWS_SIZE: usize = 50;

const MV_CONC: f32 = 50.0; // Monovalent cation concentration (mM)
const DV_CONC: f32 = 3.0; // Divalent cation concentration (mM)
const DNTP_CONC: f32 = 0.0; // dNTP concentration (mM)
const DNA_CONC: f32 = 250.0; // Primer concentration (nM)
const ANNEALING_TEMP: f32 = 25.0; // Annealing temperature (°C)

const PRIMER_MIN_TM: f32 = 30.0;
const PRIMER_MAX_TM: f32 = 60.0;
const PRIMER_MAX_SELF_ANY_TH: f32 = DEFAULT_MIN_TM - 10.0;
const PRIMER_MAX_SELF_END_TH: f32 = DEFAULT_MIN_TM - 10.0;
const PRIMER_MAX_HAIRPIN_TH: f32 = DEFAULT_MIN_TM - 10.0;

struct SequenceRecord {
    name: String,
    sequence: String,
}

const SEQ_DIR_FWD: u8 = 0x00;
const SEQ_DIR_REV: u8 = 0x01;

#[derive(Clone)]
struct KmerRecord {
    word: String,
    direction: u8,
}

impl PartialEq for KmerRecord {
    fn eq(&self, other: &Self) -> bool {
        self.word == other.word && self.direction == other.direction
    }
}

impl Hash for KmerRecord {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.word.hash(state);
        self.direction.hash(state);
    }
}

impl Eq for KmerRecord {}

#[derive(Clone)]
struct KmerFrequency<'a> {
    kmer: &'a KmerRecord,
    frequency: usize,
}

impl PartialEq for KmerFrequency<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.kmer.word == other.kmer.word
    }
}

impl Eq for KmerFrequency<'_> {}

impl Hash for KmerFrequency<'_> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.kmer.word.hash(state);
    }
}

#[derive(Clone)]
struct KmerStat {
    word: String,
    direction: u8,
    gc_percent: f32,
    mean: f32,
    std: f32,
    tm: f32,
    tm_ok: bool,
    self_any_th: f32,
    self_end_th: f32,
    hairpin_th: f32,
    runs: bool,
}

struct Segment<'a> {
    sequence: &'a SequenceRecord,
    partition_no: u16,
    index: usize,
    kmers: [Vec<KmerRecord>; 2],
}

struct SegmentManager<'a> {
    segments: Vec<Segment<'a>>,
}

impl Hash for Segment<'_> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.name.hash(state);
        self.partition_no.hash(state);
    }
}

impl PartialEq for Segment<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.sequence.name == other.sequence.name && self.partition_no == other.partition_no
    }
}

impl Eq for Segment<'_> {}

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
            "--op",
            "1.53",
            "--ep",
            "0.123",
            "--jtt",
            "200",
            &filepath.clone(),
        ])
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
        .filter(|kmer| kmer.iter().all(|c| "ATCGU".contains(*c)))
        .map(|kmer| kmer.iter().collect())
        .unique()
        .collect()
}

fn partitioning_sequence(sequence: &str, size: usize, overlap_size: usize) -> Vec<String> {
    sequence
        .chars()
        .collect::<Vec<char>>()
        .windows(size)
        .step_by(overlap_size)
        .map(|window| window.iter().collect())
        .collect()
}

fn get_sequence_on_search_windows(
    sequence: &str,
    search_windows_size: usize,
) -> (String, String) {
    let first = &sequence[..search_windows_size];
    let second = &sequence[sequence.len() - search_windows_size..];
    (first.to_string(), second.to_string())
}

struct PartitioningOption {
    segment_size: usize,
    overlap_size: usize,
    window_size: usize,
    kmer_size: usize,
}

fn get_segment_manager(records: &[SequenceRecord], opt: PartitioningOption) -> SegmentManager {
    let mut manager = SegmentManager {
        segments: Vec::new(),
    };

    if opt.overlap_size < opt.window_size {
        panic!("Overlap windows size must be greater or equal than search windows size");
    }

    for record in records.iter() {
        let partitions =
            partitioning_sequence(&record.sequence, opt.segment_size, opt.overlap_size);
        for (j, partition) in partitions.iter().enumerate() {
            let (start, end) = get_sequence_on_search_windows(partition, opt.window_size);
            let start_kmers = find_kmers(&start, opt.kmer_size);
            let end_kmers = find_kmers(&end, opt.kmer_size);
            let mut kmers: [Vec<KmerRecord>; 2] = [Vec::new(), Vec::new()];
            for kmer in start_kmers.iter() {
                kmers[0].push(KmerRecord {
                    word: kmer.clone(),
                    direction: SEQ_DIR_FWD,
                });
            }
            for kmer in end_kmers.iter() {
                kmers[1].push(KmerRecord {
                    word: reverse_complement(kmer),
                    direction: SEQ_DIR_REV,
                });
            }
            manager.segments.push(Segment {
                sequence: record,
                partition_no: j as u16,
                index: manager.segments.len(),
                kmers,
            });
        }
    }

    manager
}

fn make_kmer_segments_windows_mapping<'a>(
    segments: &'a Vec<Segment<'a>>,
) -> HashMap<&'a KmerRecord, Vec<u32>> {
    let mut kmer_segments_mapping: HashMap<&'a KmerRecord, Vec<u32>> = HashMap::new();
    for segment in segments.iter() {
        for (direction, kmers) in segment.kmers.iter().enumerate() {
            for kmer in kmers.iter() {
                if kmer.direction != direction as u8 {
                    continue;
                }
                kmer_segments_mapping
                    .entry(kmer)
                    .or_default()
                    .push(segment.index as u32);
            }
        }
    }
    kmer_segments_mapping
}

fn find_most_freq_kmer<'a>(
    segments: &'a Vec<Segment>,
    direction: u8,
    ignored_segments_windows: HashSet<u32>,
) -> Option<KmerFrequency<'a>> {
    let mut kmer_freq_map: HashMap<&KmerRecord, usize> = HashMap::new();

    for (idx, segment) in segments.iter().enumerate() {
        for (window_direction, kmers) in segment.kmers.iter().enumerate() {
            if window_direction != direction as usize {
                continue;
            }
            let key = idx as u32;
            if ignored_segments_windows.contains(&key) {
                continue;
            }
            for kmer in kmers.iter() {
                kmer_freq_map
                    .entry(kmer)
                    .and_modify(|f| *f += 1)
                    .or_insert(1);
            }
        }
    }

    kmer_freq_map
        .iter()
        .max_by_key(|(_, &v)| v)
        .map(|(k, &f)| KmerFrequency {
            kmer: k,
            frequency: f,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let sequence = "ATCGAA";
        assert_eq!(reverse_complement(sequence), "TTCGAT");
    }

    #[test]
    fn test_get_search_windows() {
        let sequence = "AACCTTGGAACCTTG-".to_string();
        let (first, second) = get_sequence_on_search_windows(&sequence, 5);
        assert_eq!(first, "AACCT");
        assert_eq!(second, "CTTG-");
    }

    #[test]
    fn test_find_kmers() {
        let sequence = "AACCTTGGAACCTTG-";
        let kmers = find_kmers(sequence, 5);
        assert_eq!(kmers.len(), 8);
        assert!(kmers.contains(&"AACCT".to_string()));
        assert!(kmers.contains(&"ACCTT".to_string()));
        assert!(kmers.contains(&"CCTTG".to_string()));
        assert!(kmers.contains(&"CTTGG".to_string()));
        assert!(kmers.contains(&"TTGGA".to_string()));
        assert!(kmers.contains(&"TGGAA".to_string()));
        assert!(kmers.contains(&"GGAAC".to_string()));
        assert!(kmers.contains(&"GAACC".to_string()));
    }

    #[test]
    fn test_get_segments() {
        // search windows = (AACCT)(TGGAA)
        // AAC, freq=2
        // ACC, freq=3
        // CCT, freq=3
        // reverse complement = AAGGT -> TTCCA
        // TTC, freq=3
        // TCC, freq=3
        // CCA, freq=3
        let records = vec![
            SequenceRecord {
                name: "seq1".to_string(),
                sequence: "AACCTTGGAACCTTGG".to_string(),
            },
            SequenceRecord {
                name: "seq2".to_string(),
                sequence: "AACCTTGGAACCTTG-".to_string(),
            },
            SequenceRecord {
                name: "seq3".to_string(),
                sequence: "-ACCTTGGAACCTT-G".to_string(),
            },
        ];
        for rec in records.iter() {
            println!("seq={} ---", rec.name);
            for p in partitioning_sequence(&rec.sequence, 10, 5) {
                println!("-- partition={}", p);
                let (i, j) = get_sequence_on_search_windows(&p, 5);
                println!("---- search_windows(start)={}", i);
                println!("---- search_windows(end)={}", j);
                let rev = reverse_complement(&j);
                println!("---- search_windows(end/rev)={}", rev);
            }
        }
        let opt = PartitioningOption {
            segment_size: 10,
            overlap_size: 5,
            window_size: 5,
            kmer_size: 3,
        };
        let manager = get_segment_manager(&records, opt);
        assert_eq!(manager.segments.len(), 6);
        let first_segment = manager.segments.first().unwrap();
        for kmer in first_segment.kmers.iter() {
            for k in kmer.iter() {
                println!("kmer={}", k.word);
            }
        }
        assert_eq!(manager.segments.first().unwrap().kmers[0].len(), 3);
        assert_eq!(manager.segments.get(1).unwrap().kmers[1].len(), 3);
    }

    #[test]
    fn test_partitioning_sequence() {
        let sequence = "AACCTTGGAACCTTGG";
        let partitions = partitioning_sequence(sequence, 10, 5);
        assert_eq!(partitions.len(), 2);
        assert_eq!(partitions[0], "AACCTTGGAA");
        assert_eq!(partitions[1], "TGGAACCTTG");
    }

    #[test]
    fn test_make_kmer_segments_mapping() {
        let seq_1 = &SequenceRecord {
            name: "seq1".to_string(),
            sequence: "ACTGAGGTATTA".to_string(),
        };
        let seq_2 = &SequenceRecord {
            name: "seq2".to_string(),
            sequence: "ACTGAGGTGGAA".to_string(),
        };
        let manager: SegmentManager = SegmentManager {
            segments: vec![
                Segment {
                    sequence: seq_1,
                    partition_no: 0,
                    index: 0,
                    kmers: [
                        vec![
                            KmerRecord {
                                word: "ACT".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "CTG".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "TGA".to_string(),
                                direction: 0,
                            },
                        ],
                        vec![
                            KmerRecord {
                                word: "TAA".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "AAT".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "ATA".to_string(),
                                direction: 1,
                            },
                        ],
                    ],
                },
                Segment {
                    sequence: seq_2,
                    partition_no: 1,
                    index: 1,
                    kmers: [
                        vec![
                            KmerRecord {
                                word: "ACT".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "CTG".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "TGA".to_string(),
                                direction: 0,
                            },
                        ],
                        vec![
                            KmerRecord {
                                word: "TTC".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "TCC".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "CCA".to_string(),
                                direction: 1,
                            },
                        ],
                    ],
                },
            ],
        };

        let kmer_segments_mapping = make_kmer_segments_windows_mapping(&manager.segments);
        assert_eq!(kmer_segments_mapping.len(), 9);
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "ACT".to_string(),
                    direction: 0
                })
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "CTG".to_string(),
                    direction: 0
                })
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "TGA".to_string(),
                    direction: 0
                })
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "TAA".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "AAT".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "ATA".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "TTC".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "TCC".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            kmer_segments_mapping
                .get(&KmerRecord {
                    word: "CCA".to_string(),
                    direction: 1
                })
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn test_find_most_freq_kmer() {
        let seq_1 = SequenceRecord {
            name: "seq1".to_string(),
            sequence: "ACTGAGGTATTA".to_string(),
        };
        let seq_2 = SequenceRecord {
            name: "seq2".to_string(),
            sequence: "ACAGGGGTGGAA".to_string(),
        };
        let manager: SegmentManager = SegmentManager {
            segments: vec![
                Segment {
                    sequence: &seq_1,
                    partition_no: 0,
                    index: 0,
                    kmers: [
                        vec![
                            KmerRecord {
                                word: "ACT".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "CTG".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "TGA".to_string(),
                                direction: 0,
                            },
                        ],
                        vec![
                            KmerRecord {
                                word: "TAA".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "AAT".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "ATA".to_string(),
                                direction: 1,
                            },
                        ],
                    ],
                },
                Segment {
                    sequence: &seq_2,
                    partition_no: 1,
                    index: 1,
                    kmers: [
                        vec![
                            KmerRecord {
                                word: "ACT".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "CAG".to_string(),
                                direction: 0,
                            },
                            KmerRecord {
                                word: "TGG".to_string(),
                                direction: 0,
                            },
                        ],
                        vec![
                            KmerRecord {
                                word: "TTC".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "TCC".to_string(),
                                direction: 1,
                            },
                            KmerRecord {
                                word: "CCA".to_string(),
                                direction: 1,
                            },
                        ],
                    ],
                },
            ],
        };

        let result = find_most_freq_kmer(&manager.segments, 0, HashSet::new());
        assert!(result.is_some());
        let kmer_freq = result.unwrap();
        assert_eq!(kmer_freq.kmer.word, "ACT");
        assert_eq!(kmer_freq.frequency, 2);
    }
}

fn find_candidates_kmers<'a>(
    segment_manager: &'a SegmentManager,
    direction: u8,
) -> Option<Vec<KmerFrequency<'a>>> {
    let mut candidate_kmers: Vec<KmerFrequency> = Vec::new();
    let kmer_segments_windows_mappings =
        make_kmer_segments_windows_mapping(&segment_manager.segments);
    let mut ignored_segments_windows: HashSet<u32> = HashSet::new();

    for iter_no in 0..MAX_ITERATIONS {
        log::trace!("Iteration: {}", iter_no + 1);
        let kmer_freq = match find_most_freq_kmer(
            &segment_manager.segments,
            direction,
            ignored_segments_windows.clone(),
        ) {
            Some(k) => {
                if k.frequency == 1 {
                    log::trace!(
                        "Iteration: {}, only 1 shared window found, stop ...",
                        iter_no
                    );
                    break;
                }
                k
            }
            None => {
                log::trace!("Iteration: {}, no k-mers found, stop ...", iter_no);
                break;
            }
        };
        candidate_kmers.push(kmer_freq.clone());

        // update ignored segments
        let mut count = 0;
        for idx in kmer_segments_windows_mappings.get(&kmer_freq.kmer).unwrap() {
            count += 1;
            ignored_segments_windows.insert(*idx);
        }
        log::debug!(
            "Iteration: {}, direction: {} winner: {}, windows removed: {}, total removed: {}",
            iter_no,
            direction,
            kmer_freq.kmer.word,
            count,
            ignored_segments_windows.len()
        );

        if kmer_freq.frequency < MAX_MISMATCH_SEGMENTS {
            log::info!("Max mismatch segments reached, exiting...");
            break;
        }
    }

    if candidate_kmers.is_empty() {
        return None;
    }

    Some(
        candidate_kmers
            .into_iter()
            .map(|k| KmerFrequency {
                kmer: k.kmer,
                frequency: k.frequency,
            })
            .collect(),
    )
}

fn get_kmer_stats(kmer_records: Vec<KmerFrequency>) -> Vec<KmerStat> {
    // first, finding the threshold for Tm
    let primers: Vec<String> = kmer_records.iter().map(|k| k.kmer.word.clone()).collect();

    let params = CheckPrimerParams {
        min_tm: DEFAULT_MIN_TM,
        max_tm: DEFAULT_MAX_TM,
    };
    let check_primers_result = check_primers(&primers, params);
    if check_primers_result.is_err() {
        panic!("Error while checking primers");
    }

    let primer_info_list = check_primers_result.unwrap();
    let mut primer_info_map = HashMap::new();
    for info in &primer_info_list {
        primer_info_map.entry(info.id).or_insert(info);
    }
    let (mean, std) = get_tm_stat(&primer_info_list);

    kmer_records
        .iter()
        .map(|kmer_freq| {
            let primer_info = match primer_info_map.get(&kmer_freq.kmer.word.as_str()) {
                Some(info) => *info,
                _ => &PrimerInfo::new(),
            };
            KmerStat {
                word: kmer_freq.kmer.word.clone(),
                direction: kmer_freq.kmer.direction,
                mean,
                std,
                gc_percent: primer_info.gc,
                tm: primer_info.tm,
                tm_ok: tm_in_threshold(primer_info.tm, mean, std),
                self_any_th: primer_info.self_any_th,
                self_end_th: primer_info.self_end_th,
                hairpin_th: primer_info.hairpin_th,
                runs: is_run(kmer_freq.kmer.word.clone()),
            }
        })
        .collect()
}

/**
 * Get the threshold for Tm
 *
 * Calculated by find (2*sd(Tm)) + mean(Tm) of the primers
 */
fn get_tm_stat(primer_info_list: &Vec<PrimerInfo>) -> (f32, f32) {
    let tm_values: Vec<f32> = primer_info_list.iter().map(|info| info.tm).collect();
    let mean = tm_values.iter().sum::<f32>() / tm_values.len() as f32;
    let std = standard_deviation(&tm_values);
    (mean, std.standard_deviation)
}

fn tm_in_threshold(tm: f32, mean: f32, std: f32) -> bool {
    (tm - mean).abs() <= (2.0 * std)
}

/**
 * Check if the kmer  has >=5nt di-nucleotide repeats
 *
 * For example, ATATATATATGG is too many AT repeats, then return `true`.
 */
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

fn filter_kmers(stats: Vec<KmerStat>) -> Vec<KmerStat> {
    stats
        .iter()
        .filter(|kmer_stat| {
            (kmer_stat.tm > PRIMER_MIN_TM && kmer_stat.tm < PRIMER_MAX_TM)
                && kmer_stat.self_any_th < PRIMER_MAX_SELF_ANY_TH
                && kmer_stat.self_end_th < PRIMER_MAX_SELF_END_TH
                && kmer_stat.hairpin_th < PRIMER_MAX_HAIRPIN_TH
                && kmer_stat.tm_ok
                && !kmer_stat.runs
        }).cloned()
        .collect()
}
#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long)]
    output: String,
}

fn main() -> io::Result<()> {
    env_logger::init();

    let args = Args::parse();
    let filename = args.input.to_string();

    // 1. Align sequences
    log::info!("Aligning sequences...");
    let records = match align_sequences(filename) {
        Ok(records) => to_records(records)?,
        Err(e) => {
            panic!("Error aligning sequences: {}", e);
        }
    };
    log::info!("Done aligning sequences");
    if records.iter().len() == 0 {
        panic!("No sequences found in the input file");
    }

    // 2. Extracting n-grams from each sequence segments
    log::info!("Extracting n-grams from each sequence segments...");
    let options = PartitioningOption {
        segment_size: WINDOW_SIZE,
        overlap_size: OVERLAP_SIZE,
        window_size: SEARCH_WINDOWS_SIZE,
        kmer_size: KMER_SIZE,
    };
    let segment_manager = get_segment_manager(&records, options);
    let total_partitions = segment_manager
        .segments
        .iter()
        .max_by_key(|s| s.partition_no)
        .unwrap()
        .partition_no;
    log::info!(
        "Done, total partitions: {}, total segments: {}",
        total_partitions,
        segment_manager.segments.len()
    );

    // 3. Calculate frequencies of n-grams for each segment both forward/reverse
    log::info!("Calculating frequencies of k-mer for all segments...");
    log::debug!("Total segments: {}", segment_manager.segments.len());
    let candidate_kmers_fwd =
        find_candidates_kmers(&segment_manager, SEQ_DIR_FWD).unwrap_or_default();
    let candidate_kmers_rev =
        find_candidates_kmers(&segment_manager, SEQ_DIR_REV).unwrap_or_default();
    log::info!(
        "Done calculating, Total candidate k-mers: fwd: {}, rev: {}",
        candidate_kmers_fwd.len(),
        candidate_kmers_rev.len()
    );

    // 4. Filtering out unmatched criteria
    log::info!("Filtering out unmatched criteria (Tm and >5nt repeats, runs...)");
    let kmer_stats_fwd = get_kmer_stats(candidate_kmers_fwd);
    let kmer_stats_rev = get_kmer_stats(candidate_kmers_rev);

    let candidate_primers_fwd: Vec<KmerStat> = filter_kmers(kmer_stats_fwd);
    let candidate_primers_rev: Vec<KmerStat> = filter_kmers(kmer_stats_rev);

    let primers: Vec<String> = candidate_primers_fwd
        .iter()
        .chain(&candidate_primers_rev)
        .map(|s| s.word.clone())
        .collect();
    let opts = NtthalOptions {
        mv: MV_CONC,
        dv: DV_CONC,
        dntp: DNTP_CONC,
        conc: DNA_CONC,
        t: ANNEALING_TEMP,
    };
    let graph = run_ntthal(primers.clone(), opts)?;
    let mut candidate_unusable_edges: Vec<&Edge> = Vec::new();
    let mut primers_total_low_dg: HashMap<String, i32> = HashMap::new();
    // find nodes with dG < -9.0kmol-1
    for primer in primers.clone() {
        let edges = graph.get_edges_for_node(&primer);
        for edge in edges {
            let dg = edge.get_dg();
            if dg < -9000.0 {
                candidate_unusable_edges.push(edge);
                let (a, b) = graph.get_edge_nodes(edge);
                log::debug!("Edge: {} -> {} dg={}", a.id, b.id, dg);
                *primers_total_low_dg.entry(a.id.clone()).or_insert(0) += 1;
                *primers_total_low_dg.entry(b.id.clone()).or_insert(0) += 1;
            }
        }
    }
    let mut deleted_primers: HashSet<String> = HashSet::new();
    for edge in candidate_unusable_edges {
        let (a, b) = graph.get_edge_nodes(edge);
        let total_n_a = primers_total_low_dg.get(&a.id).unwrap();
        let total_n_b = primers_total_low_dg.get(&b.id).unwrap();
        if total_n_a > &1 {
            deleted_primers.insert(a.id.clone());
        }
        if total_n_b > &1 {
            deleted_primers.insert(b.id.clone());
        }
        if total_n_a == &1 && !deleted_primers.contains(b.id.as_str()) {
            deleted_primers.insert(a.id.clone());
        }
        if total_n_b == &1 && !deleted_primers.contains(a.id.as_str()) {
            deleted_primers.insert(b.id.clone());
        }
    }
    log::debug!("Deleted primers: {:?}", deleted_primers);
    let good_delta_g_fwd_primers: Vec<KmerStat> = candidate_primers_fwd
        .iter()
        .filter(|p| !deleted_primers.contains(p.word.as_str())).cloned()
        .collect();
    let good_delta_g_rev_primers: Vec<KmerStat> = candidate_primers_rev
        .iter()
        .filter(|p| !deleted_primers.contains(p.word.as_str())).cloned()
        .collect();
    log::info!(
        "Done filtering out unmatched, primers left fwd={}, rev={}",
        good_delta_g_fwd_primers.len(),
        good_delta_g_rev_primers.len()
    );

    // 5. Output the primers
    log::info!("Outputting primers...");
    let output_file = args.output.to_string();
    let mut writer = csv::Writer::from_path(output_file)?;
    writer.write_record(["direction", "name", "primers", "gc", "avg", "std", "tm"])?;
    let candidate_primers = vec![good_delta_g_fwd_primers, good_delta_g_rev_primers];
    for candidates in candidate_primers {
        for (idx, primer) in candidates.iter().enumerate() {
            let direction = if primer.direction == SEQ_DIR_FWD {
                "F"
            } else {
                "R"
            };
            writer.write_record([
                direction,
                &format!("Primer_{}_{}", idx, direction),
                &*primer.word,
                &format!("{:.2}", primer.gc_percent / 100.0),
                &format!("{:.2}", primer.mean),
                &format!("{:.2}", primer.std),
                &format!("{:.2}", primer.tm),
            ])?;
        }
    }
    writer.flush().expect("Error writing output to file");
    log::info!("Done outputting primers");

    Ok(())
}
