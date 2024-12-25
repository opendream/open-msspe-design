mod align;
mod delta_g;
mod types;

use clap::Parser;
use itertools::Itertools;
use ngrams::Ngram;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::{self};
use std_dev::standard_deviation;
use types::{Args, Config, PartitioningOption, SequenceRecord};

const SEQ_DIR_FWD: u8 = 0x00;
const SEQ_DIR_REV: u8 = 0x01;

#[derive(Hash, Clone)]
struct KmerRecord {
    word: String,
    direction: u8,
}

impl PartialEq for KmerRecord {
    fn eq(&self, other: &Self) -> bool {
        self.word == other.word && self.direction == other.direction
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

struct KmerStat {
    word: String,
    direction: u8,
    gc_percent: f32,
    mean: f32,
    std: f32,
    tm: f32,
    tm_ok: bool,
    repeats: bool,
    runs: bool,
    delta_g: f32,
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
    sequence: &String,
    search_windows_size: usize,
) -> (String, String) {
    let first = &sequence[..search_windows_size];
    let second = &sequence[sequence.len() - search_windows_size..];
    (first.to_string(), second.to_string())
}

fn get_segment_manager(records: &Vec<SequenceRecord>, opt: PartitioningOption) -> SegmentManager {
    let mut manager = SegmentManager {
        segments: Vec::new(),
    };

    for (_, record) in records.iter().enumerate() {
        let partitions =
            partitioning_sequence(&record.sequence, opt.segment_size(), opt.overlap_size());
        for (j, partition) in partitions.iter().enumerate() {
            let (start, end) = get_sequence_on_search_windows(&partition, opt.window_size());
            let start_kmers = find_kmers(&start, opt.kmer_size());
            let end_kmers = find_kmers(&end, opt.kmer_size());
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
                    .or_insert(Vec::new())
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
            kmer: *k,
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
        let opt = PartitioningOption::new_with(10, 5, 5, 3);
        let manager = get_segment_manager(&records, opt);
        assert_eq!(manager.segments.len(), 6);
        let first_segment = manager.segments.get(0).unwrap();
        for kmer in first_segment.kmers.iter() {
            for k in kmer.iter() {
                println!("kmer={}", k.word);
            }
        }
        assert_eq!(manager.segments.get(0).unwrap().kmers[0].len(), 3);
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
        let seq_2 = SequenceRecord {
            name: "seq2".to_string(),
            sequence: "ACTGAGGTGGAA".to_string(),
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
        assert_eq!(result.is_some(), true);
        let kmer_freq = result.unwrap();
        assert_eq!(kmer_freq.kmer.word, "ACT");
        assert_eq!(kmer_freq.frequency, 2);
    }
}

fn find_candidates_kmers<'a>(
    config: &Config,
    segment_manager: &'a SegmentManager,
    direction: u8,
) -> Option<Vec<KmerFrequency<'a>>> {
    let mut candidate_kmers: Vec<KmerFrequency> = Vec::new();
    let kmer_segments_windows_mappings =
        make_kmer_segments_windows_mapping(&segment_manager.segments);
    let mut ignored_segments_windows: HashSet<u32> = HashSet::new();

    for iter_no in 0..config.max_iterations() {
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

        if kmer_freq.frequency < config.max_mismatch_segments() {
            log::info!("Max mismatch segments reached, exiting...");
            break;
        }
    }

    if candidate_kmers.len() == 0 {
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

fn get_kmer_stats(config: &Config, kmer_records: Vec<KmerFrequency>) -> Vec<KmerStat> {
    // first, finding the threshold for Tm
    let primers: Vec<String> = kmer_records.iter().map(|k| k.kmer.word.clone()).collect();
    let (mean, std) = get_tm_stat(primers);

    kmer_records
        .iter()
        .map(|kmer_freq| {
            let tm = get_tm(kmer_freq.kmer.word.clone());
            let delta_g = delta_g::calculate_delta_g(
                &kmer_freq.kmer.word.clone(),
                config.mv_conc(),
                config.dv_conc(),
                config.dntp_conc(),
                config.dna_conc(),
                config.annealing_temp(),
            )
            .unwrap_or_else(|| 0.0);
            KmerStat {
                word: kmer_freq.kmer.word.clone(),
                direction: kmer_freq.kmer.direction,
                gc_percent: get_gc_percent(kmer_freq.kmer.word.clone()),
                mean,
                std,
                tm,
                tm_ok: in_tm_threshold(kmer_freq.kmer.word.clone(), mean, std),
                repeats: is_repeats(kmer_freq.kmer.word.clone()),
                runs: is_run(kmer_freq.kmer.word.clone()),
                delta_g: delta_g as f32,
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
fn get_tm_stat(primers: Vec<String>) -> (f32, f32) {
    let mut tm_values: Vec<f32> = Vec::new();
    for primer in primers {
        tm_values.push(get_tm(primer));
    }
    let mean_tm = tm_values.iter().sum::<f32>() / tm_values.len() as f32;
    let sd_tm = standard_deviation(&tm_values);
    (mean_tm, sd_tm.standard_deviation)
}

fn in_tm_threshold(kmer: String, mean: f32, std: f32) -> bool {
    let tm_upper_bound = mean + (2.0 * std);
    let tm_lower_bound = mean - (2.0 * std);

    let tm = get_tm(kmer);

    tm <= tm_upper_bound && tm >= tm_lower_bound
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

fn filter_kmers(stats: Vec<KmerStat>) -> Vec<KmerStat> {
    stats
        .iter()
        .filter(|kmer_stat| kmer_stat.tm_ok && kmer_stat.delta_g > -9.0 && !kmer_stat.runs)
        .map(|kmer_stat| KmerStat {
            word: kmer_stat.word.clone(),
            direction: kmer_stat.direction,
            gc_percent: kmer_stat.gc_percent,
            mean: kmer_stat.mean,
            std: kmer_stat.std,
            tm: kmer_stat.tm,
            tm_ok: kmer_stat.tm_ok,
            repeats: kmer_stat.repeats,
            runs: kmer_stat.runs,
            delta_g: kmer_stat.delta_g,
        })
        .collect()
}

fn main() -> io::Result<()> {
    env_logger::init();

    let args = Args::parse();
    let config: Config = (&args).into();

    // 1. Align sequences
    log::info!("Aligning sequences...");
    let records = align::align_sequences(args.input.to_string());
    log::info!("Done aligning sequences");

    // 2. Extracting n-grams from each sequence segments
    log::info!("Extracting n-grams from each sequence segments...");
    let segment_manager = get_segment_manager(&records, PartitioningOption::new(&config));
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
        find_candidates_kmers(&config, &segment_manager, SEQ_DIR_FWD).unwrap_or_else(|| Vec::new());
    let candidate_kmers_rev =
        find_candidates_kmers(&config, &segment_manager, SEQ_DIR_REV).unwrap_or_else(|| Vec::new());
    log::info!(
        "Done calculating, Total candidate k-mers: fwd: {}, rev: {}",
        candidate_kmers_fwd.len(),
        candidate_kmers_rev.len()
    );

    // 4. Filtering out unmatched criteria
    log::info!("Filtering out unmatched criteria (Tm and >5nt repeats, runs...)");
    let kmer_stats_fwd = get_kmer_stats(&config, candidate_kmers_fwd);
    let kmer_stats_rev = get_kmer_stats(&config, candidate_kmers_rev);
    let candidate_primers_fwd: Vec<KmerStat> = filter_kmers(kmer_stats_fwd);
    let candidate_primers_rev: Vec<KmerStat> = filter_kmers(kmer_stats_rev);
    log::info!(
        "Done filtering out unmatched, primers left fwd={}, rev={}",
        candidate_primers_fwd.len(),
        candidate_primers_rev.len()
    );

    // 5. Output the primers
    log::info!("Outputting primers...");
    let output_file = args.output.to_string();
    let mut writer = csv::Writer::from_path(output_file)?;
    writer.write_record(&["direction", "name", "primers", "gc", "avg", "std", "tm"])?;
    let candidate_primers = vec![candidate_primers_fwd, candidate_primers_rev];
    for candidates in candidate_primers {
        for (idx, primer) in candidates.iter().enumerate() {
            let direction = if primer.direction == SEQ_DIR_FWD {
                "F"
            } else {
                "R"
            };
            writer.write_record(&[
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
