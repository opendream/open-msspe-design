use itertools::Itertools;
use ngrams::Ngram;
use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
};

use crate::{
    config::{Config, PartitioningOption},
    sequence::SequenceRecord,
};

pub const SEQ_DIR_FWD: u8 = 0x00;
pub const SEQ_DIR_REV: u8 = 0x01;

pub struct Segment<'a> {
    pub sequence: &'a SequenceRecord,
    pub partition_no: u16,
    pub index: usize,
    pub kmers: [Vec<KmerRecord>; 2],
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

#[derive(Hash, Clone)]
pub struct KmerRecord {
    pub word: String,
    pub direction: u8,
}

impl PartialEq for KmerRecord {
    fn eq(&self, other: &Self) -> bool {
        self.word == other.word && self.direction == other.direction
    }
}

impl Eq for KmerRecord {}

#[derive(Clone)]
pub struct KmerFrequency<'a> {
    pub kmer: &'a KmerRecord,
    pub frequency: usize,
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

pub struct KmerStat {
    pub word: String,
    pub direction: u8,
    pub gc_percent: f32,
    pub mean: f32,
    pub std: f32,
    pub tm: f32,
    pub tm_ok: bool,
    pub repeats: bool,
    pub runs: bool,
    pub delta_g: f32,
}

pub struct SegmentManager<'a> {
    pub config: &'a Config,
    pub segments: Vec<Segment<'a>>,
}

impl<'a> SegmentManager<'a> {
    pub fn new(records: &'a Vec<SequenceRecord>, config: &'a Config) -> Self {
        SegmentManager {
            config,
            segments: get_segments(records, PartitioningOption::new(&config)),
        }
    }

    // get total partitions
    pub fn number_of_partitions(&self) -> u16 {
        let total = self
            .segments
            .iter()
            .max_by_key(|s| s.partition_no)
            .unwrap()
            .partition_no;
        total
    }

    pub fn number_of_segments(&self) -> usize {
        self.segments.len()
    }

    pub fn find_candidates_kmers(&self, direction: u8) -> Option<Vec<KmerFrequency<'a>>> {
        let mut candidate_kmers: Vec<KmerFrequency> = Vec::new();
        let kmer_segments_windows_mappings = make_kmer_segments_windows_mapping(&self.segments);
        let mut ignored_segments_windows: HashSet<u32> = HashSet::new();
        for iter_no in 0..self.config.max_iterations() {
            log::trace!("Iteration: {}", iter_no + 1);
            let kmer_freq = match find_most_freq_kmer(
                &self.segments,
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

            if kmer_freq.frequency < self.config.max_mismatch_segments() {
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
}

fn find_most_freq_kmer<'a>(
    segments: &Vec<Segment<'a>>,
    direction: u8,
    ignored_segments_windows: HashSet<u32>,
) -> Option<KmerFrequency<'a>> {
    let mut kmer_freq_map: HashMap<KmerRecord, usize> = HashMap::new();

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
                    .entry(kmer.clone())
                    .and_modify(|f| *f += 1)
                    .or_insert(1);
            }
        }
    }

    kmer_freq_map
        .into_iter()
        .max_by_key(|(_, v)| *v)
        .map(|(k, f)| KmerFrequency {
            kmer: Box::leak(Box::new(k)),
            frequency: f,
        })
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

fn get_segments(records: &Vec<SequenceRecord>, opt: PartitioningOption) -> Vec<Segment> {
    let mut segments = Vec::new();
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
            segments.push(Segment {
                sequence: record,
                partition_no: j as u16,
                index: segments.len(),
                kmers,
            });
        }
    }
    segments
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

fn find_kmers(sequence: &str, kmer_size: usize) -> Vec<String> {
    sequence
        .chars()
        .ngrams(kmer_size)
        .filter(|kmer| kmer.iter().all(|c| "ATCGU".contains(*c)))
        .map(|kmer| kmer.iter().collect())
        .unique()
        .collect()
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

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_partitioning_sequence() {
        let sequence = "AACCTTGGAACCTTGG";
        let partitions = partitioning_sequence(sequence, 10, 5);
        assert_eq!(partitions.len(), 2);
        assert_eq!(partitions[0], "AACCTTGGAA");
        assert_eq!(partitions[1], "TGGAACCTTG");
    }

    #[test]
    fn test_reverse_complement() {
        let sequence = "ATCGAA";
        assert_eq!(reverse_complement(sequence), "TTCGAT");
    }

    #[test]
    fn test_get_segments() {
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
        let config = Config::new_for_test(10, 5, 5, 3);
        let manager = SegmentManager::new(&records, &config);
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
            config: &Config::default(),
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
            config: &Config::default(),
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
