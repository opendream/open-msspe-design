use ngrams::Ngram;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

struct FastaSequence {
    name: String,
    sequence: String,
}

fn read_fasta_file(filename: &str) -> io::Result<Vec<FastaSequence>> {
    let path = Path::new(filename);
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_name.is_empty() {
                sequences.push(FastaSequence {
                    name: current_name.clone(),
                    sequence: current_sequence.clone(),
                });
                current_sequence.clear();
            }
            current_name = line.replace(">", "").trim().to_string();
        } else {
            current_sequence.push_str(&line);
        }
    }

    if !current_name.is_empty() {
        sequences.push(FastaSequence {
            name: current_name,
            sequence: current_sequence,
        });
    }

    Ok(sequences)
}

const KMER_SIZE: usize = 13;

fn main() -> io::Result<()> {
    let filename = "samples.fasta";
    let sequences = read_fasta_file(filename)?;

    let mut seq_grams = std::collections::HashMap::new();
    for sequence in sequences.iter() {
        let filtered_grams: Vec<_> = sequence.sequence.chars().ngrams(KMER_SIZE)
            .filter(|gram| {
                gram.iter()
                    .all(|&c| c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U')
            })
            .collect();
        log::debug!(
            "{}: orig-length: {}, filtered: {}",
            sequence.name,
            sequence.sequence.len(),
            filtered_grams.len()
        );
        seq_grams.insert(sequence.name.clone(), filtered_grams);
    }

    // Find common n-grams between sequences ordered by frequency
    let mut freqs = std::collections::HashMap::new();
    for (_name, grams) in seq_grams.iter() {
        for gram in grams {
            let entry = freqs.entry(gram).or_insert(0);
            *entry += 1;
        }
    }
    println!("Total n-grams: {}", freqs.len());

    let mut common_grams = std::collections::HashMap::new();
    for sequence in sequences.iter() {
        let mut local_freqs = freqs.clone();
        let grams = seq_grams.get(&sequence.name).unwrap();
        for gram in grams {
            if let Some(count) = local_freqs.get_mut(gram) {
                *count -= 1;
                if *count == 0 {
                    local_freqs.remove(gram);
                }
            }
        }
        for (gram, count) in local_freqs {
            if count > 0 {
                let entry = common_grams.entry(gram).or_insert(0);
                *entry += 1;
            }
        }
    }

    // Print common grams ordered by frequency
    println!("Common n-grams: {}", common_grams.len());
    let mut common_grams: Vec<_> = common_grams.iter().collect();
    common_grams.truncate(100);
    common_grams.sort_by(|a, b| b.1.cmp(a.1));
    for (gram, count) in common_grams {
        println!("{}: {}", gram.iter().collect::<String>(), count);
    }

    Ok(())
}
