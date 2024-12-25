use ngrams::Ngram;
use std_dev::standard_deviation;

use crate::{
    config::Config,
    delta_g,
    segment::{KmerFrequency, KmerStat},
};

pub fn filter(config: &Config, candidate: Vec<KmerFrequency>) -> Vec<KmerStat> {
    let stats = get_kmer_stats(config, candidate);
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
