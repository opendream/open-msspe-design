mod config;
mod delta_g;
mod filter;
mod segment;
mod sequence;

use clap::Parser;
use config::{Args, Config};
use filter::filter;
use segment::{SegmentManager, SEQ_DIR_FWD, SEQ_DIR_REV};
use std::io::{self};

fn main() -> io::Result<()> {
    env_logger::init();

    let args = Args::parse();
    let config: Config = (&args).into();

    // 1. Align sequences
    log::info!("Aligning sequences...");
    let records = sequence::align_sequences(args.input.to_string());
    log::info!("Done aligning sequences");

    // 2. Extracting n-grams from each sequence segments
    log::info!("Extracting n-grams from each sequence segments...");
    let segment_manager = SegmentManager::new(&records, &config);
    log::info!(
        "Done, total partitions: {}, total segments: {}",
        segment_manager.number_of_partitions(),
        segment_manager.number_of_segments()
    );

    // 3. Calculate frequencies of n-grams for each segment both forward/reverse
    log::info!("Calculating frequencies of k-mer for all segments...");
    log::debug!("Total segments: {}", segment_manager.number_of_segments());
    let candidate_kmers_fwd = segment_manager
        .find_candidates_kmers(SEQ_DIR_FWD)
        .unwrap_or_else(|| Vec::new());
    let candidate_kmers_rev = segment_manager
        .find_candidates_kmers(SEQ_DIR_REV)
        .unwrap_or_else(|| Vec::new());
    log::info!(
        "Done calculating, Total candidate k-mers: fwd: {}, rev: {}",
        candidate_kmers_fwd.len(),
        candidate_kmers_rev.len()
    );

    // 4. Filtering out unmatched criteria
    log::info!("Filtering out unmatched criteria (Tm and >5nt repeats, runs...)");
    let candidate_primers_fwd = filter(&config, candidate_kmers_fwd.clone());
    let candidate_primers_rev = filter(&config, candidate_kmers_rev.clone());
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
