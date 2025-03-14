pub(crate) use crate::constants::{
    ANNEALING_TEMP, DEFAULT_NTTHAL_PATH, DEFAULT_PRIMER3_PATH, DELTA_G_THRESHOLD, DNA_CONC,
    DNTP_CONC, DV_CONC, KMER_SIZE, MAX_ITERATIONS, MAX_MISMATCH_SEGMENTS, MV_CONC, OVERLAP_SIZE,
    PRIMER_MAX_HAIRPIN_TH, PRIMER_MAX_SELF_ANY_TH, PRIMER_MAX_SELF_END_TH, PRIMER_MAX_TM,
    PRIMER_MIN_TM, SEARCH_WINDOWS_SIZE, WINDOW_SIZE,
};
use std::path::Path;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
pub struct Args {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long)]
    pub output: String,

    #[arg(long, env = "KMER_SIZE", default_value_t = KMER_SIZE)]
    pub kmer_size: usize,

    #[arg(long, env = "WINDOW_SIZE", default_value_t = WINDOW_SIZE)]
    pub window_size: usize,

    #[arg(long, env = "OVERLAP_SIZE", default_value_t = OVERLAP_SIZE)]
    pub overlap_size: usize,
    #[arg(long, env = "MAX_MISMATCH_SEGMENTS", default_value_t = MAX_MISMATCH_SEGMENTS)]
    pub max_mismatch_segments: usize,
    #[arg(long, env = "MAX_ITERATIONS", default_value_t = MAX_ITERATIONS)]
    pub max_iterations: usize,
    #[arg(long, env = "SEARCH_WINDOWS_SIZE", default_value_t = SEARCH_WINDOWS_SIZE)]
    pub search_windows_size: usize,
    #[arg(long, env = "MV_CONC", default_value_t = MV_CONC)]
    pub mv_conc: f32,
    #[arg(long, env = "DV_CONC", default_value_t = DV_CONC)]
    pub dv_conc: f32,
    #[arg(long, env = "DNTP_CONC", default_value_t = DNTP_CONC)]
    pub dntp_conc: f32,
    #[arg(long, env = "DNA_CONC", default_value_t = DNA_CONC)]
    pub dna_conc: f32,
    #[arg(long, env = "ANNEALING_TEMP", default_value_t = ANNEALING_TEMP)]
    pub annealing_temp: f32,
    #[arg(long, env = "MIN_TM", default_value_t = PRIMER_MIN_TM)]
    pub min_tm: f32,
    #[arg(long, env = "MAX_TM", default_value_t = PRIMER_MAX_TM)]
    pub max_tm: f32,
    #[arg(long, env = "MAX_SELF_DIMER_ANY_TM", default_value_t = PRIMER_MAX_SELF_ANY_TH)]
    pub max_self_dimer_any_tm: f32,
    #[arg(long, env = "MAX_SELF_DIMER_END_TM", default_value_t = PRIMER_MAX_SELF_END_TH)]
    pub max_self_dimer_end_tm: f32,
    #[arg(long, env = "MAX_HAIRPIN_TM", default_value_t = PRIMER_MAX_HAIRPIN_TH)]
    pub max_hairpin_tm: f32,
    #[arg(long, env = "DELTA_G_THRESHOLD", default_value_t = DELTA_G_THRESHOLD, help = "Threshold for dG, default is -9000.0 J/mol")]
    pub delta_g_threshold: f32,

    // logic based config
    #[arg(
        group = "flag",
        long,
        env = "KEEP_ALL",
        default_value = "false",
        value_parser = ["true", "false"],
        help = "Ignores all filtering and does NOT remove any primers."
    )]
    pub keep_all: String,

    #[arg(
        group = "flag",
        long,
        env = "CHECK_CROSS_DIMERS",
        default_value = "true",
        value_parser = ["true", "false"],
        help = "\
            Calculates nearest neighbor thermodynamic model for delta G values for every primer pair \
            across entire pool of primers, and removes primer with most secondary structures, or lowest Tm."
    )]
    pub check_cross_dimers: String,

    #[arg(
        group = "flag",
        long,
        env = "CHECK_SELF_DIMERS",
        default_value = "true",
        value_parser = ["true", "false"],
        help = "Calculates Tm for a self-dimer of an individual primer sequence."
    )]
    pub check_self_dimers: String,

    #[arg(
        group = "flag",
        long,
        env = "CHECK_HAIRPIN",
        default_value = "true",
        value_parser = ["true", "false"],
        help = "Calculates Tm for a hairpin of an individual primer sequence."
    )]
    pub check_hairpin: String,

    #[arg(
        group = "flag",
        long,
        env = "TM_STDDEV",
        default_value_t = 2.0,
        help = "Set the number of standard deviations away from the mean of the tm values."
    )]
    pub tm_stddev: f32,

    #[arg(
        group = "flag",
        long,
        env = "DISABLE_TM_STDDEV",
        default_value = "false",
        value_parser = ["true", "false"],
        help = "\
            Turns off tm-stddev config. Use if you do not want strictly similar tm values \
            across all primers."
    )]
    pub disable_tm_stddev: String,

    #[arg(
        group = "flag",
        long,
        env = "DO_ALIGN",
        default_value = "true",
        value_parser = ["true", "false"],
        help = "Does MAFFT multiple sequence alignment."
    )]
    pub do_align: String,

    // vendor binary path
    #[arg(long, env ="NTTHAL", default_value_t = DEFAULT_NTTHAL_PATH.to_string())]
    pub ntthal: String,

    #[arg(long, env = "PRIMER3", default_value_t = DEFAULT_PRIMER3_PATH.to_string())]
    pub primer3: String,
}

#[derive(Clone)]
pub struct PrimerConfig {
    pub kmer_size: usize,
    pub min_tm: f32,
    pub max_tm: f32,
    pub max_self_dimer_any_tm: f32,
    pub max_self_dimer_end_tm: f32,
    pub max_hairpin_tm: f32,
}

#[derive(Clone)]
pub struct ProgramConfig {
    pub ntthal_path: String,
    pub primer3_path: String,

    pub max_iterations: usize,
    pub max_mismatch_segments: usize,

    pub keep_all: bool,
    pub check_cross_dimers: bool,
    pub check_self_dimers: bool,
    pub check_hairpin: bool,
    pub tm_stddev: f32,
    pub disable_tm_stddev: bool,
    pub do_align: bool,

    pub(crate) primer_config: PrimerConfig,
}

pub fn find_executable(name: &str, exact: bool) -> Option<String> {
    // use provided path if it exists
    let path = Path::new(name);
    if std::fs::metadata(path).is_ok() {
        return Some(path.to_str().unwrap().to_string());
    }

    if exact {
        return None;
    }

    // concat paths with local bin in "./bin" directory
    let paths = format!("{}:./bin", std::env::var("PATH").unwrap());
    log::debug!("paths: {}", paths);
    for path in paths.split(':') {
        let path = Path::new(path).join(name);
        log::debug!("path: {:?}", path);
        if std::fs::metadata(&path).is_ok() {
            return Some(path.to_str().unwrap().to_string());
        }
    }

    None
}
