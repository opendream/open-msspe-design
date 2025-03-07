use crate::constants::{
    ANNEALING_TEMP, DEFAULT_NTTHAL_PATH, DEFAULT_PRIMER3_PATH, DELTA_G_THRESHOLD, DNA_CONC,
    DNTP_CONC, DV_CONC, KMER_SIZE, MAX_ITERATIONS, MAX_MISMATCH_SEGMENTS, MV_CONC, OVERLAP_SIZE,
    PRIMER_MAX_HAIRPIN_TH, PRIMER_MAX_SELF_ANY_TH, PRIMER_MAX_SELF_END_TH, PRIMER_MAX_TM,
    PRIMER_MIN_TM, SEARCH_WINDOWS_SIZE, WINDOW_SIZE,
};

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
pub struct Args {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long)]
    pub output: String,

    #[arg(long, default_value_t = KMER_SIZE)]
    pub kmer_size: usize,

    #[arg(long, default_value_t = WINDOW_SIZE)]
    pub window_size: usize,

    #[arg(long, default_value_t = OVERLAP_SIZE)]
    pub overlap_size: usize,
    #[arg(long, default_value_t = MAX_MISMATCH_SEGMENTS)]
    pub max_mismatch_segments: usize,
    #[arg(long, default_value_t = MAX_ITERATIONS)]
    pub max_iterations: usize,
    #[arg(long, default_value_t = SEARCH_WINDOWS_SIZE)]
    pub search_windows_size: usize,
    #[arg(long, default_value_t = MV_CONC)]
    pub mv_conc: f32,
    #[arg(long, default_value_t = DV_CONC)]
    pub dv_conc: f32,
    #[arg(long, default_value_t = DNTP_CONC)]
    pub dntp_conc: f32,
    #[arg(long, default_value_t = DNA_CONC)]
    pub dna_conc: f32,
    #[arg(long, default_value_t = ANNEALING_TEMP)]
    pub annealing_temp: f32,
    #[arg(long, default_value_t = PRIMER_MIN_TM)]
    pub min_tm: f32,
    #[arg(long, default_value_t = PRIMER_MAX_TM)]
    pub max_tm: f32,
    #[arg(long, default_value_t = PRIMER_MAX_SELF_ANY_TH)]
    pub max_self_dimer_any_tm: f32,
    #[arg(long, default_value_t = PRIMER_MAX_SELF_END_TH)]
    pub max_self_dimer_end_tm: f32,
    #[arg(long, default_value_t = PRIMER_MAX_HAIRPIN_TH)]
    pub max_hairpin_tm: f32,
    #[arg(long, default_value_t = DELTA_G_THRESHOLD, help = "Threshold for dG, default is -9000.0 J/mol")]
    pub delta_g_threshold: f32,

    // logic based config
    #[arg(
        group = "flag",
        long,
        default_value = "no",
        value_parser = ["yes", "no"],
        help = "Ignores all filtering and does NOT remove any primers."
    )]
    pub keep_all: String,

    #[arg(
        group = "flag",
        long,
        default_value = "yes",
        value_parser = ["yes", "no"],
        help = "\
        Calculates nearest neighbor thermodynamic model for delta G values for every primer pair \
        across entire pool of primers, and removes primer with most secondary structures, or lowest Tm."
    )]
    pub check_cross_dimers: String,

    #[arg(
        group = "flag",
        long,
        default_value = "yes",
        value_parser = ["yes", "no"],
        help = "Calculates Tm for a self-dimer of an individual primer sequence."
    )]
    pub check_self_dimers: String,

    #[arg(
        group = "flag",
        long,
        default_value = "yes",
        value_parser = ["yes", "no"],
        help = "Calculates Tm for a hairpin of an individual primer sequence."
    )]
    pub check_hairpin: String,

    #[arg(
        group = "flag",
        long,
        default_value = "yes",
        value_parser = ["yes", "no"],
        help = "Removes primers with melting temperature greater than 2 standard deviations from mean of the Tm values for all primers in set. Default is true."
    )]
    pub strict_tm_range: String,
    #[arg(
        group = "flag",
        long,
        default_value = "yes",
        value_parser = ["yes", "no"],
        help = "Does MAFFT multiple sequence alignment."
    )]
    pub do_align: String,

    // vendor binary path
    #[arg(long, default_value_t = DEFAULT_NTTHAL_PATH.to_string())]
    pub ntthal: String,

    #[arg(long, default_value_t = DEFAULT_PRIMER3_PATH.to_string())]
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
    pub max_iterations: usize,
    pub max_mismatch_segments: usize,

    pub keep_all: bool,
    pub check_cross_dimers: bool,
    pub check_self_dimers: bool,
    pub check_hairpin: bool,
    pub strict_tm_range: bool,
    pub do_align: bool,

    pub(crate) primer_config: PrimerConfig,
}
