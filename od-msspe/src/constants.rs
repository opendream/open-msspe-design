pub const KMER_SIZE: usize = 13;
pub const WINDOW_SIZE: usize = 500;
pub const OVERLAP_SIZE: usize = 250;
pub const MAX_MISMATCH_SEGMENTS: usize = 1;
pub const MAX_ITERATIONS: usize = 1000;
pub const SEARCH_WINDOWS_SIZE: usize = 50;
// Monovalent cation concentration (mM)
pub const MV_CONC: f32 = 50.0;
// Divalent cation concentration (mM)
pub const DV_CONC: f32 = 3.0;
// dNTP concentration (mM)
pub const DNTP_CONC: f32 = 0.0;
// Primer concentration (nM)
pub const DNA_CONC: f32 = 250.0;
// Annealing temperature (°C)
pub const ANNEALING_TEMP: f32 = 25.0;
pub const PRIMER_MIN_TM: f32 = 30.0;
pub const PRIMER_MAX_TM: f32 = 60.0;
pub const PRIMER_MAX_SELF_ANY_TH: f32 = PRIMER_MIN_TM - 10.0;
pub const PRIMER_MAX_SELF_END_TH: f32 = PRIMER_MIN_TM - 10.0;
pub const PRIMER_MAX_HAIRPIN_TH: f32 = PRIMER_MIN_TM - 10.0;
pub const DELTA_G_THRESHOLD: f32 = -9000.0;
pub const SEQ_DIR_FWD: u8 = 0x00;
pub const SEQ_DIR_REV: u8 = 0x01;

pub const DEFAULT_NTTHAL_PATH: &str = "ntthal";
pub const DEFAULT_PRIMER3_PATH: &str = "primer3_core";
