use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None)]
pub struct Args {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long)]
    pub output: String,

    #[arg(long, default_value_t = 13)]
    pub kmer_size: usize,

    #[arg(long, default_value_t = 500)]
    pub window_size: usize,

    #[arg(long, default_value_t = 250)]
    pub overlap_size: usize,

    #[arg(long, default_value_t = 1)]
    pub max_mismatch_segments: usize,

    #[arg(long, default_value_t = 1000)]
    pub max_iterations: usize,

    #[arg(long, default_value_t = 50)]
    pub search_windows_size: usize,

    #[arg(long, default_value_t = 50.0)]
    pub mv_conc: f64,

    #[arg(long, default_value_t = 3.0)]
    pub dv_conc: f64,

    #[arg(long, default_value_t = 0.0)]
    pub dntp_conc: f64,

    #[arg(long, default_value_t = 250.0)]
    pub dna_conc: f64,

    #[arg(long, default_value_t = 25.0)]
    pub annealing_temp: f64,
}

pub struct PartitioningOption {
    segment_size: usize,
    overlap_size: usize,
    window_size: usize,
    kmer_size: usize,
}

impl PartitioningOption {
    pub fn new(config: &Config) -> Self {
        if config.overlap_size() > config.window_size() {
            panic!("overlap_size must be lesser than or equal to window_size");
        }
        PartitioningOption {
            segment_size: config.window_size(),
            overlap_size: config.overlap_size(),
            window_size: config.search_windows_size(),
            kmer_size: config.kmer_size(),
        }
    }

    // constructor for testing
    #[allow(dead_code)]
    pub fn new_with(
        segment_size: usize,
        overlap_size: usize,
        window_size: usize,
        kmer_size: usize,
    ) -> Self {
        PartitioningOption {
            segment_size,
            overlap_size,
            window_size,
            kmer_size,
        }
    }

    pub fn segment_size(&self) -> usize {
        self.segment_size
    }

    pub fn overlap_size(&self) -> usize {
        self.overlap_size
    }

    pub fn window_size(&self) -> usize {
        self.window_size
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }
}

/*
 * private properties, provide only getters
 * to avoid direct access to the properties
 * from outside the module
 */
pub struct Config {
    window_size: usize,
    kmer_size: usize,
    overlap_size: usize,
    max_mismatch_segments: usize,
    max_iterations: usize,
    search_windows_size: usize,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    annealing_temp: f64,
}

impl Config {
    // default constructor
    #[allow(dead_code)]
    pub fn default() -> Self {
        Config {
            window_size: 500,
            kmer_size: 13,
            overlap_size: 250,
            max_mismatch_segments: 1,
            max_iterations: 1000,
            search_windows_size: 50,
            mv_conc: 50.0,
            dv_conc: 3.0,
            dntp_conc: 0.0,
            dna_conc: 250.0,
            annealing_temp: 25.0,
        }
    }

    // constructor for testing
    #[allow(dead_code)]
    pub fn new_for_test(
        segment_size: usize,
        overlap_size: usize,
        window_size: usize,
        kmer_size: usize,
    ) -> Self {
        Config {
            window_size: segment_size,
            kmer_size: kmer_size,
            overlap_size: overlap_size,
            max_mismatch_segments: 1,
            max_iterations: 1000,
            search_windows_size: window_size,
            mv_conc: 50.0,
            dv_conc: 3.0,
            dntp_conc: 0.0,
            dna_conc: 250.0,
            annealing_temp: 25.0,
        }
    }

    pub fn window_size(&self) -> usize {
        self.window_size
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    pub fn overlap_size(&self) -> usize {
        self.overlap_size
    }

    pub fn max_mismatch_segments(&self) -> usize {
        self.max_mismatch_segments
    }

    pub fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    pub fn search_windows_size(&self) -> usize {
        self.search_windows_size
    }

    pub fn mv_conc(&self) -> f64 {
        self.mv_conc
    }

    pub fn dv_conc(&self) -> f64 {
        self.dv_conc
    }

    pub fn dntp_conc(&self) -> f64 {
        self.dntp_conc
    }

    pub fn dna_conc(&self) -> f64 {
        self.dna_conc
    }

    pub fn annealing_temp(&self) -> f64 {
        self.annealing_temp
    }
}

impl From<&Args> for Config {
    fn from(args: &Args) -> Self {
        Config {
            kmer_size: args.kmer_size,
            window_size: args.window_size,
            overlap_size: args.overlap_size,
            max_mismatch_segments: args.max_mismatch_segments,
            max_iterations: args.max_iterations,
            search_windows_size: args.search_windows_size,
            mv_conc: args.mv_conc,
            dv_conc: args.dv_conc,
            dntp_conc: args.dntp_conc,
            dna_conc: args.dna_conc,
            annealing_temp: args.annealing_temp,
        }
    }
}
