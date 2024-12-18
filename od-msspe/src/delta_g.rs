use std::env::current_dir;
use std::process::Command;

/// Function to call ntthal and calculate ΔG for given sequences
pub fn calculate_delta_g(seq_a: &str, seq_b: &str, mv: f64, dv: f64, dntp: f64, conc: f64) -> Option<f64> {
    // Execute ntthal with the given sequences and conditions
    log::trace!("Calculating ΔG for sequences: {} and {}", seq_a, seq_b);
    let path = format!("{}/primer3_config/", current_dir().unwrap().display());
    let mut cmd = Command::new("./bin/ntthal");
    cmd.args(&[
            "-a", "HAIRPIN",
            "-mv", format!("{:.2}", mv).as_str(),
            "-dv", format!("{:.2}", dv).as_str(),
            "-d", format!("{:.2}", dntp).as_str(),
            "-n", format!("{:.2}", conc).as_str(),
            "-path", &*path,
        ]);
    if seq_a == seq_b {
        cmd.arg("-s1").arg(seq_a);
    } else {
        cmd.arg("-s1").arg(seq_a).arg("-s2").arg(seq_b);
    }
    let output = cmd.output().expect("Failed to execute ntthal");

    // Check if command execution was successful
    if !output.status.success() {
        log::error!("ntthal failed with status: {}, msg: {}", output.status, &String::from_utf8_lossy(&output.stderr)[..100]);
        return None;
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let re = regex::Regex::new(r"dG = (-?\d+(\.\d+)?)").unwrap();

    let caps = stdout.lines()
        .find(|line| re.captures(line).is_some())
        .map(|line| re.captures(line).unwrap()[1].to_string());
    match caps {
        Some(caps) => {
            let d_g = caps.parse::<f64>().unwrap();
            log::debug!("ΔG: {}", d_g);
            Some(d_g)
        },
        None => {
            None
        }
    }
}
