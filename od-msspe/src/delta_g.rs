use std::process::Command;

/// Function to call ntthal and calculate ΔG for given sequences
pub fn calculate_delta_g(seq_a: &str, seq_b: &str, mv: u32, dv: u32, dntp: u32, conc: u32) -> Option<f64> {
    // Execute ntthal with the given sequences and conditions
    let output = Command::new("./bin/ntthal")
        .args(&[
            "-a", seq_a,
            "-b", seq_b,
            "-mv", &mv.to_string(),
            "-dv", &dv.to_string(),
            "-d", &dntp.to_string(),
            "-n", &conc.to_string(),
            "-path", "./primer3_config/",
        ])
        .output();

    // Check if command execution was successful
    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout);
            for line in stdout.lines() {
                if line.contains("dG =") {
                    if let Some(delta_g) = line.split_whitespace().last() {
                        return delta_g.parse::<f64>().ok();
                    }
                }
            }
            None
        }
        Err(e) => {
            eprintln!("Error executing ntthal: {}", e);
            None
        }
    }
}
