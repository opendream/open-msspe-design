use crate::config::ProgramConfig;
use crate::graphdb::{Edge, GraphDB};
use crate::reverse_complement;
use std::collections::HashMap;
use std::env::current_dir;
use std::io::Write;
use std::process::{Command, Stdio};

impl Edge {
    pub fn get_dg(&self) -> f32 {
        match self.attributes.get("dg") {
            Some(dg) => dg.parse::<f32>().unwrap(),
            None => 0.0,
        }
    }
}

pub struct NtthalOptions {
    pub mv: f32,
    pub dv: f32,
    pub dntp: f32,
    pub conc: f32,
    pub t: f32,
    pub dg: f32,
}

pub fn parse_ntthal_output(input: &str, output: String, delta_g_threshold: f32) -> GraphDB {
    let mut graph = GraphDB::new();

    let mut output_lines = output.lines();
    for input_line in input.lines() {
        if let Some(output_line) = output_lines.next() {
            match output_line.split_whitespace().nth(13) {
                Some(dg_txt) => {
                    let dg = dg_txt.parse::<f32>().unwrap();
                    if dg < delta_g_threshold {
                        let primers = input_line.split(",").collect::<Vec<&str>>();
                        let primer_a = primers[0].to_string();
                        let primer_b = primers[1].to_string();
                        graph.add_node(primer_a.clone());
                        graph.add_node(primer_b.clone());

                        // saving primers
                        let mut attrs = HashMap::new();
                        attrs.insert("dg".to_string(), format!("{:.2}", dg));
                        graph.add_edge(&primer_a, &primer_b, attrs);
                    }
                }
                None => {
                    log::debug!("cannot parse output line: {}", output_line);
                }
            }
        }
        // flush to next group
        output_lines.nth(3);
    }

    graph
}

pub fn format_ntthal_input(primers: &[String], program_config: ProgramConfig) -> String {
    let primer_config = &program_config.primer_config;
    let mut output = "".to_string();
    for a in primers {
        for b in primers {
            // Skip if the primers are the same and program_config.check_self_dimers is false
            if !program_config.check_self_dimers && (a == b || reverse_complement(b) == *a) {
                continue;
            }
            // Skip if program_config.check_cross_dimers is false
            if !program_config.check_cross_dimers {
                continue;
            }
            if a.len() != primer_config.kmer_size || b.len() != primer_config.kmer_size {
                log::debug!("primers not valid: a:{}, b:{}", a, b);
            }
            output.push_str(&format!("{},{}\n", a, b));
        }
    }
    output.trim().to_string()
}

pub fn run_ntthal(
    primers: Vec<String>,
    opts: NtthalOptions,
    program_config: ProgramConfig,
) -> Result<GraphDB, std::io::Error> {
    // Execute ntthal with the given sequences and conditions
    log::trace!("Calculating ΔG for {} sequences", primers.len());
    let path = format!("{}/primer3_config/", current_dir()?.display());

    // Check if macOS, use default primer3_core, else use primer3_core from system.
    let mut cmd = Command::new(program_config.clone().ntthal_path)
        .args([
            "-a",
            "ANY",
            "-mv",
            format!("{:.2}", opts.mv).as_str(),
            "-dv",
            format!("{:.2}", opts.dv).as_str(),
            "-n",
            format!("{:.2}", opts.dntp).as_str(),
            "-d",
            format!("{:.2}", opts.conc).as_str(),
            "-t",
            format!("{:.2}", opts.t).as_str(),
            "-path",
            &*path,
            "-i",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    log::debug!(
        "Executing ntthal with args: \
        -a ANY \
        -mv {:.2} \
        -dv {:.2} \
        -n {:.2} \
        -d {:.2} \
        -t {:.2} \
        -path {}",
        opts.mv,
        opts.dv,
        opts.dntp,
        opts.conc,
        opts.t,
        path
    );

    // Write the input sequences to the stdin of the process
    let input = format_ntthal_input(&primers, program_config.clone());
    let input_clone = input.clone();
    if let Some(mut stdin) = cmd.stdin.take() {
        std::thread::spawn(move || {
            stdin
                .write_all(input_clone.as_bytes())
                .expect("failed to write to stdin");
        });
    }
    log::debug!("Done writing ntthal input");

    // Read the output from the stdout of the process
    let output = cmd.wait_with_output()?;
    if !output.status.success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "ntthal process failed",
        ));
    }

    // Process the output as needed
    let result = String::from_utf8_lossy(&output.stdout);
    Ok(parse_ntthal_output(&input, result.to_string(), opts.dg))
}

#[cfg(test)]
mod tests {
    use crate::config::{PrimerConfig, ProgramConfig};
    use crate::delta_g::{format_ntthal_input, parse_ntthal_output};
    use crate::graphdb::get_edge_id;

    #[test]
    pub fn test_format_ntthal_input() {
        let primers = vec!["GAAGCAGTATTTT".to_string(), "AATATAGAGGCTG".to_string()];
        let program_config = ProgramConfig {
            ntthal_path: "".to_string(),
            primer3_path: "".to_string(),
            max_iterations: 0,
            max_mismatch_segments: 0,
            keep_all: false,
            check_cross_dimers: true,
            check_self_dimers: true,
            check_hairpin: false,
            tm_stddev: 2.0,
            disable_tm_stddev: false,
            do_align: false,
            primer_config: PrimerConfig {
                kmer_size: 13,
                min_tm: 30.0,
                max_tm: 60.0,
                max_self_dimer_any_tm: 20.0,
                max_self_dimer_end_tm: 20.0,
                max_hairpin_tm: 20.0,
            },
        };
        let result = format_ntthal_input(&primers, program_config.clone());
        let expected = "\
            GAAGCAGTATTTT,GAAGCAGTATTTT\n\
            GAAGCAGTATTTT,AATATAGAGGCTG\n\
            AATATAGAGGCTG,GAAGCAGTATTTT\n\
            AATATAGAGGCTG,AATATAGAGGCTG"
            .to_string();
        assert_eq!(result, expected);
    }

    #[test]
    pub fn test_parse_ntthal_output() {
        let input = "\
            AGGCCTATATCCA,GAAGCAGTATTTT\n\
            GCACTTGATGTGA,GAAGCAGTATTTT\n\
            CTGAAGCAGTATT,GCATCTTTCCCTT\n\
            CTGAAGCAGTATT,AATTGTGTGGATT\n\
            AGTCCTGCGTGAT,TGGCCTACATCAG\n\
            "
        .to_string();
        let output = "\
            Calculated thermodynamical parameters for dimer:        dS = -75.3988   dH = -25700     dG = -2315.07   t = -35.9834\n\
            SEQ           AG  CTATATCCA\n\
            SEQ             GC\n\
            STR             CG\n\
            STR     TTTTATGA  AAG------\n\
            Calculated thermodynamical parameters for dimer:        dS = -65.3976   dH = -22500     dG = -2216.94   t = -44.4018\n\
            SEQ           GCA   GATGTGA\n\
            SEQ              CTT\n\
            STR              GAA\n\
            STR     TTTTATGAC   G------\n\
            Calculated thermodynamical parameters for dimer:  dS = -101.596   dH = -33500     dG = -3209.05   t = -24.1908\n\
            SEQ       CT  A  AGTATT\n\
            SEQ         GA GC\n\
            STR         CT CG\n\
            STR TTCCCTTT  A  ------\n\
            Calculated thermodynamical parameters for dimer:  dS = -54.2976   dH = -17700     dG = -1511.18   t = -70.3113\n\
            SEQ CTG  G AGTATT---\n\
            SEQ    AA C\n\
            STR    TT G\n\
            STR      A GTGTGTTAA\n\
            Calculated thermodynamical parameters for dimer:  dS = -141.872   dH = -45800     dG = -3500.74   t = -11.1906\n\
            SEQ AGTC   CG  AT-----\n\
            SEQ     CTG  TG\n\
            STR     GAC  AC\n\
            STR        T-  ATCCGGT".to_string();
        let graph = parse_ntthal_output(&input, output, 100000.00);

        let id_a = &"AGGCCTATATCCA".to_string();
        let id_b = &"GAAGCAGTATTTT".to_string();
        let id_c = &"GCACTTGATGTGA".to_string();
        let id_d = &"CTGAAGCAGTATT".to_string();
        let id_e = &"GCATCTTTCCCTT".to_string();
        let id_f = &"TGGCCTACATCAG".to_string();
        let node_a = graph.get_node(id_a);
        let node_b = graph.get_node(id_b);
        let node_c = graph.get_node(id_c);
        let node_d = graph.get_node(id_d);
        let node_e = graph.get_node(id_e);
        let node_f = graph.get_node(id_f);
        assert!(node_a.is_some());
        assert!(node_b.is_some());
        assert!(node_c.is_some());
        assert!(node_d.is_some());
        assert!(node_e.is_some());
        assert!(node_f.is_some());

        let edges_a = graph.get_edges_for_node(id_a);
        let edges_b = graph.get_edges_for_node(id_b);
        assert_eq!(edges_a.len(), 1);
        assert_eq!(edges_b.len(), 2);

        let edge_ab_id = get_edge_id(id_a, id_b);
        let edge_ab = graph.get_edge(&edge_ab_id);
        assert_eq!(edge_ab.unwrap().get_dg(), -2315.07);

        let edge_de_id = get_edge_id(id_d, id_e);
        let edge_de = graph.get_edge(&edge_de_id);
        assert_eq!(edge_de.unwrap().get_dg(), -3209.05);
    }
}
