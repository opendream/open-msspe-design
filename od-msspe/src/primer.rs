use std::io::Write;
use std::process::{Command, Stdio};
use std::string::ParseError;

/// Find primer related values like Tm, Hairpin, Dimers using Primer3

#[derive(Clone)]
pub struct PrimerInfo<'a> {
    pub id: &'a str,
    pub tm: f32,
    pub gc: f32,
    pub self_any_th: f32,
    pub self_end_th: f32,
    pub hairpin_th: f32,
}

impl<'a> PrimerInfo<'a> {
    pub(crate) fn new() -> PrimerInfo<'a> {
        PrimerInfo {
            id: "",
            tm: 0.0,
            gc: 0.0,
            self_any_th: 0.0,
            self_end_th: 0.0,
            hairpin_th: 0.0,
        }
    }
    pub(crate) fn reset(&mut self) {
        self.id = "";
        self.tm = 0.0;
        self.gc = 0.0;
        self.self_any_th = 0.0;
        self.self_end_th = 0.0;
        self.hairpin_th = 0.0;
    }
}

pub struct CheckPrimerParams {
    pub min_tm: f32,
    pub max_tm: f32,
    pub primer3_path: String,
}

/// Try to parse Primer3 output as PrimerInfo
/// Example: primer3_core output
/// SEQUENCE_ID=example1
/// SEQUENCE_PRIMER=AGCCCGTGTAAAC
/// PRIMER_TASK=check_primers
/// PRIMER_MIN_SIZE=13
/// PRIMER_MIN_TM=30.0
/// PRIMER_MAX_TM=60.0
/// PRIMER_LEFT_NUM_RETURNED=1
/// PRIMER_RIGHT_NUM_RETURNED=0
/// PRIMER_INTERNAL_NUM_RETURNED=0
/// PRIMER_PAIR_NUM_RETURNED=0
/// PRIMER_LEFT_0_PENALTY=23.273240
/// PRIMER_LEFT_0_SEQUENCE=AGCCCGTGTAAAC
/// PRIMER_LEFT_0=0,13
/// PRIMER_LEFT_0_TM=43.727
/// PRIMER_LEFT_0_GC_PERCENT=53.846
/// PRIMER_LEFT_0_SELF_ANY_TH=0.00
/// PRIMER_LEFT_0_SELF_END_TH=0.00
/// PRIMER_LEFT_0_HAIRPIN_TH=0.00
/// PRIMER_LEFT_0_END_STABILITY=2.0100
/// =
/// SEQUENCE_ID=example2
pub fn parse_primer3_output(result: &str) -> Result<Vec<PrimerInfo<'static>>, ParseError> {
    let mut output: Vec<PrimerInfo> = Vec::new();

    // split lines with \n and group each output when met `=` at line char idx 0
    let mut primer_info = PrimerInfo {
        id: "",
        tm: 0.0,
        gc: 0.0,
        self_any_th: 0.0,
        self_end_th: 0.0,
        hairpin_th: 0.0,
    };
    for line in result.lines() {
        match line {
            "=" => {
                if primer_info.id.is_empty() {
                    break;
                }
                output.push(primer_info.clone());
                primer_info.reset();
            }
            _ => {
                let key_value: Vec<&str> = line.split("=").collect();
                match key_value[0] {
                    "SEQUENCE_ID" => {
                        primer_info.id = Box::leak(key_value[1].to_string().into_boxed_str())
                    }
                    "PRIMER_LEFT_0_TM" => primer_info.tm = key_value[1].parse::<f32>().unwrap(),
                    "PRIMER_LEFT_0_GC_PERCENT" => {
                        primer_info.gc = key_value[1].parse::<f32>().unwrap()
                    }
                    "PRIMER_LEFT_0_SELF_ANY_TH" => {
                        primer_info.self_any_th = key_value[1].parse::<f32>().unwrap()
                    }
                    "PRIMER_LEFT_0_SELF_END_TH" => {
                        primer_info.self_end_th = key_value[1].parse::<f32>().unwrap()
                    }
                    "PRIMER_LEFT_0_HAIRPIN_TH" => {
                        primer_info.hairpin_th = key_value[1].parse::<f32>().unwrap()
                    }
                    _ => (),
                }
            }
        }
    }

    Ok(output)
}

/// Format k-mer string as primer3_core input
/// Example:
/// SEQUENCE_ID=example1
/// SEQUENCE_PRIMER=AGCCCGTGTAAAC
/// PRIMER_TASK=check_primers
/// PRIMER_MIN_SIZE=13
/// PRIMER_MIN_TM=30.0
/// PRIMER_MAX_TM=60.0
/// =
pub fn format_primer3_input(primers: &[String], params: &CheckPrimerParams) -> String {
    let mut input = String::new();
    for primer in primers {
        input.push_str(&format!("SEQUENCE_ID={}\n", primer));
        input.push_str(&format!("SEQUENCE_PRIMER={}\n", primer));
        input.push_str("PRIMER_TASK=check_primers\n");
        input.push_str("PRIMER_MIN_SIZE=13\n");
        input.push_str(&format!("PRIMER_MIN_TM={:.2}\n", params.min_tm));
        input.push_str(&format!("PRIMER_MAX_TM={:.2}\n", params.max_tm));
        // @see https://primer3.org/manual#PRIMER_OPT_TM
        input.push_str(&format!("PRIMER_OPT_TM={:.2}\n", params.max_tm));
        input.push_str("PRIMER_PICK_ANYWAY=1\n");
        input.push_str("=\n");
    }
    input
}

/// Call ./bin/primer3_core and provide input as pipeline.
pub fn check_primers(
    primers: &[String],
    params: CheckPrimerParams,
) -> Result<Vec<PrimerInfo<'static>>, std::io::Error> {
    let mut primer_info_list = Vec::new();

    let input = format_primer3_input(primers, &params);
    log::debug!("using primer3 binary: {:?}", params.primer3_path);
    let mut cmd = Command::new(params.primer3_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    let mut stdin = cmd.stdin.take().expect("Failed to open stdin");
    stdin.write_all(input.as_bytes())?;
    drop(stdin);

    let output = cmd.wait_with_output()?;
    let stdout = String::from_utf8(output.stdout).unwrap();
    let parsed_result = parse_primer3_output(&stdout).unwrap();
    primer_info_list.extend(parsed_result);

    Ok(primer_info_list)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::find_executable;

    fn get_test_primer3_params() -> CheckPrimerParams {
        let primer3_path = find_executable("primer3_core", false).unwrap();
        CheckPrimerParams {
            min_tm: 29.0,
            max_tm: 59.0,
            primer3_path,
        }
    }

    #[test]
    fn test_parse_primer3_output() {
        let result = "\
            SEQUENCE_ID=example1\n\
            SEQUENCE_PRIMER=AGCCCGTGTAAAC\n\
            PRIMER_TASK=check_primers\n\
            PRIMER_MIN_SIZE=13\n\
            PRIMER_MIN_TM=30.0\n\
            PRIMER_MAX_TM=60.0\n\
            PRIMER_LEFT_NUM_RETURNED=1\n\
            PRIMER_RIGHT_NUM_RETURNED=0\n\
            PRIMER_INTERNAL_NUM_RETURNED=0\n\
            PRIMER_PAIR_NUM_RETURNED=0\n\
            PRIMER_LEFT_0_PENALTY=23.273240\n\
            PRIMER_LEFT_0_SEQUENCE=AGCCCGTGTAAAC\n\
            PRIMER_LEFT_0=0,13\n\
            PRIMER_LEFT_0_TM=43.727\n\
            PRIMER_LEFT_0_GC_PERCENT=53.846\n\
            PRIMER_LEFT_0_SELF_ANY_TH=0.00\n\
            PRIMER_LEFT_0_SELF_END_TH=0.00\n\
            PRIMER_LEFT_0_HAIRPIN_TH=0.00\n\
            PRIMER_LEFT_0_END_STABILITY=2.0100\n\
            =";
        let output = parse_primer3_output(result);

        assert!(output.is_ok());
        let primer_info_list = output.unwrap();
        assert_eq!(primer_info_list[0].id, "example1");
        assert_eq!(primer_info_list[0].tm, 43.727);
        assert_eq!(primer_info_list[0].gc, 53.846);
        assert_eq!(primer_info_list[0].self_any_th, 0.00);
        assert_eq!(primer_info_list[0].self_end_th, 0.00);
        assert_eq!(primer_info_list[0].hairpin_th, 0.00);
    }

    #[test]
    fn test_format_primer3_input() {
        let primers = vec!["AGCCCGTGTAAAC".to_string()];
        let params = get_test_primer3_params();
        let result = format_primer3_input(&primers, &params);
        assert_eq!(
            result,
            "\
            SEQUENCE_ID=AGCCCGTGTAAAC\n\
            SEQUENCE_PRIMER=AGCCCGTGTAAAC\n\
            PRIMER_TASK=check_primers\n\
            PRIMER_MIN_SIZE=13\n\
            PRIMER_MIN_TM=29.00\n\
            PRIMER_MAX_TM=59.00\n\
            PRIMER_OPT_TM=59.00\n\
            PRIMER_PICK_ANYWAY=1\n\
            =\n"
        )
    }

    #[test]
    fn test_check_primers() {
        let primers = vec!["AGCCCGTGTAAAC".to_string()];
        let params = get_test_primer3_params();
        let result = check_primers(&primers, params);
        assert!(result.is_ok());
        let primer_info_list = result.unwrap();
        assert_eq!(primer_info_list[0].id, "AGCCCGTGTAAAC");
        assert_eq!(primer_info_list[0].tm, 43.727);
        assert_eq!(primer_info_list[0].gc, 53.846);
        assert_eq!(primer_info_list[0].self_any_th, 0.00);
        assert_eq!(primer_info_list[0].self_end_th, 0.00);
        assert_eq!(primer_info_list[0].hairpin_th, 0.00);
    }
}
