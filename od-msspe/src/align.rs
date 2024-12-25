use crate::types::SequenceRecord;
use seq_io::fasta::{Reader, Record};
use std::io::BufReader;

/**
 * Aligns sequences using MAFFT
 */
pub fn align_sequences(filepath: String) -> Vec<SequenceRecord> {
    let output = std::process::Command::new("mafft")
        .args([
            "--auto",
            "--quiet",
            "--thread",
            "-1",
            "--op",
            "1.53",
            "--ep",
            "0.123",
            "--jtt",
            "200",
            &filepath.clone(),
        ])
        .output()
        .expect("failed to execute MAFFT");

    return to_records(output.stdout);
}

fn to_records(src: Vec<u8>) -> Vec<SequenceRecord> {
    let mut reader = Reader::new(BufReader::new(src.as_slice()));
    let mut records = Vec::new();

    while let Some(result) = reader.next() {
        let record = result.unwrap();
        let name = record.id().unwrap().to_string();
        let sequence = String::from_utf8(record.full_seq().to_vec())
            .unwrap()
            .to_uppercase()
            .replace("U", "T");
        records.push(SequenceRecord { name, sequence });
    }

    if records.is_empty() {
        panic!("No sequences found in the input file");
    }

    return records;
}

#[cfg(test)]
mod tests {
    use super::*;

    use proptest::proptest;

    #[test]
    fn test_to_records() {
        let fasta_data = b">seq1\nAACCTTGGAACCTTG\n>seq2\nAACCTTGGAACCTTG-\n";
        let records = to_records(fasta_data.to_vec());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "seq1");
        assert_eq!(records[0].sequence, "AACCTTGGAACCTTG");
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "AACCTTGGAACCTTG-");
    }

    #[test]
    // test lower case
    fn test_to_records_lowercase() {
        let fasta_data = b">seq1\naaccttggaaccttg\n>seq2\naaccttggaaccttg-\n";
        let records = to_records(fasta_data.to_vec());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "seq1");
        assert_eq!(records[0].sequence, "AACCTTGGAACCTTG");
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "AACCTTGGAACCTTG-");
    }

    #[test]
    // test U to T conversion
    fn test_to_records_uracil() {
        let fasta_data = b">seq1\nAACCTTGGUUACCTTG\n>seq2\nAACCTTGGUUACCTTG-\n";
        let records = to_records(fasta_data.to_vec());
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "seq1");
        assert_eq!(records[0].sequence, "AACCTTGGTTACCTTG");
        assert_eq!(records[1].name, "seq2");
        assert_eq!(records[1].sequence, "AACCTTGGTTACCTTG-");
    }

    proptest! {
        #[test]
        fn test_to_records_with_random_data(name in "[a-zA-Z0-9]{1,10}", seq in "[ACGTUacgtu-]{1,100}") {
            let fasta_data = format!(">{}\n{}\n", name, seq).into_bytes();
            let records = to_records(fasta_data);
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].name, name);
            let expected_sequence = seq.to_uppercase().replace("U", "T");
            assert_eq!(records[0].sequence, expected_sequence);
        }
    }
}
