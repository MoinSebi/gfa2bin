use crate::core::bfile::count_lines;
use crate::core::core::MatrixWrapper;
use bitvec::vec::BitVec;
use clap::ArgMatches;
use log::info;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/// # View main function
///
/// Convert a bed file to VCF file
/// Vcf is bim + bed
pub fn view_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin view'");
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    // Read the bed file

    info!("Plink file: {}", plink_file);
    info!("Output file: {}", output_prefix);

    info!("Writing output (vcf)");
    write_vcf(plink_file, output_prefix)?;
    Ok(())
}

/// Read the plink and write a vcf-like file
///
/// All-in-one function
/// Does represent haplotypes
/// Not sure of version vcf format 4.2
pub fn write_vcf(
    filename_prefix: &str,
    output_prefix: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Write the header
    let file = File::create(output_prefix).unwrap();
    let mut writer = std::io::BufWriter::new(file);
    writeln!(writer, "##fileformat=VCFv4.2").expect("Error writing to file");
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .expect("Error writing to file");

    // Read the bim file
    let fam = read_first_column_from_tsv(&format!("{}{}", filename_prefix, ".fam"))?;
    write!(
        writer,
        "{}",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".to_string()
            + &fam.join("\t")
            + "\n"
    )
    .expect("Error writing to file");

    let mut mw = MatrixWrapper::new();
    let bim_count = count_lines(&format!("{}{}", filename_prefix, ".bim"))?;
    info!("bim count: {}", bim_count);
    let fam_count = count_lines(&format!("{}{}", filename_prefix, ".fam"))?;
    mw.read_bed(
        &format!("{}{}", filename_prefix, ".bed"),
        fam_count,
        bim_count,
    )
    .expect("Error reading bed file");

    let file_bim = File::open(format!("{}{}", filename_prefix, ".bim"))?;
    let reader_bim = BufReader::new(file_bim);
    let lines_bim = reader_bim.lines();

    for (index, line_bim) in lines_bim.enumerate() {
        let line1 = line_bim?;
        let b = &mw.matrix_bit[index];
        writeln!(
            writer,
            "graph\t{}\t.\t{}\t-\t{}\tPASS\t{}\tGT\t{}",
            line1.split_whitespace().nth(3).unwrap(),
            line1.split_whitespace().nth(3).unwrap(),
            mw.feature.to_string1(),
            "GT=".to_string() + &(b.len() / 2).to_string(),
            bitvec2vcf_string(b)
        )?;
    }

    Ok(())
}

/// # Read the first column from a tsv file
///
/// - Samples of a fam file
fn read_first_column_from_tsv(file_path: &str) -> Result<Vec<String>, std::io::Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut first_columns = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let columns: Vec<&str> = line.split_whitespace().collect();

        if let Some(first_column) = columns.first() {
            first_columns.push((*first_column).to_string());
        }
    }

    Ok(first_columns)
}

impl MatrixWrapper {
    /// # Write a VCF file with dumb header
    pub fn write_vcf(&self, filename: &str) {
        let file = File::create(filename).unwrap();
        let mut writer = std::io::BufWriter::new(file);
        writeln!(writer, "##fileformat=VCFv4.2").expect("Error writing to file");
        writeln!(
            writer,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .expect("Error writing to file");

        let fam_entries = self
            .fam_entries
            .iter()
            .map(|x| x.split_whitespace().next().unwrap().to_string())
            .collect::<Vec<String>>();
        let bims = self
            .bim_entries
            .iter()
            .map(|x| x.split_whitespace().nth(3).unwrap().to_string())
            .collect::<Vec<String>>();

        write!(
            writer,
            "{}",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".to_string()
                + &fam_entries.join("\t")
                + "\n"
        )
        .expect("Test");

        for (x, y) in self.matrix_bit.iter().zip(bims.iter()) {
            writeln!(
                writer,
                "graph\t{}\t.\t{}\t-\t{}\tPASS\t{}\tGT\t{}",
                y,
                y,
                self.feature.to_string1(),
                "GT=".to_string() + &(x.len() / 2).to_string(),
                bitvec2vcf_string(x)
            )
            .expect("Error writing to file");
        }
    }
}

/// # Convert bitvector to string to 0 and 1 string
///
/// - Using chunks in pairs of two
// Convert bitvec to vcf string
pub fn bitvec2vcf_string(bitvec: &BitVec<u8>) -> String {
    let mut s = String::new();
    let chunks = bitvec.chunks(2);
    for x in chunks {
        if !x[0] && !x[1] {
            s.push_str("0/0");
        } else if !x[0] && x[1] {
            s.push_str("0/1");
        } else if x[0] && !x[1] {
            s.push_str("1/0");
        } else if x[0] && x[1] {
            s.push_str("1/1");
        }
        s.push('\t');
    }
    s
}
