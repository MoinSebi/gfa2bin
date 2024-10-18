use crate::core::bfile::count_lines;

use crate::remove::remove_main::copy_file;

use clap::ArgMatches;
use log::info;

use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::io::{BufReader, BufWriter};

/// Window function
///
/// Reading a ped file return "genotypes" which reflect windows over the entries
/// We assume that the entries that in variation graphs we have some kind of pan-genomic order in the order of the entries which reflect haplotypes
pub fn split_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let plink_file = matches.value_of("plink").unwrap();
    let out_file = matches.value_of("output").unwrap();
    let number_splits = matches.value_of("splits").unwrap().parse::<usize>().expect("Error parsing splits");

    info!("Splitting bim file: {}", format!("{}.bim", plink_file));

    let lines = count_lines(&format!("{}.bim", plink_file))?;
    let fam_lines = count_lines(&format!("{}.fam", plink_file))?;
    info!("Number of samples: {}", fam_lines);
    info!("Number of variants: {}", lines);
    // split bim
    info!("Splitting bim file: {}", format!("{}.bim", plink_file));
    split_file(
        &format!("{}.bim", plink_file),
        out_file,
        number_splits,
        lines,
    )?;
    info!("Splitting fam file: {}", format!("{}.fam", plink_file));
    split_fam(&format!("{}.fam", plink_file), number_splits, out_file)?;
    info!("Splitting bed file: {}", format!("{}.bed", plink_file));
    split_bed(
        &format!("{}.bed", plink_file),
        number_splits,
        fam_lines,
        lines,
        out_file,
    )?;
    info!("Done");
    Ok(())
}


/// Split a file into n files
/// Plain text - same length every file
fn split_file(
    filename: &str,
    output_prefix: &str,
    n: usize,
    number_lines: usize,
) -> io::Result<()> {
    let input_file = File::open(filename)?;
    let reader = BufReader::new(input_file);
    // Create the output files
    let mut output_files = Vec::with_capacity(n);
    for i in 0..n {
        let file_name = format!("{}.{}.bim", output_prefix, i + 1);
        let output_file = BufWriter::new(File::create(file_name)?);
        output_files.push(output_file);
    }

    // Write lines to the appropriate output file
    let mut current_file_index = 0;
    let mut line_count = 0;
    let lines_per_file = (number_lines + n - 1) / n; // Ceiling division

    for line in reader.lines() {
        let line = line?;
        writeln!(output_files[current_file_index], "{}", line)?;

        line_count += 1;
        if line_count >= lines_per_file {
            line_count = 0;
            current_file_index += 1;
            if current_file_index >= n {
                break;
            }
        }
    }

    Ok(())
}


/// "Split" fam file
///
/// Comment: We don't actually split the file, we just copy it n times
pub fn split_fam(fam_file: &str, splits: usize, output_prefix: &str) -> io::Result<()> {
    for x in 1..=splits {
        let file_name = format!("{}.{}.fam", output_prefix, x);
        copy_file(fam_file, &file_name)?;
    }
    Ok(())
}


/// Split a bed file into n files
///
/// The bed file is split by the number of samples
fn split_bed(
    filename: &str,
    n: usize,
    sample_size: usize,
    var_size: usize,
    output_prefix: &str,
) -> io::Result<()> {
    // Open the input file
    let mut input_file = File::open(filename)?;

    // Create output files and write the initial 3 zero bytes
    let mut output_files = Vec::with_capacity(n);
    for i in 0..n {
        let file_name = format!("{}.{}.bed", output_prefix, i + 1);
        let mut output_file = File::create(file_name)?;
        // Write 3 zero bytes at the beginning of each file
        output_file.write_all(&[108, 27, 1])?;
        output_files.push(output_file);
    }

    // Prepare a buffer for reading from the input file
    let sample_size = (sample_size as f64 / 4_f64).ceil() as usize;
    let mut bytes_per_file = sample_size * (var_size as f64 / n as f64).ceil() as usize;
    info!("Bytes per new file: {}", bytes_per_file);
    let mut start = 0;
    let total_len = File::open(filename)?.metadata()?.len() - 3;
    // Read from the input file and write to the output files
    for x in 0..n {
        if (total_len - start) < bytes_per_file as u64 {
            bytes_per_file = (total_len - start) as usize;
        }
        let mut buffer = vec![0; bytes_per_file];
        input_file.read_exact(&mut buffer)?;
        let mut output_file = &output_files[x];
        output_file.write_all(&buffer)?;
        start += bytes_per_file as u64;
    }

    Ok(())
}
