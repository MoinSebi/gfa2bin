use crate::core::bfile::count_lines;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use clap::ArgMatches;
use log::{info, warn};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Error, Write};

/// Function for "gfa2bin remove"
///
/// This function removed entries (SNPs) or path by name or index.
///
/// Input is a single plink (bed, bim, fam) file.
pub fn remove_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Input parameters
    let plink_file = matches.value_of("plink").unwrap();

    // Output parameters
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.remove");

    let mut mw = MatrixWrapper::new();
    let bim_count = count_lines(&format!("{}{}", plink_file, ".bim"))?;
    let fam_count = count_lines(&format!("{}{}", plink_file, ".fam"))?;

    // Read the plink file
    mw.read_bed(&format!("{}{}", plink_file, ".bed"), fam_count, bim_count)?;

    if !(matches.is_present("genotypes")
        || matches.is_present("samples")
        || matches.is_present("index"))
    {
        warn!("One of the required parameters is present.");
    }

    // Feature based removement
    if matches.is_present("gindex") || matches.is_present("genotypes") {
        if matches.is_present("gindex") && matches.is_present("genotypes") {
            panic!("You can't use both index and feature at the same time.");
        } else if matches.is_present("gindex") {
            info!("Removing by index.");
            let index = read_file_to_vector(matches.value_of("gindex").unwrap()).unwrap();
            let aa = index.iter().map(|x| x.parse::<usize>().unwrap()).collect();
            process_file2(
                &format!("{}{}", plink_file, ".bim"),
                &format!("{}{}", output_prefix, ".bim"),
                &aa,
            )?;
            mw.remove_by_index(&aa);
        } else if matches.is_present("genotypes") {
            info!("Removing by genotype name.");
            let index = read_file_to_vector(matches.value_of("genotypes").unwrap())?;
            let f = process_file(
                &format!("{}{}", plink_file, ".bim"),
                &format!("{}{}", output_prefix, ".bim"),
                &index.iter().cloned().collect(),
                3,
            )?;
            mw.remove_by_index(&f);
        }

        copy_file(
            &format!("{}{}", plink_file, ".fam"),
            &format!("{}{}", output_prefix, ".fam"),
        )?;
    }

    // If samples are to be removed
    if matches.is_present("samples") || matches.is_present("sindex") {
        if matches.is_present("samples") && matches.is_present("sindex") {
            panic!("You can't use both paths and pindex at the same time.");
        } else if matches.is_present("samples") {
            let index = read_file_to_vector(matches.value_of("samples").unwrap())?;
            let f = process_file(
                &format!("{}{}", plink_file, ".fam"),
                &format!("{}{}", output_prefix, ".fam"),
                &index.iter().cloned().collect(),
                0,
            )?;
            mw.remove_index_samples(&f);
        } else if matches.is_present("sindex") {
            let index = read_file_to_vector(matches.value_of("sindex").unwrap()).unwrap();
            let aa = index.iter().map(|x| x.parse::<usize>().unwrap()).collect();
            process_file2(
                &format!("{}{}", plink_file, ".fam"),
                &format!("{}{}", output_prefix, ".fam"),
                &aa,
            )?;
            mw.remove_by_index(&aa);
        }

        if !matches.is_present("genotypes") && !matches.is_present("gindex") {
            copy_file(
                &format!("{}{}", plink_file, ".bim"),
                &format!("{}{}", output_prefix, ".bim"),
            )?;
        }
    }

    mw.write_bed(0, output_prefix, Feature::Node, 1);
    Ok(())
}
use std::fs;
use std::io;

/// Copy a file from filename1 to filename2
pub fn copy_file(filename1: &str, filename2: &str) -> io::Result<()> {
    // Read the contents of filename1
    let contents = fs::read(filename1)?;

    // Write the contents to filename2
    fs::write(filename2, contents)?;

    Ok(())
}

/// Read a file to a vector of strings
fn read_file_to_vector(file_path: &str) -> Result<Vec<String>, Error> {
    let file = File::open(file_path)?;
    let lines: Result<Vec<String>, _> = BufReader::new(file).lines().collect();
    lines
}

pub fn process_file(
    input_file: &str,
    output_file: &str,
    hs: &HashSet<String>,
    index_check: usize,
) -> Result<Vec<usize>, std::io::Error> {
    // Open input file for reading
    let file_in = File::open(input_file)?;
    let reader = BufReader::new(file_in);

    // Open output file for writing
    let file_out = File::create(output_file)?;
    let mut writer = BufWriter::new(file_out);

    let mut index1 = Vec::new();
    // Iterate over lines in input file
    for (index, line) in reader.lines().enumerate() {
        let line = line?; // unwrap the line or propagate error
        let line_split = line.split_whitespace().collect::<Vec<&str>>();
        // 3 for bim, 0 for fam
        if !hs.contains(line_split[index_check]) {
            writeln!(writer, "{}", line)?;
        } else {
            index1.push(index);
        }
    }
    info!("{:?} entries have been removed.\n", index1);

    // Flush the buffer to ensure all data is written to file
    writer.flush()?;

    Ok(index1)
}

/// Remove entries by index from the matrix
pub fn process_file2(
    input_file: &str,
    output_file: &str,
    index_vec: &Vec<usize>,
) -> Result<(), std::io::Error> {
    // Open input file for reading
    let file_in = File::open(input_file)?;
    let reader = BufReader::new(file_in);

    // Open output file for writing
    let file_out = File::create(output_file)?;
    let mut writer = BufWriter::new(file_out);

    let mut i = 0;
    let mut lines = reader.lines();
    let mut current_line_index = 0;

    // Iterate over lines in input file
    // Read lines from the reader
    for line_result in lines.by_ref() {
        let line = line_result?;

        // Check if we need to skip this line based on index_vec
        if i < index_vec.len() && current_line_index == index_vec[i] {
            i += 1;
            current_line_index += 1; // Increment line index for the next line
            continue; // Skip writing this line
        }

        // Write the line to the writer
        writeln!(writer, "{}", line)?;

        // Increment the line index after writing
        current_line_index += 1;

        // Check if we've processed all indices in index_vec
        if i == index_vec.len() {
            break; // No need to process further if we've handled all indices
        }
    }

    // Write any remaining lines from the reader
    for line_result in lines {
        let line = line_result?;
        writeln!(writer, "{}", line)?;
    }

    // Flush the buffer to ensure all data is written to file
    writer.flush()?;

    Ok(())
}

impl MatrixWrapper {
    //---------------------------------------------------------------------------------------------
    // Remove samples

    /// Remove samples by name

    /// Remove samples by index
    pub fn remove_index_samples(&mut self, to_be_removed: &Vec<usize>) {
        let mut i = 0;
        while i < self.matrix_bit.len() {
            let a = &mut self.matrix_bit[i];
            let mut f = 0;
            for x in to_be_removed.iter() {
                a.remove(*x * 2 - f * 2 + 1);
                a.remove(*x * 2 - f * 2);
                f += 1;
            }
            i += 1;
        }
    }

    //---------------------------------------------------------------------------------------------
    // Remove entries

    /// Remove entries by index from the matrix
    ///
    /// Removed in:
    /// - genome names
    /// - matrix bit
    /// - bim entries
    pub fn remove_by_index(&mut self, index: &Vec<usize>) {
        if index.is_empty() {
            return;
        }

        // Sort indices in ascending order to ensure we remove from end to start correctly
        let mut sorted_index = index.clone();
        sorted_index.sort_unstable();

        // Remove elements from self.matrix_bit based on indices in sorted_index
        let mut removal_count = 0;
        let mut rmi = self.matrix_bit.len() - 1; // Start from the end of the vector

        for &x in sorted_index.iter().rev() {
            let remove_index = rmi - x; // Calculate the actual index to remove
            if remove_index < self.matrix_bit.len() {
                self.matrix_bit.remove(remove_index);
                removal_count += 1;
            }
            rmi -= 1;
        }

        // Optionally, clear the entire vector if all elements are to be removed
        if removal_count == self.matrix_bit.len() {
            self.matrix_bit.clear();
        }
    }
}
