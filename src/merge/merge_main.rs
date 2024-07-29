use crate::core::bfile::count_lines;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::remove::remove_main::copy_file;
use crate::window::window_main::{iterate_test, read_write_bim};
use clap::ArgMatches;
use log::info;
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};

/// Window function
///
/// Reading a ped file return "genotypes" which reflect windows over the entries
/// We assume that the entries that in variation graphs we have some kind of pan-genomic order in the order of the entries which reflect haplotypes
pub fn merge_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let plink_list = matches.value_of("plink").unwrap();
    let out_file = matches.value_of("output").unwrap();

    let input_list = read_list(plink_list)?;
    let names = clear_names(input_list)?;
    let fams = check_fams(&names)?;
    if fams == false {
        panic!("Fam files are not the same");
    }
    copy_file(
        &format!("{}{}", names[0], ".fam"),
        &format!("{}{}", out_file, ".fam"),
    )?;
    merge_bim(&names, &(out_file.to_string() + ".bim"))?;
    bed_merge(&names, &(out_file.to_string() + ".bed"))?;

    Ok(())
}

pub fn read_list(file: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut list: Vec<String> = Vec::new();
    let file = fs::File::open(file)?;
    let reader = io::BufReader::new(file);
    for line in reader.lines() {
        list.push(line?);
    }
    Ok(list)
}

pub fn clear_names(names: Vec<String>) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut new_names: Vec<String> = Vec::new();
    for x in names.iter() {
        if x.ends_with(".bed") {
            let a = x.split(".").collect::<Vec<&str>>();
            new_names.push(a[0..a.len() - 2].join("."));
        } else {
            new_names.push(x.clone());
        }
    }
    Ok(new_names)
}

pub fn check_fams(fams: &Vec<String>) -> Result<bool, Box<dyn std::error::Error>> {
    if fams.len() == 0 {
        return Ok(false);
    }
    let a = fs::read_to_string(fams[0].to_string() + ".fam")?;

    for x in fams.iter().skip(1) {
        let b = fs::read_to_string(x.to_string() + ".fam")?;
        if a != b {
            return Ok(false);
        }
    }
    Ok(true)
}

pub fn merge_bim(fams: &Vec<String>, output: &str) -> io::Result<()> {
    let output_file = "output.txt";

    // Create or truncate the output file
    let mut output = fs::File::create(output)?;

    // Read each file and write its content to the output file
    for file in fams {
        // Read the content of the current file
        let content = fs::read_to_string(file.to_string() + ".bim")?;

        // Write the content to the output file
        writeln!(output, "{}", content)?;
    }

    info!("Files have been concatenated into {}", output_file);

    Ok(())
}

pub fn bed_merge(files: &Vec<String>, output_file: &str) -> io::Result<()> {
    // Open the output file
    let mut output = File::create(output_file)?;

    // Copy the first file as-is
    if let Some(first_file) = files.first() {
        let mut input = File::open(first_file)?;
        io::copy(&mut input, &mut output)?;
    }

    // Process and append the remaining files
    for file in files.iter().skip(1) {
        let mut input = File::open(file)?;
        let mut buffer = Vec::new();
        input.read_to_end(&mut buffer)?;

        // Skip the first 3 bytes
        if buffer.len() > 3 {
            let remaining_bytes = &buffer[3..];
            output.write_all(remaining_bytes)?;
        }
    }

    println!("Files have been concatenated into {}", output_file);

    Ok(())
}
