use crate::remove::remove_main::copy_file;



use clap::ArgMatches;
use log::info;
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};

/// # Merge main
///
/// Merge multiple PLINK files togther. This includes BED, BIM and FAM files
///
/// Comment: Fam files are only checked if they contain the same content, bim files are simply concatenated, and BED files are trimmed ([3:]) and concatenated
pub fn merge_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin merge'");

    let plink_list = matches.value_of("bed-list").unwrap();
    let out_file = matches.value_of("output").unwrap();

    info!("BED file list: {}", plink_list);
    info!("Output prefix: {}\n", out_file);

    let input_list = read_list(plink_list)?;
    let names = clear_names(input_list)?;

    let fams = check_fams(&names)?;

    if !fams {
        panic!("Fam files are not the same");
    }

    info!("Merging FAM files");
    copy_file(
        &format!("{}{}", names[0], ".fam"),
        &format!("{}{}", out_file, ".fam"),
    )?;

    info!("Merge BIM files");
    merge_bim(&names, &(out_file.to_string() + ".bim"))?;

    info!("Merge BED files");
    merge_bed(&names, &(out_file.to_string() + ".bed"))?;

    Ok(())
}

/// # Read a file (line by line
///
/// Each line one path
pub fn read_list(file: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut list: Vec<String> = Vec::new();
    let file = fs::File::open(file).expect("Could not open file");
    let reader = io::BufReader::new(file);
    for line in reader.lines() {
        list.push(line?);
    }
    Ok(list)
}

/// # Remove suffix from the string
///
/// Comment: If bed files are the input, remove those to get the prefix name
pub fn clear_names(names: Vec<String>) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut new_names: Vec<String> = Vec::new();
    for x in names.iter() {
        if x.ends_with(".bed") {
            let a = x.split('.').collect::<Vec<&str>>();
            new_names.push(a[0..a.len() - 1].join("."));
        }
    }
    Ok(new_names)
}

/// # FAM file checker
pub fn check_fams(fams: &Vec<String>) -> Result<bool, Box<dyn std::error::Error>> {
    if fams.is_empty() {
        return Ok(false);
    }
    let mut fam_content_first = "".to_string();
    if !fams.is_empty() {
        fam_content_first += &fs::read_to_string(fams[0].to_string() + ".fam")
            .unwrap_or_else(|_| panic!("Could not read FAM file {}", fams[0].to_string() + ".fam"));
    }
    for fam_entry in fams.iter().skip(1) {
        let fam_content = fs::read_to_string(fam_entry.to_string() + ".fam")?;
        if fam_content_first != fam_content {
            return Ok(false);
        }
    }
    Ok(true)
}

/// # Merge multiple plain-text files
///
/// Here - Merge bim (PLINK) file
pub fn merge_bim(fams: &Vec<String>, output_file: &str) -> io::Result<()> {
    // Create or truncate the output file
    let mut output = fs::File::create(output_file).expect("Failed to create output file");

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

/// # Merge BED files
///
/// Start with the first entry as a base, then add the others as well. Adding the rest is done bz removing the three header bytes.
pub fn merge_bed(files: &Vec<String>, output_file: &str) -> io::Result<()> {
    // Open the output file
    let mut output = File::create(output_file).expect("Failed to create output file");

    // Copy the first file as-is
    if let Some(first_file) = files.first() {
        let mut input = File::open(first_file.to_string() + ".bed")?;
        io::copy(&mut input, &mut output)?;
    }

    // Process and append the remaining files
    for file in files.iter().skip(1) {
        let mut input = File::open(file.to_string() + ".bed")?;
        let mut buffer = Vec::new();
        input.read_to_end(&mut buffer)?;

        // Skip the first 3 bytes
        if buffer.len() > 3 {
            let remaining_bytes = &buffer[3..];
            output
                .write_all(remaining_bytes)
                .expect("Failed to write to output file");
        }
    }

    Ok(())
}
