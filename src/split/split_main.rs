use crate::core::bfile::count_lines;

use crate::remove::remove_main::copy_file;

use clap::ArgMatches;
use log::info;

use rayon::iter::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, Write};
use std::io::{BufReader, BufWriter};

/// # Window function
///
/// Reading a ped file return "genotypes" which reflect windows over the entries
/// We assume that the entries that in variation graphs we have some kind of pan-genomic order in the order of the entries which reflect haplotypes
pub fn split_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin split'");

    let plink_file = matches.value_of("plink").unwrap();
    let out_file = matches.value_of("output").unwrap();
    let number_splits = matches
        .value_of("splits")
        .unwrap()
        .parse::<usize>()
        .expect("Error parsing splits");
    let threads = matches
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .expect("Error parsing threads");

    info!("Splitting file: {}", plink_file);
    info!("Number of splits: {}", number_splits);
    info!("Output prefix: {}\n", out_file);

    let lines = count_lines(&format!("{}.bim", plink_file))?;
    let fam_lines = count_lines(&format!("{}.fam", plink_file))?;
    info!("Number of samples: {}", fam_lines);
    info!("Number of variants: {}", lines);

    info!(
        "Splitting PLINK BIM file: {}",
        format!("{}.bim", plink_file)
    );
    split_file(
        &format!("{}.bim", plink_file),
        out_file,
        number_splits,
        "bim",
        threads,
    )?;
    info!(
        "Splitting PLINK FAM file: {}",
        format!("{}.fam", plink_file)
    );
    split_fam(&format!("{}.fam", plink_file), number_splits, out_file)
        .expect("Error splitting FAM file");
    info!(
        "Splitting PLINK BED file: {}",
        format!("{}.bed", plink_file)
    );
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

/// # Split a plain-text file
///
/// Split a file into n files
/// Plain text - same length every file
fn split_file(
    filename: &str,
    output_prefix: &str,
    splits: usize,
    output_suffix: &str,
    threads: usize,
) -> io::Result<()> {
    let index = index_file(filename, splits);

    // rayon number of threads pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads) // Limit to 4 threads
        .build()
        .unwrap();

    pool.install(|| {
        index.par_iter().enumerate().for_each(|(i, x)| {
            let file_name = format!("{}.{}.{}", output_prefix, i + 1, output_suffix);
            let mut output_file =
                BufWriter::new(File::create(file_name).expect("Error creating file"));

            let input_file = File::open(filename).expect("Error opening file");
            let mut buffer_input = BufReader::new(input_file);
            // Seek to the correct position in the input file
            buffer_input
                .seek(std::io::SeekFrom::Start(x[0] as u64))
                .expect("Error seeking file");

            let mut pos = x[0];
            for line in buffer_input.lines() {
                let l = line.unwrap();
                pos += l.len() + 1;
                if pos > x[1] {
                    break;
                }
                writeln!(output_file, "{}", l).expect("Error writing to file");
            }
        });
    });

    Ok(())
}

/// # Split plink fam file
///
/// **Comment: We don't actually split the file, we just copy it n times**
pub fn split_fam(fam_file: &str, splits: usize, output_prefix: &str) -> io::Result<()> {
    for split_number in 1..=splits {
        let file_name = format!("{}.{}.fam", output_prefix, split_number);
        copy_file(fam_file, &file_name)?;
    }
    Ok(())
}

/// # Split a bed file
///
/// The bed file is split by n
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
    let mut start = 0;
    let total_len = File::open(filename)
        .expect("Not able to open file")
        .metadata()
        .expect("Metadata now working")
        .len()
        - 3;
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

/// # Index a file in equal part
///
/// By lines
/// In Bytes
pub fn index_file(file: &str, number: usize) -> Vec<[usize; 2]> {
    let file = File::open(file).expect("Error opening file");
    let buffreader = BufReader::new(file);

    let mut result = vec![0];
    let mut pos = 0;
    for line in buffreader.lines() {
        pos += line.unwrap().len() + 1;
        result.push(pos);
    }
    let step_size = result.len() / number;

    let oo = (0..number + 1)
        .map(|x| result[x * step_size])
        .collect::<Vec<usize>>();
    

    (1..oo.len())
        .map(|x| [oo[x - 1], oo[x]])
        .collect::<Vec<[usize; 2]>>()
}
