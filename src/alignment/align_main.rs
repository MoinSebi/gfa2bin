use crate::alignment::pack::matrix_pack_wrapper;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use clap::ArgMatches;
use log::{info, warn};

use crate::core::helper::Feature::Alignment;
use crate::graph::parser::weight_add;
use packing_lib::core::core::PackCompact;
use packing_lib::core::reader::{read_index, unpack_zstd_to_byte, wrapper_reader};
use packing_lib::normalize::convert_helper::Method;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub fn align_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // You have either a list of packs (plain-text) or a compressed pack (cat or list), but you need to provide an index
    if matches.is_present("pack")
        | (matches.is_present("pack compressed") && matches.is_present("index"))
        | (matches.is_present("bfile") && matches.is_present("index"))
    {
        info!("Aligning");
    } else {
        if matches.is_present("pack compressed") {
            panic!("You need to provide an index file");
        }
        if matches.is_present("bfile") {
            panic!("You need to provide an index file");
        } else {
            panic!("You need to provide a pack file");
        }
    }

    // Read input files
    let pack_list = matches.value_of("pack");
    let cpack = matches.value_of("pack compressed");
    let cpacklist = matches.value_of("bpacklist");
    let output_prefix = matches.value_of("output").unwrap();

    // Normalize the rows
    let absolute_thresh = matches
        .value_of("absolute-threshold")
        .unwrap_or("0")
        .parse::<u32>()
        .unwrap();
    let fraction = matches
        .value_of("fraction")
        .unwrap_or("1.0")
        .parse::<f32>()
        .unwrap();
    let method = Method::from_str(matches.value_of("method").unwrap_or("nothing"));
    let std = matches
        .value_of("std")
        .unwrap_or("0.0")
        .parse::<f32>()
        .unwrap();
    let keep_zeros = matches.is_present("keep-zeros");

    // Output modification
    let bimbam = matches.is_present("bimbam");
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    let mut pheno = f64::MAX;
    if matches.is_present("pheno") {
        pheno = matches.value_of("pheno").unwrap().parse()?;
    }
    info!("Input parameters");
    info!("Feature: {}", Alignment.to_string1());
    info!("Absolute threshold: {}", absolute_thresh);
    info!("Method: {}", method.to_string());
    info!("Relative threshold: {:?}", fraction);
    info!("Standard deviation: {}", std);
    info!("Keep zeros: {}", keep_zeros);
    info!("Split: {}", split);
    info!("Output format: {}", if bimbam { "bimbam" } else { "plink" });
    info!("Output prefix: {}", output_prefix);
    info!(
        "Dummy-Pheno: {}",
        if pheno == f64::MAX {
            "NA".to_string()
        } else {
            pheno.to_string()
        }
    );

    // Initialize the matrix wrapper
    let mut mw = MatrixWrapper::new();
    mw.feature = Feature::Alignment;
    // "Normal" pack file

    info!("Reading the input");
    if matches.is_present("pack") {
        let files_list = read_file_lines(pack_list.unwrap()).unwrap();
        let mut pcs = Vec::new();
        for file in files_list {
            let pc = PackCompact::parse_pack(&file);
            pcs.push(pc);
        }
        matrix_pack_wrapper(&mut mw, &pcs, &pcs[0].node_index);

        // Compressed back (bin/u16, seq/node)
    } else {
        // Index of the file
        let index = read_index(matches.value_of("index").unwrap());
        // Compressed pack
        if matches.is_present("pack compressed") {
            let file_pack = cpack.unwrap();
            let bytes = unpack_zstd_to_byte(file_pack);
            let pc_vec: Vec<PackCompact> = wrapper_reader(&bytes);
            matrix_pack_wrapper(&mut mw, &pc_vec, &index);
        }
        // Compressed pack list
        if matches.is_present("cpacklist") {
            let cpack_list = read_file_lines(cpacklist.unwrap()).unwrap();
            let mut pc_vec = Vec::new();
            for x in cpack_list {
                let bytes = PackCompact::read_wrapper(&x);
                pc_vec.push(bytes);
            }
            matrix_pack_wrapper(&mut mw, &pc_vec, &index);
        }
    }

    mw.sample_index_u16 = mw
        .sample_names
        .iter()
        .enumerate()
        .map(|x| [x.0, x.0])
        .collect();

    let feature_enum = Feature::Alignment;
    let mut thresh = Vec::new();
    if mw.matrix_f32.is_empty() {
        for x in mw.matrix_u16.iter() {
            let mut b = x.clone();
            thresh.push(PackCompact::threshold(
                &mut b,
                keep_zeros,
                fraction,
                std,
                method,
            ));
        }
    } else {
        for x in mw.matrix_f32.iter() {
            let mut b = x.clone();
            thresh.push(PackCompact::threshold(
                &mut b,
                keep_zeros,
                fraction,
                std,
                method,
            ));
        }
    }
    mw.write_wrapper(
        bimbam,
        split,
        output_prefix,
        thresh,
        feature_enum,
        pheno,
        !keep_zeros,
    );
    Ok(())
}

/// Read a file and return a vector of lines
///
/// They contains files
fn read_file_lines(file_path: &str) -> io::Result<Vec<String>> {
    // Open the file
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    // Create a vector to store the entries
    let mut entries: Vec<String> = Vec::new();

    // Iterate over each line in the file
    for line in reader.lines() {
        // Add the line to the vector
        if let Ok(entry) = line {
            if Path::is_file(Path::new(&entry)) {
                entries.push(entry);
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("{} is not a file", entry),
                ));
            }
        }
    }

    // Return the vector of entries
    Ok(entries)
}
