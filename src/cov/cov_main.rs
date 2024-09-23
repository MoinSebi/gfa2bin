use crate::cov::pack::{bin2bin, f32_to_bin, f32_to_f32, init_geno_names, init_matrix, matrick_pack_wrapper, matrix_pack_wrapper, read_pack1, remove_duplicates, wrapper_reader123};
use crate::core::core::MatrixWrapper;
use crate::core::helper::{index2node_seq, Feature};
use clap::ArgMatches;
use log::{info, warn};

use crate::core::helper::Feature::Alignment;

use bitvec::order::Msb0;
use bitvec::prelude::BitVec;
use packing_lib::core::core::{DataType, PackCompact};
use packing_lib::core::reader::{read_index, unpack_zstd_to_byte, wrapper_reader};
use packing_lib::normalize::convert_helper::Method;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub fn cov_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // You have either a list of packs (plain-text) or a compressed pack (cat or list), but you need to provide an index
    if matches.is_present("pack")
        | (matches.is_present("pack compressed") && matches.is_present("index"))
        | (matches.is_present("pc-list") && matches.is_present("index"))
    {
        info!("Using packlist");
    } else {
        if matches.is_present("pack compressed") {
            panic!("You need to provide an index file");
        }
        if matches.is_present("pc-list") {
            panic!("You need to provide an index file");
        } else {
            panic!("You need to provide a pack file");
        }
    }

    // Read input files
    let pack_list = matches.value_of("pack");
    let cpack = matches.value_of("pack compressed");
    let cpacklist = matches.value_of("pc-list");
    let output_prefix = matches.value_of("output").unwrap();

    // Normalize the rows
    let absolute_thresh = matches
        .value_of("absolute-threshold")
        .unwrap_or("0")
        .parse::<u32>()
        .unwrap();
    let mut fraction = matches
        .value_of("fraction")
        .unwrap_or("0.0")
        .parse::<f32>()
        .unwrap();
    let mut method = Method::from_str(matches.value_of("method").unwrap_or("nothing"));
    let keep_zeros = matches.is_present("keep-zeros");
    let want_node = matches.is_present("node");
    // Output modification
    let bimbam = matches.is_present("bimbam");

    if !matches.is_present("fraction")
        && method == Method::Nothing
        && !matches.is_present("absolute-threshold")
    {
        info!("No method or fraction given, using default");
        method = Method::Percentile;
        fraction = 0.1;
    }

    if absolute_thresh > 0 {
        method = Method::Absolute;
        fraction = 0.0;
        // if not, give method and fraction
    } else if method == Method::Nothing && fraction != 0.0 {
        warn!("No method or fraction given");
        panic!("Exiting");
    } else if fraction == 0.0 {
        warn!("Relative threshold is 0");
        panic!("Exiting");
    } else if method == Method::Nothing {
        warn!("No method or fraction given");
        panic!("Exiting");
    }


    let mut pheno = f64::MAX;
    if matches.is_present("pheno") {
        pheno = matches.value_of("pheno").unwrap().parse()?;
    }

    info!("Input parameters");
    info!("Feature: {}", Alignment.to_string1());
    info!("Absolute threshold: {}", absolute_thresh);
    info!("Method: {}", method.to_string());
    info!("Relative threshold: {:?}", fraction);
    info!("Standard deviation: {}", 0.0);
    info!("Keep zeros: {}", keep_zeros);
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

    info!("Reading the input");
    if matches.is_present("pack") {
        // Read the first data
        let files_list = read_file_lines(pack_list.unwrap()).unwrap();
        let mut first_file = read_pack1(true, &files_list[0]);
        init_geno_names(&mut mw, &mut first_file, want_node, &vec![]);
        init_matrix(&mut mw, &mut first_file, want_node, bimbam, files_list.len());

        for (index, x) in files_list.iter().enumerate() {
            let mut pc = PackCompact::parse_pack(&x);
            matrick_pack_wrapper(
                &mut mw,
                &mut pc,
                want_node,
                keep_zeros,
                fraction,
                method,
                bimbam,
                &vec![],
                index,
            );
        }
        // Compressed back (bin/u16, seq/node)
    } else {
        // Index of the file
        let index_file = read_index(matches.value_of("index").unwrap());
        // Compressed pack list
        if matches.is_present("pc-list") {
            let cpack_list = read_file_lines(cpacklist.unwrap()).unwrap();
            let mut first_file = read_pack1(false, &cpack_list[0]);
            init_geno_names(&mut mw, &mut first_file, want_node, &index_file);
            init_matrix(&mut mw, &mut first_file, want_node, bimbam, cpack_list.len());
            for (index, x) in cpack_list.iter().enumerate() {
                let mut pc = PackCompact::read_wrapper(&x);
                matrick_pack_wrapper(
                    &mut mw,
                    &mut pc,
                    want_node,
                    keep_zeros,
                    fraction,
                    method,
                    bimbam,
                    &index_file,
                    index,
                );
            }
        } else if matches.is_present("pack compressed") {
            println!("dashkdjash");
            let file_pack = cpack.unwrap();
            let buffer = unpack_zstd_to_byte(file_pack);
            let (_kind, _include_all, _bin, _method, _relative, _std, _thresh, bytes, _length, _name) =
                PackCompact::get_meta(&buffer);
            let mut chunks = buffer.chunks((bytes + 86) as usize);
            let number_chunks = chunks.len();
            let mut first_file = wrapper_reader123(&chunks.nth(0).unwrap());
            init_geno_names(&mut mw, &mut first_file, want_node, &index_file);
            init_matrix(&mut mw, &mut first_file, want_node, bimbam, number_chunks);
            let mut chunks = buffer.chunks((bytes + 86) as usize);

            for (index, chunk) in chunks.into_iter().enumerate() {
                println!("dasjkhdjkas");
                let mut pc = wrapper_reader123(&chunk);
                println!("joö {}", pc.normalized_coverage.len());
                matrick_pack_wrapper(
                    &mut mw,
                    &mut pc,
                    want_node,
                    keep_zeros,
                    fraction,
                    method,
                    bimbam,
                    &index_file,
                    index,
                );
            }
        }
    }

    mw.sample_index_u16 = mw
        .sample_names
        .iter()
        .enumerate()
        .map(|x| [x.0, x.0])
        .collect();

    let feature_enum = Feature::Alignment;

    let thresh = mw
        .matrix_f32
        .iter()
        .map(|x| x.iter().cloned().fold(0.0, f32::max))
        .collect();


    mw.write_wrapper(
        bimbam,
        1,
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
