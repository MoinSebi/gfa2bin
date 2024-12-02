use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::cov::pack::{
    init_geno_names, init_matrix, matrick_pack_wrapper, read_pack_wrapper, wrapper_reader123,
};
use clap::ArgMatches;
use log::{info, warn};

use crate::core::helper::Feature::Alignment;

use packing_lib::core::core::PackCompact;
use packing_lib::core::reader::{read_index, unpack_zstd_to_byte};
use packing_lib::normalize::convert_helper::Method;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub fn cov_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin cov'");
    // You have either a list of packs (plain-text) or a compressed pack (cat or list), but you need to provide an index
    if matches.is_present("pack")
        || (matches.is_present("pack compressed") && matches.is_present("index"))
        || (matches.is_present("pc-list") && matches.is_present("index"))
    {
        info!("Files provided");
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
        .unwrap()
        .parse::<f32>()
        .expect("Could not parse fraction");
    let mut method = Method::from_str(matches.value_of("method").unwrap());

    let keep_zeros = matches.is_present("keep-zeros");
    let want_node = !matches.is_present("sequence");
    // Output modification

    // Output
    let bimbam = matches.is_present("bimbam");

    let mut dyna = true;
    if absolute_thresh > 0 {
        dyna = false;
        // if not, give method and fraction
    } else if method == Method::Nothing {
        panic!("You need to provide a method");
    } else if fraction == 0.0 || fraction < 0.0 {
        panic!("You need to provide a fraction");
    }

    let mut pheno = f64::MAX;
    if matches.is_present("pheno") {
        pheno = matches.value_of("pheno").unwrap().parse()?;
    }

    info!("Feature: {}", Alignment.to_string1());
    info!(
        "Absolute threshold: {}",
        if absolute_thresh == 0 {
            "None".to_string()
        } else {
            absolute_thresh.to_string()
        }
    );
    info!("Method: {}", method.to_string());
    info!("Fraction: {}", fraction);
    info!("Keep zeros: {}", keep_zeros);
    info!(
        "Dummy-Pheno: {}",
        if pheno == f64::MAX {
            "NA".to_string()
        } else {
            pheno.to_string()
        }
    );
    info!("Type: {}", if want_node { "Node" } else { "Sequence" });
    info!("Output format: {}", if bimbam { "bimbam" } else { "PLINK" });
    info!("Output prefix: {}\n", output_prefix);

    // Initialize the matrix wrapper
    let mut mw = MatrixWrapper::new();
    mw.feature = Feature::Alignment;

    info!("Reading the input");
    if matches.is_present("pack") {
        info!("Reading plain-text pack");
        // Read the first data
        let files_list = read_file_lines(pack_list.unwrap()).unwrap();
        let mut pack_first = read_pack_wrapper(true, &files_list[0][1]);
        init_geno_names(&mut mw, &mut pack_first, want_node, &vec![]);
        init_matrix(
            &mut mw,
            &mut pack_first,
            want_node,
            bimbam,
            files_list.len(),
        );



        for (index, x) in files_list.iter().enumerate() {
            let mut pc = PackCompact::parse_pack(&x[1]);

            if pc.node_index != pack_first.node_index {
                panic!("The pack files are not the same");
            }

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
                &x[0],
                absolute_thresh,
            );
        }
        // Compressed back (bin/u16, seq/node)
    } else {
        // Index of the file
        let index_file = read_index(matches.value_of("index").unwrap());
        // Compressed pack list
        if matches.is_present("pc-list") {
            info!("Reading pc-list");
            // Read the samples and path
            let cpack_list = read_file_lines(cpacklist.unwrap()).unwrap();

            // Read the first file
            let mut pack_first = read_pack_wrapper(false, &cpack_list[0][1]);

            if !pack_first.is_sequence && !want_node {
                panic!("The first file is not a sequence, but you want a node");
            }
            // Geno names based on the index
            init_geno_names(&mut mw, &mut pack_first, want_node, &index_file);

            //
            init_matrix(
                &mut mw,
                &mut pack_first,
                want_node,
                bimbam,
                cpack_list.len(),
            );
            let mut pack_first = read_pack_wrapper(false, &cpack_list[0][1]);

            for (index, x) in cpack_list.iter().enumerate() {
                let mut pc = PackCompact::read_wrapper(&x[1]);
                if pc.node_index != pack_first.node_index {
                    panic!("The pack files are not the same");
                }
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
                    &x[0],
                    absolute_thresh
                );
            }


            // Concatenated compressed pack
        } else if matches.is_present("pack compressed") {
            info!("Reading 'pack compressed'");
            let file_pack = cpack.unwrap();
            let buffer = unpack_zstd_to_byte(file_pack);
            let (
                _kind,
                _include_all,
                _bin,
                _method,
                _relative,
                _std,
                _thresh,
                bytes,
                _length,
                _name,
            ) = PackCompact::get_meta(&buffer);

            // Chunks
            let mut chunks = buffer.chunks((bytes + 86) as usize);
            let number_chunks = chunks.len();


            let mut pack_first = wrapper_reader123(chunks.next().unwrap());
            init_geno_names(&mut mw, &mut pack_first, want_node, &index_file);
            init_matrix(&mut mw, &mut pack_first, want_node, bimbam, number_chunks);
            let chunks = buffer.chunks((bytes + 86) as usize);
            pack_first.node_index = vec![];

            for (index, chunk) in chunks.enumerate() {
                let mut pc = wrapper_reader123(chunk);
                if pc.node_index != pack_first.node_index {
                    panic!("The pack files are not the same");
                }

                let name = pc.name.clone();
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
                    &name,
                    absolute_thresh
                );
            }

        }
    }

    // Sample index
    mw.sample_index_u16 = mw
        .sample_names
        .iter()
        .enumerate()
        .map(|x| [x.0, x.0])
        .collect();

    // Set feature style
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

/// Read a file and return each line in a vector
///
/// Check if entry is a path
fn read_file_lines(file_path: &str) -> io::Result<Vec<[String; 2]>> {
    // Open the file
    let file = File::open(file_path).expect("Can not open file");
    let reader = io::BufReader::new(file);

    // Create a vector to store the entries
    let mut entries: Vec<[String; 2]> = Vec::new();

    // Iterate over each line in the file
    for line in reader.lines() {
        // Add the line to the vector
        if let Ok(entry) = line {
            let entry_split = entry.split_whitespace().collect::<Vec<&str>>();
            if Path::is_file(Path::new(&entry_split[1])) {
                entries.push([entry_split[0].to_string(), entry_split[1].to_string()]);
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
