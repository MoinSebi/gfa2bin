use crate::alignment::pack::matrix_pack_wrapper;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use clap::ArgMatches;
use log::info;

use packing_lib::core::core::PackCompact;
use packing_lib::core::reader::{
    get_meta, read_index, unpack_zstd_to_byte, wrapper_bool, wrapper_u16,
};
use std::fs::File;
use std::io::{self, BufRead};


pub fn align_main(matches: &ArgMatches) {
    // You have either a pack or a compressed pack (cat or list), but you need to provide a index
    if matches.is_present("pack")
        | (matches.is_present("pack compressed") && matches.is_present("index"))
        | (matches.is_present("bfile") && matches.is_present("index"))
    {
        info!("Aligning");
    } else {
        info!("No input file");
        return;
    }

    // Read input files
    let pack_list = matches.value_of("pack");
    let cpack = matches.value_of("pack compressed");
    let cpacklist = matches.value_of("bpacklist");
    let output_prefix = matches.value_of("output").unwrap();
    let bimbam = matches.is_present("bimbam");
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    // Initialize the matrix wrapper
    let mut mw = MatrixWrapper::new();

    // "Normal" pack file
    if matches.is_present("pack") {
        let files_list = read_file_lines(pack_list.unwrap()).unwrap();
        let mut pcs = Vec::new();
        for file in files_list {
            let pc = PackCompact::parse_pack(&file);
            pcs.push(pc);
        }
        matrix_pack_wrapper(&mut mw, &pcs, &pcs[0].node);

        // Compressed back (bin/u16, seq/node)
    } else {
        // Index of the file
        let index = read_index(matches.value_of("index").unwrap());

        // Compressed pack
        if matches.is_present("pack compressed") {
            let file_pack = cpack.unwrap();
            let bytes = unpack_zstd_to_byte(file_pack);
            let meta_data = get_meta(&bytes);
            let pc_vec;
            if meta_data.1 {
                info!("Reading bool pack");
                pc_vec = wrapper_bool(&bytes);
            } else {
                info!("Reading u16 pack");
                pc_vec = wrapper_u16(&bytes);
            }
            matrix_pack_wrapper(&mut mw, &pc_vec, &index);
        }
        // Compressed pack list
        if matches.is_present("cpacklist") {
            let cpack_list = read_file_lines(cpacklist.unwrap()).unwrap();
            let mut pc_vec = Vec::new();
            for x in cpack_list {
                let bytes = PackCompact::wrapp(&x);
                pc_vec.push(bytes);
            }
            matrix_pack_wrapper(&mut mw, &pc_vec, &index);
        }
    }

    //mw.remove_non_info();

    let feature_enum = Feature::Alignment;
    if bimbam {
        info!("Writing the bimbam");
        let chunk_size = (mw.matrix_u16.len() / split) + 1;
        let chunks = mw.matrix_u16.chunks(chunk_size);
        let len = chunks.len();

        for (index, _y) in chunks.enumerate() {
            mw.write_bimbam(index, output_prefix, &feature_enum, len, &vec![0.0; 0]);
            mw.write_phenotype_bimbam(index, output_prefix, len);
        }
    } else {
        // Output
        info!("Writing the output");
        let chunk_size = (mw.matrix_bit.len() / split) + 1;
        let chunks = mw.matrix_bit.chunks(chunk_size);

        let len = chunks.len();
        for (index, _y) in chunks.enumerate() {
            //write_bed2(y, output_prefix, feature, index, len);
            mw.write_fam(index, output_prefix, feature_enum, len);
            mw.write_bed(index, output_prefix, feature_enum, len);
            mw.write_bim(index, output_prefix, &feature_enum, len);
        }
    }
}

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
            entries.push(entry);
        }
    }

    // Return the vector of entries
    Ok(entries)
}
