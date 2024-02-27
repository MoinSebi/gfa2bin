use crate::alignment::pack::matrix_pack_wrapper;
use crate::core::core::MatrixWrapper;
use clap::ArgMatches;
use log::info;
use packing_lib::convert::convert_helper::OutputType::Pack;
use packing_lib::core::core::PackCompact;
use packing_lib::core::reader::{get_meta, read_index, unpack_zstd_to_byte, wrapper_bool, wrapper_u16, zstd_decode};
use std::fs::File;
use std::io::{self, BufRead};
use std::ops::Index;

pub fn align_main(matches: &ArgMatches) {
    if matches.is_present("pack") | (matches.is_present("pack compressed") && matches.is_present("index")) | (matches.is_present("bfile")  && matches.is_present("index")){
        info!("Aligning");
    } else {
        info!("No input file");
        return;
    }

    let pack_list = matches.value_of("pack");
    let cpack = matches.value_of("pack compressed");
    let cpacklist = matches.value_of("bpacklist");

    let mut mw = MatrixWrapper::new();
    let mut index_normal: Vec<u32> = Vec::new();

    if matches.is_present("pack") {
        let file_pack = pack_list.unwrap();
        let files_list = read_file_lines(file_pack).unwrap();
        let mut pcs = Vec::new();
        for file in files_list {
            let pc = PackCompact::parse_pack(&file);
            pcs.push(pc);
        }
        matrix_pack_wrapper(&mut mw, &pcs, &pcs[0].node);

    } else {
        let index = read_index(matches.value_of("index").unwrap());


        if matches.is_present("pack compressed") {
            let file_pack = cpack.unwrap();
            let bytes = unpack_zstd_to_byte(file_pack);
            let j = get_meta(&bytes);
            let mut aa;
            if j.1 {
                info!("Reading bool pack");
                aa = wrapper_bool(&bytes);
            } else {
                info!("Reading u16 pack");
                aa = wrapper_u16(&bytes);
            }
            matrix_pack_wrapper(&mut mw, &aa, &index);
        }
        if matches.is_present("cpacklist") {
            let file_pack = cpacklist.unwrap();
            let a = read_file_lines(file_pack).unwrap();
            let mut aa = Vec::new();
            for x in a {
                let bytes = PackCompact::wrapp(&x);
                aa.push(bytes);
            }
            matrix_pack_wrapper(&mut mw, &aa, &index);
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
