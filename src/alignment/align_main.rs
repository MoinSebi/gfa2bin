use crate::alignment::pack::{matrix_pack_bit_v2, matrix_pack_u16_v2};
use crate::core::core::MatrixWrapper;
use crate::helper::get_thresh;
use clap::ArgMatches;
use log::info;

pub fn align_main(matches: &ArgMatches) {
    if matches.is_present("pack") | matches.is_present("bpack") | matches.is_present("bfile") {
        info!("Aligning");
    } else {
        info!("No input file");
        return;
    }
    let _pack = matches.value_of("pack");
    let _bpack = matches.value_of("bpack");
    let _bpacklist = matches.value_of("bpacklist");

    let mut mw = MatrixWrapper::new();
    let mut index_normal: Vec<u32> = Vec::new();

    if matches.is_present("pack") {
        let file_pack = matches.value_of("pack").unwrap();
        let j = get_thresh(file_pack);
        if j == 0 {
            info!("Reading u16 pack");
            matrix_pack_u16_v2(file_pack, &mut mw, &mut index_normal);
        } else {
            info!("Reading bit pack");
            matrix_pack_bit_v2(file_pack, &mut mw, &mut index_normal);
        }
    }


    // // We only transpose once!
    // if matrix.matrix_bin.is_empty(){
    //     let k: Vec<String> = matrix.core_names.values().cloned().collect();
    //     let mut hit = 0;
    //
    //
    //
    //     matrix.matrix_core = transpose_generic(&matrix.matrix_core);
    //     matrix.shape = (matrix.matrix_core.len(), matrix.matrix_core[0].len());
    // } else {
    //     if ! matrix.transposed {
    //         matrix.matrix_bin = transpose_bitvec(&matrix.matrix_bin);
    //         matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
    //     }
    // };
    //
    // info!("Shape is {:?}", matrix.shape);
    //
    //
    //
    //
    //
    //
    //
    // if matches.is_present("names") {
    //     matrix.write_names(_output);
    // }
    //
    //
    //
    // // Make dummy binary
    // if matrix.matrix_bin.is_empty() {
    //     info!("Make binary");
    //     matrix.make_binary(1);
    // }
}
