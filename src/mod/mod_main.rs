use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::r#mod::input_data::{read_paths, FileData};
use clap::ArgMatches;
use log::info;

pub fn mod_main(matches: &ArgMatches) {
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.mod");
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let mut mw = MatrixWrapper::new();
    mw.bfile_wrapper(plink_file);
    if matches.is_present("feature") || matches.is_present("paths") {
        info!("All good");
    }

    if matches.is_present("features") {
        let feature_file = matches.value_of("features").unwrap();
        let data = FileData::from_file(feature_file);

        if mw.feature != data.feature {
            panic!("Feature is not the same");
        }
        mw.remove_feature(&data);
    }
    if matches.is_present("paths") {
        let paths = matches.value_of("paths").unwrap();
        let paths = read_paths(paths);
        mw.remove_paths(&paths);
    }
    let feature = mw.feature;
    info!("Writing the output");

    let chunk_size = (mw.matrix_bin.len() / split) + 1;
    let chunks = mw.matrix_bin.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        //write_bed2(y, output_prefix, feature, index, len);
        mw.write_fam(index, output_prefix, feature, len);
        mw.write_bed(index, output_prefix, feature, len);
        mw.write_bim(index, output_prefix, &feature, len);
    }
}
