use clap::ArgMatches;
use log::info;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::r#mod::input_data::{FileData, read_paths};

pub fn mod_main(matches: &ArgMatches) {
    let plink_file = matches.value_of("plink");

    let mut mw = MatrixWrapper::new();
    mw.bfile_wrapper(plink_file.unwrap());

    if (matches.is_present("feature") || matches.is_present("paths")) || matches.is_present("graph"){
        info!("All good");
    }

    if matches.is_present("feature"){
        let feature = Feature::from_str(matches.value_of("feature").unwrap());
        let feature_file = matches.value_of("feature").unwrap();
        let mut data = FileData::from_file(feature_file);

        if mw.feature != feature {
            panic!("Feature is not the same");
        }

        //mw.remove_feature(&data)
    }
    if matches.is_present("paths"){
        let paths = matches.value_of("paths").unwrap();
        let paths = read_paths(paths);
        mw.remove_paths(&paths);
    }





}
