use crate::core::core::MatrixWrapper;

use crate::r#mod::input_data::{read_paths, FileData};
use clap::ArgMatches;
use log::{debug, info};


/// Function for "gfa2bin mod"
///
/// This function removed entries (SNPs) or path by name or index.
///
/// Input is a plink bed
pub fn remove_main(matches: &ArgMatches) {
    // Input parameters
    let plink_file = matches.value_of("plink").unwrap();

    // Output parameters
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.mod");
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();


    let mut mw = MatrixWrapper::new();
    mw.bfile_wrapper(plink_file);
    if matches.is_present("feature") || matches.is_present("paths") || matches.is_present("index") {
        debug!("One of the required parameters is present.");
    }

    if matches.is_present("index") {
        let index = FileData::from_file(matches.value_of("index").unwrap());
        let aa = index.data.iter().map(|x| *x as usize).collect();
        mw.remove_by_index(aa);
    }

    if matches.is_present("features") {
        let feature_file = matches.value_of("features").unwrap();
        let data = FileData::from_file(feature_file);
        if mw.feature != data.feature {
            panic!("Feature is not the same");
        }
        info!("Data: {:?}", data.data);
        info!("Matrix (before): {:?}", mw.matrix_bit.len());
        mw.remove_feature(&data);
        info!("Matrix (after): {:?}", mw.matrix_bit.len());
    }
    if matches.is_present("paths") {
        let paths = matches.value_of("paths").unwrap();
        let paths = read_paths(paths);
        mw.remove_samples(&paths);
    }

    if matches.is_present("non-info") {
        mw.remove_non_info();
    }

    let feature = mw.feature;
    info!("Writing the output");
    println!("{}", mw.bim_entries.len());

    let chunk_size = (mw.matrix_bit.len() / split) + 1;
    let chunks = mw.matrix_bit.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        //write_bed2(y, output_prefix, feature, index, len);
        mw.write_fam(index, output_prefix, feature, len);
        mw.write_bed(index, output_prefix, feature, len);
        mw.write_bim(index, output_prefix, &feature, len);
    }
}

use std::cmp::Ord;

pub fn remove_feature<T, U, V>(mut a: T, mut d: U) -> Vec<usize>
where
    T: Iterator<Item = V> + Clone,
    U: Iterator<Item = V>,
    V: Ord + Clone,
{
    let mut uu = Vec::new();
    let mut a_idx = 0;

    let mut a_next = a.next();
    let mut d_next = d.next();

    while let (Some(ai), Some(di)) = (a_next.clone(), d_next.clone()) {
        if ai < di {
            a_next = a.next();
            a_idx += 1;
        } else if ai > di {
            d_next = d.next();
        } else {
            // Element found in both iterators
            uu.push(a_idx);
            a_next = a.next();
            d_next = d.next();
            a_idx += 1;
        }
    }
    uu
}
