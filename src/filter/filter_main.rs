use crate::block::block_main::{blocks_node, node_size, wrapper_blocks};
use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, percentile, Feature};
use crate::r#mod::input_data::FileData;
use crate::subpath::subpath_main::{function1, gfa_index, subpath_wrapper};
use crate::window::window_main::iterate_test;
use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use hashbrown::HashMap;
use log::info;
use std::collections::HashSet;
use std::ffi::c_ushort;

/// Block main function
///
/// Easy block function
/// Extract the subpath from a graph for each node
pub fn filter_main(matches: &ArgMatches) {
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap();

    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    // Read the bed file
    let mut maf = 0.0;
    let mut MAF = 1.0;
    if matches.is_present("maf") {
        maf = matches.value_of("maf").unwrap().parse::<f64>().unwrap();
    }
    if matches.is_present("MAF") {
        MAF = matches.value_of("MAF").unwrap().parse::<f64>().unwrap();
    }
    info!("MAF: {}", maf);
    info!("MAF: {}", MAF);
    let mut mw = MatrixWrapper::new();
    let feature = mw.feature;
    mw.bfile_wrapper(plink_file);

    println!("Matrix size: {}", mw.matrix_bit.len());
    mw.filet_maf(maf, MAF);

    println!("Matrix size: {}", mw.matrix_bit.len());

    println!("Matrix size: {}", mw.matrix_bit[0].len());

    let mut a1 = 1;
    let mut a2 = usize::MAX;
    if matches.is_present("entry-max"){
        a2 = matches.value_of("entry-max").unwrap().parse::<usize>().unwrap();
    }
    if matches.is_present("entry-min"){
        a1 = matches.value_of("entry-min").unwrap().parse::<usize>().unwrap();
    }
    mw.filter_path(a1, a2);
    println!("Matrix size: {}", mw.matrix_bit[0].len());

    mw.write_chunks(split, output_prefix, mw.feature);
}

impl MatrixWrapper {
    pub fn filet_maf(&mut self, maf: f64, MAF: f64) {
        let mut b = Vec::new();
        for (i, (x, o)) in self
            .matrix_bit
            .iter()
            .zip(self.geno_names.iter_mut())
            .enumerate()
        {
            let mut count = 0;
            for y in x.iter() {
                if y == true {
                    count += 1;
                }
            }
            let maf1 = count as f64 / x.len() as f64;
            if maf1 < maf || maf1 > MAF {
                b.push(i);
            }
        }
        self.remove_by_index(b);
    }

    pub fn filter_path(&mut self, a1: usize, a2: usize) {
        let mut b = Vec::new();
        for i in 0..self.matrix_bit[0].len()/2 {
            let mut c = 0;
            for x in self.matrix_bit.iter() {
                if x[i * 2] == true || x[i * 2 + 1] == false {
                    c += 1;
                }
            }
            if c <= a1 || c >= a2 {
                b.push(i);
            }
        }


        self.remove_samples2(&b);
    }
}
