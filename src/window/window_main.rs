use crate::core::core::MatrixWrapper;
use crate::core::helper::{Feature, merge_u32_to_u64};
use crate::r#mod::input_data::{read_paths, FileData};
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use clap::ArgMatches;
use log::info;
use std::net::UdpSocket;

/// Window function
///
/// Reading a ped file return "genotypes" which reflect windows over the entries
/// We assume that the entries that in variation graphs we have some kind of pan-genomic order in the order of the entries which reflect haplotypes
pub fn window_main(matches: &ArgMatches) {
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin2.mod");
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let mut mw = MatrixWrapper::new();
    let feature = mw.feature;
    mw.bfile_wrapper(plink_file);
    let mut mw = iterate_test(&mw, window);
    mw.feature = Feature::MWindow;
    mw.window_size = window;
    mw.make_counter();

    let chunk_size = (mw.matrix_bit.len() / split) + 1;
    let chunks = mw.matrix_bit.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        mw.write_fam(index, output_prefix, feature, len);
        mw.write_bed(index, output_prefix, feature, len);
        mw.write_bim(index, output_prefix, &feature, len);
    }
}

/// Wrapper around the matrix in sliding window
pub fn iterate_test(mw: &MatrixWrapper, window: usize) -> MatrixWrapper {
    let num_path = mw.matrix_bit[0].len() / 2;
    let mut mw_new = MatrixWrapper::new();
    for x in window..mw.matrix_bit.len() - window {
        let mut bv2 = Vec::new();
        for y in 0..num_path {
            let mut bv: Vec<[bool; 2]> = Vec::new();

            for z in x - window..x + window {
                bv.push([mw.matrix_bit[z][y * 2], mw.matrix_bit[z][y * 2 + 1]]);
            }
            bv2.push((y, bv));
        }
        bv2.sort_by(|a, b| a.1.cmp(&b.1));
        let f = get_index(&bv2);
        let flen = f.len();
        mw_new.matrix_bit.extend(f);
        mw_new
            .geno_names
            .extend((0..flen).map(|_| merge_u32_to_u64(x as u32, flen as u32)));
    }
    mw_new.shape = (mw_new.matrix_bit.len(), mw_new.matrix_bit[0].len());
    mw_new.fam_entries = mw.fam_entries.clone();
    return mw_new;
}

pub fn get_index(vv: &Vec<(usize, Vec<[bool; 2]>)>) -> Vec<BitVec<u8>> {
    let mut pp = Vec::new();
    let mut last = &vv[0].1;
    pp.push(vec![vv[0].0]);
    for x in vv.iter().skip(1) {
        if &x.1 != last {
            pp.push(vec![x.0]);
            last = &x.1;
        } else {
            pp.last_mut().unwrap().push(x.0);
        }
    }
    getbv(&pp, vv.len())
}

pub fn getbv(vv: &Vec<Vec<usize>>, len: usize) -> Vec<BitVec<u8>> {
    let mut gg = Vec::new();
    for x in vv {
        let mut oo: BitVec<u8, Lsb0> = BitVec::<u8, Lsb0>::repeat(false, len * 2);

        for y in x {
            oo.set(*y * 2, true);
            oo.set(*y * 2 + 1, true);
        }
        gg.push(oo);
    }
    gg
}

impl MatrixWrapper {
    pub fn make_counter(&mut self) {
        let mut last = &self.geno_names[0];
        let mut counter = 0;
        self.window_number.push(counter);

        for x in self.geno_names.iter().skip(1) {
            if x == last {
                self.window_number.push(counter);
                counter += 1;
            } else {
                counter = 0;
                self.window_number.push(counter);
                last = x;
            }
        }
    }
}
