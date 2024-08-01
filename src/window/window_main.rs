use crate::core::bfile::count_lines;
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::remove::remove_main::copy_file;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use clap::ArgMatches;


use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

/// Window function
///
/// Reading a bed file and return "genotypes" which reflect windows over multiple entries
/// We assume that the entries that are present in variation graphs have some kind of pan-genomic order. Otherwise, this makes not that much sense.
pub fn window_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let plink_file = matches.value_of("plink").unwrap();
    let out_file = matches.value_of("output").unwrap();
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let mut block = None;
    if let Some(blocks_path) = matches.value_of("blocks") {
        let file = File::create(blocks_path)?;
        block = Some(BufWriter::new(file));
    }

    // Read the bed file
    let mut mw = MatrixWrapper::new();
    let bim_count = count_lines(&format!("{}{}", plink_file, ".bim"))?;
    let fam_count = count_lines(&format!("{}{}", plink_file, ".fam"))?;
    mw.read_bed(&format!("{}{}", plink_file, ".bed"), fam_count, bim_count)?;
    let (mw, index) = iterate_test(&mw, window, &mut block)?;

    mw.write_bed(0, out_file, Feature::Node, 1);
    read_write_bim(
        &mw,
        &index,
        &format!("{}{}", plink_file, ".bim"),
        &format!("{}{}", out_file, ".bim"),
        window,
    )
    .unwrap();
    copy_file(
        &format!("{}{}", plink_file, ".fam"),
        &format!("{}{}", out_file, ".fam"),
    )
    .unwrap();

    Ok(())
}

/// Wrapper around the matrix in sliding window
///
/// Iterate over the matrix
/// No bim entries
pub fn iterate_test(
    mw: &MatrixWrapper,
    window: usize,
    blocks: &mut Option<BufWriter<File>>,
) -> Result<(MatrixWrapper, Vec<usize>), std::io::Error> {
    let num_path = mw.matrix_bit[0].len() / 2;
    let mut index = Vec::new();

    let mut mw_new = MatrixWrapper::new();
    if mw.matrix_bit.len() < window {
        panic!("Window size is larger than the number of entries");
    }
    for x in window..mw.matrix_bit.len() - window {
        let mut bv2 = Vec::new();
        for y in 0..num_path {
            let mut bv: Vec<[bool; 2]> = Vec::new();

            for z in x - window..x + window {
                bv.push([mw.matrix_bit[z][y * 2], mw.matrix_bit[z][y * 2 + 1]]);
            }
            bv2.push((y, bv));
        }
        // sort by the bitvec
        bv2.sort_by(|a, b| a.1.cmp(&b.1));
        let f = get_index(&bv2, blocks, x)?;
        let flen = f.len();
        mw_new.matrix_bit.extend(f);
        index.push(flen);
    }
    mw_new.shape = (mw_new.matrix_bit.len(), mw_new.matrix_bit[0].len());
    Ok((mw_new, index))
}

///
pub fn get_index(
    vv: &Vec<(usize, Vec<[bool; 2]>)>,
    blocks: &mut Option<BufWriter<File>>,
    index: usize,
) -> Result<Vec<BitVec<u8>>, std::io::Error> {
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
    if let Some(b) = blocks {
        writeln!(b, "{}\t{:?}", index, pp)?;
    }
    Ok(getbv(&pp, vv.len()))
}

/// Create a bitvector
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

pub fn read_write_bim(
    _mw: &MatrixWrapper,
    index: &Vec<usize>,
    input: &str,
    output: &str,
    window: usize,
) -> Result<(), std::io::Error> {
    let input_file = std::fs::File::open(input).unwrap();
    let input_reader = std::io::BufReader::new(input_file);

    let file_out = File::create(output)?;
    let mut output_reader = std::io::BufWriter::new(file_out);

    let mut i = 0;

    for line in input_reader.lines().skip(window).enumerate() {
        let line = line.1?; // unwrap the line or propagate error
        let line_split = line.split_whitespace().collect::<Vec<&str>>();
        for x in 0..index[i] {
            writeln!(
                output_reader,
                "graph\t.\t0\tW{}:{}_{}\tA\tT",
                line_split[3], window, x
            )?;
        }
        i += 1;
        if i == index.len() {
            break;
        }
    }

    // Flush the buffer to ensure all data is written to file
    output_reader.flush()?;

    Ok(())
}
