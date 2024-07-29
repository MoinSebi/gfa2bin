use crate::core::core::MatrixWrapper;

use clap::ArgMatches;

use crate::core::bfile::count_lines;
use crate::core::helper::Feature;
use crate::remove::remove_main::process_file2;
use log::info;

/// Block main function
///
/// Easy block function
/// Extract the subpath from a graph for each node
pub fn filter_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Read the arguments from the command line
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

    info!("Input file: {}", plink_file);
    info!("maf: {}", maf);
    info!("MAF: {}", MAF);
    info!("Split size: {}", split);
    info!("Output prefix: {}", output_prefix);

    let mut mw = MatrixWrapper::new();
    let bim_count = count_lines(&format!("{}{}", plink_file, ".bim"))?;
    let fam_count = count_lines(&format!("{}{}", plink_file, ".fam"))?;
    mw.read_bed(&format!("{}{}", plink_file, ".bed"), fam_count, bim_count)?;

    info!(
        "Matrix size (SNPs X Samples): {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit[0].len()
    );

    info!(
        "Matrix size (SNPs X Samples): {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit[0].len()
    );

    let mut a1 = 0.0;
    let mut a2 = 1.0;
    if matches.is_present("entry-max") {
        a2 = matches
            .value_of("entry-max")
            .unwrap()
            .parse::<f64>()
            .unwrap();
    }
    if matches.is_present("entry-min") {
        a1 = matches
            .value_of("entry-min")
            .unwrap()
            .parse::<f64>()
            .unwrap();
    }
    let o2 = mw.filter_path(a1, a2);
    info!(
        "Matrix size (SNPs X Samples) - after path: {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit.first()
            .ok_or("Matrix is now empty after path removal")?
            .len()
    );
    let o = mw.filet_maf(maf, MAF);
    info!(
        "Matrix size (SNPs X Samples) - after maf: {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit.first()
            .ok_or("Matrix is now empty after MAF/maf removal ")?
            .len()
    );

    process_file2(
        &format!("{}{}", plink_file, ".bim"),
        &format!("{}{}", output_prefix, ".bim"),
        &o,
    )?;
    process_file2(
        &format!("{}{}", plink_file, ".fam"),
        &format!("{}{}", output_prefix, ".fam"),
        &o2,
    )?;

    mw.write_bed(0, output_prefix, Feature::Node, 1);

    Ok(())
}

impl MatrixWrapper {
    /// Filter bit matrix by MAF/maf
    pub fn filet_maf(&mut self, maf: f64, MAF: f64) -> Vec<usize> {
        let mut remove_index_vec = Vec::new();
        for (index, bitvec) in self.matrix_bit.iter().enumerate() {
            let mut count = 0;
            for y in bitvec.iter() {
                if y == true {
                    count += 1;
                }
            }
            let maf1 = count as f64 / bitvec.len() as f64;
            if maf1 < maf || maf1 > MAF {
                remove_index_vec.push(index);
            }
        }
        self.remove_by_index(&remove_index_vec);
        remove_index_vec
    }

    /// Filter matrix by path MAF/maf
    pub fn filter_path(&mut self, lower_limit: f64, upper_limit: f64) -> Vec<usize> {
        let mut remove_index_vec = Vec::new();
        println!("{:?}", self.matrix_bit[0].len() / 2);
        for i in 0..self.matrix_bit[0].len() / 2 {
            let mut c = 0;

            for x in self.matrix_bit.iter() {
                if x[i * 2] || x[i * 2 + 1] {
                    c += 1;
                }
            }
            let a = c as f64 / self.matrix_bit.len() as f64;
            if a < lower_limit || a >= upper_limit {
                remove_index_vec.push(i);
            }
        }

        self.remove_index_samples(&remove_index_vec);
        remove_index_vec
    }
}
