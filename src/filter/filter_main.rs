use crate::core::core::MatrixWrapper;

use clap::ArgMatches;

use crate::core::bfile::count_lines;
use crate::core::helper::Feature;
use crate::remove::remove_main::read_write_filter_index;
use log::info;

/// # Filter main function
///
/// Filter plink by MAF/maf and path MAF/maf
///
pub fn filter_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Read the arguments from the command line
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap();

    let mut maf = matches
        .value_of("maf")
        .unwrap()
        .parse::<f64>()
        .expect("Error parsing maf");

    let mut MAF = matches
        .value_of("MAF")
        .unwrap()
        .parse::<f64>()
        .expect("Error parsing MAF");
    let mut missing_rate = matches
        .value_of("missing-rate")
        .unwrap()
        .parse::<f64>()
        .expect("Error parsing missing_rate");

    info!("Input file: {}", plink_file);
    info!("maf: {}", maf);
    info!("MAF: {}", MAF);
    info!(
        "mac: {}",
        if matches.is_present("mac") {
            matches
                .value_of("mac")
                .unwrap()
                .parse::<usize>()
                .expect("mac is not a number")
                .to_string()
        } else {
            "None".to_string()
        }
    );
    info!(
        "MAC: {}",
        if matches.is_present("MAC") {
            matches
                .value_of("MAC")
                .unwrap()
                .parse::<usize>()
                .expect("MAC is not a number")
                .to_string()
        } else {
            "None".to_string()
        }
    );
    info!("Missing-rate: {}", missing_rate);
    info!(
        "Missing-count: {}",
        if matches.is_present("missing_count") {
            matches
                .value_of("missing_count")
                .unwrap()
                .parse::<usize>()
                .expect("missing_count is not a number")
                .to_string()
        } else {
            "None".to_string()
        }
    );
    info!("Output prefix: {}", output_prefix);

    let mut mw = MatrixWrapper::new();
    let bim_count = count_lines(&format!("{}{}", plink_file, ".bim"))?;
    let fam_count = count_lines(&format!("{}{}", plink_file, ".fam"))?;
    mw.read_bed(&format!("{}{}", plink_file, ".bed"), fam_count, bim_count)?;

    if matches.is_present("mac") {
        maf = matches
            .value_of("mac")
            .unwrap()
            .parse::<f64>()
            .expect("Error parsing mac")
            / fam_count as f64;
    }
    if matches.is_present("MAC") {
        MAF = matches
            .value_of("MAC")
            .unwrap()
            .parse::<f64>()
            .expect("Error parsing MAC")
            / fam_count as f64;
    }

    if matches.is_present("missing-count") {
        missing_rate = matches
            .value_of("missing_count")
            .unwrap()
            .parse::<f64>()
            .expect("Error parsing missing_count")
            / bim_count as f64;
    }

    info!(
        "Matrix size (SNPs X Samples): {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit[0].len()
    );

    let mut remove_index_samples = Vec::new();
    if missing_rate != 0.0 {
        info!(
            "Filtering by missing rate or missing count: {}",
            missing_rate
        );
        remove_index_samples = mw.filter_path(missing_rate, 1.0);
    }

    info!(
        "Matrix size (SNPs X Samples) - After sample filtering: {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit
            .first()
            .ok_or("Matrix is now empty after path removal")?
            .len()
    );

    let mut remove_index_genotypes = Vec::new();
    if maf != 0.0 || MAF != 1.0 {
        info!("Filtering by MAF/maf: {} {}", maf, MAF);
        remove_index_genotypes = mw.filet_maf(maf, MAF);
    }

    info!(
        "Matrix size (SNPs X Samples) - after maf: {} {}",
        mw.matrix_bit.len(),
        mw.matrix_bit
            .first()
            .ok_or("Matrix is now empty after MAF/maf removal ")?
            .len()
    );

    info!("Writing BED file");
    mw.write_bed(0, output_prefix, Feature::Node, 1);

    info!("Writing BIM file");
    read_write_filter_index(
        &format!("{}{}", plink_file, ".bim"),
        &format!("{}{}", output_prefix, ".bim"),
        &remove_index_genotypes,
    )
    .expect("Error writing BIM file");

    info!("Writing FAM file");
    read_write_filter_index(
        &format!("{}{}", plink_file, ".fam"),
        &format!("{}{}", output_prefix, ".fam"),
        &remove_index_samples,
    )
    .expect("Error writing FAM file");

    Ok(())
}

impl MatrixWrapper {
    /// # Filter bit matrix by MAF/maf
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

    /// # Filter matrix by path MAF/maf
    pub fn filter_path(&mut self, lower_limit: f64, upper_limit: f64) -> Vec<usize> {
        let mut remove_index_vec = Vec::new();

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
