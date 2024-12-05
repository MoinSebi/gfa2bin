use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use std::fs::File;
use std::io::Read;

use crate::graph::parser::{diploid_adder, gfa_reader};

use clap::ArgMatches;
use gfa_reader::{Gfa, Pansn};
use log::{info, warn};
use packing_lib::core::core::PackCompact;
use packing_lib::normalize::convert_helper::Method;
use std::path::Path;
use std::process;

/// # Main function for 'gfa2bin graph'
pub fn graph_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin graph'");

    // Check graph file
    let graph_file: &str = matches.value_of("gfa").unwrap();

    // Input parameters
    let feature1 = matches.value_of("feature").unwrap_or("node");
    let output_feature = if ["node", "dirnode", "edge"].contains(&feature1) {
        feature1
    } else {
        warn!("Feature {} is not supported", feature1);
        warn!("Only node, dirnode and edge are supported");
        process::exit(1);
    };

    let feature_enum = Feature::from_str(output_feature);
    let mut sep = matches.value_of("PanSN").unwrap();

    let bimbam_output = matches.is_present("bimbam");
    let output_prefix = matches.value_of("output").unwrap();

    // Threshold
    let mut absolute_thresh = matches
        .value_of("absolute-threshold")
        .unwrap()
        .parse::<u32>()
        .expect("Error: Absolute threshold is not a number");

    let method = Method::from_str(
        matches
            .value_of("method").unwrap_or("nothing")
    );
    let keep_zeros = matches.is_present("keep-zeros");

    // Bin is for faster computation
    let mut bin = false;
    if absolute_thresh == 1 {
        bin = true;
    }

    let max_scale = matches.is_present("max-scale");

    // Dummy phenotype
    let mut pheno = f64::MAX;
    if matches.is_present("pheno") {
        pheno = matches.value_of("pheno").unwrap().parse()?;
    }

    let threads = matches
        .value_of("threads")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    info!("Input parameters");
    info!("Graph file: {}", graph_file);
    info!("Feature: {} -> {}", feature1, output_feature);
    info!(
        "Pan-SN: {}",
        if sep == "\n" {
            "None".to_string()
        } else {
            format!("{:?}", sep)
        }
    );
    info!("Absolute threshold: {}", absolute_thresh);
    info!("Method: {}", method.to_string());
    info!(
        "Fraction: {}",
        if matches.is_present("fraction") {
            matches.value_of("fraction").unwrap()
        } else {
            "None"
        }
    );
    info!("Binary: {}", bin);
    info!("Keep zeros: {}", keep_zeros);
    info!("Max value scaling (only bimbam): {}", max_scale);
    info!("Threads: {}", threads);
    info!(
        "Dummy-Pheno: {}",
        if pheno == f64::MAX {
            "NA".to_string()
        } else {
            pheno.to_string()
        }
    );
    info!(
        "Output format: {}",
        if bimbam_output { "bimbam" } else { "PLINK" }
    );
    info!("Output prefix: {}\n", output_prefix);

    let mut fraction = 0.0;
    if matches.is_present("fraction") {
        fraction = matches
            .value_of("fraction")
            .unwrap()
            .parse::<f32>()
            .expect("Error: Fraction is not a number");
    }
    let mut dynamic = false;
    if fraction > 1.0 || fraction < 0.0 {
        panic!("Fraction is not between 0 and 1");
    }
    if method != Method::Nothing && fraction == 0.0 {
        panic!("Method is given but fraction is not.");
    } else if method == Method::Nothing && fraction != 0.0 {
        panic!("Fraction is given but method is not.");
    } else if method != Method::Nothing && fraction != 0.0 {
        dynamic = true;
    }
    info!("Dynamic threshold: {}", dynamic);

    info!("Read the graph");
    // Read the graph and wrapper
    let mut graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file_multi(graph_file, threads);
    if graph.paths.is_empty() && sep == "\n" {
        sep = "#"
    }

    graph.walk_to_path(sep);

    // Wrapper on PanSN
    let wrapper: Pansn<u32, (), ()> = Pansn::from_graph(&graph.paths, sep);

    // Check diploid
    let mut is_diploid = false;
    for x in wrapper.genomes.iter() {
        if x.haplotypes.len() == 2 {
            is_diploid = true;
        }
        if x.haplotypes.len() > 2 {
            warn!("More than 2 haplotypes");
            warn!(
                "Haplotypes are {}",
                x.haplotypes
                    .iter()
                    .map(|x| x.name.clone())
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            warn!("Will only take the first 2 haplotypes")
        }
    }

    info!("Diploid: {}", is_diploid);
    info!("Number of samples: {}", wrapper.genomes.len());
    info!("Number of paths: {}", graph.paths.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the index");
    mw.create_index(&graph, feature_enum);

    info!("Read the graph into matrix");
    gfa_reader(&mut mw, &wrapper, bin, feature_enum);

    // Threshold calculation
    let mut thresh = Vec::new();

    // If max_scale is true, threshold needs to be adjusted
    if max_scale && bimbam_output {
        for x in mw.matrix_u16.iter() {
            thresh.push(*x.iter().max().ok_or("Error: Empty vector")? as f32)
        }
    } else {
        if !dynamic {
            thresh = vec![absolute_thresh as f32; mw.geno_names.len()];
        } else {
            for x in mw.matrix_u16.iter() {
                let mut count_vec = x.clone();
                diploid_adder(&mw.sample_index_u16, &mut count_vec);

                thresh.push(PackCompact::threshold(
                    &mut count_vec,
                    keep_zeros,
                    fraction,
                    0.0,
                    method,
                ));
            }
        }
    }

    mw.write_wrapper(
        bimbam_output,
        1,
        output_prefix,
        thresh,
        feature_enum,
        pheno,
        !keep_zeros,
    );
    Ok(())
}
