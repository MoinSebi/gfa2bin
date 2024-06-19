use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use std::fs::File;
use std::io::Read;

use crate::graph::parser::gfa_reader;

use clap::ArgMatches;
use gfa_reader::{Gfa, Pansn};
use log::{info, warn};
use std::path::Path;
use std::process;
use std::thread::sleep;
use packing_lib::core::core::PackCompact;
use packing_lib::normalize::convert_helper::Method;

pub fn graph_main(matches: &ArgMatches) {


    // Check graph file
    let graph_file: &str = if Path::new(matches.value_of("gfa").unwrap()).exists() {
        matches.value_of("gfa").unwrap()
    } else {
        warn!("No file with such name");
        process::exit(0x0100);
    };

    // Input parameters
    let feature1 = matches.value_of("feature").unwrap_or("node");
    let feature = if ["node", "dirnode", "edge"].contains(&feature1) {
        feature1
    } else {
        "node"
    };
    let feature_enum = Feature::from_str(feature);
    let sep = matches.value_of("pansn").unwrap_or(" ");
    let need_edges = feature != "node";                             // If you something else than a node, you need to parse the edges


    // Output
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let bimbam_output = matches.is_present("bimbam");
    let output_prefix = matches.value_of("output").unwrap();



    // Threshold
    let mut absolute_thresh = matches
        .value_of("absolute-threshold")
        .unwrap_or("0")
        .parse::<u32>()
        .unwrap();
    let mut fraction = matches
        .value_of("fraction")
        .unwrap_or("1.0")
        .parse::<f32>()
        .unwrap();
    let mut method = Method::from_str(matches.value_of("method").unwrap_or("nothing"));
    let mut std = matches.value_of("std").unwrap_or("0.0").parse::<f32>().unwrap();
    let include_all = matches.is_present("non-covered");

    // Bin is for faster computation
    let mut bin = false;
    if matches.is_present("absolute-threshold") {
        if absolute_thresh == 1 && !bimbam_output {
            bin = true;
        }
    }

    if !matches.is_present("absolute-threshold")
        && !matches.is_present("method")
    {
        if !bimbam_output {
            bin = true;
            absolute_thresh = 1;
        }
    }





    info!("Input parameters");
    info!("Graph file: {}", graph_file);
    info!("Feature: {} -> {}", feature1, feature);
    info!("Absolute threshold: {}", absolute_thresh);
    info!("Method: {:?}", method.to_string());
    info!("Relative threshold: {:?}", fraction);
    info!("Separator: {:?}", sep);
    info!("Binary: {}", bin);
    info!("Exclude non-covered: {}", !include_all);
    info!("Split: {}", split);
    info!(
        "Output format: {}",
        if bimbam_output { "bimbam" } else { "plink" }
    );
    info!("Output prefix: {}", output_prefix);



    info!("Reading the graph");
    // Read the graph and wrapper
    let mut graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file(graph_file);
    graph.walk_to_path(sep);
    if matches.is_present("paths") {
        let mut file = File::open(matches.value_of("paths").unwrap()).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        // Remove paths
        graph.paths.retain(|x| !contents.contains(&x.name));
    }

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
            warn!("Haplotypes are {}", x.haplotypes.iter().map(|x| x.name.clone()).collect::<Vec<String>>().join(", "));
            warn!("Will only take the first 2 haplotypes")
        }
    }

    info!("Diploid: {}", is_diploid);
    info!("Number of samples: {}", wrapper.genomes.len());
    info!("Number of paths: {}", graph.paths.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the index");
    sleep(std::time::Duration::from_secs(10));
    mw.make_index(&graph, feature_enum);
    sleep(std::time::Duration::from_secs(10));
    info!("Index size: {}", mw.geno_names.len());
    info!("Index size: {}", mw.geno_names.capacity());

    info!("Create matrix");
    gfa_reader(&mut mw, &wrapper, bin, feature_enum);
    info!("Matrix size: {}, {}", mw.matrix_bit.len(), mw.matrix_bit[0].len());
    let mut thresh = Vec::new();
    for x in mw.matrix_u16.iter(){
        let mut b = x.clone();
        thresh.push(PackCompact::threshold(&mut b, include_all, fraction, std, method));
    }

    if !bin {
        if !bimbam_output {
            mw.make_bin_row(&thresh);
            info!("Number of entries: {}", mw.matrix_bit.len());
        } else {
            info!("Number of entries: {}", mw.matrix_u16.len());
        }
    }

    if !mw.matrix_bit.is_empty() {
        mw.remove_non_info();
        info!("Number of entries (after remove): {}", mw.matrix_bit.len());
    }

    mw.g(bimbam_output, split, output_prefix, thresh, feature_enum);
}
