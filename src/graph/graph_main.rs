use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use std::fs::File;
use std::io::Read;

use crate::graph::parser::gfa_reader;

use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use log::{info, warn};
use std::path::Path;
use std::process;

pub fn graph_main(matches: &ArgMatches) {
    // Check graph file
    let _graph_file: &str;
    let graph_file: &str = if Path::new(matches.value_of("gfa").unwrap()).exists() {
        matches.value_of("gfa").unwrap()
    } else {
        warn!("No file with such name");
        process::exit(0x0100);
    };

    // Input parameters
    let feature1 = matches.value_of("feature").unwrap_or("node");

    let threshold = matches
        .value_of("threshold")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let sep = matches.value_of("pansn").unwrap_or("#");

    let mut bin = matches.is_present("bin");

    // Output
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let bimbam_output = matches.is_present("bimbam");
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.graph");

    // Modifications
    bin = bin || ((!bimbam_output) && threshold == 1);
    let feature = if ["node", "dirnode", "edge"].contains(&feature1) {
        feature1
    } else {
        "node"
    };
    let feature_enum = Feature::from_str(feature);
    let need_edges = feature != "node";

    info!("Input parameters");
    info!("Graph file: {}", graph_file);
    info!("Feature: {} -> {}", feature1, feature);
    info!("Threshold: {}", threshold);
    info!("Separator: {:?}", sep);
    info!("Binary: {}", bin);
    info!("Split: {}", split);
    info!(
        "Output format: {}",
        if bimbam_output { "bimbam" } else { "plink" }
    );
    info!("Output prefix: {}", output_prefix);

    info!("Reading the graph");
    // Read the graph and wrapper
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file(graph_file, need_edges);
    if matches.is_present("paths") {
        let mut file = File::open(matches.value_of("paths").unwrap()).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();
        graph.paths.retain(|x| !contents.contains(&x.name));
    }

    let wrapper: Pansn<NCPath> = Pansn::from_graph(&graph.paths, sep);

    let mut is_diploid = false;
    for x in wrapper.genomes.iter() {
        if x.haplotypes.len() == 2 {
            is_diploid = true;
        }
        if x.haplotypes.len() > 2 {
            warn!("More than 2 haplotypes");
            process::exit(0x0100);
        }
    }

    info!("Diploid: {}", is_diploid);

    info!("Number of samples: {}", wrapper.genomes.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the index");
    mw.make_index(&graph, feature_enum);

    info!("Create matrix");
    gfa_reader(&mut mw, &wrapper, bin, feature_enum);

    if !mw.matrix_bit.is_empty() {
        info!(
            "Shape is {:?} - {}",
            mw.matrix_bit.len(),
            mw.matrix_bit[0].len()
        );
        mw.remove_non_info();
        info!(
            "Shape (after remove) is {:?} - {}",
            mw.matrix_bit.len(),
            mw.matrix_bit[0].len()
        );
    }

    if bimbam_output {
        info!("Writing the bimbam");
        let chunk_size = (mw.matrix_u16.len() / split) + 1;
        let chunks = mw.matrix_u16.chunks(chunk_size);
        let len = chunks.len();

        for (index, _y) in chunks.enumerate() {
            mw.write_bimbam(index, output_prefix, &feature_enum, len, 1);
            mw.write_phenotype_bimbam(index, output_prefix, len);
        }
    } else {
        // Output
        info!("Writing the output");
        let chunk_size = (mw.matrix_bit.len() / split) + 1;
        let chunks = mw.matrix_bit.chunks(chunk_size);

        let len = chunks.len();
        for (index, _y) in chunks.enumerate() {
            //write_bed2(y, output_prefix, feature, index, len);
            mw.write_fam(index, output_prefix, feature_enum, len);
            mw.write_bed(index, output_prefix, feature_enum, len);
            mw.write_bim(index, output_prefix, &feature_enum, len);
        }
    }
}
