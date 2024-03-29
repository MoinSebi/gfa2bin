use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;

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
    let feature1 = matches.value_of("Feature").unwrap_or("node");
    let threshold = matches
        .value_of("threshold")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let sep = matches.value_of("sep").unwrap_or("#");

    let mut bin = matches.is_present("bin");
    let diploid = matches.is_present("diploid");

    // Output
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let output_format = matches.value_of("format").unwrap();
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.graph");

    // Modifications
    bin = bin || (output_format == "plink" && threshold == 1);
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
    info!("Diploid: {}", diploid);
    info!("Split: {}", split);
    info!("Output format: {}", output_format);
    info!("Output prefix: {}", output_prefix);

    info!("Reading the graph");
    // Read the graph and wrapper
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file(graph_file, need_edges);
    let wrapper: Pansn<NCPath> = Pansn::from_graph(&graph.paths, sep);

    info!("Number of samples: {}", wrapper.genomes.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the index");
    mw.make_index(&graph, feature_enum);

    info!("Create matrix");
    gfa_reader(&mut mw, &wrapper, bin, feature_enum);

    info!(
        "Shape is {:?} - {}",
        mw.matrix_bin.len(),
        mw.matrix_bin[0].len()
    );
    // Filter the matrix
    mw.remove_non_info();
    info!(
        "Shape is {:?} - {}",
        mw.matrix_bin.len(),
        mw.matrix_bin[0].len()
    );
    // Output
    info!("Writing the output");
    let chunk_size = (mw.matrix_bin.len() / split) + 1;
    let chunks = mw.matrix_bin.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        //write_bed2(y, output_prefix, feature, index, len);
        mw.write_fam(index, output_prefix, feature_enum, len);
        mw.write_bed(index, output_prefix, feature_enum, len);
        mw.write_bim(index, output_prefix, &feature_enum, len);
    }
}
