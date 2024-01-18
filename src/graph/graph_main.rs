use std::path::Path;
use std::process;
use bitvec::macros::internal::funty::Fundamental;
use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use log::{info, warn};
use crate::core::core::{MatrixWrapper};
use crate::core::helper::Feature;
use crate::core::writer::{write_bed2, write_bim_dirnode, write_bim_edges, write_bim_nodes};
use crate::graph::parser::{gfa_nodes_reader2};
use crate::helper::make_dir_name;

pub fn graph_main(matches: &ArgMatches){

    // Check graph file
    let graph_file: &str;
    if Path::new(matches.value_of("gfa").unwrap()).exists() {
        graph_file = matches.value_of("gfa").unwrap();
    } else {
        warn!("No file with such name");
        process::exit(0x0100);
    }

    // Input parameters
    let mut feature1 = matches.value_of("Feature").unwrap_or("node");
    let threshold = matches.value_of("threshold").unwrap_or("1").parse::<usize>().unwrap();
    let sep = matches.value_of("sep").unwrap_or("#");

    let mut bin = matches.is_present("bin");
    let diploid = matches.is_present("diploid");



    // Output
    let split = matches.value_of("split").unwrap_or("1").parse::<usize>().unwrap();
    let output_format = matches.value_of("format").unwrap();
    let output_prefix = matches.value_of("output").unwrap_or("gfa2bin.graph");

    // Modifications
    bin = bin || (output_format == "plink" && threshold == 1);
    let feature = if ["node", "dirnode", "edge"].contains(&feature1) { feature1 } else { "node" };
    let feature_enum = Feature::from_str(feature);

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
    graph.parse_gfa_file_direct(graph_file, true);
    let mut wrapper: Pansn<NCPath> = Pansn::from_graph(&graph.paths, sep);


    info!("Number of samples: {}", wrapper.genomes.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the matrix");
    mw.make_index(&graph, feature_enum);

    if feature == "node"{
        //dir_nodes2(&graph, &haplo, &wrapper, &mut mw, false);

        gfa_nodes_reader2(&mut mw, &wrapper, &graph, bin, &feature_enum)
    }


    info!("Shape is {:?} - {}", mw.matrix_bin.len(), mw.matrix_bin[0].len());

    info!("Filtering matrix");
    //if feature == "node"{
    //    mw.filter_shared2();
    //}


    info!("Write files");
    // Output
    if output_format == "plink" {
        let chunk_size = (mw.matrix_bin.len() / split) + 1;
        let chunks = mw.matrix_bin.chunks(chunk_size);

        let len = chunks.len();
        for (index, y) in chunks.enumerate(){
            //write_bed2(y, output_prefix, feature, index, len);
            mw.write_fam(index, output_prefix, feature, len);
            mw.write_bed(index, output_prefix, feature, len);
            mw.write_bim(index, output_prefix, &feature_enum, len);
        }

    }





}