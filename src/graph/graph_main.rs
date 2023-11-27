use std::path::Path;
use std::process;
use bitvec::macros::internal::funty::Fundamental;
use clap::ArgMatches;
use gfa_reader::{GraphWrapper, NCGfa, NCPath};
use log::{info, warn};
use crate::core::core::MatrixWrapper;
use crate::core::writer::{write_bed2, write_bim_dirnode, write_bim_edges, write_bim_nodes};
use crate::graph::parser::{dir_nodes2, Haplotype, gfa_nodes_reader, egdes1};
use crate::helper::make_dir_name;

pub fn graph_main(matches: &ArgMatches){

    // Input
    let graph_file: &str;
    if Path::new(matches.value_of("gfa").unwrap()).exists() {
        graph_file = matches.value_of("gfa").unwrap();
    } else {
        warn!("No file with such name");
        process::exit(0x0100);
    }

    // Input parameters
    let mut feature1 = matches.value_of("feature").unwrap_or("node");
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
    // Start the tool - read the graph and wrapper
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file_direct(graph_file, true);
    let mut wrapper: GraphWrapper<NCPath> = GraphWrapper::new();
    wrapper.from_gfa(&graph.paths, &sep);

    // Get the haplotypes
    let mut haplo = Haplotype::from_wrapper(&wrapper, &' '.to_string());
    if diploid{
        haplo = Haplotype::from_wrapper(&wrapper, sep);
    }

    info!("Number of samples: {}", haplo.haplotype.len());

    // This is the matrix
    let mut mw = MatrixWrapper::new();

    info!("Create the matrix"); 
    // These are the names
    let mut node_name: Vec<usize> = Vec::with_capacity(graph.nodes.len());
    let mut dir_node_name: Vec<(usize, bool)> = Vec::new();
    let mut edges_names: Vec<(u32, bool, u32, bool)> = Vec::new();
    if feature == "node"{
        node_name = (1..graph.nodes.len()+1).collect();
        //dir_nodes2(&graph, &haplo, &wrapper, &mut mw, false);

        gfa_nodes_reader(&graph, &haplo, &wrapper, &mut mw, bin)
    } else if feature == "dirnode" {
        dir_node_name = make_dir_name(&graph.nodes.len());

        dir_nodes2(&graph, &haplo, &wrapper, &mut mw, bin);
    } else {
        //edges_names = (0..graph.edges.as_ref().unwrap().len()).collect();
        egdes1( &graph, &wrapper, &haplo, &mut mw, bin, &mut edges_names);
    }
    mw.transposed = true;



    if output_format == "plink" && mw.matrix_bin.len() == 0{
        mw.make_binary2(threshold as u32, &haplo);
    }
    info!("Shape is {:?} - {}", mw.matrix_bin.len(), mw.matrix_bin[0].len());

    info!("Filtering matrix");
    if feature == "node"{
        mw.filter_shared1(&mut node_name);
    } else if feature == "dirnode" {
        mw.filter_shared1(&mut dir_node_name);
    } else {
        mw.filter_shared1(&mut edges_names);
    }
    info!("New shape is {:?} - {:?}", mw.matrix_bin.len(), mw.matrix_bin[0].len());



    info!("Write files");
    // Output
    if output_format == "plink" {
        let chunk_size = (mw.matrix_bin.len() / split) + 1;
        let chunks = mw.matrix_bin.chunks(chunk_size);
        let mut chunk1 = node_name.chunks(chunk_size);
        let mut chunk2 = dir_node_name.chunks(chunk_size);
        let mut chunk3 = edges_names.chunks(chunk_size);

        let len = chunks.len();
        for (index, y) in chunks.enumerate(){
            write_bed2(y, output_prefix, feature, index, len);
            mw.write_fam(index, output_prefix, feature, len);
            if ! node_name.is_empty(){
                write_bim_nodes(chunk1.next().unwrap(), output_prefix, feature, index, len);
            } else if ! dir_node_name.is_empty(){
                write_bim_dirnode(chunk2.next().unwrap(), output_prefix, feature, &index, len);

            } else {
                write_bim_edges(chunk3.next().unwrap(), output_prefix, feature, &index, len);
            }
        }

    }





}