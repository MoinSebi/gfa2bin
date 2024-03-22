use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use std::fs::File;
use std::io::Read;

use crate::graph::parser::gfa_reader;

use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use log::{info, warn};
use packing_lib::convert::convert_helper::Method;
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
    let sep = matches.value_of("pansn").unwrap_or(" ");

    // Output
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();
    let bimbam_output = matches.is_present("bimbam");
    let output_prefix = matches.value_of("output").unwrap();

    // Modifications
    let feature = if ["node", "dirnode", "edge"].contains(&feature1) {
        feature1
    } else {
        "node"
    };
    let feature_enum = Feature::from_str(feature);
    let need_edges = feature != "node";
    let mut absolute_thresh = matches
        .value_of("absolute-threshold")
        .unwrap_or("0")
        .parse::<u32>()
        .unwrap();
    let mut relative_thresh = matches
        .value_of("relative-threshold")
        .unwrap_or("0")
        .parse::<u32>()
        .unwrap();
    let mut method = Method::from_str(matches.value_of("method").unwrap_or("nothing"));

    let mut bin = false;
    if matches.is_present("absolute-threshold") {
        if absolute_thresh == 1 && !bimbam_output {
            bin = true;
        }
    }
    if !matches.is_present("absolute-threshold")
        && !matches.is_present("relative-threshold")
        && !matches.is_present("method")
    {
        if !bimbam_output {
            bin = true;
            absolute_thresh = 1;
        }
        if bimbam_output {
            method = Method::Percentile;
            relative_thresh = 100;
        }
    }

    if matches.is_present("method") && !matches.is_present("relative-threshold") {
        relative_thresh = 100;
    }

    info!("Input parameters");
    info!("Graph file: {}", graph_file);
    info!("Feature: {} -> {}", feature1, feature);
    info!("Absolute threshold: {}", absolute_thresh);
    info!("Method: {:?}", method.to_string());
    info!("Relative threshold: {:?}", relative_thresh);
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

    let thresh = mw.make_thresh(absolute_thresh, relative_thresh as u16, method);

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

    if bimbam_output {
        info!("Writing the bimbam");
        let chunk_size = (mw.matrix_u16.len() / split) + 1;
        let chunks = mw.matrix_u16.chunks(chunk_size);
        let len = chunks.len();

        for (index, _y) in chunks.enumerate() {
            mw.write_bimbam(index, output_prefix, &feature_enum, len, &thresh);
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
