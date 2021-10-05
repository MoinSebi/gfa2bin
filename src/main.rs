mod matrix;
mod counting;
mod matrix_wrapper;
mod helper;

use clap::{App, Arg};
use gfaR_wrapper::{NGfa, GraphWrapper};
use crate::matrix_wrapper::{matrix_node, matrix_edge, matrix_dir_node, MatrixWrapper};





fn main() {
    let matches = App::new("panSV")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        .arg(Arg::new("gfa")
            .short('g')
            .long("gfa")
            .about("Sets the input file to use")
            .takes_value(true)
            .default_value("/home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .about("Output prefix")
            .takes_value(true)
            .default_value("gfa2bin.default"))
        .arg(Arg::new("type")
            .short('t')
            .long("type")
            .about("Type of the matrix (when on graph)")
            .takes_value(true))
        .arg(Arg::new("bed")
            .long("bed")
            .about("Output bed + bim file"))
        .arg(Arg::new("bimbam")
            .long("bimbam")
            .about("Output bimbam format"))
        .arg(Arg::new("traversal")
            .long("traversal")
            .about("Additional traversaloutput file")
            .takes_value(true))
        .arg(Arg::new("delimter")
            .short('d')
            .long("delimter")
            .about("Delimter between genome and chromosome number")
            .takes_value(true))
        .arg(Arg::new("threshhold")
            .long("Copy number threshold")
            .about("Normalize to this number")
            .takes_value(true))
        .arg(Arg::new("decopy")
            .short('c')
            .long("decopy this")
            .takes_value(true))
        .arg(Arg::new("bubble input")
            .short('b')
            .long("bubble")
            .takes_value(true))

        .get_matches();


    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa


    // removed .default stuff
    let _input = matches.value_of("gfa").unwrap();
    let _output: &str = matches.value_of("output").unwrap();
    // Read the graph
    let mut graph = NGfa::new();
    graph.from_graph(_input);

    // Make graph, wrapper
    let mut gwrapper: GraphWrapper = GraphWrapper::new();
    gwrapper.fromNGfa(&graph, "_");


    println!("{} genomes and {} paths", gwrapper.genomes.len(), graph.paths.len());
    let type_out;

    if matches.is_present("bed") & matches.is_present("bimbam"){
        type_out = "all"
    }
    else if matches.is_present("bed"){
        type_out = "bed"

    } else if matches.is_present("bimbam") {
        type_out = "bimbam"
    } else {
        type_out = "bed"
    }

    // Make matrix
    //let h = test1(&gwrapper, &graph);
    let mat_node: MatrixWrapper<u32>;
    let mat_dir: MatrixWrapper<(u32, bool)>;
    let mat_edge: MatrixWrapper<(u32, bool, u32, bool)>;

    if matches.is_present("type"){
        let values: &str = matches.value_of("type").unwrap();
        if values.contains('n'){
            mat_node = matrix_node(&gwrapper, &graph);
            mat_node.write(type_out, _output, "node");
        }
        if values.contains('e'){
            mat_edge = matrix_edge(&gwrapper, &graph);
            mat_edge.write(type_out, _output, "edge");

        }
        if values.contains('d'){
            mat_dir = matrix_dir_node(&gwrapper, &graph);
            mat_dir.write(type_out,  _output, "dir");
        }
    } else {
        mat_node = matrix_node(&gwrapper, &graph);
        mat_node.write(type_out, _output, "node");
    }






}
