mod matrix;
mod counting;
use gfaR::{Gfa};
use crate::matrix::{matrix_wrapper, test1, test22};
use std::collections::HashMap;
use clap::{App, Arg};
use gfaR_wrapper::{NGfa, GraphWrapper};
use std::mem::size_of_val;
use std::{thread, time::Duration};


fn all_nodes(graph: &Gfa) -> Vec<String>{
    let mut g1: Vec<String> = graph.nodes.keys().cloned().collect();
    g1.sort();
    g1
}



fn main() {
    let matches = App::new("panSV")
        .version("0.1.0")
        .author("Sebastian V")
        .about("Bubble detection")
        .arg(Arg::new("gfa2")
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
            .default_value("panSV.output"))
        .arg(Arg::new("traversal")
            .short('t')
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

    // removed .default stuff
    let _input = matches.value_of("gfa2").unwrap();
    let _output: &str = matches.value_of("output").unwrap();


    println!("We read a graph");
    // Change this for read in
    let g1 = "/home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa";
    //let g1 = "/home/svorbrugg_local/Rust/data/sixRef.noC.gfa";

    //let g1 = "/home/svorbrugg_local/Rust/data/chr1.wfmash.n20.a90.s10000.p1,19,39,3,81,1.seqwish.sort.smooth.sort.noC.gfa";

    // Read the graph
    let mut graph = NGfa::new();
    graph.from_graph(g1);

    // Make graph, wrapper
    let mut gwrapper: GraphWrapper = GraphWrapper::new();
    gwrapper.fromNGfa(&graph, "_");


    println!("{}", gwrapper.genomes.len());

    // Make matrix
    println!("matrix1");
    //let h = test1(&gwrapper, &graph);
    let h = test22(&gwrapper, &graph);

    // Write matrix
    h.matrix.write_bed(_output);
    h.bim_helper();
    h.matrix.write_bimbam();





}
