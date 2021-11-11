mod matrix;
mod counting;
mod matrix_wrapper;
mod helper;
mod writer;

use clap::{App, Arg};
use gfaR_wrapper::{NGfa, GraphWrapper};
use crate::matrix_wrapper::{matrix_node, matrix_edge, matrix_dir_node, MatrixWrapper, matrix_node_coverage, matrix_node_coverage2, matrix_node10, MatrixWrapper2, write_matrix, matrix_edge2, matrix_dir_node2, matrix_pack_bit, matrix_pack_u16, remove_bimap};
use std::process;
use crate::matrix::Matrix;
use crate::helper::{write_genomes, get_thresh};
use bimap::BiMap;
use crate::writer::{write_reduce, write_bimap};


fn main() {
    let matches = App::new("panSV")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        // Input
        .arg(Arg::new("gfa")
            .short('g')
            .long("gfa")
            .about("Sets the input file to use")
            .takes_value(true))
        .arg(Arg::new("pack")
            .short('p')
            .long("pack")
            .about("if the input is coverage")
            .takes_value(true))
        .arg(Arg::new("delimter")
            .short('d')
            .long("delimter")
            .about("Delimter between genome and chromosome number")
            .takes_value(true))
        .arg(Arg::new("removeNames")
            .about("Only this"))

        .arg(Arg::new("type")
            .short('t')
            .long("type")
            .about("Type of the matrix (when on graph)")
            .takes_value(true))


        // Computational
        .arg(Arg::new("threshhold")
            .long("Copy number threshold")
            .about("Normalize to this number")
            .takes_value(true))

        // Output
        .arg(Arg::new("genomes")
            .about("Output just the genomes"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .about("Output prefix")
            .takes_value(true)
            .default_value("gfa2bin.default"))
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
        .arg(Arg::new("smart")
            .long("smart")
            .about("Reduce the number of tests"))


        .get_matches();


    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa


    // Check if input
    if !(matches.is_present("gfa") | matches.is_present("pack")){
        eprintln!("No input");
        process::exit(0x0100);
    }


    // This is to decide which output
    let type_out;
    let _output: &str = matches.value_of("output").unwrap();

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

    // Check if del is set or not
    let del: &str;
    if matches.is_present("delimter"){
        del = matches.value_of("delimter").unwrap();
    } else {
        del = " ";
    }


    let mut matrix = MatrixWrapper2::new();
    let mut index_normal: BiMap<u32, usize> = BiMap::new();
    let mut index_dir: BiMap<(u32, bool), usize> = BiMap::new();
    let mut index_edge: BiMap<(u32, bool, u32, bool), usize> = BiMap::new();
    if matches.is_present("gfa"){
        let _input = matches.value_of("gfa").unwrap();
        let _output: &str = matches.value_of("output").unwrap();
        // Read the graph
        let mut graph = NGfa::new();
        graph.from_graph(_input);

        // Make graph, wrapper
        let mut gwrapper: GraphWrapper = GraphWrapper::new();
        gwrapper.fromNGfa(&graph, del);
        if matches.is_present("type"){
            let values: &str = matches.value_of("type").unwrap();
            if values.contains('n'){
                matrix_node10(&gwrapper, &graph, & mut matrix, & mut index_normal);
            }
            if values.contains('e'){
                matrix_edge2(&gwrapper, &graph, & mut matrix, & mut index_edge);

            }
            if values.contains('d'){
                matrix_dir_node2(&gwrapper, &graph, & mut matrix, & mut index_dir);
            }
        } else {
            matrix_node10(&gwrapper, &graph, & mut matrix, & mut index_normal);
        }
    } else {
        if matches.is_present("pack") {
            let file_pack = matches.value_of("pack").unwrap();
            let j = get_thresh(&file_pack);
            if j == 0{
                matrix_pack_u16(file_pack, & mut matrix, & mut index_normal);
            } else {
                matrix_pack_bit(file_pack, & mut matrix, & mut index_normal);
            }
        }
    }

    if matches.is_present("names"){
        matrix.write_names();
    }

    if matches.is_present("removeNames"){
        matrix.remove_genomes(matches.value_of("removeNames").unwrap())
    }


    matrix.make_binary(1);
    let mut remove_this: Vec<u32> = Vec::new();
    if matches.is_present("filter"){
        remove_this = matrix.filter();
    }

    if matches.is_present("reduce"){
        let k = matrix.reduce_comb1();
        write_reduce(&k.0, &k.1);
    }




    //Output
    if matches.is_present("type"){
        let values: &str = matches.value_of("type").unwrap();
        if values.contains('n'){
            write_matrix(& mut matrix, type_out,  _output, "node");
        }
        if values.contains('e'){
            write_matrix(& mut matrix, type_out,  _output, "edge");

        }
        if values.contains('d'){
            write_matrix(& mut matrix, type_out,  _output, "dir");
        }
    } else {
        write_matrix(& mut matrix, type_out,  _output, "gfa2bin");
    }


    // THEN FILTER ROWS (BIMAP)
    if !index_normal.is_empty(){
        remove_bimap(& mut index_normal, remove_this);
        write_bimap(& index_normal);
    } else if !index_dir.is_empty() {
        remove_bimap(& mut index_dir, remove_this);
        write_bimap(& index_normal);
    } else if !index_edge.is_empty() {
        remove_bimap(& mut index_edge, remove_this);
        write_bimap(& index_normal);
    }








    // Comment: Because different wrapper are used, we need to check everything by itself
    let mut i: & mut Matrix;
    // Read from a pack file
    if matches.is_present("pack") {
        let ii = matches.value_of("pack").unwrap();
        let j =0;
        let mut gw: MatrixWrapper<u32>;
        if j == 0{
            gw = matrix_node_coverage2(ii);
        } else {
            gw = matrix_node_coverage(ii);
            gw.matrix.filter();

        }
        gw.matrix.filter();

        gw.write(type_out, _output , "node");

    } else if matches.is_present("gfa") {
        let _input = matches.value_of("gfa").unwrap();
        let _output: &str = matches.value_of("output").unwrap();
        // Read the graph
        let mut graph = NGfa::new();
        graph.from_graph(_input);

        // Make graph, wrapper
        let mut gwrapper: GraphWrapper = GraphWrapper::new();
        gwrapper.fromNGfa(&graph, del);

        write_genomes(&gwrapper);

        eprintln!("{} genomes and {} paths", gwrapper.genomes.len(), graph.paths.len());

        // Make matrix
        //let h = test1(&gwrapper, &graph);
        let mut mat_node: MatrixWrapper<u32>;
        let mut mat_dir: MatrixWrapper<(u32, bool)>;
        let mut mat_edge: MatrixWrapper<(u32, bool, u32, bool)>;

        if matches.is_present("type"){
            let values: &str = matches.value_of("type").unwrap();
            if values.contains('n'){
                mat_node = matrix_node(&gwrapper, &graph);
                i = & mut mat_node.matrix;
                let rem = i.filter();
                mat_node.remove_stuff(rem);
                mat_node.write(type_out, _output, "node");
            }
            if values.contains('e'){
                mat_edge = matrix_edge(&gwrapper, &graph);
                i = & mut mat_edge.matrix;
                let rem = i.filter();
                mat_edge.remove_stuff(rem);
                mat_edge.write(type_out, _output, "edge");

            }
            if values.contains('d'){
                mat_dir = matrix_dir_node(&gwrapper, &graph);
                i = & mut mat_dir.matrix;
                let rem = i.filter();
                mat_dir.remove_stuff(rem);
                mat_dir.write(type_out,  _output, "dir");
            }
        } else {
            mat_node = matrix_node(&gwrapper, &graph);
            i = & mut mat_node.matrix;
            let rem = i.filter();
            mat_node.remove_stuff(rem);
            mat_node.write(type_out, _output, "node");
        }
    }


    else {
        eprintln!("No input!");
        process::exit(0x0100);
    }



}

#[cfg(test)]
mod main {
    use crate::matrix_wrapper::{matrix_node_coverage2, matrix_node, matrix_dir_node, matrix_edge, matrix_node_coverage, matrix_node10, MatrixWrapper2};
    use gfaR_wrapper::{GraphWrapper, NGfa};
    use bimap::BiMap;
    use crate::writer::write_reduce;

    #[test]
    fn names() {
        let mut index_normal: BiMap<u32, usize> = BiMap::new();
        let mut matrix = MatrixWrapper2::new();
        let mut graph = NGfa::new();

        graph.from_graph("/home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa");
        let mut gwrapper: GraphWrapper = GraphWrapper::new();
        gwrapper.fromNGfa(&graph, " ");
        if index_normal.is_empty(){
            println!("HOLDDSADAS");
        }
        matrix_node10(&gwrapper, &graph, & mut matrix, & mut index_normal);
        if index_normal.is_empty(){
            println!("dajksdjakda");
        } else {
            println!("daksljdklasjdas");
        }
        matrix.write_names();
        matrix.remove_genomes("holyshit12");

        matrix.make_binary(1);
        let k = matrix.reduce_comb1();
        println!("HOLY {}", matrix.matrix_bin.len());
        write_reduce(&k.0, &k.1);

    }


}
