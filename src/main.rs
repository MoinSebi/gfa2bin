mod matrix;
mod counting;
mod matrix_wrapper;
mod helper;
use clap::{App, Arg};
use gfaR_wrapper::{NGfa, GraphWrapper};
use crate::matrix_wrapper::{matrix_node, matrix_edge, matrix_dir_node, MatrixWrapper, matrix_node_coverage, matrix_node_coverage2};
use packing_lib::helper::u8_u16;
use packing_lib::reader::get_file_as_byte_vec;
use std::process;
use crate::matrix::Matrix;


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




    // Comment: Because different wrapper are used, we need to check everything by itself
    let mut i: & mut Matrix;
    // Read from a pack file
    if matches.is_present("pack") {
        let ii = matches.value_of("pack").unwrap();
        let j = u8_u16(&mut & get_file_as_byte_vec(ii)[7..9]);
        let mut gw: MatrixWrapper<u32>;
        if j == 0{
            gw = matrix_node_coverage2(ii);
        } else {
            gw = matrix_node_coverage(ii);
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
    use packing_lib::helper::u8_u16;
    use packing_lib::reader::get_file_as_byte_vec;
    use crate::matrix_wrapper::{MatrixWrapper, matrix_node_coverage2, matrix_node, matrix_dir_node, matrix_edge};
    use gfaR_wrapper::{GraphWrapper, NGfa};

    #[test]
    fn exploration() {
        let j = u8_u16(&mut & get_file_as_byte_vec("/home/svorbrugg_local/Rust/gfa2bin/pack.test.bin")[7..9]);
        eprintln!("Number {}", j);
        let mut gw: MatrixWrapper<u32> = matrix_node_coverage2("/home/svorbrugg_local/Rust/gfa2bin/pack.test.bin");

        gw.matrix.filter();
    }

    #[test]
    fn normal_run() {
        let input = "/home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa";
        let mut graph = NGfa::new();
        graph.from_graph(input);

        // Make graph, wrapper
        let mut gwrapper: GraphWrapper = GraphWrapper::new();
        gwrapper.fromNGfa(&graph, "_");
        eprintln!("LONG {}", gwrapper.genomes.len());
        //gwrapper.fromNGfa(&graph, " ");
        eprintln!("LONG {}", gwrapper.genomes.len());
        eprintln!("LONG {}", graph.paths.len());
        let mat_node = matrix_edge(&gwrapper, &graph);
        let o = mat_node.matrix.filter();
        mat_node.matrix.reduce_comb();

    }

}
