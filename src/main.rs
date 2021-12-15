mod matrix;
mod counting;
mod matrix_wrapper;
mod helper;
mod writer;
mod find;

use clap::{App, Arg};
use gfaR_wrapper::{NGfa, GraphWrapper};
use std::process;
use crate::helper::{get_thresh};
use bimap::BiMap;
use crate::writer::{write_reduce, write_bimap, write_bed};
use crate::matrix_wrapper::{MatrixWrapper2, matrix_edge2, matrix_node10, matrix_dir_node2, matrix_pack_u16, matrix_pack_bit, write_matrix, remove_bimap, write_bim, write_bimhelper};
use std::env::args;


fn main() {
    let matches = App::new("panSV")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        // Input
        .subcommand(App::new("convert")
            .about("Convert GFA or PACK format to binary data matrix")
            .version("0.1.0")
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
            .long("rgenomes")
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
        .arg(Arg::new("reduce")
            .long("reduce")
            .about("Reduce to minimum"))
        .arg(Arg::new("filter")
            .long("filter")
            .about("Filter this"))

        // Output
        .arg(Arg::new("genomes")
            .long("genomes")
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
            .short('j')
            .about("Reduce the number of tests"))
            .arg(Arg::new("number")
                .short('n')
                .about("Number of output")
                .takes_value(true)))
        .subcommand(App::new("find")
            .version("1.0.1")
            .about("Find significant hits")
            .arg(Arg::new("GWAS Output")
                .short('g')
                .long("GEMMA")
                .about("GEMMA ouptu file")
                .takes_value(true))
            .arg(Arg::new("Helper")
                .about("helper file")
                .takes_value(true))
            .arg(Arg::new("Number")
                .about(""))
            .arg(Arg::new("One hit"))
            .arg(Arg::new("bed"))
            .arg(Arg::new("sequence"))
            .arg(Arg::new("ACC")))


        .get_matches();


    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa


    if let Some(ref matches) = matches.subcommand_matches("convert") {
        // Check if input
        if !(matches.is_present("gfa") | matches.is_present("pack")) {
            eprintln!("No input");
            process::exit(0x0100);
        }


        // This is to decide which output
        let type_out;
        let _output: &str = matches.value_of("output").unwrap();

        if matches.is_present("bed") & matches.is_present("bimbam") {
            type_out = "all"
        } else if matches.is_present("bed") {
            type_out = "bed"
        } else if matches.is_present("bimbam") {
            type_out = "bimbam"
        } else {
            type_out = "bed"
        }

        // Check if del is set or not
        let del: &str;
        if matches.is_present("delimter") {
            del = matches.value_of("delimter").unwrap();
        } else {
            del = " ";
        }


        let mut matrix = MatrixWrapper2::new();
        let mut index_normal: BiMap<u32, usize> = BiMap::new();
        let mut index_dir: BiMap<(u32, bool), usize> = BiMap::new();
        let mut index_edge: BiMap<(u32, bool, u32, bool), usize> = BiMap::new();
        if matches.is_present("gfa") {
            let _input = matches.value_of("gfa").unwrap();
            let _output: &str = matches.value_of("output").unwrap();
            // Read the graph
            let mut graph = NGfa::new();
            graph.from_graph(_input);

            // Make graph, wrapper
            let mut gwrapper: GraphWrapper = GraphWrapper::new();
            gwrapper.fromNGfa(&graph, del);
            if matches.is_present("type") {
                let values: &str = matches.value_of("type").unwrap();
                if values.contains('n') {
                    matrix_node10(&gwrapper, &graph, &mut matrix, &mut index_normal);
                }
                if values.contains('e') {
                    matrix_edge2(&gwrapper, &graph, &mut matrix, &mut index_edge);
                }
                if values.contains('d') {
                    matrix_dir_node2(&gwrapper, &graph, &mut matrix, &mut index_dir);
                }
            } else {
                matrix_node10(&gwrapper, &graph, &mut matrix, &mut index_normal);
            }
        } else {
            if matches.is_present("pack") {
                let file_pack = matches.value_of("pack").unwrap();
                let j = get_thresh(&file_pack);
                if j == 0 {
                    eprintln!("Reading u16 pack");
                    matrix_pack_u16(file_pack, &mut matrix, &mut index_normal);
                } else {
                    eprintln!("Reading bit pack");
                    matrix_pack_bit(file_pack, &mut matrix, &mut index_normal);
                }
            }
        }

        if matches.is_present("names") {
            matrix.write_names(_output);
        }

        if matches.is_present("removeNames") {
            matrix.remove_genomes(matches.value_of("removeNames").unwrap())
        }


        if matrix.matrix_bin.is_empty() {
            eprintln!("Make binary");
            matrix.make_binary(1);
        }
        let mut remove_this: Vec<u32> = Vec::new();
        if matches.is_present("filter") {
            eprintln!("Filtering");
            remove_this = matrix.filter();
        }

        if matches.is_present("reduce") {
            eprintln!("Reducing combinations");
            let k = matrix.reduce_combinations_test();
            write_reduce(&k.0, &k.1, _output, "gfa2bin");
        }


        if matches.is_present("split") {
            eprintln!("Splitting matrix");
            let o = matrix.split_bin(10);
            for (index, x) in o.enumerate() {
                write_bed(x, "testsplit", &*index.to_string())
            }
        }

        //Output
        if matches.is_present("type") {
            let values: &str = matches.value_of("type").unwrap();
            if values.contains('n') {
                eprintln!("Node type");
                write_matrix(&mut matrix, type_out, _output, "node");
            }
            if values.contains('e') {
                eprintln!("Edge type");
                write_matrix(&mut matrix, type_out, _output, "edge");
            }
            if values.contains('d') {
                eprintln!("Directed ");
                write_matrix(&mut matrix, type_out, _output, "dir");
            }
        } else {
            write_matrix(&mut matrix, type_out, _output, "gfa2bin");
        }

        eprintln!("Remove and writing");
        // THEN FILTER ROWS (BIMAP)
        if !index_normal.is_empty() {
            remove_bimap(&mut index_normal, remove_this);
            //write_bim(& index_normal,_output, "gfa2bin");
            write_bimhelper(&index_normal, _output, "test");
        } else if !index_dir.is_empty() {
            remove_bimap(&mut index_dir, remove_this);
            //write_bim(& index_dir,_output, "gfa2bin");
            write_bimhelper(&index_dir, _output, "test");
        } else if !index_edge.is_empty() {
            remove_bimap(&mut index_edge, remove_this);
            //write_bim(& index_edge,_output, "gfa2bin");
            write_bimhelper(&index_edge, _output, "test");
        }
    }



}

#[cfg(test)]
mod main {
    use gfaR_wrapper::{GraphWrapper, NGfa};
    use bimap::BiMap;
    use crate::writer::write_reduce;
    use crate::matrix_wrapper::{MatrixWrapper2, matrix_node10};

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
        matrix.write_names("test_data/test");
        matrix.remove_genomes("holyshit12");

        matrix.make_binary(1);
        let k = matrix.reduce_combinations_test();
        println!("LL {} ", k.1.len());
        let k2 = matrix.reduce_combinations();
        // println!("LL dasd a{} ", matrix.matrix_bin.len());
        // println!("LL {} {}", k.1.len(), k2.1.len());
        // println!("LL {} {}", k.1[3671], k2.1[3671]);
        // println!("LL {} {}", k.0[3671], k2.0[3671]);
        // println!("HOLY {}", matrix.matrix_bin.len());
        write_reduce(&k2.0, &k2.1, "test_data/test", "gfa2bin");

    }


}
