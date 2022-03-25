mod counting;
mod helper;
mod convert2;
mod find;


use clap::{App, Arg};
use gfaR_wrapper::{GraphWrapper, NGfa};
use std::process;
use crate::helper::get_thresh;
use bimap::BiMap;
use std::env::args;
use chrono::Local;
use env_logger::{Builder, Target};
use log::{info, LevelFilter, warn};
use std::io::Write;
use crate::convert2::core::{MatrixWrapper, remove_bimap};
use crate::convert2::gfa::{matrix_dir_node, matrix_edge, matrix_edge2, matrix_node, matrix_node_wrapper2};
use crate::convert2::pack::{matrix_pack_bit, matrix_pack_u16};
use crate::convert2::writer::{write_bed_split, write_bimhelper, write_matrix, write_reduce};


fn main() {
    let matches = App::new("panSV")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .about("verbose "))
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
            .arg(Arg::new("delimiter")
                .short('d')
                .long("delimiter")
                .about("Delimiter between genome and chromosome number")
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
        .arg(Arg::new("threshold")
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
            .arg(Arg::new("split")
                .long("split")
                .about("Split the output bed and bim"))
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
                    .takes_value(true))
                        .arg(Arg::new("threads")
                            .short('t')
                            .takes_value(true)
                            .about("Number of threads if multithreading")))
        .subcommand(App::new("find")
            .version("1.0.1")
            .about("Find significant hits")
            .arg(Arg::new("GWAS Output")
                .short('g')
                .long("GEMMA")
                .about("GEMMA ouptut file")
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

    // Checking verbose
    if matches.is_present("verbose"){
        Builder::new()
            .format(|buf, record| {
                writeln!(buf,
                         "{} [{}] - {}",
                         Local::now().format("%Y-%m-%dT%H:%M:%S"),
                         record.level(),
                         record.args()
                )
            })
            .filter(None, LevelFilter::Trace)
            .target(Target::Stderr)
            .init();
    } else {
        Builder::new()
            .format(|buf, record| {
                writeln!(buf,
                         "{} [{}] - {}",
                         Local::now().format("%Y-%m-%dT%H:%M:%S"),
                         record.level(),
                         record.args()
                )
            })
            .filter(None, LevelFilter::Info)
            .target(Target::Stderr)
            .init();
    }

    /// You want to convert stuff
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
        if matches.is_present("delimiter") {
            del = matches.value_of("delimiter").unwrap();
        } else {
            del = " ";
        }


        let mut matrix = MatrixWrapper::new();
        let mut index_normal: BiMap<u32, usize> = BiMap::new();
        let mut index_dir: BiMap<(u32, bool), usize> = BiMap::new();
        let mut index_edge: BiMap<(u32, bool, u32, bool), usize> = BiMap::new();

        // Check if gfa or coverage
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
                    matrix_node_wrapper2(&gwrapper, &graph, &mut matrix, &mut index_normal, &2);
                }
                if values.contains('e') {
                    matrix_edge2(&gwrapper, &graph, &mut matrix, &mut index_edge, &2);
                }
                if values.contains('d') {
                    matrix_dir_node(&gwrapper, &graph, &mut matrix, &mut index_dir, &2);
                }
            } else {
                matrix_node_wrapper2(&gwrapper, &graph, &mut matrix, &mut index_normal, &(10 as usize));
            }
        } else {
            if matches.is_present("pack") {
                let file_pack = matches.value_of("pack").unwrap();
                let j = get_thresh(&file_pack);
                if j == 0 {
                    info!("Reading u16 pack");
                    matrix_pack_u16(file_pack, &mut matrix, &mut index_normal);
                } else {
                    info!("Reading bit pack");
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
            info!("Make binary");
            matrix.make_binary(1);
        }
        let mut remove_this: Vec<u32> = Vec::new();
        if matches.is_present("filter") {
            info!("Filtering");
            remove_this = matrix.filter();
        }

        if matches.is_present("reduce") {
            info!("Reducing combinations");
            let k = matrix.reduce_combinations_test();
            write_reduce(&k.0, &k.1, _output, "gfa2bin");
        }


        if matches.is_present("split") {
            info!("Splitting matrix");
            let o = matrix.split_bin(10);
            for (index, x) in o.enumerate() {
                write_bed_split(x, "testsplit", &*index.to_string())
            }
        }

        //Output
        if matches.is_present("type") {
            let values: &str = matches.value_of("type").unwrap();
            if values.contains('n') {
                info!("Node type");
                write_matrix(&mut matrix, type_out, _output, "node");
            }
            if values.contains('e') {
                info!("Edge type");
                write_matrix(&mut matrix, type_out, _output, "edge");
            }
            if values.contains('d') {
                info!("Directed ");
                write_matrix(&mut matrix, type_out, _output, "dir");
            }
        } else {
            write_matrix(&mut matrix, type_out, _output, "gfa2bin");
        }

        info!("Remove and writing");
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


