mod counting;
mod helper;
mod convert;
mod find;


use clap::{App, Arg};
use gfaR_wrapper::{GraphWrapper, NGfa};
use std::process;
use crate::helper::{get_thresh, transpose_bitvec, transpose_generic};
use bimap::BiMap;
use chrono::Local;
use env_logger::{Builder, Target};
use log::{info, LevelFilter, warn};
use std::io::Write;
use std::path::Path;
use crate::convert::core::{MatrixWrapper, remove_bimap, write_bim2};
use crate::convert::fam::Fam;
use crate::convert::gfa::{matrix_dir_node, matrix_edge, matrix_node_wrapper};
use crate::convert::pack::{matrix_pack_bit, matrix_pack_u16};
use crate::convert::writer::{write_bed_split, write_bimhelper, write_matrix, write_reduce};


fn main() {
    let matches = App::new("gfa2bin")
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
            .arg(Arg::new("fam")
                .short('f')
                .long("fam")
                .about("For filtering and removing information from the matrix")
                .takes_value(true))
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
            .takes_value(true)
            .hidden(true))
        .arg(Arg::new("smart")
            .long("smart")
            .short('j')
            .about("Reduce the number of tests"))
        .arg(Arg::new("threads")
            .long("threads")
            .about("Number of threads if multithreading")
            .takes_value(true)
            .default_value("1"))
        .arg(Arg::new("number")
            .short('n')
            .about("Number of output")
            .takes_value(true))
        )



        // Will work on this later
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



    // Check subcommands
    if let Some(ref matches) = matches.subcommand_matches("convert") {
        let threads: usize = matches.value_of("threads").unwrap().parse().unwrap();
        // Check if input exists
        if !(matches.is_present("gfa") | matches.is_present("pack")) {
            eprintln!("No input");
            process::exit(0x0100);
        }
        // Check additional fam file
        let mut _fam = Fam::new();
        if matches.is_present("fam"){
            _fam = Fam::from_file(matches.value_of("fam").unwrap())

        }


        // This is to decide which output
        // bed is binary, bimbam is normalized f16 (or something), all is both
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
            let mut _input: &str = "not relevant";
            if Path::new(matches.value_of("gfa").unwrap()).exists() {
                _input = matches.value_of("gfa").unwrap();
            } else {
                warn!("No file with such name");
                process::exit(0x0100);
            }
            let _output: &str = matches.value_of("output").unwrap();
            // Read the graph
            let mut graph = NGfa::new();
            graph.from_graph(_input);

            // Make graph, wrapper
            let mut gwrapper: GraphWrapper = GraphWrapper::new();
            gwrapper.from_ngfa(&graph, del);
            info!("test {:?}", gwrapper.genomes.iter().map(|x| x.0.clone()));

            // Nodes, edges, or directed nodes
            if matches.is_present("type") {
                let values: &str = matches.value_of("type").unwrap();
                if values.contains('n') {
                    matrix_node_wrapper(&gwrapper, &graph, &mut matrix, &mut index_normal, &threads);
                }
                if values.contains('e') {
                    matrix_edge(&gwrapper, &graph, &mut matrix, &mut index_edge, &threads);
                }
                if values.contains('d') {
                    matrix_dir_node(&gwrapper, &graph, &mut matrix, &mut index_dir, &threads);
                }
            } else {
                matrix_node_wrapper(&gwrapper, &graph, &mut matrix, &mut index_normal, &(threads));
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
        //--------------------------------------------------------------------------------------------------------------------------





        // We only transpose once!
        if matrix.matrix_bin.is_empty(){
            let k: Vec<String> = matrix.column_name.values().cloned().collect();
            let mut hit = 0;

            if !_fam.family_id.is_empty(){
                for (key, value) in k.iter().enumerate(){
                    if !_fam.family_id.contains(value){
                        matrix.remove_col((key as u32) - hit);
                        hit += 1;
                    }
                }
            }

            matrix.matrix_core = transpose_generic(&matrix.matrix_core);
            matrix.shape = (matrix.matrix_core.len(), matrix.matrix_core[0].len());
        } else {
            matrix.matrix_bin = transpose_bitvec(&matrix.matrix_bin);
            matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());

        };

        info!("Shape is {:?}", matrix.shape);







        if matches.is_present("names") {
            matrix.write_names(_output);
        }

        // This is not working atm
        if matches.is_present("removeNames") {
            matrix.remove_genomes(matches.value_of("removeNames").unwrap())
        }


        // Make dummy binary
        if matrix.matrix_bin.is_empty() {
            info!("Make binary");
            matrix.make_binary(1);
        }

        // Apply a filter
        let mut remove_this: Vec<u32> = Vec::new();
        if matches.is_present("filter") {
            info!("Filtering");
            remove_this = matrix.filter_shared();
        }

        if matches.is_present("reduce") {
            info!("Reducing combinations");
            let k = matrix.reduce_combinations_test();
            write_reduce(&k.0, &k.1, _output, "gfa2bin");
        }


        if matches.is_present("split") {
            info!("Splitting matrix");
            let o = matrix.split_bin(10);
            println!("long {}", matrix.matrix_bin.len());
            println!("long {}", matrix.matrix_core.len());
            let size = ((matrix.matrix_bin.len() as f64)/(10 as f64)).ceil() as usize;
            let mut count = 0;
            for (index, x) in o.enumerate() {
                println!("index {}", index);
                write_bed_split(x, _output, &*index.to_string());
                if count+size > matrix.matrix_bin.len(){
                    count = matrix.matrix_bin.len()-size;
                }
                write_bim2(count, count+size, _output, &*index.to_string());
                count += size;
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


