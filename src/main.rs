mod find;
mod helper;
mod logging;
mod plink;

mod alignment;
mod core;
mod graph;

use crate::alignment::align_main::align_main;
use crate::core::core::MatrixWrapper;
use crate::graph::graph_main::graph_main;
use crate::logging::newbuilder;
use clap::{App, Arg};

fn main() {
    let matches = App::new("gfa2bin")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .about("verbose "),
        )
        .subcommand(
            App::new("graph")
                .about("Convert from a graph directly (GFA1 format)")
                .version("0.1.0")
                .help_heading("Input parameters")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("Sets the input file to use")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("Feature")
                        .short('f')
                        .long("Feature")
                        .about("Feature to check")
                        .takes_value(true)
                        .default_value("node"),
                )
                .arg(
                    Arg::new("sep")
                        .long("sep")
                        .about("Separator between genome and chromosome number")
                        .takes_value(true),
                )
                .help_heading("Thresholds")
                .arg(
                    Arg::new("threshold")
                        .short('t')
                        .long("threshold")
                        .about("Set a absolute threshold")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("diploid")
                        .short('d')
                        .long("diploid")
                        .about("Diploid dataset"),
                )
                .help_heading("Output parameter")
                .arg(
                    Arg::new("split")
                        .short('s')
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("Output")
                        .takes_value(true)
                        .about("Output prefix"),
                )
                .arg(
                    Arg::new("format")
                        .long("format")
                        .about("Plink or bimbamOutput plink format (bed, bim, fam)")
                        .default_value("plink"),
                ),
        )
        .subcommand(
            App::new("align")
                .about("Convert from a alignment (pack or bpack)")
                .version("0.1.0")
                .help_heading("Input parameters")
                .arg(
                    Arg::new("pack")
                        .short('p')
                        .long("packlist")
                        .about("List of pack files")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("bpack")
                        .short('b')
                        .long("bpack")
                        .about("concatenated bpack file")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("bpack-list")
                        .long("pback-list")
                        .about("List with all bpack files")
                        .takes_value(true),
                )
                .help_heading("Thresholds")
                .arg(
                    Arg::new("Absolute threshold")
                        .short('t')
                        .long("threshold")
                        .about("Set a absolute threshold")
                        .takes_value(true),
                )
                .help_heading("Output parameter")
                .arg(
                    Arg::new("Plink")
                        .short('p')
                        .long("plink")
                        .about("Output plink format"),
                )
                .arg(
                    Arg::new("bimbam")
                        .short('b')
                        .long("bimbam")
                        .about("Output bimbam format"),
                ),
        )
        // Input
        .subcommand(
            App::new("plink")
                .about("Convert GFA or PACK format to binary data matrix")
                .version("0.1.0")
                .help_heading("Input parameters")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("Sets the input file to use")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("pack")
                        .short('p')
                        .long("pack")
                        .about("if the input is coverage")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("bfile")
                        .long("bfile")
                        .about("bfile")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("delimiter")
                        .short('d')
                        .long("delimiter")
                        .about("Delimiter between genome and chromosome number")
                        .takes_value(true),
                )
                .arg(Arg::new("removeNames").long("rgenomes").about("Only this"))
                .arg(
                    Arg::new("fam")
                        .short('f')
                        .long("fam")
                        .about("For filtering and removing information from the matrix")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("type")
                        .short('t')
                        .long("type")
                        .about("Type of the matrix (when on graph)")
                        .takes_value(true),
                )
                .help_heading("Modification")
                .arg(
                    Arg::new("threshold")
                        .long("Copy number threshold")
                        .about("Normalize to this number")
                        .takes_value(true),
                )
                .arg(Arg::new("reduce").long("reduce").about("Reduce to minimum"))
                .arg(Arg::new("filter").long("filter").about("Filter this"))
                .help_heading("Output parameters")
                .arg(
                    Arg::new("genomes")
                        .long("genomes")
                        .about("Output just the genomes"),
                )
                .arg(
                    Arg::new("split")
                        .long("split")
                        .about("Split the output bed and bim")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix")
                        .takes_value(true)
                        .default_value("gfa2bin.default"),
                )
                .arg(Arg::new("bed").long("bed").about("Output bed + bim file"))
                .arg(
                    Arg::new("bimbam")
                        .long("bimbam")
                        .about("Output bimbam format"),
                )
                .arg(
                    Arg::new("traversal")
                        .long("traversal")
                        .about("Additional traversaloutput file")
                        .takes_value(true)
                        .hidden(true),
                )
                .arg(
                    Arg::new("smart")
                        .long("smart")
                        .short('j')
                        .about("Reduce the number of tests"),
                )
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .about("Number of threads if multithreading")
                        .takes_value(true)
                        .default_value("1"),
                ),
        )
        // Will work on this later
        .subcommand(
            App::new("find")
                .version("1.0.1")
                .about("Find significant hits")
                .arg(
                    Arg::new("GWAS Output")
                        .short('g')
                        .long("GEMMA")
                        .about("GEMMA ouptut file")
                        .takes_value(true),
                )
                .arg(Arg::new("Helper").about("helper file").takes_value(true))
                .arg(Arg::new("Number").about(""))
                .arg(Arg::new("One hit"))
                .arg(Arg::new("bed"))
                .arg(Arg::new("sequence"))
                .arg(Arg::new("ACC")),
        )
        .get_matches();

    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa

    // Checking verbose
    newbuilder(&matches);

    if let Some(matches) = matches.subcommand_matches("graph") {
        graph_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("align") {
        align_main(matches);
    } else if let Some(_matches) = matches.subcommand_matches("mod") {
        println!("dasjkdjas");
    }
}
