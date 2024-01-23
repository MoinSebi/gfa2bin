mod find;
mod helper;
mod logging;
mod alignment;
mod core;
mod graph;
mod r#mod;

use crate::alignment::align_main::align_main;
use crate::graph::graph_main::graph_main;
use crate::logging::newbuilder;
use clap::{App, Arg};
use crate::r#mod::mod_main::mod_main;

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

        .subcommand(
            App::new("mod")
                .about("Modify a plink file made by a graph")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("graph")
                        .short('g')
                        .long("graph")
                        .about("Graph file")
                        .takes_value(true),
                )

                .arg(
                    Arg::new("features")
                        .short('f')
                        .long("feature")
                        .about("Feature list to remove (one per line)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("paths")
                        .short('p')
                        .long("paths")
                        .about("Path to remove (one per line)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true),
                )
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
        mod_main(&matches);
    }
}
