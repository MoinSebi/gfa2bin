mod alignment;
mod block;
mod core;
mod filter;
mod find;
mod graph;
mod helper;
mod logging;
mod r#mod;
mod subpath;
mod view;
mod window;

use crate::alignment::align_main::align_main;
use crate::block::block_main::block_main;
use crate::filter::filter_main::filter_main;
use crate::find::find_main::find_main;
use crate::graph::graph_main::graph_main;
use crate::logging::newbuilder;
use crate::r#mod::mod_main::mod_main;
use crate::subpath::subpath_main::subpath_main;
use crate::view::view_main::view_main;
use crate::window::window_main::window_main;
use clap::{App, AppSettings, Arg};

fn main() {
    let matches = App::new("gfa2bin")
        .version("0.1.0")
        .author("Sebastian V")
        .about("gfa2bin")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .about("verbose "),
        )
        .subcommand(
            App::new("graph")
                .about("Conversion from a graph directly (GFA1 format)")
                .version("0.1.0")
                .setting(AppSettings::ArgRequiredElseHelp)


                .help_heading("Input parameters")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("Sets the input file to use")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("feature")
                        .short('f')
                        .long("feature")
                        .about("Specific the observed feature (node, dirnode, edge)")
                        .takes_value(true)
                        .default_value("node"),
                )
                .arg(
                    Arg::new("pansn")
                        .long("pansn")
                        .about("PanSN-spec separator")
                        .takes_value(true),
                )

                .help_heading("Modifier")
                .arg(
                    Arg::new("paths")
                        .long("paths")
                        .about("Ignore these paths (one per line)")
                        .takes_value(true),
                )
                .help_heading("Thresholds")
                .arg(Arg::new("absolute-threshold")
                    .short('a')
                    .long("absolute-threshold")
                    .about("Set a absolute threshold")
                    .takes_value(true)
                    .display_order(0))
                // Modification
                .arg(Arg::new("method")
                    .short('m')
                    .long("method")
                    .about("Normalization method (mean|median|percentile|nothing) [default: nothing]")
                    .takes_value(true)
                    .display_order(1)
                )

                .arg(Arg::new("fraction")
                    .long("fraction")
                    .about("Fraction")
                    .takes_value(true)
                    .display_order(2)
                )
                .arg(Arg::new("standard-deviation")
                    .short('s')
                    .long("std")
                    .about("Adjust your threshold by decreasing if by X * standard deviation")
                    .takes_value(true)
                    .display_order(2))
                .arg(Arg::new("non-covered")
                    .long("non-covered")
                    .about("Include non-covered entries (nodes or sequences) for dynamic threshold calculations (e.g mean)")
                    .display_order(4)
                )



                .help_heading("Output parameter")
                .arg(
                    Arg::new("split")
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .about("Output prefix for all files")
                        .required(true),
                )
                .arg(
                    Arg::new("bimbam")
                        .long("bimbam")
                        .about("Output in BIMBAM format")),
        )
        .subcommand(
            App::new("align")
                .about("Conversion from a alignment (pack or bpack)")
                .version("0.1.0")
                .setting(AppSettings::ArgRequiredElseHelp)


                .help_heading("Input parameters")
                .arg(
                    Arg::new("pack")
                        .short('p')
                        .long("packlist")
                        .about("List of pack files")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("pack compressed")
                        .short('c')
                        .long("packcompressed")
                        .about("concatenated bpack file")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("pc list")
                        .long("pc-list")
                        .about("File with pc files (one per line)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("index")
                        .short('i')
                        .long("index")
                        .about("Index file")
                        .takes_value(true),
                )


                .help_heading("Entry thresholds")
                .arg(Arg::new("absolute-threshold")
                    .long("eabsolute-threshold")
                    .about("Set a absolute threshold")
                    .takes_value(true)
                    .display_order(0))
                // Modification
                .arg(Arg::new("emethod")
                    .long("emethod")
                    .about("Normalization method (mean|median|percentile|nothing) [default: nothing]")
                    .takes_value(true)
                    .display_order(1)
                )

                .arg(Arg::new("efraction")
                    .long("efraction")
                    .about("Fraction")
                    .takes_value(true)
                    .display_order(2)
                )
                .arg(Arg::new("estandard-deviation")
                    .long("estd")
                    .about("Adjust your threshold by decreasing if by X * standard deviation")
                    .takes_value(true)
                    .display_order(2))
                .arg(Arg::new("enon-covered")
                    .long("enon-covered")
                    .about("Include non-covered entries (nodes or sequences) for dynamic threshold calculations (e.g mean)")
                    .display_order(4)
                )





                .help_heading("Output parameter")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("split")
                        .short('s')
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                )
                .arg(
                    Arg::new("bimbam")
                        .long("bimbam")
                        .about("Output bimbam format [default: plink]"),
                ),
        )




        .subcommand(
            App::new("mod")
                .about("Remove samples from you graph/alignment-based plink file")
                .version("0.1.0")
                .help_heading("Input parameters")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
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
                        .long("paths")
                        .about("Path to remove (one per line)")
                        .takes_value(true),
                )
                .help_heading("Intrinsic Modification")
                .arg(
                    Arg::new("non-info")
                        .long("non-info")
                        .about(
                            "Remove all entries which hold no information (all true or all false)",
                        )
                        .takes_value(true),
                )
                .help_heading("Output parameter")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("split")
                        .short('s')
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                )
                .arg(
                    Arg::new("index")
                        .short('i')
                        .long("index")
                        .takes_value(true)
                        .about("Remove the entries on this specific index (0-based)"),
                ),
        )
        // Will work on this later
        .subcommand(
            App::new("find")
                .version("1.0.1")
                .about("Find features in the graph and return a BED file for further analysis")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("Sets the input file to use")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("features")
                        .short('f')
                        .long("feature")
                        .about("Feature list to remove (one per line)")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .about("Length of the feature")
                        .takes_value(true)
                        .default_value("200"),
                ),
        )
        .subcommand(
            App::new("window")
                .version("1.0.1")
                .about("Find features in the graph and return a BED file for further analysis")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .about("Length of the feature")
                        .takes_value(true)
                        .default_value("200"),
                ),
        )
        .subcommand(
            App::new("subpath")
                .version("1.0.1")
                .about("Find features in the graph and return a BED file for further analysis")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("graph")
                        .about("GFA input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .about("Length of the feature")
                        .takes_value(true)
                        .default_value("200"),
                ),
        )



        .subcommand(
            App::new("block")
                .version("1.0.1")
                .about("Genotyping by pan-genomic blocks")

                .help_heading("Input parameters")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("graph")
                        .about("GFA input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("PanSN")
                        .long("PanSN")
                        .about("PanSN-spec separator")
                        .takes_value(true),
                )


                .help_heading("Parameter")
                .arg(
                    Arg::new("window")
                        .short('w')
                        .long("window")
                        .about("Window size (in nodes)")
                        .takes_value(true)
                        .default_value("500"),
                )
                .arg(
                    Arg::new("step")
                        .short('s')
                        .long("step")
                        .about("Step")
                        .takes_value(true)
                        .default_value("100"),
                )
                .arg(Arg::new("distance")
                    .short('d')
                    .long("distance")
                    .about("Distance till breaking the block")
                    .takes_value(true)
                    .default_value("10000"))


                .help_heading("Output parameter")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("split")
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                ),
        )

        // This is fine
        .subcommand(
            App::new("view")
                .version("1.0.0")
                .about("Convert BED to VCF")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("PLINK input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                ),
        )


        .subcommand(
            App::new("filter")
                .version("1.0.0")
                .about("Filter a PLINK file")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("maf")
                        .short('m')
                        .long("maf")
                        .about("Major allele frequency")
                        .takes_value(true)
                        .default_value("0.05"),
                )
                .arg(
                    Arg::new("MAF")
                        .short('M')
                        .long("MAF")
                        .about("Minor allele frequency")
                        .takes_value(true)
                        .default_value("0.95"),
                )
                .arg(
                    Arg::new("entry min")
                        .short('e')
                        .long("entries-min")
                        .about("Minimum amount of entries")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("entry max")
                        .short('E')
                        .long("entries-max")
                        .about("Maximum amount of entries")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("split")
                        .long("split")
                        .takes_value(true)
                        .about("Split output in multiple files"),
                ),
        )
        .get_matches();

    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa

    // Checking verbose
    newbuilder(&matches);

    if let Some(matches) = matches.subcommand_matches("graph") {
        graph_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("align") {
        align_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("mod") {
        mod_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("find") {
        find_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("window") {
        window_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("subpath") {
        subpath_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("block") {
        block_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("view") {
        view_main(matches);
    } else if let Some(matches) = matches.subcommand_matches("filter") {
        filter_main(matches);
    } else {
        println!("No subcommand was used");
    }
}
