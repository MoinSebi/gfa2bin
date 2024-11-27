mod core;
mod cov;
mod filter;
mod find;
mod graph;
mod helper;
mod logging;
mod merge;
mod nearest;
mod remove;
mod split;
mod subpath;
mod view;
mod window;

use crate::filter::filter_main::filter_main;
use crate::find::find_main::find_main;
use crate::graph::graph_main::graph_main;
use crate::logging::newbuilder;
use crate::merge::merge_main::merge_main;
use crate::remove::remove_main::remove_main;
use crate::split::split_main::split_main;
use crate::subpath::subpath_main::subpath_main;
use crate::view::view_main::view_main;
use crate::window::window_main::window_main;
use clap::{App, AppSettings, Arg};

use crate::cov::cov_main::cov_main;
use crate::nearest::nearest_main::nearest_main;
use std::error::Error;


fn main() -> Result<(), Box<dyn Error>> {
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
                .about("Convert GFA file (v1) to plink (bed, bim, fam). \n \
                                Threshold modifier are used used to create a (1) presence-absence matrix or (2) scaling in bimbam format.")

                .setting(AppSettings::ArgRequiredElseHelp)


                .help_heading("Input options")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("Input GFA file")
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
                .arg(Arg::new("method")
                    .short('m')
                    .long("method")
                    .about("Normalization method (mean|median|percentile|nothing) [default: nothing]")
                    .takes_value(true)
                    .display_order(1)
                )
                .arg(Arg::new("fraction")
                    .long("fraction")
                    .about("Adjust your threshold by multiplying it by this value")
                    .takes_value(true)
                    .display_order(2)
                )
                .arg(Arg::new("standard-deviation")
                    .short('s')
                    .long("std")
                    .about("Adjust your threshold by decreasing if by X * standard deviation")
                    .takes_value(true)
                    .display_order(2))
                .arg(Arg::new("keep-zeros")
                    .long("keep-zeros")
                    .about("Include non-covered entries (nodes or sequences) for dynamic threshold calculations (e.g mean) [default: false]")
                    .display_order(4)
                )
                .arg(Arg::new("max")
                    .long("max")
                    .about("Use max value for scaling (only in BIMBAM)")
                    .display_order(5)
                )



                .help_heading("Output parameter")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .about("Output prefix for all files")
                        .required(true),
                )
                .arg(Arg::new("pheno")
                    .long("pheno")
                         .about("Phenotype value")
                         .takes_value(true)
                )
                .arg(
                    Arg::new("bimbam")
                        .long("bimbam")
                        .about("Output in BIMBAM format [default: plink]"),),
        )



        .subcommand(
            App::new("cov")
                .about("Conversion from a coverage information (plain-text or compressed)")
                .setting(AppSettings::ArgRequiredElseHelp)


                .help_heading("Input parameters")
                .arg(
                    Arg::new("pack")
                        .short('p')
                        .long("packlist")
                        .about("List of plain-text pack files (one per line). Tab separated with sample name. Example [tair10   /path/to/file]")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("pack compressed")
                        .short('c')
                        .long("cp")
                        .about("Concatenated bpack file (packing)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("pc-list")
                        .short('l')
                        .long("cp-list")
                        .about("List of compressed pack files (one per line). Tab separated with sample name. Example [tair10   /path/to/file]")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("index")
                        .short('i')
                        .long("index")
                        .about("Index file is needed for compressed pack")
                        .takes_value(true),
                )


                .help_heading("Threshold modifier")
                .arg(Arg::new("absolute-threshold")
                    .long("absolute-threshold")
                    .about("Set a absolute threshold")
                    .takes_value(true)
                    .display_order(0))
                // Modification
                .arg(Arg::new("method")
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
                .arg(Arg::new("keep-zeros")
                    .long("keep-zeros")
                    .about("Include non-covered entries (nodes or sequences) for dynamic threshold calculations (e.g mean)")
                    .display_order(4)
                )
                .arg(Arg::new("node")
                    .long("node")
                    .about("Use node instead of sequence")
                    .display_order(5)
                )
                .arg(Arg::new("no-default")
                    .long("no-default")
                    .about("Do not use default values (only works with bimbam)")
                    .display_order(6)
                )




                .help_heading("Output parameter")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(Arg::new("pheno")
                    .long("pheno")
                    .about("Phenotype value")
                    .takes_value(true)
                )
                .arg(
                    Arg::new("bimbam")
                        .long("bimbam")
                        .about("Output bimbam format [default: plink]"),
                ),
        )




        .subcommand(
            App::new("remove")
                .about("Remove samples or genotypes from you graph/coverage-based plink file")
                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file (prefix)")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Modification parameters")
                .arg(
                    Arg::new("genotypes")
                        .short('g')
                        .long("genotypes")
                        .about("List of genotypes to remove (one per line)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("gindex")
                        .long("gindex")
                        .takes_value(true)
                        .about("List of index (genotypes) to remove (one per line, (0-based))"),
                )
                .arg(
                    Arg::new("samples")
                        .long("samples")
                        .about("List of samples to remove (one per line)")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("sindex")
                        .long("sindex")
                        .takes_value(true)
                        .about("List of index (samples) to remove (one per line, (0-based))"),
                )

                .arg(
                    Arg::new("non-info")
                        .long("non-info")
                        .about(
                            "Remove all genotypes which hold no information (all true or all false)",
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
        )
        // Will work on this later
        .subcommand(
            App::new("find")
                .about("Find features in the graph and return a BED file for further analysis")
                .help_heading("Input options")
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

                .help_heading("Output options")
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
                .about("Find features in the graph and return a BED file for further analysis")

                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Window options")
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .about("Length of the feature")
                        .takes_value(true)
                        .default_value("200"),
                )
                .arg(
                    Arg::new("blocks")
                        .long("blocks")
                        .short('b')
                        .about("Output blocks [default: false]")
                        .takes_value(true),)

                .help_heading("Output options")
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
            App::new("subpath")
                .about("Find features in the graph and return a BED file for further analysis")

                .help_heading("Input options")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("graph")
                        .about("GFA input file")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Subpath options")
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .about("Length of the feature")
                        .takes_value(true)
                        .default_value("5"),
                )
                .arg(
                    Arg::new("blocks")
                        .long("blocks")
                        .short('b')
                        .about("Output blocks [default: false]")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .short('t')
                        .about("Number of threads")
                        .takes_value(true)
                        .default_value("1")
                )

                .help_heading("Output options")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
        )





        // This is fine
        .subcommand(
            App::new("view")
                .about("Convert BED to VCF")

                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("PLINK input file")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Output options")
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
                .about("Filter a PLINK file. (1) 'SNPs' by allele frequency or (2) path/samples by matrix coverage")

                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink input file")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Filter options")
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

                .help_heading("Output options")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
        )


        .subcommand(
            App::new("merge")
                .about("Merge the multiple plink files into one. Must be the same samples (fam).")

                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("List of BED files (one per line). Bim and Fam files should have the same prefix.")
                        .takes_value(true)
                        .required(true)

                )

                .help_heading("Output options")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink file")
                        .takes_value(true)
                        .required(true),
                )
        )

        .subcommand(
            App::new("split")
                .about("Split PLINK files into multiple files")
                .help_heading("Input options")
                .arg(
                    Arg::new("plink")
                        .short('p')
                        .long("plink")
                        .about("Plink prefix OR BED files")
                        .takes_value(true)
                        .required(true)

                )

                .help_heading("Split parameter")
                .arg(
                    Arg::new("splits")
                        .short('s')
                        .long("splits")
                        .about("Number of splits")
                        .takes_value(true)
                        .required(true)
                )

                .help_heading("Output options")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink files")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Performance")
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .about("Number of threads")
                        .takes_value(true)
                        .default_value("1")
                )
        )

        .subcommand(
            App::new("nearest")
                .about("Nearest node to reference node")

                .help_heading("Input options")
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .about("GFA input file")
                        .takes_value(true)
                        .required(true),
                )

                .help_heading("Reference options")
                .arg(Arg::new("prefix")
                    .short('p')
                    .long("prefix")
                    .about("Prefix of the reference nodes")
                    .takes_value(true))
                .arg(Arg::new("references")
                    .long("references")
                    .about("Reference nodes")
                    .takes_value(true))

                .help_heading("Node options ")
                .arg(Arg::new("nodes")
                    .long("nodes")
                    .about("Nodes to find the nearest reference node")
                    .takes_value(true))

                .help_heading("Output options")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .about("Output prefix for the new plink files")
                        .takes_value(true)
                        .required(true),
                )
        )
        .get_matches();

    //cargo run -- -g /home/svorbrugg_local/Rust/data/AAA_AAB.cat.gfa

    // Checking verbose
    newbuilder(&matches);

    if let Some(matches) = matches.subcommand_matches("graph") {
        graph_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("cov") {
        cov_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("remove") {
        remove_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("find") {
        find_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("window") {
        window_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("subpath") {
        subpath_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("view") {
        view_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("filter") {
        filter_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("merge") {
        merge_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("split") {
        split_main(matches)
    } else if let Some(matches) = matches.subcommand_matches("nearest") {
        nearest_main(matches)
    } else {
        println!("No subcommand was used");
        Ok(())
    }
}
