use assert_cmd::prelude::*; // Add methods on commands
use std::fs;
use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g
/// -f node
/// --pansn #
///
/// 4 Samples
fn gfa_nodes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("data/output/gfa2bin.graph.node")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.graph.node.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (8 * 2));
    assert_eq!(buffer[3], 127);
    assert_eq!(buffer[4], 0);
    //fs::remove_file("data/output/graph.node.bed")?;
    fs::remove_file("data/output/graph.node.bim")?;
    fs::remove_file("data/output/graph.node.fam")?;

    Ok(())
}

