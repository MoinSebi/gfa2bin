use assert_cmd::prelude::*; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn gfa_edges() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.edges")
        .arg("-f")
        .arg("edge")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut b = File::open("data/output/graph.edges.edge.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (11 * 2));

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f edge
fn gfa_nodes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.nodes")
        .arg("-f")
        .arg("node");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.nodes.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (8 * 2));

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f dirnode
fn gfa_dir_node() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.dirnodes")
        .arg("-f")
        .arg("dirnode");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.dirnodes.dirnode.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (8 * 2));

    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}
