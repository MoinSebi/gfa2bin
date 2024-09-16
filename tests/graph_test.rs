use assert_cmd::prelude::*; // Add methods on commands
use std::fs;
use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for "gfa2bin graph"
///
/// Nodes
fn gfa_nodes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.graph.node")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.graph.node.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();

    // Buffer should be 8 samples + header
    assert_eq!(buffer.len(), 3 + (9 * 2));

    // First "real" byte is 00
    assert_eq!(buffer[3], 127);
    // Second "real" byte is 000000000
    assert_eq!(buffer[4], 0);
    fs::remove_file("data/output/gfa2bin.graph.node.bed")?;
    fs::remove_file("data/output/gfa2bin.graph.node.bim")?;
    fs::remove_file("data/output/gfa2bin.graph.node.fam")?;

    Ok(())
}

#[test]
/// Test for "gfa2bin graph"
///
/// Dirnodes
fn gfa_dirnodes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("data/output/gfa2bin.graph.dirnode")
        .arg("-f")
        .arg("dirnode")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.graph.dirnode.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();

    // Buffer should be 8 samples + header
    assert_eq!(buffer.len(), 3 + (8 * 2));

    // First "real" byte is 00
    assert_eq!(buffer[3], 127);
    // Second "real" byte is 000000000
    assert_eq!(buffer[4], 0);
    fs::remove_file("data/output/gfa2bin.graph.dirnode.bed")?;
    fs::remove_file("data/output/gfa2bin.graph.dirnode.bim")?;
    fs::remove_file("data/output/gfa2bin.graph.dirnode.fam")?;

    Ok(())
}

#[test]
/// Test for "gfa2bin graph"
///
/// Edges
fn gfa_edges() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("data/output/gfa2bin.graph.edge")
        .arg("-f")
        .arg("edge")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.graph.edge.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();

    // Buffer should be 8 samples + header
    assert_eq!(buffer.len(), 3 + (11 * 2));

    // First "real" byte is 00
    assert_eq!(buffer[3], 63);
    // Second "real" byte is 000000000
    assert_eq!(buffer[4], 0);
    fs::remove_file("./data/output/gfa2bin.graph.edge.bed")?;
    fs::remove_file("./data/output/gfa2bin.graph.edge.bim")?;
    fs::remove_file("./data/output/gfa2bin.graph.edge.fam")?;

    Ok(())
}
