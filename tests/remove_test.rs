use assert_cmd::prelude::*;
use std::fs; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

// cargo run -- graph -g data/example_data/testGraph.gfa -f node -o data/example_data/node.remove

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn remove_genotype() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.graph.node2")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("remove")
        .arg("-p")
        .arg("./data/output/gfa2bin.graph.node2")
        .arg("-o")
        .arg("data/output/remove.feature.node")
        .arg("--genotypes")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/additional_input/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/remove.feature.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + ((9 - 2) * 2));
    assert_eq!(buffer[3], 127);
    //fs::remove_file("data/output/remove.feature.node.fam")?;
    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn remove_samples() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.remove")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("remove")
        .arg("-p")
        .arg("./data/output/gfa2bin.remove")
        .arg("-o")
        .arg("data/output/gfa2bin.remove.samples")
        .arg("--samples")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/additional_input/mod_path.txt");
    cmd.assert().success();
    fs::remove_file("./data/output/gfa2bin.remove.bed")?;
    fs::remove_file("./data/output/gfa2bin.remove.bim")?;
    fs::remove_file("./data/output/gfa2bin.remove.fam")?;

    let mut b = File::open("data/output/gfa2bin.remove.samples.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 9);
    assert_eq!(buffer[3], 31);
    //fs::remove_file("data/output/remove.feature.node.fam")?;
    fs::remove_file("data/output/gfa2bin.remove.samples.bed")?;
    fs::remove_file("data/output/gfa2bin.remove.samples.bim")?;
    fs::remove_file("data/output/gfa2bin.remove.samples.fam")?;

    Ok(())
}
