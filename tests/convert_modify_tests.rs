use assert_cmd::prelude::*;
use std::fs; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

// cargo run -- graph -g data/example_data/testGraph.gfa -f node -o data/example_data/node.mod

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn mod_feature() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("mod")
        .arg("-p")
        .arg("data/example_data/node.mod")
        .arg("-o")
        .arg("data/output/mod.feature.node")
        .arg("-f")
        .arg("data/example_data/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/mod.feature.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + ((8 - 2) * 2));
    assert_eq!(buffer[3], 127);
    fs::remove_file("data/output/mod.feature.node.bed")?;
    fs::remove_file("data/output/mod.feature.node.bim")?;
    fs::remove_file("data/output/mod.feature.node.fam")?;
    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn mod_path1() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("mod")
        .arg("-p")
        .arg("data/example_data/node.mod")
        .arg("-o")
        .arg("data/output/mod.path.node")
        .arg("--paths")
        .arg("data/example_data/mod_path.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/mod.path.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 8);
    assert_eq!(buffer[3], 31);
    fs::remove_file("data/output/mod.path.node.bed")?;
    fs::remove_file("data/output/mod.path.node.bim")?;
    fs::remove_file("data/output/mod.path.node.fam")?;

    Ok(())
}
