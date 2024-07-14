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
fn remove_feature() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("remove")
        .arg("-p")
        .arg("./jo1231")
        .arg("-o")
        .arg("data/output/remove.feature.node")
        .arg("-f")
        .arg("data/example_data/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/remove.feature.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + ((9 - 2) * 2));
    assert_eq!(buffer[3], 127);
    fs::remove_file("data/output/remove.feature.node.bed")?;
    fs::remove_file("data/output/remove.feature.node.bim")?;
    //fs::remove_file("data/output/remove.feature.node.fam")?;
    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn remove_path() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("remove")
        .arg("-p")
        .arg("./jo1231")
        .arg("-o")
        .arg("data/output/remove.path.node")
        .arg("--paths")
        .arg("data/example_data/mod_path.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/remove.path.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 9);
    assert_eq!(buffer[3], 31);
    fs::remove_file("data/output/remove.path.node.bed")?;
    fs::remove_file("data/output/remove.path.node.bim")?;
    fs::remove_file("data/output/remove.path.node.fam")?;

    Ok(())
}
