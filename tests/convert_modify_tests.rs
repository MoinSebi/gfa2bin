use assert_cmd::prelude::*; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

// cargo run -- graph -g data/example_data/testGraph2.gfa -f node -o data/example_data/node.mod

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
        .arg("data/output/tt3.node")
        .arg("-f")
        .arg("data/example_data/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/tt3.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + ((8-2) * 2));
    assert_eq!(buffer[3], 127);

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
        .arg("data/output/tt10.node")
        .arg("--paths")
        .arg("data/example_data/mod_path.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/tt10.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 8);
    assert_eq!(buffer[3], 31);

    Ok(())
}
