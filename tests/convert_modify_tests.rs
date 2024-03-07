use assert_cmd::prelude::*; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn mod1() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("mod")
        .arg("-p")
        .arg("data/output/tt2.node")
        .arg("-o")
        .arg("data/output/tt3.node")
        .arg("-f")
        .arg("data/example_data/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/tt3.node.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + ((9-2) * 2));

    Ok(())
}
