use assert_cmd::prelude::*;
use std::fs; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn nearest() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("nearest")
        .arg("-g")
        .arg("data/example_data/gfa/testGraph.gfa")
        .arg("-p")
        .arg("a")
        .arg("-o")
        .arg("data/example_data/nearest_output.txt")
        .arg("-t")
        .arg("4");
    cmd.assert().success();

    let mut b = File::open("data/example_data/nearest_output.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    println!("Content: {}", content);

    assert!(content.contains("a#1#Chr1"));

    // Remove files
    fs::remove_file("data/example_data/nearest_output.txt")?;
    Ok(())
}

#[test]
#[should_panic]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn nearest_fail() {
    let mut cmd = Command::cargo_bin("gfa2bin").expect("Error executing binary");
    cmd.arg("nearest")
        .arg("-g")
        .arg("data/example_data/gfa/testGraph.gfa")
        .arg("-p")
        .arg("a")
        .arg("-o")
        .arg("data/example_data/nearest_output1.txt")
        .arg("-t")
        .arg("4");
    cmd.assert().success();

    let mut b = File::open("data/example_data/nearest_output1.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert!(content.contains("GCF_020526725"));

    // Remove files
    fs::remove_file("data/example_data/nearest_output.txt").expect("Error removing file");
}
