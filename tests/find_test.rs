use assert_cmd::prelude::*;
use std::fs; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn find_1() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("find")
        .arg("-g")
        .arg("data/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/find.nodes1.txt")
        .arg("-f")
        .arg("data/example_data/nodes.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/find.nodes1.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("a#1#Chr1"), true);
    fs::remove_file("data/output/find.nodes1.txt")?;
    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn find_2() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("find")
        .arg("-g")
        .arg("data/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/find.nodes2.txt")
        .arg("-f")
        .arg("data/example_data/nodes.txt")
        .arg("-l")
        .arg("0");
    cmd.assert().success();

    let mut b = File::open("data/output/find.nodes2.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("a#1#Chr1"), true);
    fs::remove_file("data/output/find.nodes2.txt")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn find_3() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("find")
        .arg("-g")
        .arg("data/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/find.edges1.txt")
        .arg("-f")
        .arg("data/example_data/edges.txt");
    cmd.assert().success();

    let mut b = File::open("data/output/find.edges1.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("b#1#Chr1"), true);
    fs::remove_file("data/output/find.edges1.txt")?;
    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn find_4() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("find")
        .arg("-g")
        .arg("data/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/find.edges2.txt")
        .arg("-f")
        .arg("data/example_data/edges.txt")
        .arg("-l")
        .arg("0");
    cmd.assert().success();

    let mut b = File::open("data/output/find.edges2.txt").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("b#1#Chr1"), true);
    fs::remove_file("data/output/find.edges2.txt")?;

    Ok(())
}
