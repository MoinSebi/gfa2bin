use assert_cmd::prelude::*; // Add methods on commands

use std::process::Command;

#[test]
/// Convert plink to VCF
///
/// 4 Samples
fn split_test() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.split")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("split")
        .arg("-p")
        .arg("data/output/gfa2bin.split")
        .arg("-s")
        .arg("3")
        .arg("-o")
        .arg("data/output/gfa2bin.split.split1");
    cmd.assert().success();

    Ok(())
}

#[test]
/// Merge a split file
fn merge_test() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.split2")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("split")
        .arg("-p")
        .arg("data/output/gfa2bin.split2")
        .arg("-s")
        .arg("3")
        .arg("-o")
        .arg("data/output/gfa2bin.split2.split");
    cmd.assert().success();

    let status = Command::new("sh")
        .arg("-c")
        .arg("realpath ./data/output/gfa2bin.split2.split*.bed > data/output/gfa2bin.split2.split.list")
        .status()?;

    if status.success() {
        println!("Command executed successfully.");
    } else {
        eprintln!("Command failed with status: {}", status);
    }

    println!("dashjdhas");
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("merge")
        .arg("--bed-list")
        .arg("data/output/gfa2bin.split2.split.list")
        .arg("-o")
        .arg("data/output/gfa2bin.merge");
    cmd.assert().success();

    Ok(())
}
