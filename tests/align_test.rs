use assert_cmd::prelude::*;
use std::fs; // Add methods on commands

use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for cov subcommand
///
/// Input: pt
/// Number of samples: 2
/// Output: bed
fn cov_v1() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("cov")
        .arg("-c")
        .arg("data/example_data/packs/9986.1k.a1.a2.merge.pt")
        .arg("-i")
        .arg("data/example_data/packs/9986.1k.pi")
        .arg("-o")
        .arg("data/output/gfa2bin.cov.merge.pt");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.cov.merge.pt.bed").unwrap();
    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 67);


    let content = fs::read_to_string("data/output/gfa2bin.cov.merge.pt.bim")?;
    assert_eq!(content.lines().count(), 67);

    let content = fs::read_to_string("data/output/gfa2bin.cov.merge.pt.fam")?;
    assert_eq!(content.lines().count(), 2);
    // Buffer should be 8 samples + header
    //fs::remove_file("data/output/gfa2bin.cov.merge.pt.bed")?;
    //fs::remove_file("data/output/gfa2bin.cov.merge.pt.bim")?;
    //fs::remove_file("data/output/gfa2bin.cov.merge.pt.fam")?;

    Ok(())
}

#[test]
/// Test cov pn
fn cov_pn() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("cov")
        .arg("-c")
        .arg("data/example_data/packs/9986.1k.seq.merge.pn")
        .arg("-i")
        .arg("data/example_data/packs/9986.1k.pi")
        .arg("-o")
        .arg("data/output/gfa2bin.cov.merge.pn");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.cov.merge.pn.bed").unwrap();
    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 999);


    let content = fs::read_to_string("data/output/gfa2bin.cov.merge.pn.bim")?;
    assert_eq!(content.lines().count(), 999);

    let content = fs::read_to_string("data/output/gfa2bin.cov.merge.pn.fam")?;
    assert_eq!(content.lines().count(), 2);
    // Buffer should be 8 samples + header
    Ok(())
}
