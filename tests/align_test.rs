use assert_cmd::prelude::*; // Add methods on commands
use std::fs;
use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g
/// -f node
/// --pansn #
///
/// 4 Samples
fn align_v1() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("align")
        .arg("-c")
        .arg("data/example_data/packs/9986.1k.a2.copy.pt")
        .arg("-i")
        .arg("data/example_data/packs/9986.1k.pi")
        .arg("-o")
        .arg("data/output/gfa2bin.align.pt");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.align.pt.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();

    // Buffer should be 8 samples + header
    assert_eq!(buffer.len(), 1002);

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g
/// -f node
/// --pansn #
///
/// 4 Samples
fn align_pn() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("align")
        .arg("-c")
        .arg("data/example_data/packs/9986.1k.pn")
        .arg("-i")
        .arg("data/example_data/packs/9986.1k.pi")
        .arg("-o")
        .arg("data/output/gfa2bin.align.pn");
    cmd.assert().success();
    let mut b = File::open("data/output/gfa2bin.align.pn.bed").unwrap();

    // Read the buffer
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();

    // Buffer should be 8 samples + header
    assert_eq!(buffer.len(), 1002);

    // First "real" byte is 00
    //assert_eq!(buffer[3], 63);
    // Second "real" byte is 000000000
    //assert_eq!(buffer[4], 0);

    Ok(())
}
