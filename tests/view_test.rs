use assert_cmd::prelude::*; // Add methods on commands
use std::fs;
use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Convert plink to VCF
///
/// 4 Samples
fn view_test() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.view")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("view")
        .arg("-p")
        .arg("data/output/gfa2bin.view")
        .arg("-o")
        .arg("data/output/gfa2bin.view.vcf");
    cmd.assert().success();

    fs::remove_file("data/output/gfa2bin.view.bed")?;
    fs::remove_file("data/output/gfa2bin.view.bim")?;
    fs::remove_file("data/output/gfa2bin.view.fam")?;

    Ok(())
}
