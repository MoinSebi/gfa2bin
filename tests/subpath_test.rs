use assert_cmd::prelude::*; // Add methods on commands
use std::fs;
use std::process::Command;

#[test]
/// Convert plink to VCF
///
/// 4 Samples
fn subpath_test() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("subpath")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.subpath");
    cmd.assert().success();

    fs::remove_file("data/output/gfa2bin.subpath.bed")?;
    fs::remove_file("data/output/gfa2bin.subpath.bim")?;
    fs::remove_file("data/output/gfa2bin.subpath.fam")?;

    Ok(())
}
