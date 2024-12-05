use assert_cmd::prelude::*;
 // Add methods on commands



use std::process::Command;

// cargo run -- graph -g data/example_data/testGraph.gfa -f node -o data/example_data/node.remove

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn filter_test() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("./data/example_data/gfa/testGraph.gfa")
        .arg("-o")
        .arg("./data/output/gfa2bin.filter")
        .arg("-f")
        .arg("node")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("filter")
        .arg("-p")
        .arg("./data/output/gfa2bin.filter")
        .arg("-m")
        .arg("0.5")
        .arg("-o")
        .arg("data/output/gfa2bin.filter.0.1");
    cmd.assert().success();

    //fs::remove_file("data/output/remove.feature.node.fam")?;
    Ok(())
}
