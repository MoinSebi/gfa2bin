use assert_cmd::prelude::*; // Add methods on commands
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// default is node
fn gfa_node() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("plink")
        .arg("-g")
        .arg("/home/svorbrugg_local/Rust/gSV/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/tt1");
    cmd.assert().success();
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn gfa_edges() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/tt2");
    cmd.assert().success();
    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t d (directed node)
fn gfa_dirnode() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("plink")
        .arg("-g")
        .arg("/home/svorbrugg_local/Rust/gSV/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/tt3")
        .arg("-t")
        .arg("d");
    cmd.assert().success();
    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t d (directed node)
fn gfa_dirnode2() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("plink")
        .arg("-g")
        .arg("data/example_data/testGraph.gfa")
        .arg("-o")
        .arg("data/output/tt5")
        .arg("-t")
        .arg("d")
        .arg("-f")
        .arg("data/example_data/fam_testGraph.fam");
    cmd.assert().success();
    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}
