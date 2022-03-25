use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command;

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("convert")
        .arg("-g")
        .arg("dasdsadasd");
    cmd.assert().stderr(predicate::str::contains("No file with such name"));

    Ok(())
}
