use std::process::Command;

#[test]
fn demo_k12b_matches_upstream_fit() {
    let output = Command::new(env!("CARGO_BIN_EXE_admixrs"))
        .args(["-m", "K12b", "--sort"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = String::from_utf8(output.stdout).unwrap();
    let mut lines = report.trim().lines();
    assert_eq!(lines.next(), Some("K12b"));
    // Public demo estimates from the upstream Python implementation (tol=1e-6).
    let expected = [
        ("East Asian", 62.74),
        ("Southeast Asian", 33.38),
        ("Siberian", 3.62),
    ];
    for (line, (label, reference)) in lines.by_ref().take(3).zip(expected) {
        let (name, percent) = line.split_once(": ").unwrap();
        assert_eq!(name, label);
        let value: f64 = percent.strip_suffix('%').unwrap().parse().unwrap();
        assert!((value - reference).abs() < 0.03);
    }
    assert_eq!(lines.count(), 9);
}

#[test]
fn upstream_options_work_with_short_and_long_names() {
    let temp = tempfile::tempdir().unwrap();
    let data = temp.path();
    // A uniform frequency row has an exact and fast solution in either model.
    for (model, components) in [("K7b", 7), ("K12b", 12)] {
        std::fs::write(data.join(format!("{model}.alleles")), "rs1 A G\n").unwrap();
        std::fs::write(
            data.join(format!("{model}.{components}.F")),
            "0.5 ".repeat(components),
        )
        .unwrap();
    }
    let demo = "rs1\t1\t1\tAG\n";
    std::fs::write(data.join("demo_genome_23andme.txt"), demo).unwrap();
    std::fs::write(data.join("sample.txt"), demo).unwrap();
    for (options, filename) in [
        (
            vec![
                "-f",
                "sample.txt",
                "-v",
                "wegene",
                "-m",
                "K7b",
                "K12b",
                "-o",
                "-z",
                "-t",
                "1e-4",
            ],
            "admixture_results.txt",
        ),
        (
            vec![
                "--file=sample.txt",
                "--vendor",
                "wegene",
                "--models",
                "unknown",
                "--models",
                "K7b",
                "K12b",
                "--output",
                "report.txt",
                "--zhongwen",
                "--tolerance",
                "1e-4",
            ],
            "report.txt",
        ),
        (
            vec!["-f", "-m", "K7b", "K12b", "-t", "-z", "--output=demo.txt"],
            "demo.txt",
        ),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_admixrs"))
            .current_dir(data)
            .args(options)
            .args(["--sort", "--ignore-zeros", "--data-dir", "."])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(std::fs::read(data.join(filename)).unwrap(), result.stdout);
        let report = String::from_utf8(result.stdout).unwrap();
        assert!(report.starts_with("K7b\n"));
        assert!(report.contains("\n\nK12b\n"));
        assert!(report.contains("东亚"));
    }
}
