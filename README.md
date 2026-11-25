# AdmixRS

AdmixRS estimates ancestry component proportions from consumer DNA genotypes on your own computer. It is a Rust port of Steven Liu’s Admix, with the same 27 active reference models, six input formats, and English and Chinese population labels. [Upstream Admix](https://github.com/stevenliuyi/admix)

## Installation

Download the `.tar.gz` asset for Windows x86-64, Linux x86-64, macOS Apple Silicon, or macOS Intel. Extract it with `tar -xzf <archive.tar.gz>`. Inside the extracted folder, the executable (`admixrs.exe` on Windows, `admixrs` elsewhere) reads the reference models from the adjacent `data/` directory. You can move the models to another location and pass `--data-dir DIR`. [GitHub Releases](https://github.com/pzweuj/AdmixRS/releases)

Run `./admixrs -m K12b` from the extracted directory, or `.\admixrs.exe -m K12b` in PowerShell. To use the commands below from any directory, add the extracted folder to your `PATH`.

To build from source, install stable Rust and run:

```sh
cargo build --release --locked
cargo run --release --locked -- -m K12b
```

## Using existing Admix commands

Replace `admix` with `admixrs` in an existing command. For example, to analyze a WeGene file with two models, display the largest components first, and use Chinese labels:

```sh
admixrs -f private/sample.txt -v wegene -m K7b K36 -z --sort --ignore-zeros
```

To run every model on a 23andMe file and save the report:

```sh
admixrs --file private/sample.txt --vendor 23andme --output admixture_results.txt
```

| Upstream option | Behavior |
| --- | --- |
| `-f [FILE]`, `--file [FILE]` | Read this genotype file. With no file, run the public 23andMe demo. |
| `-v VENDOR`, `--vendor VENDOR` | Accepts `23andme` (default), `ancestry`, `ftdna`, `ftdna2`, `wegene`, and `myheritage`. |
| `-m MODEL ...`, `--models MODEL ...` | Run one or more space-separated models in the supplied order. Omit to run all 27. Model names are case-sensitive. |
| `-o [FILE]`, `--output [FILE]` | Write a copy of the report. A bare `-o` writes `admixture_results.txt` in the current directory. |
| `-z`, `--zhongwen` | Display population names in Chinese. The default is English. |
| `-t [VALUE]`, `--tolerance [VALUE]` | Set the optimizer tolerance. The default is `1e-3`. |
| `--sort` | Sort components by decreasing proportion within each model. |
| `--ignore-zeros` | Hide components below 0.01%. |
| `-h`, `--help` | Show command help. |

AdmixRS also accepts `--list-models` to print the model names and `--data-dir DIR` to locate a separate reference directory. Reports go to standard output. Standard error shows the number of matched SNPs and the total reference SNP count for each model.

## Analysis and memory use

AdmixRS indexes the distinct marker IDs used by the selected models, reads the genotype file once, and retains calls only for those markers. Common numeric rsIDs use integer keys. Each model loads its matched frequency rows and releases them after fitting, so multiple models can run in sequence.

The likelihood uses the upstream allele orientation. Column three of each `.alleles` row corresponds to `F`, and column two corresponds to `1 − F`. The projected-gradient solver fits nonnegative weights that sum to one. Its tolerance bounds the remaining improvement in total log likelihood. The upstream SciPy solver uses a different stopping rule, so the last reported digits can vary.

Percentages describe the reference components of the selected model. Compare components within a model and use their labels in that context. [Model descriptions and sources](data/README.md)

## Privacy and licensing

Keep personal inputs in `private/` and saved reports in `results/`. Git ignores both directories and the default `admixture_results.txt` report. Tests use synthetic calls and the upstream public demo. Release archives include only the executable, documentation, license, public demo, and registered reference files.

AdmixRS code is licensed under GPL-3.0-only. Reference calculator files belong to their respective authors and retain their own terms. [GPLv3 license](LICENSE) · [Reference data credits](data/README.md)

## Publishing a release

Push a version tag such as `v0.1.0`, or choose **Actions → Build and release → Run workflow**. A manual run on a branch creates a release for that commit using the version in `Cargo.toml` prefixed with `v`. The optional tag field overrides that default. To repair an older release, select its tag as the workflow ref. Publishing a release through the GitHub UI also triggers the build.

Each runner tests the native executable and runs the public demo from its extracted archive. Once all four builds pass, the workflow attaches the archives to the corresponding GitHub Release. A rerun fills in missing assets.
