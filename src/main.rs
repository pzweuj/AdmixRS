mod fit;
mod genome;
mod model;

use clap::Parser;
use genome::Vendor;
use std::error::Error;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    version,
    args_override_self = true,
    infer_long_args = true,
    about = "Estimate admixture proportions from SNP raw data"
)]
struct Args {
    /// Genotype file (default: the upstream 23andMe demo)
    #[arg(short, long, num_args = 0..=1)]
    file: Option<Option<PathBuf>>,
    /// Input file format
    #[arg(short, long, value_enum, default_value = "23andme")]
    vendor: Vendor,
    /// Models in report order (default: all)
    #[arg(short, long, num_args = 1.., action = clap::ArgAction::Set)]
    models: Vec<String>,
    /// Save a copy of the report
    #[arg(short, long, num_args = 0..=1, default_missing_value = "admixture_results.txt")]
    output: Option<PathBuf>,
    /// Display Chinese population labels
    #[arg(short, long)]
    zhongwen: bool,
    /// Absolute tolerance on the total log-likelihood optimality bound
    #[arg(short, long, num_args = 0..=1, default_value_t = 1e-3, default_missing_value = "1e-3")]
    tolerance: f64,
    /// Sort by decreasing proportion
    #[arg(long)]
    sort: bool,
    /// Hide proportions below 0.01%
    #[arg(long)]
    ignore_zeros: bool,
    /// Directory containing the allele and frequency files
    #[arg(long)]
    data_dir: Option<PathBuf>,
    /// Print the available models and exit
    #[arg(long)]
    list_models: bool,
}

fn default_data() -> PathBuf {
    if let Ok(exe) = std::env::current_exe()
        && let Some(dir) = exe.parent().map(|parent| parent.join("data"))
        && dir.is_dir()
    {
        return dir;
    }
    Path::new(env!("CARGO_MANIFEST_DIR")).join("data")
}

fn run(args: Args) -> Result<(), Box<dyn Error>> {
    let catalog = model::catalog();
    if args.list_models {
        for model in catalog {
            println!("{}", model.name);
        }
        return Ok(());
    }
    if !args.tolerance.is_finite() || args.tolerance <= 0.0 {
        return Err("tolerance must be finite and positive".into());
    }
    let models: Vec<_> = if args.models.is_empty() {
        catalog.iter().collect()
    } else {
        args.models
            .iter()
            .map(|name| {
                catalog
                    .iter()
                    .find(|model| model.name == name)
                    .ok_or_else(|| format!("unknown model '{name}'; see --list-models"))
            })
            .collect::<Result<_, _>>()?
    };
    let data = args.data_dir.unwrap_or_else(default_data);
    let (input, vendor) = match args.file.flatten().filter(|p| !p.as_os_str().is_empty()) {
        Some(path) => (path, args.vendor),
        None => {
            eprintln!("Using the bundled 23andMe demo.");
            (
                data.join("demo_genome_23andme.txt"),
                Vendor::TwentyThreeAndMe,
            )
        }
    };
    let file = File::open(&input).map_err(|e| format!("{}: {e}", input.display()))?;
    if let Some(path) = &args.output
        && path.exists()
        && path.canonicalize()? == input.canonicalize()?
    {
        return Err("output path points to the genotype input".into());
    }
    let mut genome = genome::Genome::default();
    for model in &models {
        model.index_markers(&data, &mut genome)?;
    }
    genome::read(BufReader::new(file), vendor, &mut genome)?;
    let mut report = String::new();
    for model in models {
        let obs = model
            .load(&data, &genome)
            .map_err(|e| format!("{}: {e}", model.name))?;
        let q = fit::fit(&obs, args.tolerance).map_err(|e| format!("{}: {e}", model.name))?;
        eprintln!(
            "{}: {} / {} model SNPs matched",
            model.name,
            obs.counts.len(),
            obs.total_markers
        );
        writeln!(report, "{}", model.name)?;
        let mut rows: Vec<_> = model.labels.iter().zip(q).collect();
        if args.sort {
            rows.sort_by(|a, b| b.1.total_cmp(&a.1));
        }
        for (&(en, zh), fraction) in rows {
            if args.ignore_zeros && fraction < 1e-4 {
                continue;
            }
            let label = if args.zhongwen { zh } else { en };
            writeln!(report, "{label}: {:.2}%", 100.0 * fraction)?;
        }
        writeln!(report)?;
    }
    if let Some(path) = args.output {
        let dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut tmp = tempfile::NamedTempFile::new_in(dir)?;
        tmp.write_all(report.as_bytes())?;
        tmp.persist(&path)?;
    }
    io::stdout().lock().write_all(report.as_bytes())?;
    Ok(())
}

fn main() -> ExitCode {
    match run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("admixrs: {error}");
            ExitCode::FAILURE
        }
    }
}
