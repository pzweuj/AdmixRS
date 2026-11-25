"""Package a native build with the reference models enabled by the CLI."""

import argparse
from collections import Counter
from pathlib import Path
import tarfile
import zipfile


def release_files(binary):
    yield binary, binary.name
    for name in ("LICENSE", "README.md", "data/README.md", "data/demo_genome_23andme.txt"):
        yield Path(name), name

    rows = Path("src/populations.tsv").read_text(encoding="utf-8").splitlines()
    components = Counter(row.split("\t", 1)[0] for row in rows)
    for model, count in components.items():
        for name in (f"data/{model}.alleles", f"data/{model}.{count}.F"):
            yield Path(name), name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    args.archive.parent.mkdir(parents=True, exist_ok=True)
    if args.archive.suffix == ".zip":
        folder = args.archive.stem
        with zipfile.ZipFile(args.archive, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for source, name in release_files(args.binary):
                archive.write(source, f"{folder}/{name}")
    elif args.archive.name.endswith(".tar.gz"):
        folder = args.archive.name.removesuffix(".tar.gz")
        with tarfile.open(args.archive, "w:gz", compresslevel=6) as archive:
            for source, name in release_files(args.binary):
                archive.add(source, arcname=f"{folder}/{name}", recursive=False)
    else:
        parser.error("expected a .zip or .tar.gz archive")


if __name__ == "__main__":
    main()
