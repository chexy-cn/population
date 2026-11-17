"""Record local dataset sizes without copying or compressing event data."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for directory in sorted(args.root.iterdir()):
        if not directory.is_dir() or not directory.name.startswith(("catalog_", "bbh_catalogs_")):
            continue
        files = [path for path in directory.rglob("*") if path.is_file()]
        catalogs = [path for path in files if path.suffix == ".csv" and path.stem[8:].isdigit()
                    and path.stem.startswith("catalog_")]
        sizes = [path.stat().st_size for path in catalogs]
        rows.append({
            "directory": directory.name, "file_count": len(files), "catalog_count": len(catalogs),
            "bytes": sum(path.stat().st_size for path in files),
            "mean_catalog_bytes": sum(sizes) / len(sizes) if sizes else None,
            "largest_catalog_bytes": max(sizes) if sizes else None,
        })
    args.output.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"recorded {len(rows)} local directories in {args.output}")


if __name__ == "__main__":
    main()
