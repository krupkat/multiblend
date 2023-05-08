#!/usr/bin/env python3

import argparse
import yaml
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Deduplicate clang-tidy diagnostics"
    )
    parser.add_argument("fixes_path", help="Path to fixes exported from clang-tidy")
    args = parser.parse_args()

    with open(args.fixes_path, "r") as file:
        fixes = yaml.safe_load(file)

    diagnostics = fixes.get("Diagnostics", [])
    dedup = {
        "Diagnostics": list(),
        "MainSourceFile": ""
    }

    seen_diagnostics = set()

    for diagnostic in diagnostics:
        diagnostic_hash = yaml.safe_dump(diagnostic)
        if diagnostic_hash not in seen_diagnostics:
            seen_diagnostics.add(diagnostic_hash)
            dedup["Diagnostics"].append(diagnostic)

    backup_filename = args.fixes_path + ".original"
    if os.path.exists(backup_filename):
        os.remove(backup_filename)

    os.rename(args.fixes_path, backup_filename)

    with open(args.fixes_path, "w") as output_file:
        yaml.safe_dump(dedup, output_file)
