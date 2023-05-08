#!/usr/bin/env python3

import argparse
import yaml
import os
import collections

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check clang-tidy fixes statistics, deduplicating"
    )
    parser.add_argument("fixes_path", help="Path to fixes exported from clang-tidy")
    args = parser.parse_args()

    with open(args.fixes_path, "r") as file:
        fixes = yaml.safe_load(file)

    diagnostics = fixes.get("Diagnostics", [])
    stats = collections.defaultdict(set)

    for diagnostic in diagnostics:
        name = diagnostic["DiagnosticName"]
        location = "{}:{}".format(
            os.path.normpath(diagnostic["DiagnosticMessage"]["FilePath"]),
            diagnostic["DiagnosticMessage"]["FileOffset"],
        )
        stats[name].add(location)

    sorted_stats = sorted(stats.items(), key=lambda x: (-len(x[1]), x[0]))
    for name, locations in sorted_stats:
        print("{}: {}".format(name, len(locations)))
