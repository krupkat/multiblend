#!/usr/bin/env bash

if [ -z "$1" ]; then
    target="build/Release"
else
    target="$1"
fi

run-clang-tidy -p "$target" -export-fixes fixes.yaml
./misc/python/clang_tidy_deduplicate.py fixes.yaml
clang-apply-replacements .

./misc/format.sh
