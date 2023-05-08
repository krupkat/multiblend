#!/usr/bin/env bash

run-clang-tidy -p build/Release -export-fixes fixes.yaml
./misc/python/clang_tidy_deduplicate.py fixes.yaml
clang-apply-replacements .

./misc/format.sh
