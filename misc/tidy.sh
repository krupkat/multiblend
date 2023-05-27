#!/usr/bin/env bash

if [ -z "$1" ]; then
    target="build/Release"
else
    target="$1"
fi

if [ -z "$2" ]; then
    result="fixes.yaml"
else
    result="$2"
fi

run-clang-tidy `pwd`/src `pwd`/tests "-p=$target" -quiet "-export-fixes=$result"
