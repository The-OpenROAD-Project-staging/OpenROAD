#!/usr/bin/env bash

## SPDX-License-Identifier: BSD-3-Clause
## Copyright (c) 2024-2026, The OpenROAD Authors

# Bazel-mode helper used by save_ok / save_defok / save_guideok.
# Each test's artifact lives under bazel-testlogs/<package>/<name>-<lang>_test/
# either as test.outputs/results/<name>-<lang>.<ext> (unzipped, current
# bazel default) or test.outputs/outputs.zip (zipped, older bazel).
#
# Usage: bazel_save.sh <dest_ext> <src_ext> <test_name>...
#   dest_ext   Extension to write next to the test (e.g. ok, defok, guideok).
#   src_ext    Artifact extension under results/ (log, def, guide).
#   test_name  One or more test stems. Targets are derived as
#              //<package>:<name>-tcl_test then //<package>:<name>-py_test.

set -e

if [ $# -lt 3 ]; then
    echo "usage: $0 <dest_ext> <src_ext> <test_name>..." >&2
    exit 2
fi

dest_ext=$1
src_ext=$2
shift 2

if ! command -v bazel >/dev/null 2>&1; then
    echo "bazel not on PATH; cannot extract .${dest_ext} from bazel-testlogs" >&2
    exit 1
fi

workspace=$(bazel info workspace 2>/dev/null || true)
testlogs=$(bazel info bazel-testlogs 2>/dev/null || true)
if [ -z "$workspace" ] || [ -z "$testlogs" ]; then
    echo "not inside a bazel workspace" >&2
    exit 1
fi

pkg=$(realpath --relative-to="$workspace" "$PWD")

for test_name in "$@"; do
    saved=0
    for lang_ext in tcl py; do
        target="//${pkg}:${test_name}-${lang_ext}_test"
        out_dir="${testlogs}/${pkg}/${test_name}-${lang_ext}_test/test.outputs"
        artifact="results/${test_name}-${lang_ext}.${src_ext}"
        # Run the test so the cache reflects the current code; failure
        # is expected (mismatched .ok is why we're here). Build errors
        # and target-not-found stay visible on stderr.
        bazel test "$target" --test_summary=terse >/dev/null || true
        if [ -f "${out_dir}/${artifact}" ]; then
            cp "${out_dir}/${artifact}" "${test_name}.${dest_ext}"
            echo "${test_name}"
            saved=1
            break
        fi
        zip="${out_dir}/outputs.zip"
        if [ -f "$zip" ]; then
            tmp="${test_name}.${dest_ext}.tmp"
            if unzip -p "$zip" "$artifact" > "$tmp" 2>/dev/null \
                    && [ -s "$tmp" ]; then
                mv "$tmp" "${test_name}.${dest_ext}"
                echo "${test_name}"
                saved=1
                break
            fi
            rm -f "$tmp"
        fi
    done
    if [ "$saved" -eq 0 ]; then
        echo "\"${test_name}\" ${src_ext} file not found in bazel-testlogs" >&2
    fi
done
