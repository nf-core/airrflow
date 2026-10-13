#!/usr/bin/env bash
# Written by Gisela Gabernet and released under the MIT license (2020).
# Reworked to fetch and build with sourcerer instead of the fetch/ref2igblast
# shell scripts. sourcerer's source names match the database types (imgt,
# airrc-imgt), so the type is passed through directly.
set -euo pipefail

DATABASE_TYPE=imgt

while getopts "d:" OPT; do
    case "$OPT" in
    d)  DATABASE_TYPE=$OPTARG
        ;;
    esac
done

echo "Fetching databases with immvoke (${DATABASE_TYPE})..."

# 'all' is every species immvoke supports for this source, in one call, with
# provenance merged across them.
immvoke "${DATABASE_TYPE}" download all --outdir immvoke_out
mv immvoke_out/reference_base reference_base

# Build the IgBLAST databases from the reference tree: cleans the FASTAs, runs
# makeblastdb, and mirrors the NCBI internal_data / optional_file trees.
immvoke reference build reference_base --out igblast_base

# A fetch that silently returns less than the source holds would only surface as a
# bad annotation, so check the databases the source is expected to carry. OGRDB
# curates immunoglobulin sets only, so TR is not expected there.
expected="human_ig_v human_ig_d human_ig_j"
if [ "${DATABASE_TYPE}" != "ogrdb" ]; then
    expected="${expected} human_tr_v human_tr_d human_tr_j"
fi
missing=""
for db in ${expected}; do
    [ -f "igblast_base/database/${db}.ndb" ] || missing="${missing} ${db}"
done
if [ -n "${missing}" ]; then
    echo "The ${DATABASE_TYPE} fetch built no database for:${missing}" >&2
    exit 1
fi

echo "FetchDBs process finished."
