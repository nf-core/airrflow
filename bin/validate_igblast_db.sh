#!/usr/bin/env bash
# Written by Ayelet Peres and released under the MIT license.

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR="igblast_base"
REFERENCE_FASTA_DIR=""
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    echo "Usage: $(basename "$0") -i <igblast_dir> -r <reference_fasta_dir> [-o <output_dir>]" >&2
}

while getopts "i:o:r:h" OPT; do
    case "$OPT" in
        i) INPUT_DIR=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        r) REFERENCE_FASTA_DIR=$OPTARG ;;
        h)
            usage
            exit 0
            ;;
        \?)
            usage
            exit 1
            ;;
    esac
done

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

canonical_stem() {
    local stem=$1

    if [[ $stem =~ ^(human|mouse)_(ig|tr)_([vdjc])$ ]]; then
        echo "${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}"
        return 0
    fi

    if [[ $stem =~ ^aa_(human|mouse)_(ig|tr)_v$ ]]; then
        echo "${BASH_REMATCH[0]}"
        return 0
    fi

    if [[ $stem =~ ^(imgt|airrc-imgt|airrc)_(human|mouse)_(ig|tr)_([vdjc])$ ]]; then
        echo "${BASH_REMATCH[2]}_${BASH_REMATCH[3]}_${BASH_REMATCH[4]}"
        return 0
    fi

    if [[ $stem =~ ^(imgt|airrc-imgt|airrc)_aa_(human|mouse)_(ig|tr)_v$ ]]; then
        echo "aa_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}_v"
        return 0
    fi

    return 1
}

rename_if_needed() {
    local source=$1
    local target=$2

    if [[ $source == "$target" ]]; then
        return 0
    fi

    if [[ -e $target ]]; then
        fail "Cannot white-label '$source' because '$target' already exists."
    fi

    mv "$source" "$target"
}

normalize_fasta_dir() {
    local fasta_dir=$1
    local path name stem canonical

    shopt -s nullglob
    for path in "$fasta_dir"/*.fasta; do
        name=$(basename "$path")
        stem=${name%.fasta}
        if canonical=$(canonical_stem "$stem"); then
            rename_if_needed "$path" "$fasta_dir/${canonical}.fasta"
        fi
    done
    shopt -u nullglob
}

clear_rebuildable_database_files() {
    local db_dir=$1
    local path name stem canonical

    shopt -s nullglob
    for path in "$db_dir"/*; do
        [[ -f $path ]] || continue
        name=$(basename "$path")
        stem=${name%%.*}
        if canonical=$(canonical_stem "$stem"); then
            rm -f "$path"
        fi
    done
    shopt -u nullglob
}

require_dir() {
    local path=$1
    [[ -d $path ]] || fail "Missing required directory '$path' in supplied IgBLAST cache."
}

compare_canonical_fasta_dirs() {
    local observed_dir=$1
    local expected_dir=$2
    local observed expected basename
    local observed_files=()
    local expected_files=()

    shopt -s nullglob
    for observed in "$observed_dir"/*.fasta; do
        basename=$(basename "$observed")
        [[ $basename =~ ^(human|mouse)_(ig|tr)_[vdjc]\.fasta$ || $basename =~ ^aa_(human|mouse)_(ig|tr)_v\.fasta$ ]] || continue
        observed_files+=("$basename")
    done
    for expected in "$expected_dir"/*.fasta; do
        basename=$(basename "$expected")
        [[ $basename =~ ^(human|mouse)_(ig|tr)_[vdjc]\.fasta$ || $basename =~ ^aa_(human|mouse)_(ig|tr)_v\.fasta$ ]] || continue
        expected_files+=("$basename")
    done
    shopt -u nullglob

    mapfile -t observed_files < <(printf '%s\n' "${observed_files[@]}" | sort -u)
    mapfile -t expected_files < <(printf '%s\n' "${expected_files[@]}" | sort -u)

    [[ ${#observed_files[@]} -eq ${#expected_files[@]} ]] || return 1
    for ((i = 0; i < ${#observed_files[@]}; i++)); do
        [[ ${observed_files[$i]} == "${expected_files[$i]}" ]] || return 1
        cmp -s "${observed_dir}/${observed_files[$i]}" "${expected_dir}/${expected_files[$i]}" || return 1
    done
}

reference_fastas_match_igblast_fastas() {
    local igblast_fasta_dir=$1
    local expected_root

    for database_type in imgt airrc-imgt; do
        expected_root=$(mktemp -d)
        "${SCRIPT_DIR}/ref2igblast.sh" -i "$REFERENCE_FASTA_DIR" -o "$expected_root" -d "$database_type" -f
        if compare_canonical_fasta_dirs "$igblast_fasta_dir" "$expected_root/fasta"; then
            rm -rf "$expected_root"
            return 0
        fi
        rm -rf "$expected_root"
    done

    return 1
}

rebuild_database_from_fasta_dir() {
    local fasta_dir=$1
    local temp_rebuild_input

    temp_rebuild_input=$(mktemp -d)
    cp "$fasta_dir"/*.fasta "$temp_rebuild_input"/
    "${SCRIPT_DIR}/ref2igblast.sh" -c "$temp_rebuild_input" -o "$OUTPUT_DIR"
    rm -rf "$temp_rebuild_input"
}

validate_canonical_assets() {
    local fasta_dir=$1
    local db_dir=$2
    local path name stem
    local -A fasta_stems=()
    local -A db_stems=()

    shopt -s nullglob

    for path in "$fasta_dir"/*.fasta; do
        name=$(basename "$path")
        stem=${name%.fasta}
        if [[ $stem =~ ^(human|mouse)_(ig|tr)_[vdjc]$ || $stem =~ ^aa_(human|mouse)_(ig|tr)_v$ ]]; then
            fasta_stems["$stem"]=1
        fi
    done

    for path in "$db_dir"/*; do
        [[ -f $path ]] || continue
        name=$(basename "$path")
        stem=${name%%.*}
        if [[ $stem =~ ^(human|mouse)_(ig|tr)_[vdjc]$ || $stem =~ ^aa_(human|mouse)_(ig|tr)_v$ ]]; then
            db_stems["$stem"]=1
        fi
    done

    shopt -u nullglob

    [[ ${#fasta_stems[@]} -gt 0 ]] || fail "No canonical FASTA assets were found in '$fasta_dir'."
    [[ ${#db_stems[@]} -gt 0 ]] || fail "No canonical BLAST database assets were found in '$db_dir'."

    for stem in "${!fasta_stems[@]}"; do
        [[ -n ${db_stems["$stem"]+x} ]] || fail "Missing database files for canonical stem '$stem'."
    done
}

[[ -n $INPUT_DIR ]] || fail "You must provide an input directory with -i."
[[ -d $INPUT_DIR ]] || fail "Input IgBLAST cache directory not found: $INPUT_DIR"
[[ -n $REFERENCE_FASTA_DIR ]] || fail "You must provide a reference FASTA directory with -r."
[[ -d $REFERENCE_FASTA_DIR ]] || fail "Reference FASTA directory not found: $REFERENCE_FASTA_DIR"
if [[ -e $OUTPUT_DIR ]] && [[ $(realpath "$INPUT_DIR") == $(realpath "$OUTPUT_DIR") ]]; then
    fail "Output directory must differ from input directory: $OUTPUT_DIR"
fi

REFERENCE_FASTA_DIR=$(realpath "$REFERENCE_FASTA_DIR")
mkdir -p "$OUTPUT_DIR"
cp -a "$INPUT_DIR"/. "$OUTPUT_DIR"/

require_dir "$OUTPUT_DIR/database"
require_dir "$OUTPUT_DIR/fasta"
require_dir "$OUTPUT_DIR/internal_data"
require_dir "$OUTPUT_DIR/optional_file"

normalize_fasta_dir "$OUTPUT_DIR/fasta"
reference_fastas_match_igblast_fastas "$OUTPUT_DIR/fasta" || \
    fail "Canonical FASTAs in reference_igblast/fasta do not match the FASTAs derivable from reference_fasta."
clear_rebuildable_database_files "$OUTPUT_DIR/database"
rebuild_database_from_fasta_dir "$OUTPUT_DIR/fasta"
validate_canonical_assets "$OUTPUT_DIR/fasta" "$OUTPUT_DIR/database"
