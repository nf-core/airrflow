#!/usr/bin/env bash
# Written by Ayelet Peres and released under the MIT license.

set -euo pipefail

OUTDIR="."
DATABASE_TYPE="imgt"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FASTA_ONLY=false
REFERENCE_DIR=""
CANONICAL_FASTA_DIR=""

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo "  -i  Input directory containing fetched references."
    echo "  -c  Input directory containing canonical white-labeled FASTAs."
    echo "  -o  Output directory for the built database."
    echo "  -d  Database type: imgt or airrc-imgt. Defaults to imgt."
    echo "  -f  Build canonical FASTA output only."
    echo "  -h  This message."
}

while getopts "i:c:o:d:fh" OPT; do
    case "$OPT" in
        i) REFERENCE_DIR=$OPTARG ;;
        c) CANONICAL_FASTA_DIR=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        d) DATABASE_TYPE=$OPTARG ;;
        f) FASTA_ONLY=true ;;
        h) usage; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
    esac
done

if [[ -n $REFERENCE_DIR && -n $CANONICAL_FASTA_DIR ]]; then
    fail "Use either -i or -c, not both."
fi

if [[ -z $REFERENCE_DIR && -z $CANONICAL_FASTA_DIR ]]; then
    fail "You must specify an input directory using -i or -c."
fi

mkdir -p "${OUTDIR}"
OUTDIR=$(realpath "${OUTDIR}")
mkdir -p "${OUTDIR}/fasta" "${OUTDIR}/database"

concat_to_file() {
    local dest=$1
    local pattern file
    local files=()
    shift
    shopt -s nullglob
    for pattern in "$@"; do
        for file in $pattern; do
            files+=("$file")
        done
    done
    shopt -u nullglob
    if [ ${#files[@]} -gt 0 ]; then
        cat "${files[@]}" > "$dest"
    else
        : > "$dest"
    fi
}

build_matching_files() {
    local input_dir=$1
    local pattern=$2
    local dbtype=$3
    local file basename
    local built_any=false

    shopt -s nullglob
    for file in "${input_dir}"/*.fasta; do
        basename=$(basename "$file")
        [[ $basename =~ $pattern ]] || continue
        built_any=true
        "${SCRIPT_DIR}/clean_imgtdb.py" "$file" "${OUTDIR}/fasta/${basename}"
        if [[ $FASTA_ONLY == false ]]; then
            makeblastdb -parse_seqids -dbtype "$dbtype" -in "${OUTDIR}/fasta/${basename}" \
                -out "${OUTDIR}/database/${basename%.fasta}"
        fi
    done
    shopt -u nullglob

    $built_any
}

build_from_canonical_fastas() {
    local input_dir=$1

    [[ -d $input_dir ]] || fail "Canonical FASTA directory not found: $input_dir"
    input_dir=$(realpath "$input_dir")

    build_matching_files "$input_dir" '^(human|mouse)_(ig|tr)_[vdjc]\.fasta$' nucl || \
        fail "No canonical nucleotide FASTA files were found in '$input_dir'."
    build_matching_files "$input_dir" '^aa_(human|mouse)_(ig|tr)_v\.fasta$' prot || true
}

prepare_canonical_fastas_from_reference_tree() {
    local germdir=$1
    local tmpdir

    [[ -d $germdir ]] || fail "Reference directory not found: $germdir"
    germdir=$(realpath "$germdir")
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' RETURN

    case "$DATABASE_TYPE" in
        imgt)
            BUILD_AA=true
            ;;
        airrc-imgt)
            BUILD_AA=false
            ;;
        *)
            fail "Unknown database type: $DATABASE_TYPE"
            ;;
    esac

    for SPECIES in human mouse; do
        for CHAIN in IG TR; do
            echo "|---- ${SPECIES} ${CHAIN} ----|"
            for SEGMENT in V D J; do
                echo "|---- ${SEGMENT} ----|"
                F=$(echo "${SPECIES}_${CHAIN}_${SEGMENT}.fasta" | tr '[:upper:]' '[:lower:]')
                if [[ $DATABASE_TYPE == "imgt" ]]; then
                    concat_to_file "${tmpdir}/${F}" "${germdir}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta"
                else
                    concat_to_file "${tmpdir}/${F}" \
                        "${germdir}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta" \
                        "${germdir}/${SPECIES}/vdj/airrc_${SPECIES}_${CHAIN}?${SEGMENT}.fasta"
                fi
            done

            F=$(echo "${SPECIES}_${CHAIN}_c.fasta" | tr '[:upper:]' '[:lower:]')
            echo "|---- C ----|"
            if [[ $DATABASE_TYPE == "imgt" ]]; then
                concat_to_file "${tmpdir}/${F}" "${germdir}/${SPECIES}/constant/imgt_${SPECIES}_${CHAIN}?C.fasta"
            else
                concat_to_file "${tmpdir}/${F}" \
                    "${germdir}/${SPECIES}/constant/imgt_${SPECIES}_${CHAIN}?C.fasta" \
                    "${germdir}/${SPECIES}/constant/airrc_${SPECIES}_${CHAIN}?C.fasta"
            fi

            if [[ $BUILD_AA == true ]]; then
                F=$(echo "aa_${SPECIES}_${CHAIN}_v.fasta" | tr '[:upper:]' '[:lower:]')
                echo "|---- AA V ----|"
                concat_to_file "${tmpdir}/${F}" "${germdir}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?V.fasta"
            fi
        done
    done

    build_from_canonical_fastas "$tmpdir"
}

if [[ -n $CANONICAL_FASTA_DIR ]]; then
    build_from_canonical_fastas "$CANONICAL_FASTA_DIR"
else
    prepare_canonical_fastas_from_reference_tree "$REFERENCE_DIR"
fi
