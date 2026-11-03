#!/usr/bin/env bash
# Written by Ayelet Peres and released under the MIT license.

set -euo pipefail

OUTDIR="."
DATABASE_TYPE="imgt"
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OGRDB_PY="${HERE}/fetch_ogrdb_release_meta.py"
TMPDIR=$(mktemp -d)
DATE=$(date +"%Y.%m.%d")

usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -d  Database type: imgt or airrc-imgt. Defaults to imgt."
    echo "  -h  This message."
}

cleanup() {
    rm -rf "$TMPDIR"
}

validate_fasta() {
    local file=$1
    if [ ! -s "$file" ]; then
        echo "Error: FASTA file ${file} is empty" >&2
        exit 1
    fi
    awk '
        /^>/ {if (seen && !seq) exit 1; seen=1; seq=0; next}
        NF {if (!seen) exit 1; seq=1}
        END {exit !(seen && seq)}
    ' "$file" || {
        echo "Error: FASTA file ${file} is malformed" >&2
        exit 1
    }
}

species_url() {
    case "$1" in
        human) echo "Homo%20sapiens" ;;
        mouse) echo "Mus%20musculus" ;;
    esac
}

species_label() {
    case "$1" in
        human) echo "Homo sapiens" ;;
        mouse) echo "Mus musculus" ;;
    esac
}

species_replace() {
    case "$1" in
        human) echo "Homo sapiens" ;;
        mouse) echo "Mus musculus" ;;
    esac
}

fetch_imgt_chain() {
    local species_key=$1 chain=$2 query=$3 subdir=$4 prefix=$5
    local species="$(species_label "$species_key")"
    local url_species="$(species_url "$species_key")"
    local tmp="${TMPDIR}/${species_key}_${chain}.html"
    local dest="${OUTDIR}/${species_key}/${subdir}/${prefix}_${species_key}_${chain}.fasta"
    wget -q "https://www.imgt.org/genedb/GENElect?query=${query}+${chain}&species=${url_species}" -O "$tmp"
    awk '/<pre>/{i++}/<\/pre>/{j++} j==2{exit} i==2 && j==1 && $0 !~ /^<pre>$/ {print}' "$tmp" > "$dest"
    sed -i.bak "s/${species}/${species// /_}/g" "$dest"
    rm -f "$tmp" "${dest}.bak"
    validate_fasta "$dest"
}

fetch_imgt_aa_chain() {
    local species_key=$1 chain=$2
    local species="$(species_label "$species_key")"
    local url_species="$(species_url "$species_key")"
    local tmp="${TMPDIR}/${species_key}_${chain}_aa.html"
    local dest="${OUTDIR}/${species_key}/vdj_aa/imgt_aa_${species_key}_${chain}.fasta"
    wget -q "https://www.imgt.org/genedb/GENElect?query=7.3+${chain}&species=${url_species}" -O "$tmp"
    awk '/<pre>/{i++}/<\/pre>/{j++} j==2{exit} i==2 && j==1 && $0 !~ /^<pre>$/ {print}' "$tmp" > "$dest"
    sed -i.bak "s/${species}/${species// /_}/g" "$dest"
    rm -f "$tmp" "${dest}.bak"
    validate_fasta "$dest"
}

fetch_imgt_leader() {
    local species_key=$1 chain=$2
    local species="$(species_label "$species_key")"
    local url_species="$(species_url "$species_key")"
    local tmp="${TMPDIR}/${species_key}_${chain}_leader.html"
    local dest="${OUTDIR}/${species_key}/leader/imgt_${species_key}_${chain}L.fasta"
    wget -q "https://www.imgt.org/genedb/GENElect?query=8.1+${chain}V&species=${url_species}&IMGTlabel=L-PART1+L-PART2" -O "$tmp"
    awk '/<pre>/{i++}/<\/pre>/{j++} j==2{exit} i==2 && j==1 && $0 !~ /^<pre>$/ {print}' "$tmp" > "$dest"
    sed -i.bak "s/${species}/${species// /_}/g" "$dest"
    rm -f "$tmp" "${dest}.bak"
    validate_fasta "$dest"
}

append_airrc_yaml() {
    local species=$1 set_name=$2 version=$3 release_date=$4 doi=$5 zenodo_id=$6 zenodo_url=$7
    cat >> "${OUTDIR}/AIRRC.yaml" <<EOF
  - species: ${species}
    set: ${set_name}
    version: ${version}
    release_date: ${release_date}
    doi: ${doi}
    zenodo_record_id: ${zenodo_id}
    zenodo_url: ${zenodo_url}
EOF
}

fetch_ogrdb_set() {
    local species_key=$1 locus=$2 set_name=$3 subdir=$4
    shift 4
    local species="$(species_label "$species_key")"
    local prefix="${TMPDIR}/$(echo "${species_key}_${set_name}" | tr -cs 'A-Za-z0-9' '_')_"
    local version release_date doi zenodo_id zenodo_url chain seg src dest
    IFS=$'\t' read -r version release_date doi zenodo_id zenodo_url < <("$OGRDB_PY" "$species" "$set_name")
    "$OGRDB_PY" download "$species" "$set_name" "$version" "$prefix"
    for chain in "$@"; do
        seg=${chain: -1}
        if [ "$seg" = "V" ] && [[ "$chain" =~ ^IG[HKL]V$ ]]; then
            src="${prefix}V_gapped.fasta"
        else
            src="${prefix}${seg}.fasta"
        fi
        dest="${OUTDIR}/${species_key}/${subdir}/airrc_${species_key}_${chain}.fasta"
        validate_fasta "$src"
        mv "$src" "$dest"
    done
    append_airrc_yaml "$species_key" "$set_name" "$version" "$release_date" "$doi" "$zenodo_id" "$zenodo_url"
}

write_imgt_yaml() {
    cat > "${OUTDIR}/IMGT.yaml" <<EOF
source:  https://www.imgt.org/genedb
date:    ${DATE}
species:
    - human:Homo+sapiens
    - mouse:Mus
EOF
}

write_airrc_yaml() {
    cat > "${OUTDIR}/AIRRC.yaml" <<EOF
source: https://ogrdb.airr-community.org
date: ${DATE}
sets:
EOF
}

fetch_imgt_database() {
    local species_key chain query
    for species_key in human mouse; do
        mkdir -p "${OUTDIR}/${species_key}/vdj" "${OUTDIR}/${species_key}/vdj_aa" "${OUTDIR}/${species_key}/leader" "${OUTDIR}/${species_key}/constant"
        for chain in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ TRAV TRAJ TRBV TRBD TRBJ TRDV TRDD TRDJ TRGV TRGJ; do
            fetch_imgt_chain "$species_key" "$chain" 7.14 vdj imgt
        done
        for chain in IGHV IGKV IGLV TRAV TRBV TRDV TRGV; do
            fetch_imgt_aa_chain "$species_key" "$chain"
        done
        for chain in IGH IGK IGL TRA TRB TRG TRD; do
            fetch_imgt_leader "$species_key" "$chain"
        done
        for chain in IGHC IGKC IGLC TRAC TRBC TRGC TRDC; do
            query=14.1
            if [ "$species_key" = "mouse" ] && { [ "$chain" = "IGKC" ] || [ "$chain" = "IGLC" ]; }; then
                query=7.5
            fi
            fetch_imgt_chain "$species_key" "$chain" "$query" constant imgt
        done
    done
    write_imgt_yaml
}

fetch_airrc_imgt_database() {
    local species_key chain query
    mkdir -p \
        "$OUTDIR/human/vdj" "$OUTDIR/human/constant" \
        "$OUTDIR/mouse/vdj" "$OUTDIR/mouse/constant"

    write_imgt_yaml
    write_airrc_yaml

    fetch_ogrdb_set human IGH "IGH_VDJ" vdj IGHV IGHD IGHJ
    fetch_ogrdb_set human IGK "IGKappa_VJ" vdj IGKV IGKJ
    fetch_ogrdb_set human IGL "IGLambda_VJ" vdj IGLV IGLJ
    fetch_ogrdb_set human IGH "IGHC" constant IGHC

    fetch_ogrdb_set mouse IGH "C57BL/6 IGH" vdj IGHV IGHD IGHJ
    fetch_ogrdb_set mouse IGK "C57BL/6J IGKV" vdj IGKV
    fetch_ogrdb_set mouse IGK "IGKJ (all strains)" vdj IGKJ
    fetch_ogrdb_set mouse IGL "C57BL/6J IGLV" vdj IGLV
    fetch_ogrdb_set mouse IGL "IGLJ (all strains)" vdj IGLJ

    for species_key in human mouse; do
        for chain in TRAV TRAJ TRBV TRBD TRBJ TRDV TRDD TRDJ TRGV TRGJ; do
            fetch_imgt_chain "$species_key" "$chain" 7.14 vdj imgt
        done
        if [ "$species_key" = "human" ]; then
            for chain in IGKC IGLC; do
                fetch_imgt_chain "$species_key" "$chain" 14.1 constant imgt
            done
        else
            fetch_imgt_chain "$species_key" IGHC 14.1 constant imgt
            fetch_imgt_chain "$species_key" IGKC 7.5 constant imgt
            fetch_imgt_chain "$species_key" IGLC 7.5 constant imgt
        fi
        for chain in TRAC TRBC TRGC TRDC; do
            fetch_imgt_chain "$species_key" "$chain" 14.1 constant imgt
        done
    done
}

while getopts ":o:d:h" opt; do
    case "$opt" in
        o) OUTDIR=$OPTARG ;;
        d) DATABASE_TYPE=$OPTARG ;;
        h) usage; exit 0 ;;
        :) echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

trap cleanup EXIT

case "$DATABASE_TYPE" in
    imgt) fetch_imgt_database ;;
    airrc-imgt) fetch_airrc_imgt_database ;;
    *) echo "Unknown database type: $DATABASE_TYPE" >&2; exit 1 ;;
esac
