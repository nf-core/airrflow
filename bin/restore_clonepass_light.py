#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


HEAVY_LOCI = {"IGH", "TRB", "TRD"}
LIGHT_LOCI = {"IGK", "IGL", "TRA", "TRG"}
CLONE_META_COLS = ["clone_id", "clone_size_count", "clone_size_freq", "mu_freq"]


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def infer_locus(row):
    if row.get("locus"):
        return row["locus"]
    calls = row.get("c_call") or row.get("v_call") or ""
    match = re.match(r"^((?:IG[HKL])|(?:TR[ABDG]))", calls)
    return match.group(1) if match else ""


parser = argparse.ArgumentParser()
parser.add_argument("--input-collapse", required=True)
parser.add_argument("--clone-pass", required=True)
args = parser.parse_args()

collapse_rows = read_tsv(args.input_collapse)
clone_rows = read_tsv(args.clone_pass)

if not collapse_rows or not clone_rows:
    raise SystemExit(0)
if "cell_id" not in collapse_rows[0] or "cell_id" not in clone_rows[0]:
    raise SystemExit(0)

for row in collapse_rows:
    row["locus"] = infer_locus(row)
for row in clone_rows:
    row["locus"] = infer_locus(row)

clone_fields = list(clone_rows[0].keys())
retained_cells = {row["cell_id"] for row in clone_rows}
existing_sequences = {row["sequence_id"] for row in clone_rows if row.get("sequence_id")}

heavy_meta = {}
heavy_rows_by_cell = {}
for row in clone_rows:
    if row["locus"] not in HEAVY_LOCI:
        continue
    heavy_rows_by_cell.setdefault(row["cell_id"], []).append(row)

for cell_id, rows in heavy_rows_by_cell.items():
    clone_ids = {row.get("clone_id", "") for row in rows if row.get("clone_id", "")}
    if len(clone_ids) == 1:
        heavy_meta[cell_id] = {col: rows[0].get(col, "") for col in CLONE_META_COLS if col in rows[0]}

restored_rows = []
for row in collapse_rows:
    if row["cell_id"] not in retained_cells:
        continue
    if row["locus"] not in LIGHT_LOCI:
        continue
    if row.get("sequence_id") in existing_sequences:
        continue
    restored = {field: "" for field in clone_fields}
    for field in clone_fields:
        if field in row:
            restored[field] = row[field]
    if row["cell_id"] in heavy_meta:
        for field, value in heavy_meta[row["cell_id"]].items():
            restored[field] = value
    restored_rows.append(restored)

if not restored_rows:
    raise SystemExit(0)

combined_rows = clone_rows + restored_rows
combined_rows.sort(key=lambda row: (row.get("cell_id", ""), row.get("sequence_id", "")))
write_tsv(args.clone_pass, combined_rows, clone_fields)
