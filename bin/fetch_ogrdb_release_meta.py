#!/usr/bin/env python3
# Written by Ayelet Peres and released under the MIT license.

import html
import http.cookiejar
import json
import re
import sys
import urllib.parse
import urllib.request


API = "https://ogrdb.airr-community.org/api_v2"
SITE = "https://ogrdb.airr-community.org/germline_sets"


def die(msg):
    print(msg, file=sys.stderr)
    raise SystemExit(1)


def get_text(url, data=None, opener=None):
    req = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": "airrflow-fetch-ogrdb-local"},
    )
    open_fn = opener.open if opener else urllib.request.urlopen
    with open_fn(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


def get_json(url):
    return json.loads(get_text(url))


def norm_version(value):
    text = str(value)
    return text[:-2] if text.endswith(".0") else text


def species_id(species):
    data = get_json(f"{API}/germline/species")
    for item in data.get("species", []):
        if item.get("label") == species:
            return item["id"]
    die(f"OGRDB species not found: {species}")


def set_id(species, set_name):
    sid = species_id(species)
    data = get_json(f"{API}/germline/sets/{sid}")
    for item in data.get("germline_species", []):
        if item.get("germline_set_name") == set_name:
            return item["germline_set_id"]
    die(f"OGRDB germline set not found for {species}: {set_name}")


def latest_meta(germline_set_id):
    safe_id = urllib.parse.quote(germline_set_id, safe=".")
    data = get_json(f"{API}/germline/set/{safe_id}/latest")
    record = data["GermlineSet"][0]
    return str(record["release_version"]), record["release_date"][:10]


def doi_for(species, set_name, version, release_date):
    cj = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
    page = get_text(SITE, opener=opener)
    csrf = re.search(r'name="csrf_token" type="hidden" value="([^"]+)"', page)
    if not csrf:
        die("Failed to read OGRDB CSRF token")
    payload = urllib.parse.urlencode(
        {"csrf_token": csrf.group(1), "species": species, "submit": "Show Germline Sets"}
    ).encode()
    page = get_text(SITE, data=payload, opener=opener)
    for row in re.findall(r"<tr>.*?</tr>", page, re.S):
        cells = [
            html.unescape(re.sub(r"<[^>]+>", "", cell)).strip()
            for cell in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
        ]
        doi = re.search(r'<a href="https://doi.org/([^"]+)"', row)
        if not cells or cells[0] != set_name or not doi or release_date not in cells:
            continue
        row_version = norm_version(cells[cells.index(release_date) - 1])
        if row_version == version:
            return doi.group(1)
    die(f"Failed to resolve DOI for {species} {set_name} v{version}")


def parse_fasta(text):
    records = {}
    current = None
    seq = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records[current] = "".join(seq).upper()
            current = line[1:].split()[0]
            seq = []
        else:
            if current is None:
                die("Malformed FASTA from OGRDB")
            seq.append(line.strip())
    if current is not None:
        records[current] = "".join(seq).upper()
    return records


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")


def request_set(germline_set_id, version, species, fmt):
    safe_id = urllib.parse.quote(germline_set_id, safe=".")
    ex = "_ex" if species == "Homo sapiens" else ""
    if fmt == "gapped":
        path = f"{API}/germline/set/{safe_id}/{version}/gapped{ex}"
    else:
        path = f"{API}/germline/set/{safe_id}/{version}/ungapped{ex}"
    return get_text(path)


def download_split(species, set_name, version, prefix):
    gid = set_id(species, set_name)
    ungapped = parse_fasta(request_set(gid, version, species, "ungapped"))
    gapped = parse_fasta(request_set(gid, version, species, "gapped"))

    split_u = {
        "V": {k: v for k, v in ungapped.items() if len(k) > 3 and k[3] == "V"},
        "D": {k: v for k, v in ungapped.items() if len(k) > 3 and k[3] == "D" and len(v) < 100},
        "J": {k: v for k, v in ungapped.items() if len(k) > 3 and k[3] == "J"},
    }
    split_g_c = {
        k: v for k, v in gapped.items()
        if len(k) > 3 and (k[3] not in ["V", "D", "J"] or (k[3] == "D" and len(v) >= 100))
    }
    split_g_v = {k: v for k, v in gapped.items() if len(k) > 3 and k[3] == "V"}

    for suffix, records in split_u.items():
        if records:
            write_fasta(f"{prefix}{suffix}.fasta", records)
    if split_g_c:
        write_fasta(f"{prefix}C.fasta", split_g_c)
    if split_g_v:
        write_fasta(f"{prefix}V_gapped.fasta", split_g_v)


def main():
    if len(sys.argv) == 3:
        species, set_name = sys.argv[1], sys.argv[2]
        gid = set_id(species, set_name)
        version, release_date = latest_meta(gid)
        doi = doi_for(species, set_name, norm_version(version), release_date)
        match = re.search(r"zenodo\.(\d+)", doi)
        if not match:
            die(f"Failed to parse Zenodo record ID from DOI: {doi}")
        record_id = match.group(1)
        print("\t".join([version, release_date, doi, record_id, f"https://zenodo.org/records/{record_id}"]))
        return

    if len(sys.argv) == 6 and sys.argv[1] == "download":
        _, _, species, set_name, version, prefix = sys.argv
        download_split(species, set_name, version, prefix)
        return

    die(
        "usage:\n"
        "  fetch_ogrdb_release_meta.py 'human' 'IGH_VDJ'\n"
        "  fetch_ogrdb_release_meta.py download 'human' 'IGH_VDJ' 9.0 /tmp/prefix_"
    )


if __name__ == "__main__":
    main()
