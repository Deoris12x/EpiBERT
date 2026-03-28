#!/usr/bin/env python3
"""
Fetch human nuclear proteome from UniProt and extract IDR regions using metapredict.

Usage:
    pip install metapredict requests
    python fetch_nuclear_proteome_idrs.py

Or run in Google Colab:
    1. Go to https://colab.research.google.com
    2. Create a new notebook
    3. Paste this entire script into a cell
    4. Run it — files will appear in /content/ and can be downloaded

Steps:
    1. Query UniProt REST API for reviewed human proteins with nuclear localization
    2. Download FASTA sequences (paginated to handle large datasets)
    3. Predict IDRs using metapredict
    4. Save results as TSV and FASTA files
"""

import requests
import os
import csv
import sys
import time

# Try to install metapredict if not present (useful for Colab)
try:
    import metapredict as meta
except ImportError:
    print("Installing metapredict...")
    os.system(f"{sys.executable} -m pip install metapredict")
    import metapredict as meta

# --- Configuration ---
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nuclear_proteome_idrs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

UNIPROT_QUERY = "(organism_id:9606) AND (reviewed:true) AND (cc_scl_term:SL-0191)"
SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
PAGE_SIZE = 500

# --- Step 1: Fetch nuclear proteome from UniProt (paginated) ---
print("=" * 60)
print("Step 1: Querying UniProt for human nuclear proteins...")
print(f"Query: {UNIPROT_QUERY}")
print("=" * 60)


def fetch_all_fasta():
    """Fetch all FASTA entries using UniProt pagination."""
    all_fasta = ""
    url = SEARCH_URL
    params = {
        "query": UNIPROT_QUERY,
        "format": "fasta",
        "size": PAGE_SIZE,
    }
    page = 1

    while url:
        print(f"  Fetching page {page}...")
        for attempt in range(4):
            try:
                resp = requests.get(url, params=params, timeout=120)
                resp.raise_for_status()
                break
            except requests.RequestException as e:
                wait = 2 ** (attempt + 1)
                print(f"  Retry {attempt+1}/4 after error: {e}. Waiting {wait}s...")
                time.sleep(wait)
                if attempt == 3:
                    raise

        all_fasta += resp.text

        # Check for next page via Link header
        link = resp.headers.get("Link", "")
        if 'rel="next"' in link:
            # Extract next URL from: <URL>; rel="next"
            url = link.split(";")[0].strip("<>")
            params = {}  # params are embedded in the next URL
        else:
            url = None

        page += 1

    return all_fasta


fasta_text = fetch_all_fasta()

# Save raw FASTA
fasta_path = os.path.join(OUTPUT_DIR, "human_nuclear_proteome.fasta")
with open(fasta_path, "w") as f:
    f.write(fasta_text)

# Parse FASTA
proteins = {}
current_id = None
current_header = None
current_seq = []

for line in fasta_text.strip().split("\n"):
    if line.startswith(">"):
        if current_id:
            proteins[current_id] = {
                "header": current_header,
                "sequence": "".join(current_seq)
            }
        parts = line.split("|")
        current_id = parts[1] if len(parts) >= 2 else line[1:].split()[0]
        current_header = line
        current_seq = []
    else:
        current_seq.append(line.strip())

if current_id:
    proteins[current_id] = {
        "header": current_header,
        "sequence": "".join(current_seq)
    }

print(f"\nDownloaded {len(proteins)} human nuclear proteins from UniProt (reviewed/Swiss-Prot)")

# --- Step 2: Predict IDRs using metapredict ---
print("\n" + "=" * 60)
print("Step 2: Predicting IDRs with metapredict...")
print("=" * 60)

idr_results = []
idr_fasta_entries = []
processed = 0
errors = 0

for uniprot_id, data in proteins.items():
    seq = data["sequence"]
    if len(seq) == 0:
        continue

    try:
        idrs = meta.predict_disorder_domains(seq)
        disordered_regions = idrs.disordered_domains

        for i, idr_seq in enumerate(disordered_regions):
            start = seq.find(idr_seq)
            end = start + len(idr_seq) if start >= 0 else -1

            idr_results.append({
                "uniprot_id": uniprot_id,
                "protein_length": len(seq),
                "idr_index": i + 1,
                "idr_start": start + 1 if start >= 0 else "N/A",
                "idr_end": end if start >= 0 else "N/A",
                "idr_length": len(idr_seq),
                "idr_sequence": idr_seq
            })

            idr_fasta_entries.append(
                f">{uniprot_id}_IDR{i+1} start={start+1} end={end} len={len(idr_seq)}\n{idr_seq}"
            )

    except Exception as e:
        errors += 1
        if errors <= 5:
            print(f"  Error processing {uniprot_id}: {e}")

    processed += 1
    if processed % 500 == 0:
        print(f"  Processed {processed}/{len(proteins)} proteins...")

print(f"\nDone! Processed {processed} proteins ({errors} errors)")
print(f"Found {len(idr_results)} IDR regions total")

# --- Step 3: Save results ---
print("\n" + "=" * 60)
print("Step 3: Saving results...")
print("=" * 60)

# Save TSV
tsv_path = os.path.join(OUTPUT_DIR, "nuclear_proteome_idrs.tsv")
with open(tsv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "uniprot_id", "protein_length", "idr_index",
        "idr_start", "idr_end", "idr_length", "idr_sequence"
    ], delimiter="\t")
    writer.writeheader()
    writer.writerows(idr_results)

# Save IDR FASTA
idr_fasta_path = os.path.join(OUTPUT_DIR, "nuclear_proteome_idrs.fasta")
with open(idr_fasta_path, "w") as f:
    f.write("\n".join(idr_fasta_entries))

# Save summary
summary_path = os.path.join(OUTPUT_DIR, "summary.txt")
total_proteins_with_idrs = len(set(r["uniprot_id"] for r in idr_results))
avg_idrs = len(idr_results) / total_proteins_with_idrs if total_proteins_with_idrs > 0 else 0
avg_idr_len = sum(r["idr_length"] for r in idr_results) / len(idr_results) if idr_results else 0

with open(summary_path, "w") as f:
    f.write("Nuclear Proteome IDR Analysis Summary\n")
    f.write("=" * 40 + "\n")
    f.write(f"Total nuclear proteins queried: {len(proteins)}\n")
    f.write(f"Proteins with IDRs: {total_proteins_with_idrs}\n")
    f.write(f"Total IDR regions found: {len(idr_results)}\n")
    f.write(f"Average IDRs per protein: {avg_idrs:.2f}\n")
    f.write(f"Average IDR length: {avg_idr_len:.1f} residues\n")
    f.write(f"Processing errors: {errors}\n")

print(f"\nFiles saved to: {OUTPUT_DIR}/")
print(f"  - human_nuclear_proteome.fasta  (all nuclear protein sequences)")
print(f"  - nuclear_proteome_idrs.tsv     (IDR regions table)")
print(f"  - nuclear_proteome_idrs.fasta   (IDR sequences in FASTA)")
print(f"  - summary.txt                   (analysis summary)")
print(f"\nSummary:")
print(f"  Nuclear proteins: {len(proteins)}")
print(f"  Proteins with IDRs: {total_proteins_with_idrs}")
print(f"  Total IDRs: {len(idr_results)}")
print(f"  Avg IDR length: {avg_idr_len:.1f} residues")
