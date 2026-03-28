"""
Fetch human nuclear proteome from UniProt and extract IDR regions using metapredict.

Steps:
1. Query UniProt REST API for reviewed human proteins with nuclear localization
2. Download FASTA sequences
3. Predict IDRs using metapredict
4. Save results as TSV and FASTA files
"""

import requests
import time
import os
import csv
import metapredict as meta

OUTPUT_DIR = "/home/user/EpiBERT/nuclear_proteome_idrs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Step 1: Fetch nuclear proteome protein IDs from UniProt ---
print("=" * 60)
print("Step 1: Querying UniProt for human nuclear proteins...")
print("=" * 60)

# UniProt query: reviewed human proteins localized to nucleus
UNIPROT_QUERY = "(organism_id:9606) AND (reviewed:true) AND (cc_subcellular_location:nucleus)"
BASE_URL = "https://rest.uniprot.org/uniprotkb/stream"

params = {
    "query": UNIPROT_QUERY,
    "format": "fasta",
    "compressed": "false",
}

print(f"Query: {UNIPROT_QUERY}")
print("Downloading FASTA sequences (this may take a few minutes)...")

response = requests.get(BASE_URL, params=params, timeout=300)
response.raise_for_status()

fasta_text = response.text
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
        # Parse UniProt ID from header like >sp|P12345|NAME_HUMAN ...
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

print(f"Downloaded {len(proteins)} human nuclear proteins from UniProt (reviewed)")

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
        # Get IDR regions using metapredict
        idrs = meta.predict_disorder_domains(seq)
        disordered_regions = idrs.disordered_domains

        for i, idr_seq in enumerate(disordered_regions):
            # Find the start position of this IDR in the full sequence
            start = seq.find(idr_seq)
            end = start + len(idr_seq) if start >= 0 else -1

            idr_results.append({
                "uniprot_id": uniprot_id,
                "protein_length": len(seq),
                "idr_index": i + 1,
                "idr_start": start + 1 if start >= 0 else "N/A",  # 1-based
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

# Save summary stats
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

print(f"Saved: {fasta_path}")
print(f"Saved: {tsv_path}")
print(f"Saved: {idr_fasta_path}")
print(f"Saved: {summary_path}")
print(f"\nSummary:")
print(f"  Nuclear proteins: {len(proteins)}")
print(f"  Proteins with IDRs: {total_proteins_with_idrs}")
print(f"  Total IDRs: {len(idr_results)}")
print(f"  Avg IDR length: {avg_idr_len:.1f} residues")
