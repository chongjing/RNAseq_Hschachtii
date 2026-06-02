import csv, os

# Input TSV
INPUT = "fastq2_template_1760166804427.tsv"
# Output dirs
MANIFEST_DIR = "manifests"
CHECKSUM_FILE = "checksums.md5"

os.makedirs(MANIFEST_DIR, exist_ok=True)

checksums = []

with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        manifest = f"""STUDY\t{row['study']}
SAMPLE\t{row['sample_acc']}
NAME\t{row['library_name']}
INSTRUMENT\t{row['instrument']}
LIBRARY_NAME\t{row['library_name']}
LIBRARY_SOURCE\t{row['source']}
LIBRARY_SELECTION\t{row['selection']}
LIBRARY_STRATEGY\t{row['strategy']}
FASTQ\t{row['fq1']}
FASTQ\t{row['fq2']}
"""
        fname = os.path.join(MANIFEST_DIR, f"{row['library_name']}.manifest.txt")
        with open(fname, "w") as out:
            out.write(manifest)
        print(f"Wrote {fname}")

        # Collect checksums
        checksums.append(f"{row['fq1_md5']}  {row['fq1']}")
        checksums.append(f"{row['fq2_md5']}  {row['fq2']}")

# Write combined checksums file
with open(CHECKSUM_FILE, "w") as out:
    out.write("\n".join(checksums) + "\n")

print(f"Wrote {CHECKSUM_FILE}")
