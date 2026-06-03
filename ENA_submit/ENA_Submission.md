

## 6. ENA EMBL submission
6.1 install webin-cli https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html
```bash
cd /home/cx264/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli
wget https://github.com/enasequence/webin-cli/releases/download/9.0.1/webin-cli-9.0.1.jar
``
6.2 login Webin submission Portal to register study and samples, via https://www.ebi.ac.uk/ena/submit/webin/login

| Type       | Accession    | Unique name (alias)                  |
|------------|--------------|--------------------------------------|
| Project    | PRJEB114114  | 1b5d2a12-04d3-427d-951d-d065c9de58b7 |
| Submission | ERA36352013  | SUBMISSION-02-06-2026-11:08:22:314   |

6.3 Submission
After get sample accessions, prepare metadata table with header like

```bash
cd /home/cx264/rds/rds-scrna_spatial-6qULnBz5AIM/Chongjing_Xia/06.Sebastian/07.X204SC25060928-Z01-F001/EBI_ENA_submit
#prepare a sample file fastq2_template_1760166804427.tsv:
| sample_id   | study      | sample_acc  | instrument          | library_name | source         | selection   | strategy | layout | fq1             | fq1_md5                          | fq2             | fq2_md5                          |
|-------------|------------|-------------|---------------------|--------------|----------------|-------------|----------|--------|-----------------|----------------------------------|-----------------|----------------------------------|
| ERS30371088  | PRJEB114114   | ERS30371088  | Illumina HiSeq 2500   | D4W1    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4W1_EKRN250021315-1A_22VH2WLT4_L6_1.fq.gz    | e5af4d3f8377e9c6e36c4c38f57c667f     | D4W1_EKRN250021315-1A_22VH2WLT4_L6_2.fq.gz    | e754ae64382f11c937afd75f937954f9     |
| ERS30371089  | PRJEB114114   | ERS30371089  | Illumina HiSeq 2500   | D4W3    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4W3_EKRN250021316-1A_22VH2WLT4_L6_1.fq.gz    | 1de6f58aed243ae2c8ca7f20a735e8d2     | D4W3_EKRN250021316-1A_22VH2WLT4_L6_2.fq.gz    | 875a1d27920cc8e26d0b1c04fee1d5df     |
| ERS30371090  | PRJEB114114   | ERS30371090  | Illumina HiSeq 2500   | D4W4    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4W4_EKRN250021317-1A_22VH2WLT4_L6_1.fq.gz    | 8c1d1297b5e0263976778a47ea1fe78d     | D4W4_EKRN250021317-1A_22VH2WLT4_L6_2.fq.gz    | f8096d7943ea1a673f52bf257598760e     |
| ERS30371091  | PRJEB114114   | ERS30371091  | Illumina HiSeq 2500   | D4EF1   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4EF1_EKRN250021318-1A_22VH2WLT4_L6_1.fq.gz   | c1aff857b135c3942e8e166181d1fec0     | D4EF1_EKRN250021318-1A_22VH2WLT4_L6_2.fq.gz   | 87a0f0ef259c3174e45011443d83a18b     |
| ERS30371092  | PRJEB114114   | ERS30371092  | Illumina HiSeq 2500   | D4EF2   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4EF2_EKRN250021319-1A_22VH2WLT4_L6_1.fq.gz   | e54743b1e0a6763c670f1bd2feea424f     | D4EF2_EKRN250021319-1A_22VH2WLT4_L6_2.fq.gz   | 2e63301df1021849baf3e6756619069f     |
| ERS30371093  | PRJEB114114   | ERS30371093  | Illumina HiSeq 2500   | D4EF4   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D4EF4_EKRN250021320-1A_22VH2WLT4_L6_1.fq.gz   | ba22aeaf96b0f22e87ddffd07caef941     | D4EF4_EKRN250021320-1A_22VH2WLT4_L6_2.fq.gz   | 6affe557fbcfd6464f53d7e349774682     |
| ERS30371094  | PRJEB114114   | ERS30371094  | Illumina HiSeq 2500   | D9W1    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9W1_EKRN250021321-1A_22VH2WLT4_L6_1.fq.gz    | 731cedcab3ecea38c477f6b8fbb7ee2c     | D9W1_EKRN250021321-1A_22VH2WLT4_L6_2.fq.gz    | 4755b9d7d3c999abd18fc38e8b52e90b     |
| ERS30371095  | PRJEB114114   | ERS30371095  | Illumina HiSeq 2500   | D9W2    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9W2_EKRN250021322-1A_22VH2WLT4_L6_1.fq.gz    | 247c63782958208d63d93b58f6d040ed     | D9W2_EKRN250021322-1A_22VH2WLT4_L6_2.fq.gz    | 042ba95a7735dad192b68917872dbcca     |
| ERS30371096  | PRJEB114114   | ERS30371096  | Illumina HiSeq 2500   | D9W4    | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9W4_EKRN250021323-1A_22VH2WLT4_L6_1.fq.gz    | e6a9c5f8e7a6f0969d88bb762cc91553     | D9W4_EKRN250021323-1A_22VH2WLT4_L6_2.fq.gz    | cf623c7fdd95603f733afce2bb9d79ab     |
| ERS30371097  | PRJEB114114   | ERS30371097  | Illumina HiSeq 2500   | D9EF1   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9EF1_EKRN250021324-1A_22VH2WLT4_L6_1.fq.gz   | 56878381a97069d36128d2c2de87e6c5     | D9EF1_EKRN250021324-1A_22VH2WLT4_L6_2.fq.gz   | 5085eb08c08007e37432d16472a8f846     |
| ERS30371098  | PRJEB114114   | ERS30371098  | Illumina HiSeq 2500   | D9EF3   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9EF3_EKRN250021325-1A_22VH2WLT4_L6_1.fq.gz   | 323b535ad221315a99900cbb93485c83     | D9EF3_EKRN250021325-1A_22VH2WLT4_L6_2.fq.gz   | 577fb898292bf22a29be72a1df626732     |
| ERS30371099  | PRJEB114114   | ERS30371099  | Illumina HiSeq 2500   | D9EF4   | TRANSCRIPTOMIC  | RANDOM PCR   | RNA-Seq PAIRED  | D9EF4_EKRN250021326-1A_22VH2WLT4_L6_1.fq.gz   | 94b3718948e0a3e782a24f8eae0e407a     | D9EF4_EKRN250021326-1A_22VH2WLT4_L6_2.fq.gz   | 3d95e958c1228d24e595ce9d9e7954f9     |

# prepare reads
ln -s ../*/01.RawData/*/*.fq.gz .

# prepare manifests file
prepare `gene_manifests.py`:
```

```python
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
```

prepare submission script:

```bash
#!/usr/bin/env bash
set -euo pipefail

WEBIN_USER="Webin-69760"
WEBIN_PASS="Xia05578678409!" ###password here
JAR="$HOME/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli/webin-cli-9.0.1.jar"
OUTDIR="./webin_out"

mkdir -p "$OUTDIR"

for m in manifests/*.manifest.txt; do
  echo "Validating $m"
  /rds/project/rds-FTKWLWDeHys/programs/java/jdk-17.0.10/bin/java -jar "$JAR" -context reads -manifest "$m" \
       -userName "$WEBIN_USER" -password "$WEBIN_PASS" \
       -outputDir "$OUTDIR" -validate

  echo "Submitting $m"
  /rds/project/rds-FTKWLWDeHys/programs/java/jdk-17.0.10/bin/java -jar "$JAR" -context reads -manifest "$m" \
       -userName "$WEBIN_USER" -password "$WEBIN_PASS" \
       -outputDir "$OUTDIR" -submit
done
```

Submission.
```bash
python3 generate_manifests.py
./submit_reads.sh
```
Successful submission should look like this:
```bash
% ./submit_reads.sh
Validating manifests/arab_c1.manifest.txt
INFO : Your application version is 9.0.1
INFO : Submission(s) validated successfully.
INFO : Creating report file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/./webin_out/./webin-cli.report
Submitting manifests/arab_c1.manifest.txt
INFO : Your application version is 9.0.1
INFO : Connecting to FTP server : webin2.ebi.ac.uk
INFO : Creating report file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/./webin_out/./webin-cli.report
INFO : Uploading file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/arab_c1_1.fq.gz
INFO : Uploading file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/arab_c1_2.fq.gz
INFO : Files have been uploaded to webin2.ebi.ac.uk.
INFO : The submission has been completed successfully. The following run accession was assigned to the submission: ERR15696041
INFO : The submission has been completed successfully. The following experiment accession was assigned to the submission: ERX15100159
```