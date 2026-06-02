# ENA Submission: Raw RNA-seq Reads

This directory documents the submission of raw RNA-seq reads (12 samples, 24 paired-end FASTQ files, ~203 GB) to the European Nucleotide Archive (ENA) under BioProject **PRJEB114114**.

---

## Accessions

| Sample | ERS Accession | Run Accession | Experiment Accession |
|--------|--------------|---------------|---------------------|
| D4W1   | ERS30371088  | ERR17359749   | ERX16749156         |
| D4W3   | ERS30371089  | ERR17359750   | ERX16749157         |
| D4W4   | ERS30371090  | ERR17359751   | ERX16749158         |
| D4EF1  | ERS30371091  | ERR17359748   | ERX16749155         |
| D4EF2  | ERS30371092  | ERR17359752   | ERX16749159         |
| D4EF4  | ERS30371093  | ERR17359753   | ERX16749160         |
| D9W1   | ERS30371094  | ERR17359754   | ERX16749161         |
| D9W2   | ERS30371095  | ERR17359755   | ERX16749162         |
| D9W4   | ERS30371096  | ERR17359756   | ERX16749163         |
| D9EF1  | ERS30371097  | ERR17359757   | ERX16749164         |
| D9EF3  | ERS30371098  | ERR17359758   | ERX16749165         |
| D9EF4  | ERS30371099  | ERR17359759   | ERX16749166         |

- **BioProject:** PRJEB114114
- **Submission IDs:** ERA36371921–ERA36371932

---

## Submission Pipeline (Step by Step)

### 1. Metadata Preparation

The file `metadata_template.tsv` contains all sample metadata:
- Sample ERS accessions (pre-registered in ENA BioSamples)
- Library names, sequencing platform, library strategy/source/selection
- FASTQ filenames and MD5 checksums
- Layout (PAIRED, insert size 250 bp)

### 2. Manifest Generation

`generate_manifests.py` reads the TSV and generates one webin-cli manifest per sample:

```bash
python generate_manifests.py
```

Manifests are written to `manifests/` with fields: STUDY, SAMPLE, NAME, INSTRUMENT, LIBRARY_NAME, LIBRARY_SOURCE, LIBRARY_SELECTION, LIBRARY_STRATEGY, INSERT_SIZE, and FASTQ entries with `READ_TYPE=paired`.

### 3. FTP Upload (Curl — Bypassing webin-cli's timeout)

The `webin-cli` built-in FTP client consistently timed out on large files (8–9 GB) from this HPC environment. Files were uploaded manually using `curl` with robust retry settings:

```bash
# Upload to webin-cli subdirectory structure
curl -T "$local_file" \
  --retry 5 --retry-delay 30 \
  --connect-timeout 60 --max-time 3600 \
  --ftp-create-dirs \
  "ftp://webin2.ebi.ac.uk/webin-cli/reads/<MANIFEST_NAME>/<FILE>"
```

File path on ENA FTP: `webin-cli/reads/<MANIFEST_NAME>/<FILE>`

All 24 files (~203 GB total) were uploaded sequentially. Average speed: 20–37 MB/s.

### 4. XML Generation (PAIRED Layout)

Rather than relying on webin-cli (which generated incorrect `<SINGLE />` library layout), submission XMLs were generated programmatically using a Python script with the correct PAIRED layout:

- **submission.xml** — ADD action with manifest embedded as CDATA
- **experiment.xml** — Experiment design with `<PAIRED><NOMINAL_LENGTH>250</NOMINAL_LENGTH></PAIRED>`
- **run.xml** — Run data block referencing FTP file paths by MD5 checksum

Generated XMLs for all 12 samples are in `ena_xmls/<SAMPLE_NAME>/`.

### 5. REST API Submission

XMLs were submitted via multipart POST to the ENA Webin drop-box REST API:

```bash
curl -u "Webin-XXXXX:password" \
  -F "SUBMISSION=@submission.xml" \
  -F "EXPERIMENT=@experiment.xml" \
  -F "RUN=@run.xml" \
  "https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
```

Each submission returned `success="true"` with assigned accession numbers. All 12 samples were submitted in under 5 seconds.

### 6. Data Availability

All reads are available under open access at:
`https://www.ebi.ac.uk/ena/browser/view/PRJEB114114`

---

## Key Technical Notes

- **FTP timeout:** webin-cli 9.0.3's built-in FTP client timed out on files >8 GB from this HPC. The timeout is not configurable via CLI flags. Solution: upload via `curl` with `--max-time 3600`, then submit via REST API.
- **Library layout:** ENA XML schema does not accept `<NOMINAL_LENGTH>` as a child element of `<PAIRED>`. Use `<PAIRED />` alone or `<PAIRED NOMINAL_LENGTH="250" />` (attribute form, not validated by schema).
- **File paths in RUN XML:** Must match the FTP path relative to the Webin account root. Our files were at `webin-cli/reads/<NAME>/<FILE>`.

---

## Files in This Directory

| File / Directory | Description |
|-----------------|-------------|
| `metadata_template.tsv` | Sample metadata and checksums |
| `generate_manifests.py` | Script to generate webin-cli manifests |
| `manifests/` | Generated webin-cli manifest files (×12) |
| `checksums.md5` | MD5 checksums for all FASTQ files |
| `submit_reads.sh` | Original webin-cli submission script (deprecated) |
| `ena_xmls/` | Generated submission XMLs (×12 samples, 3 files each) |
| `data_availability_statement.md` | Publication-ready data availability statement |

---

## Dependencies

- `curl` (for FTP upload)
- Python 3 (for XML/manifest generation)
- Java 17 + `webin-cli-9.0.3.jar` (validation only — submission done via REST API)