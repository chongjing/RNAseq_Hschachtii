#!/usr/bin/env bash
set -euo pipefail

WEBIN_USER="Webin-69760"
WEBIN_PASS="Xia05578678409!"
JAR="$HOME/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli/webin-cli-9.0.3.jar"
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
