# RNAseq analysis of Heterodera schachtii (X204SC25060928-Z01)

### 1. Reference

```bash
cd /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/

# index
/data/pathology/program/STAR/bin/Linux_x86_64/STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/STAR \
  --genomeFastaFiles /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
  --sjdbGTFfile /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.gtf \
  --sjdbOverhang 149

```

### 2. Trimming
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/02.CleanData

input_dir="/data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/01.RawData"
output_dir="/data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/02.CleanData"

# Run Trimmomatic iterate through each sample
for fq1 in "$input_dir"/*/*_1.fq.gz; do
    fq2="${fq1/_1.fq.gz/_2.fq.gz}" # Get the corresponding paired-end file
    cd "${output_dir}"

    # Extract sample name
    sample_name=$(basename "$fq1" | cut -d'_' -f1-3)

    # Set directory
    mkdir -p "${sample_name}"
    echo "Directory $sample_name created"
        # Change to the newly created directory
    cd "${sample_name}"
    echo "Changed directory to $sample_name"
    output_fwd_paired="${sample_name}_1P.fq.gz"
    output_fwd_unpaired="${sample_name}_1U.fq.gz"
    output_rev_paired="${sample_name}_2P.fq.gz"
    output_rev_unpaired="${sample_name}_2U.fq.gz"
    # Run Trimmomatic
    echo "Processing ${sample_name}"
    # Need to set specific quality filtering criteria in the following
    java -jar /data/pathology/program/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -summary "${sample_name}.summary" "$fq1" "$fq2" \
    "$output_fwd_paired" "$output_fwd_unpaired" "$output_rev_paired" "$output_rev_unpaired" \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
    cd "${output_dir}"
    echo "${sample_name} finished"
done

echo "Trimmomatic process completed."
```

### 3. Aligment

```bash
#!/bin/bash
# we are using STAR for alignment
input_dir="/data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/02.CleanData"
output_dir="/data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/03.Mapping"

# Run STAR iterate through each sample
for fq1 in "$input_dir"/*/*_1P.fq.gz; do
    fq2="${fq1/_1P.fq.gz/_2P.fq.gz}" # Get the corresponding paired-end file
    cd "${output_dir}"

    # Extract sample name
    sample_name=$(basename "$fq1" | cut -d'_' -f1-1)

    # Set directory
    mkdir -p "${sample_name}"
    echo "Directory $sample_name created"
    # Change to the newly created directory
    cd "${sample_name}"
    echo "Changed directory to $sample_name"

    # Run STAR
    echo "Processing ${sample_name}"
    /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
      --runThreadN 16 \
      --genomeDir /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/STAR \
      --readFilesIn ${fq1} ${fq2} \
      --readFilesCommand zcat \
      --outFilterMultimapNmax 1 \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix ${sample_name}

    samtools index ${sample_name}Aligned.sortedByCoord.out.bam
    
    # Run mapinsights bamqc
    /data/pathology/program/mapinsights/mapinsights/mapinsights bamqc -r /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa -i ${sample_name}Aligned.sortedByCoord.out.bam -o ./
    # Run QualiMap
    /data/pathology/program/QualiMap/qualimap_v2.3/qualimap rnaseq -bam ${sample_name}Aligned.sortedByCoord.out.bam -gtf /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.gtf --java-mem-size=32G -outdir ./qualimap/ -outformat pdf

    # Get expression raw counts
    /data/pathology/program/Miniforge3/bin/htseq-count --type transcript --counts_output 3.sorted.bam.count.tsv --nprocesses 16 --max-reads-in-buffer 1000000 ${sample_name}Aligned.sortedByCoord.out.bam /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.gtf
    # Get expression FPKM TPM
    /data/pathology/program/stringtie-3.0.0.Linux_x86_64/stringtie -p 32 -G /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.gtf -e -B -A 4.sorted.FPKM.tsv ${sample_name}Aligned.sortedByCoord.out.bam

    # Master Table
    sed 's/gene\:Hsc_/Hsc_/g' 3.sorted.bam.count.tsv | sort -k1,1V > 3.sorted.bam.count.sorted.tsv
    sed 's/gene\:Hsc_/Hsc_/g' 4.sorted.FPKM.tsv | sort -k1,1V | awk '$1 ~ "Hsc_" {print $0}' > 4.sorted.rawcount.sorted.tsv
    paste -d"\t" 4.sorted.rawcount.sorted.tsv 3.sorted.bam.count.sorted.tsv | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11}' |  sed '1i GeneID\tGeneName\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\tRawCounts' > ${sample_name}.MastTable.tsv

   cd "${output_dir}"

   echo "${sample_name} finished"
done

echo "STAR alignment process completed."

```

An example of quality control results from `mapinsights` is [here](https://github.com/chongjing/RNAseq_Hschachtii/tree/main/D4W1)

### 4. Differetial Expression analysis
#### 4.1 D4EF vs D4W
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/04.DE_analysis/01.D4EF_vs_D4W
paste -d"\t" ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W3/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W4/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF2/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF4/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' | sed 's/gene://g' > 01.counts.tsv
# add header to 01.counts.tsv 
# prepare 02.sample.tsv 
SampleID        condition
D4W1    D4W
D4W3    D4W
D4W4    D4W
D4EF1   D4EF
D4EF2   D4EF
D4EF4   D4EF

conda activate R4.3.2
```
```R
library(DESeq2)
count_matrix <- read.table('01.counts.tsv', sep = "\t", header=TRUE, row.names=1)
sample_info <- read.table('02.sample.tsv', sep = "\t", header=TRUE,row.names=1)
DEseqDataSet <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)
#perform pre-filtering to keep only rows that have a count of at least 5 for a minimal number of samples (2).
keep <- rowSums(counts(DEseqDataSet) >= 5) >= 3
DEseqDataSet <- DEseqDataSet[keep,]
dim(DEseqDataSet)
# set reference level
DEseqDataSet$condition <- relevel(DEseqDataSet$condition, ref = "D4W")
dds <- DESeq(DEseqDataSet)
res <- results(dds)
res
#summary 
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
# p value cut-off as 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
# adjusted p-value < 0.05
LFC > 0 (up)       : 2637, 12%
LFC < 0 (down)     : 2487, 12%
sum(res05$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(res05), file="03.P-cutoff-0.05.csv")
```

#### 4.2 D9EF vs D9W
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/04.DE_analysis/02.D9EF_vs_D9W
paste -d"\t" ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W2/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W4/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D9EF1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D9EF3/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9EF4/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' | sed 's/gene://g' > 01.counts.tsv
# add header to 01.counts.tsv 
# prepare 02.sample.tsv 
SampleID        condition
D9W1    D9W
D9W2    D9W
D9W4    D9W
D9EF1   D9EF
D9EF3   D9EF
D9EF4   D9EF

conda activate R4.3.2
```
```R
library(DESeq2)
count_matrix <- read.table('01.counts.tsv', sep = "\t", header=TRUE, row.names=1)
sample_info <- read.table('02.sample.tsv', sep = "\t", header=TRUE,row.names=1)
DEseqDataSet <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)
#perform pre-filtering to keep only rows that have a count of at least 5 for a minimal number of samples (2).
keep <- rowSums(counts(DEseqDataSet) >= 5) >= 3
DEseqDataSet <- DEseqDataSet[keep,]
dim(DEseqDataSet)
# set reference level
DEseqDataSet$condition <- relevel(DEseqDataSet$condition, ref = "D9W")
dds <- DESeq(DEseqDataSet)
res <- results(dds)
res
#summary 
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
# p value cut-off as 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
# adjusted p-value < 0.05
LFC > 0 (up)       : 4003, 19%
LFC < 0 (down)     : 4133, 19%
sum(res05$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(res05), file="03.P-cutoff-0.05.csv")
```

#### 4.3 D9W vs D4W
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/04.DE_analysis/03.D9W_D4W
paste -d"\t" ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W2/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9W4/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W3/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4W4/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' | sed 's/gene://g' > 01.counts.tsv
# add header to 01.counts.tsv 
# prepare 02.sample.tsv 
SampleID        condition
D9W1    D9W
D9W2    D9W
D9W4    D9W
D4W1    D4W
D4W3    D4W
D4W4    D4W

conda activate R4.3.2
```
```R
library(DESeq2)
count_matrix <- read.table('01.counts.tsv', sep = "\t", header=TRUE, row.names=1)
sample_info <- read.table('02.sample.tsv', sep = "\t", header=TRUE,row.names=1)
DEseqDataSet <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)
#perform pre-filtering to keep only rows that have a count of at least 5 for a minimal number of samples (2).
keep <- rowSums(counts(DEseqDataSet) >= 5) >= 3
DEseqDataSet <- DEseqDataSet[keep,]
dim(DEseqDataSet)
# set reference level
DEseqDataSet$condition <- relevel(DEseqDataSet$condition, ref = "D4W")
dds <- DESeq(DEseqDataSet)
res <- results(dds)
res
#summary 
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
# p value cut-off as 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
# adjusted p-value < 0.05
LFC > 0 (up)       : 1456, 6.7%
LFC < 0 (down)     : 1523, 7%
sum(res05$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(res05), file="03.P-cutoff-0.05.csv")
```

#### 4.4 D9EF vs D4EF
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/04.DE_analysis/04.D9EF_vs_D4EF
paste -d"\t" ../../X204SC25060928-Z01-F001_02/03.Mapping/D9EF1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D9EF3/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_01/03.Mapping/D9EF4/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF1/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF2/3.sorted.bam.count.tsv ../../X204SC25060928-Z01-F001_02/03.Mapping/D4EF4/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' | sed 's/gene://g' > 01.counts.tsv
# add header to 01.counts.tsv 
# prepare 02.sample.tsv 
SampleID        condition
D9EF1   D9EF
D9EF3   D9EF
D9EF4   D9EF
D4EF1   D4EF
D4EF2   D4EF
D4EF4   D4EF

conda activate R4.3.2
```
```R
library(DESeq2)
count_matrix <- read.table('01.counts.tsv', sep = "\t", header=TRUE, row.names=1)
sample_info <- read.table('02.sample.tsv', sep = "\t", header=TRUE,row.names=1)
DEseqDataSet <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)
#perform pre-filtering to keep only rows that have a count of at least 5 for a minimal number of samples (2).
keep <- rowSums(counts(DEseqDataSet) >= 5) >= 3
DEseqDataSet <- DEseqDataSet[keep,]
dim(DEseqDataSet)
# set reference level
DEseqDataSet$condition <- relevel(DEseqDataSet$condition, ref = "D4EF")
dds <- DESeq(DEseqDataSet)
res <- results(dds)
res
#summary 
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
# p value cut-off as 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
# adjusted p-value < 0.05
LFC > 0 (up)       : 1862, 8.7%
LFC < 0 (down)     : 1787, 8.3%
sum(res05$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(res05), file="03.P-cutoff-0.05.csv")
```

