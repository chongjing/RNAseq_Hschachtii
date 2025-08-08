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
    # this code is working for both expression counts and SNP calling with GATK, see bellow comments
    echo "Processing ${sample_name}"
    /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
      --runThreadN 16 \
      --genomeDir /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/STAR \
      --readFilesIn ${fq1} ${fq2} \
      --readFilesCommand zcat \
      --outFilterMultimapNmax 1 \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --outSAMstrandField intronMotif \
      --outFilterIntronMotifs RemoveNoncanonical \
      --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds Yes \
      --outSAMattributes NH HI AS NM MD RG \
      --outSAMattrRGline "ID:${sample_name} PL:ILLUMINA PU:unit1 LB:lib1 SM:${sample_name}" \
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

Summarize raw counts, FPKM, TPM.
```bash
#post-processing
# check head and tail
# sed 's/gene\:Hsc_/Hsc_/g' 3.sorted.bam.count.tsv | sort -k1,1V > 3.sorted.bam.count.sorted.tsv
# sed 's/gene\:Hsc_/Hsc_/g' 4.sorted.FPKM.tsv | sort -k1,1V > 4.sorted.FPKM.sorted.tsv
# paste -d"\t" 4.sorted.FPKM.sorted.tsv 3.sorted.bam.count.sorted.tsv > 5.r1_12dpi_Co.counts.FPKM.tsv
```

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

### 5. Genetic Relationship
The idea here is to see if the D4 samples and the D9 samples are genetically different populations, from within the same pool. Here I'll apply two approaches, 1) sample clustering based on expression; 2) genetic distance based on SNPs.
#### 5.1 TomicsVis clustering
First to prepare counts and group information.
```bash
paste -d"\t" ../X204SC25060928-Z01-F001_02/03.Mapping/D4W1/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D4W3/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D4W4/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D4EF1/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D4EF2/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D4EF4/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_01/03.Mapping/D9W1/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_01/03.Mapping/D9W2/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_01/03.Mapping/D9W4/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D9EF1/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_02/03.Mapping/D9EF3/3.sorted.bam.count.tsv ../X204SC25060928-Z01-F001_01/03.Mapping/D9EF4/3.sorted.bam.count.tsv | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18"\t"$20"\t"$22"\t"$24"\t"$26"\t"$28"\t"$30"\t"$32"\t"$34"\t"$36}' | sed 's/gene://g' | sed '1i GeneID\tD4W1\tD4W3\tD4W4\tD4EF1\tD4EF2\tD4EF4\tD9W1\tD9W2\tD9W4\tD9EF1\tD9EF3\tD9EF4' > 01.counts.tsv
mamba activate R4.2.3
/data/pathology/program/Miniforge3/envs/R4.2.3/bin/R
```
```R
library(TOmicsVis)
library(readr)

counts <- read_tsv("01.counts.tsv")
group <- read_tsv("01.groups.tsv")
group <- as.data.frame(group)

# Correlation Heatmap for samples/groups based on Pearson/Spearman algorithm
pdf("02.corr_heatmap.pearson.pdf")
corr_heatmap(
  data = counts,
  corr_method = "pearson",
  cell_shape = "square",
  fill_type = "full",
#  multi_shape = FALSE,
  lable_size = 2,
  axis_angle = 45,
  axis_size = 9,
  lable_digits = 1,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_bw"
)
corr_heatmap(
  data = counts,
  corr_method = "spearman",
  cell_shape = "square",
  fill_type = "full",
#  multi_shape = FALSE,
  lable_size = 2,
  axis_angle = 45,
  axis_size = 9,
  lable_digits = 1,
  color_low = "blue",
  color_mid = "white",
  color_high = "red",
  outline_color = "white",
  ggTheme = "theme_bw"
)
dev.off()

# pca_analysis
pdf("03.cluster.pdf")
pca_plot(
  sample_gene = counts,
  group_sample = group,
  multi_shape = FALSE,
  xPC = 1,
  yPC = 2,
  point_size = 4,
  text_size = 2,
  fill_alpha = 0.10,
  border_alpha = 0.00,
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)

# tsne plot
tsne_plot(
  sample_gene = counts,
  group_sample = group,
  seed = 1,
  multi_shape = FALSE,
  point_size = 3,
  point_alpha = 0.8,
  text_size = 5,
  text_alpha = 0.60,
  fill_alpha = 0.20,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)

# UMAP, # doesn't work: Error: umap: number of neighbors must be smaller than number of items
umap_plot(
  sample_gene = counts,
  group_sample = group,
  seed = 1,
  neighbors = 10,
  multi_shape = FALSE,
  point_size = 3,
  point_alpha = 1,
  text_size = 5,
  text_alpha = 0.60,
  fill_alpha = 0.20,
  border_alpha = 0.00,
  sci_fill_color = "Sci_AAAS",
  legend_pos = "right",
  legend_dir = "vertical",
  ggTheme = "theme_bw"
)
dev.off()

# dendro_plot
pdf("04.samples_dendrop.pdf")
dendro_plot(
  data = counts,
  dist_method = "euclidean",
  hc_method = "ward.D2",
  tree_type = "rectangle",
  k_num = 5,
  palette = "npg",
  color_labels_by_k = TRUE,
  horiz = FALSE,
  label_size = 1,
  line_width = 1,
  rect = TRUE,
  rect_fill = TRUE,
  xlab = "Samples",
  ylab = "Height",
  ggTheme = "theme_light"
)
dev.off
```

#### 5.2 Genetic distance
##### 5.2.1 GATK SNP calling
There are several issue for using alignment from STAR with GATK SNP calling. First, STAR does not automatically add read group (@RG) tags, which are mandatory for GATK. (I aready updated above STAR code). Second, STAR assigns MAPQ values differently compared to DNA aligners like BWA. In RNA-seq, especially with splice junctions, STAR might assign lower MAPQ scores to reads that span introns or have multiple possible alignments. 
After solving these, do initial SNP calling for recalibration.
```bash
conda activate R4.3.2
export PATH=/data/pathology/program/Miniforge3/envs/R4.3.2/bin:$PATH

# create ref dictionary
#/data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx16g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar CreateSequenceDictionary -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa

cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/06.GeneticDistance/01.GATK
for sample in D4W1 D4W3 D4W4 D4EF1 D4EF2 D4EF4; do

## bam from STAR doesn't have RG, add it if needed
  # /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx16g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar AddOrReplaceReadGroups -I /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_02/03.Mapping/${sample}/${sample}Aligned.sortedByCoord.out.bam -O /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_02/03.Mapping/${sample}/${sample}Aligned.sortedByCoord.RG.bam --RGID ${sample} --RGLB lib1 --RGPL Illumina --RGPU unit1 --RGSM ${sample} --CREATE_INDEX true
  # #Recalculate MAPQ for RNAseq
  # /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx16g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar SplitNCigarReads -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa -I /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/03.Mapping/${sample}/${sample}Aligned.sortedByCoord.RG.bam -O /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/03.Mapping/${sample}/${sample}Aligned.sortedByCoord.MAPQ.bam
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx16g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar HaplotypeCaller -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa -I /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_02/03.Mapping/${sample}/${sample}Aligned.sortedByCoord.out.bam -O ${sample}.snps.indels.vcf.gz --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20.0
  /data/pathology/program/VCFtools/vcftools/bin/vcftools --gzvcf ${sample}.snps.indels.vcf.gz --min-alleles 2 --max-alleles 2 --minQ 5000 --recode --recode-INFO-all --out ${sample}.gold.snps.indels
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar IndexFeatureFile -I ${sample}.gold.snps.indels.recode.vcf
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar BaseRecalibrator -R /data/pathology/cxia/projects/Sebastian/SexDetermine/1.data/HS_inhertied_snps/genome/Cam_Hsc_genome1.2.fasta -I /data/pathology/cxia/projects/Sebastian/SexDetermine/1.data/HS_inhertied_snps/bams/${sample}_sorted_marked_RG.bam --known-sites ${sample}.gold.snps.indels.recode.vcf -O ${sample}.recal_data.table
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar ApplyBQSR -R /data/pathology/cxia/projects/Sebastian/SexDetermine/1.data/HS_inhertied_snps/genome/Cam_Hsc_genome1.2.fasta -I /data/pathology/cxia/projects/Sebastian/SexDetermine/1.data/HS_inhertied_snps/bams/${sample}_sorted_marked_RG.bam --bqsr-recal-file ${sample}.recal_data.table -O ${sample}.recal.bam
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar HaplotypeCaller -R /data/pathology/cxia/projects/Sebastian/SexDetermine/1.data/HS_inhertied_snps/genome/Cam_Hsc_genome1.2.fasta --emit-ref-confidence GVCF --max-mnp-distance 0 -I ${sample}.recal.bam -O ${sample}.snps.indels.g.vcf.gz

done
```