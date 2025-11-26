# RNAseq analysis of Heterodera schachtii (X204SC25060928-Z01)

### 5. Genetic Relationship
The idea here is to see if the D4 samples and the D9 samples are genetically different populations, from within the same pool. Here I'll apply three approaches, 1) sample clustering based on expression; 2) sample clustering based on genetic distance (calculated from SNPs); 3) phylogenetic relationship based on SNPs.

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

The cluster and dendrogram are:
<table>
  <tr>
    <td><img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/01.clustering_expression/03.cluster.jpeg" alt="Image 1" width="400"/></td>
    <td><img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/01.clustering_expression/04.samples_dendrop.jpeg" alt="Image 2" width="400"/></td>
  </tr>
</table>

#### 5.2 Genetic distance
##### 5.2.1 GATK SNP calling - gVCF for each sample
There are several issue for using alignment from STAR with GATK SNP calling. First, STAR does not automatically add read group (@RG) tags, which are mandatory for GATK. (I aready updated above STAR code). Second, STAR assigns MAPQ values differently compared to DNA aligners like BWA. In RNA-seq, especially with splice junctions, STAR might assign lower MAPQ scores to reads that span introns or have multiple possible alignments. 
After solving these, do initial SNP calling for recalibration.
```bash
conda activate R4.3.2
export PATH=/data/pathology/program/Miniforge3/envs/R4.3.2/bin:$PATH

# create ref dictionary
/data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx16g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar CreateSequenceDictionary -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa

cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/06.GeneticDistance/01.GATK
for sample in D4W1 D4W3 D4W4 D4EF1 D4EF2 D4EF4; do
  cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/06.GeneticDistance/01.GATK
  mkdir ${sample} && cd ${sample}
  cd ${sample}

  # STAR alignment with proper RG injection, default without RG, cannot be used for GATK
  /data/pathology/program/STAR/bin/Linux_x86_64/STAR \
    --runThreadN 16 \
    --genomeDir /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/STAR \
    --readFilesIn /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/02.CleanData/${sample}_*/${sample}_*_1P.fq.gz /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/X204SC25060928-Z01-F001_01/02.CleanData/${sample}_*/${sample}_*_2P.fq.gz \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --outSAMmultNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --alignSJDBoverhangMin 1 \
    --alignSoftClipAtReferenceEnds Yes \
    --outSAMmapqUnique 60 \
    --quantMode GeneCounts \
    --outFileNamePrefix ${sample}

  # Add RG to Star alignment
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar AddOrReplaceReadGroups \
    -I ${sample}Aligned.sortedByCoord.out.bam \
    -O ${sample}Aligned.sortedByCoord.RG.bam \
    --RGID "${sample}" \
    --RGLB "lib1" \
    --RGPL "ILLUMINA" \
    --RGPU "unit1" \
    --RGSM "${sample}" \
    --CREATE_INDEX true

  # Mark duplicates
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar MarkDuplicatesSpark \
    -I ${sample}Aligned.sortedByCoord.RG.bam \
    -O ${sample}.dedup.bam \
    --create-output-bam-index true \
    --tmp-dir ${sample}_STARtmp

  # Split ‘N’ cigar reads (This step splits reads at introns and hard-clips overhangs so HaplotypeCaller can cope with them.)
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar SplitNCigarReads \
      -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
      -I ${sample}.dedup.bam \
      -O ${sample}Aligned.sortedByCoord.MAPQ.bam \
      --tmp-dir ${sample}_STARtmp

  # Initial SNP/indel calling for known-sites
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar HaplotypeCaller \
      -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
      -I ${sample}Aligned.sortedByCoord.MAPQ.bam \
      -O ${sample}.snps.indels.vcf.gz \
      --dont-use-soft-clipped-bases \
      --standard-min-confidence-threshold-for-calling 20.0

  /data/pathology/program/VCFtools/vcftools/bin/vcftools --gzvcf ${sample}.snps.indels.vcf.gz --min-alleles 2 --max-alleles 2 --minQ 5000 --recode --recode-INFO-all --out ${sample}.gold.snps.indels
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar IndexFeatureFile -I ${sample}.gold.snps.indels.recode.vcf

  # Base quality score recalibration
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar BaseRecalibrator \
    -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
    -I ${sample}Aligned.sortedByCoord.MAPQ.bam \
    --known-sites ${sample}.gold.snps.indels.recode.vcf \
    -O ${sample}.recal_data.table

  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar ApplyBQSR \
    -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
    -I ${sample}Aligned.sortedByCoord.MAPQ.bam \
    --bqsr-recal-file ${sample}.recal_data.table \
    -O ${sample}.recal.bam

  #Final gVCF emission
  /data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx32g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar HaplotypeCaller \
    -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
    --emit-ref-confidence GVCF \
    --max-mnp-distance 0 \
    -I ${sample}.recal.bam \
    -O ${sample}.snps.indels.g.vcf.gz \
    --dont-use-soft-clipped-bases   

  cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/06.GeneticDistance/01.GATK
done
```

#### 5.2.2 GATK final SNP calling
```bash
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/06.GeneticDistance/01.GATK

/data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx64g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar CombineGVCFs \
    -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
    --variant ./D4EF1/D4EF1.snps.indels.g.vcf.gz \
    --variant ./D4EF2/D4EF2.snps.indels.g.vcf.gz \
    --variant ./D4EF4/D4EF4.snps.indels.g.vcf.gz \
    --variant ./D4W1/D4W1.snps.indels.g.vcf.gz \
    --variant ./D4W3/D4W3.snps.indels.g.vcf.gz \
    --variant ./D4W4/D4W4.snps.indels.g.vcf.gz \
    --variant ./D9EF1/D9EF1.snps.indels.g.vcf.gz \
    --variant ./D9EF3/D9EF3.snps.indels.g.vcf.gz \
    --variant ./D9EF4/D9EF4.snps.indels.g.vcf.gz \
    --variant ./D9W1/D9W1.snps.indels.g.vcf.gz \
    --variant ./D9W2/D9W2.snps.indels.g.vcf.gz \
    --variant ./D9W4/D9W4.snps.indels.g.vcf.gz \
    -O 11.AllSample.snps.indels.g.vcf.gz

/data/pathology/program/java/jdk-17.0.10+7/bin/java -Xmx64g -jar /data/pathology/program/GATK/gatk-4.5.0.0/gatk.jar GenotypeGVCFs \
    -R /data/pathology/cxia/projects/0.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.genomic.fa \
    -V 11.AllSample.snps.indels.g.vcf.gz \
    -O 11.AllSample.snps.indels.vcf.gz
```

#### 5.3.1 Genetic distance
```bash
# only keep high-quality, no missing data
/data/pathology/program/VCFtools/vcftools/bin/vcftools --gzvcf 11.AllSample.snps.indels.vcf.gz --minQ 10000 --max-missing 1.0 --maf 0.05 --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out 12.AllSample.snps.minQ10000.noMissing
/data/pathology/program/VCF2Dis-1.54/bin/VCF2Dis -InPut 11.AllSample.snps.indels.vcf.gz -OutPut 12.GeneticDis
```
The genetic distance and phylogenetic tree:
<div align="center">
  <img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/03.phylogenetic/GeneticDistance.jpg" alt="Genetic Distance" width="600"/>
</div>

<div align="center">
  <img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/03.phylogenetic/PhylogeneticTree.jpg" alt="Phylogenetic Tree" width="600"/>
</div>

#### 5.3.2 Cluster
```bash
/data/pathology/program/VCF2PCACluster/bin/VCF2PCACluster -InVCF 12.AllSample.snps.minQ10000.noMissing.recode.vcf.gz -OutPut 13.AllSample.snps.minQ10000.noMissing.cluster -InSampleGroup 13.AllSample.Group
```
<img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/02.clustering_genetic/13.AllSample.snps.minQ10000.noMissing.cluster.N.PC1_PC2.p.jpg" alt="Image 1" width="600"/>

### 6. Genetic Relationship (non-DEGs)
#### 6.1 non-DEGs in all pairwise comparisons
```bash
# extract genes with -0.5 < log2FoldChange < 0.5 && padj > 0.05
cd /data/pathology/cxia/projects/Sebastian/07.X204SC25060928-Z01-F001/04.DE_analysis
for dir in *; do echo "processing ${dir}"; 
  awk -F"," 'function abs(x) { return (x < 0) ? -x : x } abs($3)<"0.5" && $7 > "0.05" {print $1}' ${dir}/03.P-cutoff-0.05.csv | sed 's/"//g' | sort -t "_" -k3,3n > ${dir}.nonDEGs.list; 
done
```

```python
# get non-DEGs in all pairwise comparisons
def read_file_lines(path):
    with open(path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

# Replace these with your actual file paths
files = ["01.D4EF_vs_D4W.nonDEGs.list", "02.D9EF_vs_D9W.nonDEGs.list", "03.D9W_D4W.nonDEGs.list", "04.D9EF_vs_D4EF.nonDEGs.list"]

# Read and intersect
sets = [read_file_lines(f) for f in files]
common = set.intersection(*sets)

# Write to output file
with open("20.nonDEGs.list", "w") as out:
    for item in sorted(common):
        out.write(item + "\n")

print(f"✅ Intersection written to '20.nonDEGs.list' with {len(common)} items.")
```

#### 6.2 Genetic distance
```bash
cd /home/cx264/project/07.X204SC25060928-Z01-F001/07.GeneticDistance.noDEGs
while read line; do awk -F"\t" -v line=$line '$4==line {print $0}' /home/cx264/project/00.ref/1.HeteroderaSchachtii/heterodera_schachtii.PRJNA522950.WBPS19/heterodera_schachtii.PRJNA522950.WBPS19.annotations.bed >> 21.nonDGE.bed; done < 20.nonDEGs.list
sort -k1,1 -k2,2n  21.nonDGE.bed > 21.nonDGE.sort.bed

# get all SNPs within nonDEGs
vcftools --gzvcf 12.AllSample.snps.minQ10000.noMissing.recode.vcf.gz --min-alleles 2 --max-alleles 2 --minQ 10000 --bed 21.nonDGE.sort.bed --remove-indels --recode --recode-INFO-all --stdout | gzip -c > 22.AllSample.snps.noDEGs.vcf.gz
#Read 6522 BED file entries.
#After filtering, kept 65057 out of a possible 213116 Sites

/home/cx264/program/VCF2Dis-1.54/bin/VCF2Dis -InPut 22.AllSample.snps.noDEGs.vcf.gz -OutPut 23.GeneticDis.mat
```
The genetic distance and phylogenetic tree:
<div align="center">
  <img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/05.noDEGs/GeneticDistance.jpg" alt="Genetic Distance" width="600"/>
</div>

<div align="center">
  <img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/05.noDEGs/PhylogeneticTree.jpg" alt="Phylogenetic Tree" width="600"/>
</div>

#### 6.3 Cluster
```bash
/home/cx264/program/VCF2PCACluster/bin/VCF2PCACluster -InVCF 22.AllSample.snps.noDEGs.vcf.gz -OutPut 24.AllSample.snps.minQ10000.noMissing.noDEGs.cluster -InSampleGroup 22.AllSample.Group
```
<img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/05.noDEGs/24.AllSample.snps.minQ10000.noMissing.noDEGs.cluster.C.PC1_PC2.p.jpg" alt="Image 1" width="600"/>

### 7. Fst Population Differentiation (non-DEGs)
SNPs with minQ 10000 and no missing data in any sample are used for Fst calculation using VCFtools. Number of SNPs with Fst above 0.25 (highly differentiated) are summarized.

```bash
cd /home/cx264/project/07.X204SC25060928-Z01-F001/08.Fst.noDEGs
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D4W_population.txt --weir-fst-pop D4EF_population.txt --out 002.D4.Water_vs_EF
awk -F"\t" '$3>0.25 {print $0}' 002.D4.Water_vs_EF.weir.fst | less -N
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D4W_population.txt --weir-fst-pop D9W_population.txt --out 002.D4W_vs_D9W
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D4W_population.txt --weir-fst-pop D9EF_population.txt --out 002.D4W_vs_D9EF
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D4EF_population.txt --weir-fst-pop D9W_population.txt --out 002.D4EF_vs_D9W
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D4EF_population.txt --weir-fst-pop D9EF_population.txt --out 002.D4EF_vs_D9EF
vcftools --gzvcf ../07.GeneticDistance.noDEGs/22.AllSample.snps.noDEGs.vcf.gz --weir-fst-pop D9W_population.txt --weir-fst-pop D9EF_population.txt --out 002.D9W_vs_D9EF
```
<img src="https://github.com/chongjing/RNAseq_Hschachtii/blob/main/GATK/06.Fst_noDEGs/Fst_noDEGs.jpg" alt="Image 1" width="600"/>