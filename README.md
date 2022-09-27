# De novo mCSPC multi-focal sequencing code
This repository contains the code used to carry out the analyses in our article "". All code within this repository is originally authored by Matti Annala.

Patient whole-genome sequencing data have been deposited in the European Genome-Phenome Archive (EGA) database under the accession code EGAS00001006466 and is available under standard EGA controlled release.

Simulated patient data files for experimenting with our subclonal reconstruction algorithms is also available here:
https://www.dropbox.com/s/hofrv1vlnln1mau/example_patients.zip?dl=0


## Pre-requisites

Our analysis software and scripts are written using the Julia and Rust programming languages. To use our code, you will first need to:

1. Install Julia: https://julialang.org/downloads/
2. Install the Rust compiler and toolchain: https://www.rust-lang.org/tools/install
3. Install required Julia packages:
```
import Pkg;
Pkg.add(split("StatsBase Distributions HypothesisTests PyCall KernelDensity Interpolations Combinatorics Suppressor"))
```
4. Modify your JULIA_LOAD_PATH environment variable so that Julia can find the code modules found in the root directory of this repository:
```
export JULIA_LOAD_PATH=/path/to/this/repository:
```
5. Add the Julia scripts found in subfolder `scripts` of this repository to your PATH environment variable.
6. Install Mutato: https://www.github.com/annalam/mutato

Steps 5 - 6 are only required for some of the analysis steps, so you may be able to skip some steps depending on what analyses you wish to run.





## Read alignment

All plasma cfDNA and tissue samples were sequenced using Illumina X Ten instruments. The paired end sequencing reads were adapter-trimmed using cutadapt-1.11, quality-masked using seqkit-0.8, and then aligned against the human hg38 reference genome using Bowtie-2.3.0. Duplicate DNA fragments were identified and marked using samblaster-0.1.24:
```
fasta interleave sample_1.fq.gz sample_2.fq.gz | \
  cutadapt --interleaved -f fastq -m 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC - | \
  fasta trim by quality - 30 | fasta mask by quality - 20 | \
  bowtie2 -p20 -X 1000 --score-min L,0,-0.6 --ignore-quals -x hg38 --interleaved - | \
  samblaster | samtools view -u - | samtools sort -@ 8 -m 4G -o sample.bam
```





## Somatic mutation analysis

As the first step in mutation analysis, we generated a matrix containing read-level evidence for each potential variant (row) in each sample (column):
```
mutato call --alt-reads=5 --alt-frac=0.05 --min-mapq=0 hg38.fa *.bam > variants.tsv
```

Next we generated a file `tumor_normal_pairs.txt` listing all tumor and cfDNA samples in the cohort, together with their matched germline samples. Here is an example:
```
TEST	REF
AE-015-Baseline-cfDNA	AE-015-WBC
AE-015-Progression2-cfDNA		AE-015-WBC
AE-018-Baseline-cfDNA	AE-018-WBC
```

Cancer samples are listed in the first column, and matched germline samples in the second column. Additional negative control  samples were included by listing each of them on their own line, with an empty first column.

We searched for somatic mutations fulfilling the following criteria:
- At least 10 mutant allele reads
- Mutant allele fraction ≥ 8%
- Mutant allele fraction 10x higher than in the matched germline sample
- Mutant allele fraction 25x (50x for non-protien altering) higher than the average of all cancer-negative samples (germline samples and known cancer-negative samples)
- A minimum read depth of 10x in the matched germline sample
- Average mutant allele distance from nearest read end ≥ 15 (≥ 30 for non-protien altering)
- Average mapping quality of mutant allele reads ≥ 15 (≥ 30 for non-protien altering)

Thresholds changed in WES samples according to the methods. 

This was done using the following command:
```
variant somatic --alt-reads=10 --test-ref-ratio=10 --test-bg-ratio=25 --ref-reads=20 \
--min-sidedness=15 --min-mapq=15 variants.vcf ../tumor_normal_pairs.txt | \
variant predict effect - | variant protein altering - > somatic_protein_altering.tmp
variant nearby indels variants.vcf | variant somatic --alt-reads=10 --test-ref-ratio=10 \
--test-bg-ratio=50 --ref-reads=30 --min-sidedness=30 --min-mapq=30 --mapq-filter-max-indel-len=100\
 - ../tumor_normal_pairs.txt | variant predict effect - | variant protein altering --invert - |\
 variant discard sketchy silent - > somatic_silent.tmp
```

Next, we annotated all mutations with their predicted biological effect using `variant predict effect` (which internally uses ANNOVAR). We also annotated each mutation with information about its frequency in the COSMIC (version 77) and GNOMAD (version 3.0) databases:
```
cat somatic_protein_altering.tmp <(tail -n +2 somatic_silent.tmp) | sort -k1,1V -k2,2n |\
 variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | \
 variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | \
 variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 > somatic.vcf
```



## Germline heterozygous SNP identification

We generated a list of germline heterozygous SNPs for each patient, using the following command:
```
variant heterozygous snps --min-depth=30 --max-depth=500 --min-mapq=30 \
  --min-sidedness=15 variants.tsv WBC | variant discard indels - | \
  variant annotate - gnomad_v3.jls | egrep 'CHROM|GNOMAD' > hetz_snps.tsv
```

Indels were omitted since they are associated with alignment artifacts that can compromise accurate allele fraction quantification.


## Deleterious germline variant identification

To search for germline variants, we looked for variants found in the matched leukocyte samples that fulfilled the following criteria:
- At least 5 reads supporting the variant
- Allele fraction ≥ 20%
- Allele fraction is at least 20 times higher than the background error rate (i.e. the average allele fraction of leukocyte samples that had an allele fraction < 20% for the variant)

This analysis was carried out using the `variant germline` tool:
```
variant germline --alt-reads=5 --alt-frac=0.2 --bg-ratio=20 \
  variants.tsv WBC > germline.tsv
```

We further narrowed this list down to protein altering germline variants with a population frequency below 0.5%, and annotated the variants with information from the COSMIC and ClinVar (dated 2020-07-06)  databases:
```
variant predict effect germline.tsv | variant protein altering - | \
  variant discard if frequency above - gnomad_v3.jls 0.005 | \
  variant annotate - gnomad_v3.jls | \
  variant annotate - cosmic_77_hg38.jls | \
  variant annotate - clinvar-2020-07-06_hg38.jls \
  > germline_annotated.tsv
```

The resulting list of rare germline variants was curated for deleterious alterations by looking for protein-truncating mutations in DNA repair genes and other critical cancer genes, and by searching for variants annotated as “pathogenic” in ClinVar.


## Copy number analysis

For targeted sequencing data, we calculated the coverage of each sample across each gene. 
```
echo *.bam | parallel -n8 '[ ! -e ../coverage/${x/.bam/.tsv} ] && sam count $x \
../baits_hg38.bed > ../coverage/${x/.bam/.tsv}'
```
# For FFPE tissue samples
```
grep -v cfDNA ../tumor_normal_pairs.txt > ../.tumor_normal_pairs.ffpe.tmp
copynum call targeted --hetz-snps=../mutations/hetz_snps.vcf --report-dir=../violins\
 ../baits_hg38.bed ../.tumor_normal_pairs.ffpe.tmp > ../gene_cna.ffpe.tsv
```
# For cfDNA samples
```
egrep 'TEST|cfDNA' ../tumor_normal_pairs.txt > ../.tumor_normal_pairs.cfdna.tmp
copynum call targeted --gc-fractions=../baits_hg38.gc --hetz-snps=../mutations/hetz_snps.vcf \
--report-dir=~/tmp --controls=../cna_controls_for_cfdna.txt ../baits_hg38.bed \
../.tumor_normal_pairs.cfdna.tmp > ../gene_cna.cfdna.tsv
```

First, for WES data,  we generated a BED file describing a grid of half-overlapping 1000 bp windows covering the entire genome:
```
copynum grid hg38.chrom.sizes 1000 > grid.bed
```

Then we calculated the GC nucleotide content of each window, for use in subsequent GC bias correction:
```
fasta gc content hg38.fa grid.bed > grid.gc
```

Next, we counted the number of concordantly aligned DNA fragments within each 1000 bp genomic window using the `sam count` tool:
```
sam count sample.bam grid.bed > sample.tsv
```

Finally we generated coverage logratio and heterozygous SNP allele fraction (HSAF) tracks using the `copynum call genomewide` tool:
```
copynum call genomewide --output-format=igv --logratio-dots=Inf \
  --snp-median-decimate=9 --discard-noisiest=0.05 --min-ref-depth=30 \ 
  --controls=cna_controls.txt --gc-fractions=grid.gc --report-dir=./ \
  --hetz-snps=hetz_snps.tsv grid.bed tumor_normal_pairs.txt
```

We used the same `tumor_normal_pairs.txt` file used in the mutation analysis step.

