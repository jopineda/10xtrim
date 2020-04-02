# Quickstart - how to trim 10X-specific artifacts

The original purpose of 10xtrim was to improve the *tumour-only* somatic mutation calling accuracy when using 10X Genomics' linked-reads. Here we provide a step-by-step tutorial to help you get started.

**Requirements**:

* `10xtrim`
* [`samtools`](https://htslib.org)
* [`picard`](https://github.com/broadinstitute/picard)
* `python 2.7`
* [`MuTect1`](https://github.com/broadinstitute/mutect)

## Download example dataset

You can download the example dataset we will use here: ::

    wget http://s3.climb.ac.uk/nanopolish_tutorial/ecoli_2kb_region.tar.gz
    tar -xvf ecoli_2kb_region.tar.gz
    cd ecoli_2kb_region

**Details**:

* Sample: HCC1954 (Human cell line derived from a primary stage IIA, grade 3 invasive ductal carcinoma) 
* Instrument : Illumina X Ten
* Region: "chr20:7586000-7588000"
* Reference: hg19
* [Tumour phased BAM](https://support.10xgenomics.com/genome-exome/datasets/2.1.0/HCC1954T_WGS_210)
* [Normal phased BAM](https://support.10xgenomics.com/genome-exome/datasets/2.1.0/HCC1954N_WGS_210)

This is a subset of phased reads aligned to a 2kb region of hg19 reference genome. 

You should find the following files:

* ``tumour.phased.md.sorted.bam`` : tumour read alignments to reference
* ``normal.phased.md.sorted.bam`` : tumour read alignments to reference
* ``highconf.hg19.bed`` : high confidence intervals from GIAB
* ``cosmic.hg19.vcf`` : known somatic mutations from COSMIC
* ``dbsnp.hg19.vcf``  : known common variants from dbSNP

## Objective

In this tutorial, we aim to remove a 10X-specific false positive (FP) variant when calling in tumour-only mode. 
Here, we define a 10X-specific artifact as a variant NOT found when calling somatc mutations in matched tumour-normal mode.

The FP multinucleotide variant (MNV) we are removing looks like this:

<img src="chr20_7587045_pretrim.png" width="80%">

This MNV has many softclipped bases on the evidence reads, which present chimeric signatures.

For example for the evidence read with this sequence:

```
TCATAGGCCTGCTTGCCATTTATATGTCTTCTTTGGAGAAATATCTA*TT*TAGATATTTCTCCAAAGAAGACATATAAATGGCAAGCAGGCCTATGAAAAGGTGCTCAACGTTATTAATCATAGGAGAAAAGCAAATCCCCAAACTACAATG
```
The FP variant is denoted by the surrounding asterisks.

The subsections of the chimera map to nearby locations in the genome. As seen in this BLAT result:

```
   ACTIONS      QUERY   SCORE START   END QSIZE IDENTITY  CHROM           STRAND  START       END   SPAN
--------------------------------------------------------------------------------------------------------
browser details YourSeq   109    39   151   151    98.3%  chr20           +     7587036   7587148    113
browser details YourSeq    58     1    58   151   100.0%  chr20           -     7587036   7587093     58
```

And show this inverted repeat signature, that can form self-overlaps. As seen in this mfold result:

```
          10        20        30        40
.-T|                                            A
   CATAGGCCTGCTTGCCATTTATATGTCTTCTTTGGAGAAATATCT T
   GTATCCGGACGAACGGTAAATATACAGAAGAAACCTCTTTATAGA T
\ -^                                            T
       90        80        70        60        50
 
```

We use 10xtrim to identify these reads, and further softclip these reads to remove these 10X-specific artifacts.

## Pre-processing steps (already done)

### Marking duplicates, again

We recommend an additional round of marking duplicates. LongRanger provides the phased BAM file and carries out a barcode-aware markng of duplicates. Reads with missing backcodes may not be missed.

This step can occur before or after 10xtrim.

In the interest of time, we already carried out mark duplicates with the following commands on Picard:

```
# mark duplicates
java -jar [path-to-picard-tools]/MarkDuplicates.jar\
    I=tumour.phased.sorted.bam \
    O=tumour.phased.md.sorted.bam\
    M=tumour.phased.md.metrics.txt"

# index bam
samtools index tumour.phased.md.sorted.bam
```

### Identifying false positives

To identify FP, we ran MuTect1 in matched tumour-normal mode. 

A FP in this case is a variant found tumour-only mode, but not found in matched tumour-normal mode.
 
The following parameters were used:

```
java -jar /u/jpineda/tools/mutect-src/mutect/target/mutect-1.1.7.jar\
     -T MuTect -L chr20\
     -R refdata-hg19-2.1.0/fasta/genome.fa\
     -I:tumor tumour.phased.md.sorted.bam\
     -N:normal normal.phased.md.sorted.bam\
     --vcf tumour_normal.md.vcf\
     -o tumour_normal.md.out\
     --cosmic data/cosmic.hg19.vcf\
     --dbsnp data/dbsnp_138.hg19.vcf\
     --tumor_sample_name HCC1954T\
     --normal_sample_name HCC1954N\
     --normal_panel pon.hg19.mutect1.siteonly.vcf
```

To filter variants, we highly recommend using a panel of normal, high confidence filter (from GIAB), our haplotype discordant filter, and a minimum allele frequency cut off of 0.06. This can be run using our custom script `filter.py`.

```
python filter.py tumour_normal.md.vcf data/tumour.phased.md.sorted.bam data/highconf.hg19.bed HCC1954T
``` 

The output is a filterd VCF file: `tumour_normal.md.filtered.vcf`.

## Set up on OICR cluster

Load modules:

```
module load picard
module load samtools
module load java/1.7.0_21
module load java/jre1.7.0_11
```

## Compute a trimmed modified BAM file

Let's get started! First we will trim the BAM file and then sort the alignments:

```
./10xtrim -b tumour.phased.md.sorted.bam -o tumour.trimmed.stats | samtools view -Sbh | samtools sort > tumour.phased.md.trimmed.sorted.bam
samtools index tumour.phased.md.trimmed.sorted.bam
```

Then we need to fix any mate pairs where 10xtrim completely unmaps an alignment. This may cause inconsistent BAM records.


## Post-processing steps for downstream analyses:

We can use Picard's Fixmateinformation:

```
java -jar [path-to-picard-tools]/FixMateInformation.jar\
    I=tumour.phased.md.trimmed.sorted.bam\
    O=tumour.phased.md.trimmed.fixmates.bam

samtools sort tumour.phased.md.trimmed.fixmates.bam > tumour.phased.md.fixmates.sorted.bam
samtools index tumour.phased.md.trimmed.fixmates.sorted.bam
```

Call somatic mutations in tumour-only mode:
------------------------------------------------------------------------

We can now call our somatic mutations using MuTect1 in tumour-only mode:

```
java -jar /u/jpineda/tools/mutect-src/mutect/target/mutect-1.1.7.jar\
     -T MuTect -L chr20\
     -R refdata-hg19-2.1.0/fasta/genome.fa\
     -I:tumor tumour.phased.md.trimmed.fixmates.sorted.bam\
     --vcf tumour_only.md.trimmed.vcf\
     -o tumour_only.md.trimmed.out
     --cosmic data/cosmic.hg19.vcf\
     --dbsnp data/dbsnp_138.hg19.vcf\
     --tumor_sample_name HCC1954T\
     --normal_sample_name HCC1954N\
     --normal_panel data/pon.hg19.mutect1.siteonly.vcf
```

Visualize difference in IGV:
------------------------------------------------------------------------

To see how 10xtrim removed this 10X-specific artifact we can visualize the false positive in IGV with the BAMs pre and post 10xtrim. In the interest of time, I've generated the IGV screenshot of what we should expect:

<img src="chr20_7587045_posttrim.png" width="80%">

