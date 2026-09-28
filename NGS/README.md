# NGS Module 1

These scripts were used for the NGS Module 1 practical work. The main purpose is to perform basic quality checking and preprocessing of paired-end FASTQ files.

## Tools Used

* FastQC
* Trimmomatic
* Cutadapt
* Bash
* Java

## Scripts

### `ngs_module1_pipeline.sh`

This script is used when there are multiple datasets. It finds the `fastq_seq` folder for each dataset and processes them separately.

### `ngs_module1_individual.sh`

This is a simpler version used to process one dataset at a time.

## Steps

The scripts follow these main steps:

```text id="4n4x9j"
FASTQ files
    ↓
FastQC
    ↓
Trimmomatic
    ↓
FastQC
    ↓
Cutadapt
    ↓
FastQC
```

FastQC is used to check the quality of the reads before and after trimming.

Trimmomatic is used for adapter and quality trimming.

Cutadapt is used separately for adapter trimming as part of the practical.

## Input

Paired-end FASTQ files are required.

Example:

```text id="r1f6do"
sample_1.fastq
sample_2.fastq
```

For the multiple dataset script, the expected structure is:

```text id="a7c9ne"
data/
├── person1/
│   └── fastq_seq/
│       ├── sample_1.fastq
│       └── sample_2.fastq
└── person2/
    └── fastq_seq/
        ├── sample_1.fastq
        └── sample_2.fastq
```

## Trimmomatic Settings

The parameters used in the practical were:

```text id="6m3z8p"
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
SLIDINGWINDOW:4:20
MINLEN:36
```

## Output

Separate output folders are created for:

* FastQC before trimming
* Trimmomatic paired reads
* Trimmomatic unpaired reads
* FastQC after Trimmomatic
* Cutadapt trimmed reads
* FastQC after Cutadapt

The original FASTQ files and generated results are not included in this repository.

These scripts were written and used during my bioinformatics practical/project work.
