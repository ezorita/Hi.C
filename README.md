# Hi.C a Hi-C pipeline written in C

## Index
1. Download and compilation
2. Usage
3. Example

## 1. Download and compilation

Clone from this repository using `git`. Compile directly using make:
```bash
$ git clone http://github.com/ezorita/Hi.C.git
$ cd Hi.C
$ make
```

This will generate three binaries:
- `re_digest`: in-silico digestion of genomes using defined restriction enzymes.
- `parse_contacts`: reads mapped files and finds valid Hi-C contact pairs.
- `merge_contacts`: simplifies the output files of `parse_contacts`.

## 2. Usage

### 2.1. Mapping

For mapping, we use [bwa](https://github.com/lh3/bwa) with the following parameters (for hg19, Illumina 75bp PE):

```bash
$ bwa mem -P -k17 -U0 -L0,0 -T25 genome_index.fasta HiC-read1.fastq.gz HiC-read2.fastq.gz | samtools view -bS > HiC-mapped.bam
```

### 2.2. Digesting the genome

Before finding the contacts we need to precompute the fragments produced by the restriction enzyme used in Hi-C. To do so, run `re_digest` as follows:

```bash
$ re_digest [organism name] [RE name] [RE sequence] [cut fw] [cut rv]
```

All arguments are mandatory:
- **organism**: The name of the organism. You must create manually the organism in the db before performing the digestion. To do so, create a directory in the same path as `re_digest` called `db`. Inside `db` create another directory with the name of the organism (e.g. `hg`) and then place the genome of the organism in fasta format in a file called `genome.fasta` (yes, symbolic links are allowed). See the example below for more info.
- **RE name**: The name of the restriction enzyme, e.g. *MboI, DpnII, EcoRI*...
- **RE sequence**: The target sequence of the restriction enzyme, e.g. `GATC` (MboI, DpnII), GAATTC (EcoRI)...
- **cut fw/cut rv**: Cutting sites at the forward/reverse strand, from the first 5' nucleotide of the RE sequence. For instance, EcoRI produces a cut such that G'AATT,C (AATT overhang), so cut fw=1, cut rv=5.

### 2.3. Finding contacts

#### Usage

To find the contacts of your Hi-C experiment, run `parse_contacts`:

```
$ parse_contacts [organism] [RE name] [HiC-mapped.sam] [[mapq]] [[insert size]]
```

Mandatory arguments:
- **organism**: The organism as described during the digestion.
- **RE name**: The name of the restriction enzyme used in the experiment (must have been previously digested, see above).
- **HiC-mapped.sam**: The file containing the output of bwa mapping. Use `<(samtools view HiC-mapped.bam)` instead if you used .bam compression during the mapping process.

Optional arguments:
- **mapq**: The minimum mapping quality of the mapped fragments (default is 20).
- **insert size**: The maximum insert size (in bp) of the mapping technology (default 2000, Illumina).

#### Output

Running `parse_contacts` will print in the standard output all the valid contact pairs found in the mapping file using [TADbit](https://github.com/3DGenomes/TADbit) format:

| Column   | Read | Description                  |
| :-------:|:----:|------------------------------|
| 1        | both | Read name                    |
| 2        | 1    | Chromosome                   |
| 3        | 1    | Mapping locus                |
| 4        | 1    | Strand (0:forward/1:reverse) |
| 5        | 1    | Length mapped                |
| 6        | 1    | Upstream RE site position    |
| 7        | 1    | Downstream RE site position  |
| 8        | 2    | Chromosome                   |
| 9        | 2    | Mapping locus                |
| 10       | 2    | Strand (0:forward/1:reverse) |
| 11       | 2    | Length mapped                |
| 12       | 2    | Upstream RE site position    |
| 13       | 2    | Downstream RE site position  |


### 2.4. Compacting the contacts

#### Sorting the output file
In this step we will merge the contact pairs into pairs of RE fragments. However, first we need to sort the output file of `parse_contacts` with the following keys (GNU sort) `-k2,2 -k8,8 -k6,6n -k12,12n`, e.g.:

```bash
$ sort -k2,2 -k8,8 -k6,6n -k12,12n parse_contacts.out > parse_contacts_sorted.out
```

#### Usage

After sorting we can safely merge the contacts using `merge_contacts`. This scripts merges the contacts of the same restriction enzyme fragments and removes potential PCR duplicates.

```bash
$ merge_contacts [parse_contacts_sorted.out]
```

Mandatory arguments:
- **parse_contacts_sorted.out**: A file generated with `parse_contacts` and sorted with GNU sort.

#### Output

The output produced by `merge_contacts` is similar to a bed file. The first 6 columns are the restriction enzyme fragments of the two contacting sequences and the last column is the contact count:

| Column   | Description                             |
| :-------:|-----------------------------------------|
| 1        | Chromosome of fragment 1                |
| 2        | First nucleotide of fragment 1          |
| 3        | Last nucleotide of fragment 1           |
| 4        | Chromosome of fragment 2                |
| 5        | First nucleotide of fragment 2          |
| 6        | Last nucleotide of fragment 2           |
| 7        | Contact count between fragments 1 and 2 |

## 3. Example

### Introduction

> In this example we will process Hi-C data from human cells. The restriction enzyme used during the experiment is MboI. The DNA was sequenced in an Illumina platform using 75nt paired-end reads.

### Files and paths
> Assume that we have our precious paired-end reads in two files called `hic-read1.fastq.gz` and `hic-read2.fastq.gz` in the same folder as our compiled binaries `re_digest`, `parse_contacts` and `merge_contacts`.
> Also assume that we have the human genome reference file and its bwa index in `/genomes/hg/genome.fasta`.

### Processing steps
#### 1. Make a digestion reference of MboI for human genome.
```bash
# Create the digestion database.
$ mkdir -p db/hg
$ ln -s /genomes/hg/genome.fasta db/hg/genome.fasta

# Digest the human genome with MboI (GATC, 5' overhang: GATC)
$ ./re_digest hg MboI GATC 0 4
```

#### 2. Map the sequencing reads using `bwa`:
```bash
$ bwa mem -P -k17 -U0 -L0,0 -T25 /genomes/hg/genome.fasta hic-read1.fastq.gz hic-read2.fastq.gz | samtools view -bS > hic-mapped.bam
```

#### 3. Parse contacts from mapping file and sort the output:
```bash
$ ./parse_contacts hg MboI <(samtools view hic-mapped.bam) | sort -k2,2 -k8,8 -k6,6n -k12,12n > contacts_sorted.out
```

#### 4. Merge contacts:
```bash
$ ./merge_contacts contacts_sorted.out > fragment_contacts.out
```
