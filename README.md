[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


# SAPiN
---------------------------
<img src="img/IRD.png" width="300" height="100" /> <img src="img/MIVEGEC.png" width="150" height="100" />
<h2><em>S</em>ummarize <em>A</em>alignment <em>Pi</em>le by <i>N</i>ucleotide</h2> 

## Table of Contents

   * [Foreword](#foreword)
   * [Install](#install)
   * [Usage](#usage)
     * [Usage](#usage)
   * [Output](#output)     
   * [Acknowledgement](#acknowledgement)

## Foreword

This tool aims to summarize BAM read alignment by pileup or reads at each position in a tabulated way. More convenient as a mpileup format and containing extra information.

## Output

Here an example of output you would get with SAPiN

```
SEQID   POS     REF     QUAL    A       T       G       C       N       INS     DEL     IUPAC   COV     COV_ATGC    MUT_RAT APOBEC  ADAR    REGION  CODON   NUC     DESC
HPV42REF        118     C       38.28   0       16      0       1356    0       0       0       0       1374    1372    1.17    1.2     .       AATGTCAGGTA             CAG     1       gene:ID=gene-1;Name=E6@@mRNA:ID=nbis-rna-1;Parent=gene-1;Name=E6@@exon:ID=nbis-exon-1;Parent=nbis-rna-1;Name=E6@@CDS:ID=cds-1;Parent=nbis-rna-1;Name=E6
```

Here a description of the different fields

| Field | Optional | Type | Description |
| --- | --- | --- | --- |
| SEQID |  | String | The ID of the landmark used to establish the coordinate system for the current feature. |
| POS   |  | Integer | The reference position, with the 1st base having position 1 |
| REF   |  | Character | The reference base. | 
| QUAL  |  | Float | Mean Phred-scaled quality score for the sequenced position. |
| A     |  | Integer | Number of Adenine nucleotide at the position |
| T     |  | Integer | Number of Thymine nucleotide at the position |
| G     |  | Integer | Number of Guanosine nucleotide at the position |
| C     |  | Integer | Number of Cytosine nucleotide at the position |
| N     |  | Integer | Number of Unknown nucleotide at the position |
| INS     |  | Integer | Number of Insertion at the position |
| DEL     |  | Integer | Number of Deletion at the position |
| IUPAC     |  | Integer | Number of IUPAC nucleotide (minus A,T,G,C,N) at the position |
| COV     |  | Integer | Coverage at the position (including INS,DEL,IUPAC) |
| COV_ATGC    |  | Integer | Coverage at the position of A,T,G,C nucleotide only |
| MUT_RAT    |  | Float | Mutation ration (COV_ATGC/nb mutated nuc*100) |
| APOBEC    |  | Float | Mutation ration of C-to-T or G-to-A. Usefull when studying transcriptomes |
| ADAR    |  | Float | Mutation ration of A-to-G or T-to-C. Usefull when studying transcriptomes |
| REGION    |  | STRING | substring of 5 nucleotide on each side. Usefill to make pattern |
| CODON    | Only if GFF provided | STRING | substring of codon in phase/frame (/!\ do not take spliced CDS in account). |
| NUC    | Only if GFF provided | Integer | 1,2 or 3. Indicate in the CODON (previous column) which nucleotide is the one studied at the position |
| DESC    | Only if GFF provided | STRING | feature type and attributes extracted from the gff at the position |

## Install

### Prerequisite

 * python3
 * pysam
 * gffutils
 * matplotlib

 They should be automatically installed during SAPiN installation.

#### Installation with pip:

```bash
pip install git+https://github.com/Juke34/SAPiN.git
```

or if you do not have administrative rights on your machine

```bash
pip install --user git+https://github.com/Juke34/SAPiN.git
```


#### Installation with git:

Clone the repository:

```bash
git clone https://github.com/Juke34/SAPiN.git
```

Move into the folder:

```bash
cd SAPiN/
```

Install:

```bash
python setup.py install
```

or if you do not have administrative rights on your machine:

```bash
python setup.py install --user
```

#### Check installation

Executing:
```bash
sapin
```

or

```bash
sapin -h
```

will display some help.

## Update

#### Update with pip:

```bash
pip install git+https://github.com/Juke34/SAPiN.git --upgrade
```

or if you do not have administartive rights on your machine

```bash
pip install --user git+https://github.com/Juke34/SAPiN.git --upgrade
```

#### Update with git:

Move into the repository folder and execute:

```bash
git pull
cd SAPiN/
python setup.py install
```

## Uninstall

```bash
pip uninstall sapin
```

## Usage

```
sapin -a t/reference.bam -f t/reference.fasta 
```

**advanced:**
```
sapin -a t/reference.bam -f t/reference.fasta -g t/reference_agat.gff3 -cf 1000 -bqf 20 -p
```

## Parameters

| Parameter | Type | Description |
| --- | --- | --- |
|  -a, --ali     | String | Path to the BAM input file |
|  -f, --fasta     | String | Path to the reference fasta file used to align the reads against. |
|  -g, --gff     | String | Optional - Path to the reference gff |
|  -o, --output    | String | Path to the tsv output file |
|  -p, --plot    | Boolean | To plot the ratio of mutation per position (sapin_plot.svg by default. If outpout provided output.svg). |
|  -q, --quiet    | Boolean | "Decrease verbosity |
|  -v, --verbose    | Boolean | Increase verbosity |
|  -z, --gzip    | Boolean | Gzip output file |
|  -s, --shame    | Boolean | Suppress the shameless plug |
|  -cf, --cover_filter    | Integer | filter output to report only site with coverage >= <Integer> |
|  -bqf, --base_quality_filter    | Integer | filter output to report only site with base quality >= <Integer> (default 0) |
|  -mqf, --base_quality_filter    | Integer | filter output to report only site with mapping quality >= <Integer> (default 0) |
|  -mf, --mutation_filter    | Integer | filter output to report only site where the mutation ratio >= <Integer> (default 0) |
