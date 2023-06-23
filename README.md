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

This tool aims to summarize SAM/BAM read alignment by pileup or reads at each position in a tabulated way. More convenient as a mpileup format and containing extra information.

## Output

Here an example of output you would get with SAPiN

```
##COLUMN_INFO=<ID=SEQID,Type=String,Description="The ID of the landmark used to establish the coordinate system for the current feature.">
##COLUMN_INFO=<ID=POS,Type=Integer,Description="The reference position, with the 1st base having position 1.">
##COLUMN_INFO=<ID=REF,Type=String,Description="The reference base. Each base must be one of A,C,G,T,N (case insensitive).">
##COLUMN_INFO=<ID=QUAL,Type=Integer,Description="Quality: Mean Phred-scaled quality score for the sequenced position.">
##COLUMN_INFO=<ID=A,Type=Integer,Description="Number of Adenine nucleotide at the position">
##COLUMN_INFO=<ID=T,Type=Integer,Description="Number of Thymine nucleotide at the position">
##COLUMN_INFO=<ID=G,Type=Integer,Description="Number of Guanosine nucleotide at the position">
##COLUMN_INFO=<ID=C,Type=Integer,Description="Number of Cytosine nucleotide at the position">
##COLUMN_INFO=<ID=N,Type=Integer,Description="Number of Unknown nucleotide at the position">
##COLUMN_INFO=<ID=INS,Type=Integer,Description="Number of Insertion at the position">
##COLUMN_INFO=<ID=DEL,Type=Integer,Description="Number of Deletion at the position">
##COLUMN_INFO=<ID=IUPAC,Type=Integer,Description="Number of IUPAC nucleotide at the position">
##COLUMN_INFO=<ID=COV,Type=Integer,Description="Coverage at the position">
SEQID   POS REF QUAL    A   T   G   C   N   INS DEL IUPAC   COV
```
## Install

### Prerequisite

 * python3
 * pysam

## Usage

```
python3 -m SAPiN -i t/bwamem2_sorted.bam 
```

## Parameters

| Parameter | Type | Description |
| --- | --- | --- |
|  -i, --input     | String | Path to the SAM/BAM input file |
|  -o, --output    | String | Path to the tsv output file |
|  -q, --quiet    | Boolean | "Decrease verbosity |
|  -v, --verbose    | Boolean | Increase verbosity |
|  -z, --gzip    | Boolean | Gzip output file |
|  -s, --shame    | Boolean | Suppress the shameless plug |
|  -cf, --cover_filter    | Integer | filter output to report only site with coverage >= <Integer> |
|  -qf, --quality_filter    | Integer | filter output to report only site with quality >= <Integer> |
