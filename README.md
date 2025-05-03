# Transmorph

**Transmorph** is a lightweight Python tool to convert BAM files from genomic coordinates to transcriptomic coordinates using a provided transcript annotation file (GTF/GFF).

---

## Requirements

- Python 3.12
- `pandas` library (process gtf file)
- `pysam` library (for working with BAM files)

## Parameters
- `-i` `--input-bamfile` input bam file
- `-g` `--input-gtffile` input gtf file
- `-o` `--output-filename` output bam file

## Installation

From GitHub:

```bash
pip install git+https://github.com/yujiafeng8888/Transmorph.git
```
Or clone and install locally:

```bash
git clone https://github.com/yujiafeng8888/Transmorph.git
cd Transmorph
pip install .
```

## Usage

After installation, you can run it via command line:
```bash
transmorph -i input.bam -g annotation.gtf -o output_transcript.bam

```

## Author

Yujia Feng
yfeng80@jh.edu
