# Transmorph

### BAM to Transcriptomic Coordinates Converter

This tool converts a BAM file containing read alignments in genomic coordinates into transcriptomic coordinates using a given transcript annotation file (GTF format).

## Requirements

- Python 3.12
- `pandas` library (process gtf file)
- `pysam` library (for working with BAM files)

## Parameters
- `-i` `--input-bamfile` input bam file
- `-g` `--input-gtffile` input gtf file
- `-o` `--output-filename` output bam file

```python
python3 main.py -g gtffile.gtf -i bamfile.bam -o output
```
## Installation

To use the tool, you can simply download the script or clone the repository, and then run the Python script.

```bash
git clone https://github.com/yujiafeng8888/Transmorph.git
cd Transmorph
