[![Build Status](https://travis-ci.org/gear-genomics/dicey.svg?branch=master)](https://travis-ci.org/gear-genomics/dicey)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/dicey/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/dicey.svg)](https://github.com/gear-genomics/dicey/releases)
[![GitHub Issues](https://img.shields.io/github/issues/gear-genomics/dicey.svg)](https://github.com/gear-genomics/dicey/issues)

# dicey

In-silico PCR and variant primer design


## Installation from source

```bash
git clone --recursive https://github.com/gear-genomics/dicey.git
cd dicey/
make
```

This will generate the binary `src/dicey`.


## Sequence search in an indexed reference genome

Searching a large reference genome requires a pre-built index on a bgzip compressed genome.

`./src/dicey index -o hg19.fa.fm9 hg19.fa.gz`

`samtools faidx hg19.fa.gz`

The indexing step is only required once. You can then search nucleotide sequences at a user-defined edit or hamming distance.

`./src/dicey hunt -g hg19.fa.gz TCTCTGCACACACGTTGT | python scripts/json2txt.py`

You can also redirect the output in JSON format to a file.

`./src/dicey hunt -g hg19.fa.gz -o out.json.gz TCTCTGCACACACGTTGT`


