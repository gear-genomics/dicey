[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/dicey)
[![Build Status](https://travis-ci.org/gear-genomics/dicey.svg?branch=master)](https://travis-ci.org/gear-genomics/dicey)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/dicey/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/dicey.svg)](https://github.com/gear-genomics/dicey/releases)
[![GitHub Issues](https://img.shields.io/github/issues/gear-genomics/dicey.svg)](https://github.com/gear-genomics/dicey/issues)

# Installing dicey

The easiest way to get Dicey is to download a statically linked binary from the [Dicey release page](https://github.com/gear-genomics/dicey/releases) or to download Dicey from [Bioconda](https://anaconda.org/bioconda/dicey). Building from source is also possible:

`apt-get install -y build-essential g++ cmake zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev`

`git clone --recursive https://github.com/gear-genomics/dicey.git`

`cd dicey/`

`make all`

`make install`

This will generate the binary `bin/dicey`.


## Sequence search in an indexed reference genome

Searching a large reference genome requires a pre-built index on a bgzip compressed genome.

`dicey index -o hg19.fa.fm9 hg19.fa.gz`

`samtools faidx hg19.fa.gz`

The indexing step is only required once. You can then search nucleotide sequences at a user-defined edit or hamming distance.

`dicey hunt -g hg19.fa.gz TCTCTGCACACACGTTGT | python scripts/json2txt.py`

You can also redirect the output in JSON format to a file.

`dicey hunt -g hg19.fa.gz -o out.json.gz TCTCTGCACACACGTTGT`


## In-silico PCR for a set of primers

Dicey can search for multiple primer pairs, show off-target products and determine PCR amplicons.

`echo -e ">FGA_f\nGCCCCATAGGTTTTGAACTCA\n>FGA_r\nTGATTTGTCTGTAATTGCCAGC" > primers.fa`

`dicey search -g hg19.fa.gz primers.fa | python scripts/json2txt.py`

The default output is a JSON file that can also be stored in a file.

`dicey search -o out.json.gz -g hg19.fa.gz primers.fa`


## FAQ

* Dicey cannot find the primer3 config directory    
The primer3 config directory is included in the repository. Just clone the repository `git clone https://github.com/gear-genomics/dicey.git` and then use the cloned config directory `dicey search -i dicey/src/primer3_config/ -g hg19.fa.gz primers.fa`.

* The script json2txt.py is not found     
The json2txt.py python script is included in the repository. Just clone the repository `git clone https://github.com/gear-genomics/dicey.git` and then you will find the script in the `./scripts/` subdirectory.


## Questions

In case of questions feel free to send us an [email](https://www-db.embl.de/EMBLPersonGroup-PersonPicture/MailForm/?recipient=ggenomics).
