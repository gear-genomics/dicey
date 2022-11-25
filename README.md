[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/dicey/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/dicey/badges/downloads.svg)](https://anaconda.org/bioconda/dicey)
[![C/C++ CI](https://github.com/gear-genomics/dicey/workflows/C/C++%20CI/badge.svg)](https://github.com/gear-genomics/dicey/actions)
[![Docker CI](https://github.com/gear-genomics/dicey/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/geargenomics/dicey/)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/gear-genomics/dicey/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/gear-genomics/dicey.svg)](https://github.com/gear-genomics/dicey/releases)

## Installing dicey

Dicey is available as a [Bioconda package](https://anaconda.org/bioconda/dicey), as a pre-compiled statically linked binary from [Dicey's github release page](https://github.com/gear-genomics/dicey/releases), as a singularity container [SIF file](https://github.com/gear-genomics/dicey/releases) or as a minimal [Docker container](https://hub.docker.com/r/geargenomics/dicey/).

`apt-get install -y build-essential g++ cmake zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev`

`git clone --recursive https://github.com/gear-genomics/dicey.git`

`cd dicey/`

`make all`

`make install`

This will generate the binary `bin/dicey`.


## Running Dicey

`dicey -h`

## Sequence search in an indexed reference genome

Searching a large reference genome requires a pre-built index on a bgzip compressed genome.

`dicey index -o hg19.fa.fm9 hg19.fa.gz`

`samtools faidx hg19.fa.gz`

The indexing step is only required once. You can then search nucleotide sequences at a user-defined edit or hamming distance.

`dicey hunt -g hg19.fa.gz TCTCTGCACACACGTTGT | python scripts/json2txt.py`

You can also redirect the output in JSON format to a file.

`dicey hunt -g hg19.fa.gz -o out.json.gz TCTCTGCACACACGTTGT`

Pre-built genome indices for commonly used reference genomes are available for [download here](https://gear.embl.de/data/tracy/).


## In-silico PCR for a set of primers

Dicey can search for multiple primer pairs, show off-target products and determine PCR amplicons.

`echo -e ">FGA_f\nGCCCCATAGGTTTTGAACTCA\n>FGA_r\nTGATTTGTCTGTAATTGCCAGC" > primers.fa`

`dicey search -c 45 -g hg19.fa.gz primers.fa | python scripts/json2txt.py`

The default output is a JSON file that can also be stored in a file.

`dicey search -c 45 -o out.json.gz -g hg19.fa.gz primers.fa`


## Padlock probe design

Dicey can design padlock probes for imaging mRNA in single cells. You need to download an [indexed reference genome](https://gear.embl.de/data/tracy/) and a matching GTF file, e.g., for GRCh38:

`wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz`

With these files, you can then design padlock probes for a given gene using

`dicey padlock -g Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -t Homo_sapiens.GRCh38.107.gtf.gz -b data/bar.fa.gz ENSG00000136997`


## Graphical user interface

You can search primers interactively using our web application [silica](https://www.gear-genomics.com/silica/).


## FAQ

* Dicey cannot find the primer3 config directory    
The primer3 config directory is included in the repository. Just clone the repository `git clone --recursive https://github.com/gear-genomics/dicey.git` and then use the cloned config directory `dicey search -i dicey/src/primer3_config/ -g hg19.fa.gz primers.fa`.

* The script `json2txt.py` is not found     
The `json2txt.py` python script is included in the repository. Just clone the repository `git clone --recursive https://github.com/gear-genomics/dicey.git` and then you will find the script in the `./scripts/` subdirectory.


## Citation

Dicey is part of the [GEAR genomics framework](https://www.gear-genomics.com) which is described in the below publication.

Rausch, T., Fritz, M.H., Untergasser, A. and Benes, V.         
Tracy: basecalling, alignment, assembly and deconvolution of sanger chromatogram trace files.             
BMC Genomics 21, 230 (2020).             
[https://doi.org/10.1186/s12864-020-6635-8](https://doi.org/10.1186/s12864-020-6635-8)


## License

Dicey is distributed under the GPL license. Consult the accompanying [LICENSE](https://github.com/gear-genomics/dicey/blob/master/LICENSE) file for more details.


## Questions

In case of questions feel free to send us an [email](https://www-db.embl.de/EMBLPersonGroup-PersonPicture/MailForm/?recipient=ggenomics).
