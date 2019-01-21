[![Build Status](https://travis-ci.org/gear-genomics/dicey.svg?branch=master)](https://travis-ci.org/gear-genomics/dicey)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/gear-genomics/dicey/master/LICENSE)
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


## Usage

`./src/dicey hunt -g /opt/dev/saint/fm/Danio_rerio.GRCz10.dna.toplevel.fa.gz CATTACTAACATCAGT | python scripts/json2txt.py`

You can also redirect the output to a file:

`./src/dicey hunt -g /opt/dev/saint/fm/Danio_rerio.GRCz10.dna.toplevel.fa.gz -o out.json.gz CATTACTAACATCAGT`

