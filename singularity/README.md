You can build a [dicey](https://github.com/gear-genomics/dicey) singularity container (SIF file) using

`sudo singularity build dicey.sif dicey.def`

Once you have built the container you can run analysis using

`singularity exec dicey.sif dicey --version`
