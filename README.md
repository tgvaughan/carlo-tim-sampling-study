Simulation study on the effects of varying sampling rates on the
inference of multi-type birth-death model parameters.

This directory contains a nextflow pipeline for simulating 5
independent trees for each of the 16 different conditions listed in
the file `sim_params.org`, then inferring the parameters of this
simulation using both a single- and a two-type phylodynamic analysis.

Running this pipeline requires the following software:

1. R, together with the packages `tidyverse` and `bioseq`
2. OpenJDK 1.8 or later
3. Nextflow

To run the pipeline, use:
```
        nextflow run pipeline.nf
```