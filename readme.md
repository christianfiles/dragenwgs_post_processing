# DragenWGS Post Processing

## Introduction

Perform post processing of the raw variant calls from the Dragen server. Works with the DragenWGS pipeline.

Performs the following processes:

- Addional Quality Filtering of Variant Calls
- Prepares VCF for CIP Upload
- Annotates VCF with VEP
- Creates variant reports

## Install

```
git clone https://github.com/AWGL/dragenge_post_processing.git
cd dragenwgs_post_processing
conda env create -f env/dragenwgs_post_processing.yaml

source activate dragenwgs_post_processing


```

Configuration files are located in the config directory. Locations of files such as reference genome will need to be edited to get pipeline to work.

## Run

The pipeline uses nextflow as a workflow manager to make it easy to deploy the pipeline in local and cluster environments.

The starting point for running the pipeline is the results directory with the output of the DragenWGS pipeline.

To run on cluster:

```
nextflow -C config/NexteraDNAFlex/NexteraDNAFlex_pbs.config \
run \
-E \
dragenwgs_post_processing.nf \
--vcf /data/results/dragen_results/191010_D00501_0366_BH5JWHBCX3/NexteraDNAFlex/191010_D00501_0366_BH5JWHBCX3\{.vcf.gz,.vcf.gz.tbi\} \
--variables /data/results/dragen_results/191010_D00501_0366_BH5JWHBCX3/NexteraDNAFlex/\*/\*.variables \
--publish_dir /share/data/results/dragen_temp/191010_D00501_0366_BH5JWHBCX3/NexteraDNAFlex/results \
--sequencing_run 191010_D00501_0366_BH5JWHBCX3 \
-work-dir /share/data/results/dragen_temp/191010_D00501_0366_BH5JWHBCX3/NexteraDNAFlex/work
```

To run locally:
```
nextflow -C config/NexteraDNAFlex/NexteraDNAFlex_local.config run dragenwgs_post_processing.nf -resume
```

