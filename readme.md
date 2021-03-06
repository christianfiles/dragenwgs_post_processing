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

## Install UTRannotator plugin for VEP

```
Download the UTRannotation repos from GitHub (https://github.com/ImperialCardioGenetics/UTRannotator#installation)
Add this to .bash_profile: export PERL5LIB=$PERL5LIB:/export/home/$USER/UTRannotator-master

Run vep with UTRannotator:
vep --cache --dir /share/data/db/human/vep/refseq37_v100/ --fasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta --offline --cache_version 100 --species homo_sapiens --assembly GRCh37 --input_file /path/to/vcf --output_file vcfFile.vcf --vcf --refseq --plugin UTRannotator,uORF_starts_ends_GRCh37_PUBLIC.txt


```

 ### Required Reference Files

params.reference_genome - http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

params.clinvar - https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

params.vep_cache - ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache

params.vep_cache_mt - must be ensembl transcripts - ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/

params.gnotate_gnomad - https://s3.amazonaws.com/slivar/gnomad.hg37.zip

params.gnotate_spliceai - https://s3.amazonaws.com/slivar/spliceai.hg37.zip

params.mitomap_vcf - https://www.mitomap.org/cgi-bin/polymorphisms.cgi?format=vcf

params.gene_panel = config/NexteraDNAFlex/NexteraDNAFlex_virtual_panel.csv

params.whitelist = config/NexteraDNAFlex/NexteraDNAFlex_whitelist.csv

params.whitelist_mito = config/NexteraDNAFlex/NexteraDNAFlex_whitelist_mito.csv

### Variables Files

Sample specific configuration variables are stored in variables files. For an example of this see the test/variables directory.

Required variables are:

workList=test_worklist
sampleId=NA12878i
familyId=FAM001
phenotype=2
sex=1
paternalId=NA12891
maternalId=NA12892


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

To run using docker:
```
nextflow -C config/NexteraDNAFlex/NexteraDNAFlex_local_docker.config run dragenwgs_post_processing.nf -resume
```

Test

```
nextflow -E -C config/NexteraDNAFlex/NexteraDNAFlex_pbs_test.config run dragenwgs_post_processing.nf

```
