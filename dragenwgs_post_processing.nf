#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for germline variant calling
========================================================================================

Author: Joseph Halstead

========================================================================================
Define initial files
========================================================================================
*/

reference_genome = file(params.reference_genome)
gnomad_exomes = file(params.gnomad_exomes)
gnomad_genomes = file(params.gnomad_genomes)
vep_cache = file(params.vep_cache)
gene_panel = file(params.gene_panel)
whitelist = file(params.whitelist)
clinvar = file(params.clinvar)

/*
========================================================================================
Define initial channels
========================================================================================
*/

Channel
  .fromPath(params.variables)
  .ifEmpty { exit 1, "Cannot find any variables files matching: ${params.variables}" }
  .set{ variables_channel }

Channel
  .fromFilePairs(params.vcf, flat: true)
  .ifEmpty { exit 1, "Cannot find a vcf file files matching: ${params.vcf}" }
  .set { raw_vcf }

raw_vcf.into{
    raw_vcf_annotation
    raw_vcf_relatedness
    raw_vcf_cip
    raw_vcf_sample_names
}


/*
========================================================================================
Main pipeline
========================================================================================
*/

// create ped file
process create_ped {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/ped/", mode: 'copy'

    input:
    file(variables) from variables_channel.collect()

    output:
    file("${params.sequencing_run}.ped") into ped_channel

    """
    create_ped.py --variables "*.variables" > ${params.sequencing_run}.ped

    """

}

ped_channel.into{
    ped_channel_json
    ped_channel_reports
}

// create json for qiagen upload
process create_json_for_qiagen {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/qiagen_json/", mode: 'copy'

    input:
    file(ped) from ped_channel_json

    output:
    file('*.json') 

    """
    convert_ped_to_json.py --pedfile $ped --seqid $params.sequencing_run --output ./

    """

}



// Extract sample names from vcf
process get_sample_names_from_vcf{

    input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf_sample_names

    output:
    file("${params.sequencing_run}_sample_names.txt") into sample_names_from_vcf

    """
    bcftools query -l $vcf > ${params.sequencing_run}_sample_names.txt
    """

}


// Create a channel with sample names as the values
sample_names_from_vcf.splitCsv(header:['col1']).map{ row-> tuple(row.col1)}.set { samples_ch }


samples_ch.into{
    samples_ch_splitting
    samples_ch_reporting
}


// split multiallelics and normalise
process split_multiallelics_and_normalise{

    cpus params.vcf_processing_cpus

	input:
	set val(id), file(vcf), file(vcf_index) from raw_vcf_annotation 

	output:
	set val(id), file("${params.sequencing_run}.norm.vcf.gz"), file("${params.sequencing_run}.norm.vcf.gz.tbi") into normalised_vcf_channel

	"""
	zcat  < $vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}.norm.vcf
	bgzip ${params.sequencing_run}.norm.vcf
	tabix ${params.sequencing_run}.norm.vcf.gz

	"""

}


normalised_vcf_channel.into{
    for_vep_channel
    qiagen_vcf_channel
    qiagen_prefiltered_vcf_channel
}



// Annotate using VEP
process annotate_with_vep{

    cpus params.vep_cpus

    publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

    input:
    set val(id), file(normalised_vcf), file(normalised_vcf_index) from for_vep_channel

    output:
    set val(id), file("${params.sequencing_run}.norm.anno.vcf.gz"), file("${params.sequencing_run}.norm.anno.vcf.gz.tbi")  into annotated_vcf

    """
    vep \
    --verbose \
    --format vcf \
    --everything \
    --fork $params.vep_cpus \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file $normalised_vcf \
    --output_file ${params.sequencing_run}.norm.anno.vcf \
    --force_overwrite \
    --cache \
    --dir  $vep_cache \
    --fasta $reference_genome \
    --offline \
    --cache_version 94 \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq \
    --flag_pick \
    --pick_order biotype,canonical,appris,tsl,ccds,rank,length \
    --exclude_predicted \
    --custom ${gnomad_genomes},gnomADg,vcf,exact,0,AF_POPMAX \
    --custom ${gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX \
    --custom ${clinvar},clinvarSIG,vcf,exact,0,CLNSIG,CLNSIGCONF

    bgzip ${params.sequencing_run}.norm.anno.vcf
    tabix ${params.sequencing_run}.norm.anno.vcf.gz

    """
}


process create_variant_reports {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/variant_reports/", mode: 'copy'

    input:
    set val(id), file(vcf), file(vcf_index) from annotated_vcf
    each sample_names from samples_ch_reporting
    file(ped) from ped_channel_reports

    output:
    file("${params.sequencing_run}.${sample_names[0]}_variant_report.csv") optional true into variant_report_channel

    """
    germline_variant_reporter.py \
    --vcf $vcf \
    --proband_id ${sample_names[0]} \
    --ped $ped \
    --gene_list $gene_panel \
    --gnomad_ad $params.gnomad_ad \
    --gnomad_r $params.gnomad_r \
    --minqual_snps $params.snp_qual \
    --minqual_indels $params.indel_qual \
    --min_dp $params.min_dp \
    --min_gq $params.min_gq \
    --min_af_mt $params.min_mt_af \
    --output ${params.sequencing_run}.${sample_names[0]}_variant_report.csv \
    --worklist $params.worklist_id \
    --whitelist $whitelist
    """

}

// Calculate relatedness between samples
process calculate_relatedness {

    cpus params.relatedness_cpus

    publishDir "${params.publish_dir}/relatedness/", mode: 'copy'

    input:
    set val(id), file(vcf), file(vcf_index) from raw_vcf_relatedness

    output:
    file("${params.sequencing_run}.relatedness2")

    """
    vcftools --relatedness2 \
    --out $params.sequencing_run \
    --gzvcf $vcf

    """

}



process split_multisample_vcfs {

    cpus params.small_task_cpus

    input:
    set val(id), file(vcf), file(vcf_index) from qiagen_vcf_channel
    each sample_names from samples_ch_splitting

    output:
    set val(id), file("${params.sequencing_run}.${sample_names[0]}.vcf.gz"), file("${params.sequencing_run}.${sample_names[0]}.vcf.gz.tbi") into per_sample_vcf_channel

    """
    bcftools view -I -s ${sample_names[0]} $vcf > ${params.sequencing_run}.${sample_names[0]}.vcf

    bgzip ${params.sequencing_run}.${sample_names[0]}.vcf
    tabix ${params.sequencing_run}.${sample_names[0]}.vcf.gz
    """

}


per_sample_vcf_channel.into{
    whole_per_sample_vcf_channel
    filtered_per_sample_vcf_channel

}


process filter_single_sample_vcfs_for_qiagen {

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/qiagen_vcfs/", mode: 'copy'

    input:
    set val(id), file(vcf), file(vcf_index) from whole_per_sample_vcf_channel

    output:
    set file("${params.sequencing_run}.${sample_id}.qual.vcf.gz"), file("${params.sequencing_run}.${sample_id}.qual.vcf.gz.tbi")

    script:
    sample_id = vcf.name.toString().tokenize('.').get(1)

    """
    quality_filter.py --vcf $vcf \
    --min_mt_af $params.min_mt_af \
    --snp_qual $params.snp_qual \
    --indel_qual $params.indel_qual \
    --min_dp $params.min_dp \
    --min_gq $params.min_gq \
    --sampleid $sample_id > ${params.sequencing_run}.${sample_id}.qual.vcf

    bgzip ${params.sequencing_run}.${sample_id}.qual.vcf
    tabix ${params.sequencing_run}.${sample_id}.qual.vcf.gz

    """

}

process create_pre_filtered_vcfs_from_variant_reports{

    cpus params.small_task_cpus

    publishDir "${params.publish_dir}/qiagen_prefiltered_vcfs/", mode: 'copy'


    input:
    set val(id), file(normalised_vcf), file(normalised_vcf_index) from qiagen_prefiltered_vcf_channel
    each file(variant_report) from variant_report_channel


    output:
    file "${id}*filtered.vcf.gz"

    """
    sampleids=\$(get_sampleids_from_csv.py --csv $variant_report)


    for i in \$sampleids; do

        bcftools view -I -s \$i $normalised_vcf > ${params.sequencing_run}.\$i.vcf

        bgzip ${params.sequencing_run}.\$i.vcf
        tabix ${params.sequencing_run}.\$i.vcf.gz

        convert_report_to_vcf.py \
        --vcf ${params.sequencing_run}.\$i.vcf.gz \
        --csv $variant_report \
        --sample_id \$i > ${id}.\$i.filtered.vcf

        bgzip ${id}.\$i.filtered.vcf


    done

    """







}






workflow.onComplete{

    if (workflow.success){

        ran_ok = "${params.sequencing_run} success!.\n"
    }
    else{

        ran_ok = "${params.sequencing_run} fail!.\n"

    }

    def newFile = new File("${params.publish_dir}/post_processing_finished.txt")

    newFile.createNewFile()
    newFile.append(ran_ok)

}