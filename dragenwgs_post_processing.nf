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
reference_genome_index = file(params.reference_genome_index)
vep_cache = file(params.vep_cache)
gene_panel = file(params.gene_panel)
whitelist = file(params.whitelist)
whitelist_mito = file(params.whitelist_mito)
clinvar = file(params.clinvar)
vep_cache_mt = file(params.vep_cache_mt)
mitomap = file(params.mitomap_vcf)
gnotate_gnomad = file(params.gnotate_gnomad)
gnotate_spliceai = file(params.gnotate_spliceai)

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


variables_channel.into{
	variables_channel_ped
	variables_channel_variant_report
	variables_channel_mito_variant_report
}


// Chromosomes for when we do VEP in parallel
chromosome_ch = Channel.from('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' )


/*
========================================================================================
Main pipeline
========================================================================================
*/

// Create PED file
process create_ped {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/ped/", mode: 'copy'

	input:
	file(variables) from variables_channel_ped.collect()

	output:
	file("${params.sequencing_run}.ped") into ped_channel

	"""
	create_ped.py --variables "*.variables" > ${params.sequencing_run}.ped

	"""

}

ped_channel.into{
	ped_channel_json
	ped_channel_reports
	ped_channel_reports_mito
	ped_channel_upd
}

// Create json for qiagen upload
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

	cpus params.small_task_cpus

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

// Split this into multiple channels
samples_ch.into{
	samples_ch_splitting
	samples_ch_reporting
	samples_ch_reporting_mito
	samples_ch_upd
}

// Split multiallelics and normalise also annotate with spliceai + gnomad
process split_multiallelics_and_normalise{

	cpus params.vcf_processing_cpus

	input:
	set val(id), file(vcf), file(vcf_index) from raw_vcf_annotation 
	file(reference_genome_index)
	file(gnotate_spliceai)

	output:
	set val(id), file("${params.sequencing_run}.norm.vcf.gz"), file("${params.sequencing_run}.norm.vcf.gz.tbi") into normalised_vcf_channel

	"""
	zcat  < $vcf | vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}.norm.inter.vcf
	bgzip ${params.sequencing_run}.norm.inter.vcf
	tabix ${params.sequencing_run}.norm.inter.vcf.gz

	slivar expr --gnotate $gnotate_spliceai -o ${params.sequencing_run}.norm.inter2.vcf -v ${params.sequencing_run}.norm.inter.vcf.gz

	bgzip ${params.sequencing_run}.norm.inter2.vcf
	tabix ${params.sequencing_run}.norm.inter2.vcf.gz

	slivar expr --gnotate $gnotate_gnomad -o ${params.sequencing_run}.norm.vcf -v ${params.sequencing_run}.norm.inter2.vcf.gz

	bgzip ${params.sequencing_run}.norm.vcf
	tabix ${params.sequencing_run}.norm.vcf.gz

	"""

}

normalised_vcf_channel.into{
	for_splitting_channel
	qiagen_vcf_channel
	for_mito_vcf_channel
	for_upd_channel
}

// Split the normalised vcf up per chromosome for parallel processing
process split_vcf_per_chromosome{

	cpus params.vcf_processing_cpus

	input:
	set val(id), file(normalised_vcf), file(normalised_vcf_index) from for_splitting_channel
	each chromosome from chromosome_ch

	output:
	set val(chromosome), file("${params.sequencing_run}.norm.${chromosome}.vcf.gz"), file("${params.sequencing_run}.norm.${chromosome}.vcf.gz.tbi")  into for_vep_channel

	"""
	bcftools view -r $chromosome $normalised_vcf > ${params.sequencing_run}.norm.${chromosome}.vcf

	bgzip ${params.sequencing_run}.norm.${chromosome}.vcf
	tabix ${params.sequencing_run}.norm.${chromosome}.vcf.gz

	"""
}

// Annotate using VEP 
process annotate_with_vep_and_gnomad{

	cpus params.vep_cpus

	input:
	set val(chromosome), file(normalised_vcf), file(normalised_vcf_index) from for_vep_channel
	file reference_genome_index
	file vep_cache

	output:
	file("${params.sequencing_run}.norm.${chromosome}.anno.vcf") into annotated_vcf_per_chromosome

	"""
	vep \
	--verbose \
	--format vcf \
	--hgvs \
	--symbol \
	--numbers \
	--domains \
	--regulatory \
	--canonical \
	--protein \
	--biotype \
	--uniprot \
	--tsl \
	--appris \
	--gene \
	--variant_class \
	--check_existing \
	--fork $params.vep_cpus \
	--species homo_sapiens \
	--clin_sig_allele 0 \
	--assembly GRCh37 \
	--input_file $normalised_vcf \
	--output_file ${params.sequencing_run}.norm.${chromosome}.anno.vcf \
	--force_overwrite \
	--cache \
	--dir  $vep_cache \
	--fasta $reference_genome \
	--offline \
	--cache_version $params.vepversion \
	--no_escape \
	--shift_hgvs 1 \
	--vcf \
	--refseq \
	--flag_pick \
	--pick_order biotype,canonical,appris,tsl,ccds,rank,length \
	--exclude_predicted \
	--custom ${clinvar},clinvar,vcf,exact,0,CLNSIG,CLNSIGCONF

	"""
}

// Merge per chromosome vcfs back into one vcf
process merge_annotated_vcfs{

	publishDir "${params.publish_dir}/annotated_vcf/", mode: 'copy'

	input:
	file(vcfs) from annotated_vcf_per_chromosome.collect()

	output:
	set file("${params.sequencing_run}.norm.anno.vcf.gz"), file("${params.sequencing_run}.norm.anno.vcf.gz.tbi") into annotated_vcf

	"""
	bcftools concat ${vcfs.collect { "$it " }.join()} > ${params.sequencing_run}.norm.anno.unsorted.vcf

	bcftools sort -T ./ ${params.sequencing_run}.norm.anno.unsorted.vcf > ${params.sequencing_run}.norm.anno.vcf

	bgzip ${params.sequencing_run}.norm.anno.vcf
	tabix ${params.sequencing_run}.norm.anno.vcf.gz

	"""
}

// Create variant reports csvs for non mitochondrial variants
process create_variant_reports {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/variant_reports/", mode: 'copy'

	input:
	set file(vcf), file(vcf_index) from annotated_vcf
	each sample_names from samples_ch_reporting_mito
	file(ped) from ped_channel_reports
	file(variables) from variables_channel_variant_report.collect()

	output:
	file("${params.sequencing_run}.${sample_names[0]}_variant_report.csv") into variant_report_channel

	"""

	. ${sample_names[0]}.variables

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
	--max_parental_alt_ref_ratio $params.max_parental_alt_ref_ratio \
	--output ${params.sequencing_run}.${sample_names[0]}_variant_report.csv \
	--worklist \$worklistId \
	--whitelist $whitelist \
	--splice_ai $params.splice_ai_cutoff \
	--apply_panel 
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

// Split multisample VCF up per sample VCFs as QCII requires in this way
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

// Apply quality filtering script to per sample VCF
process filter_single_sample_vcfs_for_qiagen {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/qiagen_vcfs/", mode: 'copy'

	input:
	set val(id), file(vcf), file(vcf_index) from per_sample_vcf_channel

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

// Extract mitochrondrial variants and annotate with VEP - use ensembl transcripts
process get_mitochondrial_variant_and_annotate{

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/annotated_mitochrondrial_vcf/", mode: 'copy'

	input:
	set val(id), file(normalised_vcf), file(normalised_vcf_index) from for_mito_vcf_channel
	file reference_genome_index
	file vep_cache_mt


	output:
	set val(id), file("${params.sequencing_run}.mito.anno.vcf.gz"), file("${params.sequencing_run}.mito.anno.vcf.gz.tbi") into mito_vcf_channel

	"""
	bcftools view -r MT $normalised_vcf > ${params.sequencing_run}.mito.vcf
	bgzip ${params.sequencing_run}.mito.vcf
	tabix ${params.sequencing_run}.mito.vcf.gz

	vep \
	--verbose \
	--format vcf \
	--everything \
	--fork $params.vep_cpus \
	--species homo_sapiens \
	--assembly GRCh37 \
	--input_file ${params.sequencing_run}.mito.vcf.gz \
	--output_file ${params.sequencing_run}.mito.anno.vcf \
	--force_overwrite \
	--cache \
	--dir  $vep_cache_mt \
	--fasta $reference_genome \
	--offline \
	--cache_version $params.vepversion_mt \
	--no_escape \
	--shift_hgvs 1 \
	--vcf \
	--merged \
	--flag_pick \
	--pick_order biotype,canonical,appris,tsl,ccds,rank,length \
	--exclude_predicted \
	--custom ${mitomap},mitomap,vcf,exact,0,AF,AC \
	--custom ${clinvar},clinvar,vcf,exact,0,CLNSIG,CLNSIGCONF

	bgzip ${params.sequencing_run}.mito.anno.vcf
	tabix ${params.sequencing_run}.mito.anno.vcf.gz

	"""
}

// Create mitochondrial variant reports 
process create_mito_variant_reports {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/variant_reports_mitochondrial/", mode: 'copy'
	publishDir "${params.publish_dir}/variant_reports_mitochondrial/", mode: 'copy', pattern: '.command.*', saveAs: { filename -> "${sample_names[0]}_$filename" }

	input:
	set val(id), file(vcf), file(vcf_index) from mito_vcf_channel
	each sample_names from samples_ch_reporting
	file(ped) from ped_channel_reports_mito
	file(variables) from variables_channel_mito_variant_report.collect()

	output:
	file("${params.sequencing_run}.${sample_names[0]}_variant_report_mito.csv") into variant_report_mito_channel

	"""
	. ${sample_names[0]}.variables
	
	mitochondrial_variant_reporter.py \
	--vcf $vcf \
	--proband_id ${sample_names[0]} \
	--ped $ped \
	--output ${params.sequencing_run}.${sample_names[0]}_variant_report_mito.csv \
	--worklist \$worklistId \
	--whitelist $whitelist_mito \
	--min_mt_af $params.min_mt_af \
	--max_mitomap_af $params.max_mitomap_af 

	"""
}

// Create UPD plots and report
process create_upd_plots {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/upd/", mode: 'copy'

	input:
	set val(id), file(vcf), file(vcf_index) from for_upd_channel
	each sample_names from samples_ch_upd
	file(ped) from ped_channel_upd

	output:
	file("${params.sequencing_run}.${sample_names[0]}_*") optional true into upd_channel

	when:
	params.calculate_upd == true

	"""
	UPDog.py --vcf $vcf \
	--proband_id ${sample_names[0]} \
	--ped $ped \
	--output ${params.sequencing_run}.${sample_names[0]} \
	--min_dp $params.upd_min_dp \
	--min_gq $params.upd_min_gq \
	--min_qual $params.upd_min_qual \
	--p_value $params.upd_p_value \
	--block_size $params.upd_block_size \
	--min_variants_per_block $params.upd_min_variants_per_block \
	--min_blocks $params.upd_min_blocks \
	--min_proportion $params.upd_min_proportion

	"""
}


variant_report_channel.filter(  {it=~ /$params.test_sample.*/ } ).set{ variant_report_channel_testing }
variant_report_mito_channel.filter( {it=~ /$params.test_sample.*/} ).set{ variant_report_mito_channel_testing }

//  Run Tests
process test_pipeline {

	cpus params.small_task_cpus

	publishDir "${params.publish_dir}/test_results/", mode: 'copy'

	input:
	file(variant_report) from variant_report_channel_testing
	file(mt_variant_report) from variant_report_mito_channel_testing

	output:
	file 'test_report.txt'

	when:
	params.testing == true

	"""
	tests.py \
	--variant_report $variant_report \
	--mt_variant_report $mt_variant_report > test_report.txt
	"""
}

// Create marker once complete
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