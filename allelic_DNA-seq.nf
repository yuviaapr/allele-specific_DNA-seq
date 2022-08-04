#!/usr/bin/env nextflow

/*
 Analysis of DNA-seq data to identify major genomic rearrangements in cell lines.
 Authors: 
	- Yuvia A. PEREZ RICO <yuvia.perez-rico@embl.de>
*/

log.info "          DNA-seq karyotyping - version 0.0.2          "
log.info "#######################################################"
log.info "Sample description file	= ${params.sampleInfo}"
log.info "Reads per split file		= ${params.chunkSize}"
log.info "Mouse strain 1		= ${params.G1}"
log.info "Mouse strain 2		= ${params.G2}"
log.info "Effective genome size		= ${params.effectiveSize}"
log.info "Bin size used for the tracks	= ${params.binSize}"
log.info "Output directory		= ${params.outDir}"
log.info "Path to genome index		= ${params.genomeDirPath}"
log.info "Genome in 2bit format		= ${params.genome2bit}"
log.info "SNP reference file		= ${params.snpFile}"
log.info "Repeat annotations		= ${params.repeats}"
log.info "Blacklisted regions		= ${params.blacklist}"
log.info "Number of threads		= ${params.numCPUs}"
log.info "Number of threads (DeepTools)	= ${params.numCPUsDT}"
log.info "\n"

/*
 Validate input parameters
*/

if( !(params.chunkSize instanceof Number) ){
	exit 1, "Invalid chunk size = ${params.chunkSize}"
}

if( !(params.effectiveSize instanceof Number) ){
	exit 1, "Invalid effective genome size = ${params.effectiveSize}"
}

if( !(params.binSize instanceof Number) ){
	exit 1, "Invalid bin size = ${params.binSize}"
}

if( !(params.numCPUs instanceof Number) ){
	exit 1, "Invalid number of CPUs = ${params.numCPUs}"
}

if( !(params.numCPUsDT instanceof Number) ){
	exit 1, "Invalid number of CPUs for GC correction = ${params.numCPUsDT}"
}

if( !(params.G1 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd', 'PGK']) ){
	exit 1, "Invalid strain 1 name = ${params.G1}"
}

if( !(params.G2 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd', 'PGK']) ){
	exit 1, "Invalid strain 2 name = ${params.G2}"
}

/*
 Validate input files and directory to save results
*/

sdFile = file(params.sampleInfo)
if( !sdFile.exists() ){
	exit 1, "The specified sample description file does not exist = $sdFile"
}
log.info "Checking sample description file = $sdFile"

varFile = file(params.snpFile)
if( !varFile.exists() ){
	exit 1, "The specified SNP annotation file does not exist = $varFile"
}
log.info "Checking SNP annotations file = $varFile"

genome =  file(params.genome2bit)
if( !genome.exists() ){
	exit 1, "The specified genome in 2bit format does not exist = $genome"
}
log.info "Checking genome in 2bit format = $genome"

repeats = file(params.repeats)
if( !repeats.exists() ){
	exit 1, "The specified repeat annotation file does not exist = $repeats"
}
log.info "Checking repeat annotations = $repeats"

blackL = file(params.blacklist)
if( !blackL.exists() ){
	exit 1, "The specified blacklisted annotation file does not exist = $blackL"
}
log.info "Checking blacklist annotations = $blackL"

resDir = file(params.outDir)
if( !resDir.exists() && !resDir.mkdirs() ){
	exit 1, "The specified directory to save results cannot be created = $resDir\n Check file system access permission"
}
log.info "Checking results directory = $resDir"

/*
 Program versions
*/

process get_program_versions{
	publishDir "$resDir/software", mode: 'move'

	output:
	file('programs_version.txt') into programs_version

	"""
	echo nextflow ${nextflow.version} > tmp_version.txt
	fastqc -v >> tmp_version.txt
	echo trim_galore \$(trim_galore -v | awk '/version/{print\$2}') >> tmp_version.txt
	echo bowtie2 \$(bowtie2 --version | head -1 | cut -d " " -f 3) >> tmp_version.txt
	samtools --version | grep samtools >> tmp_version.txt
	echo SNPsplit \$(SNPsplit --version | awk '/Version/{print\$2}') >> tmp_version.txt
	picard SortSam --version 2>&1 | awk '{sub(/-SNAPSHOT/,"");print"picard "\$1}' >> tmp_version.txt
	deeptools --version >> tmp_version.txt
	sort tmp_version.txt > programs_version.txt
	"""

}

/*
 Create channels with the fastq files for processing and quality control
*/

// Reads1

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample, file(row.reads1)) }
	.set { samples_r1 }

// Duplicate channel

samples_r1.into {samples_quality; samples_analysis}

/*
 Step 0. Quality check
*/

process fastq_quality{
	publishDir "$resDir/qc/fastqc", mode: 'copy'

	input:
	set val(name), file(reads) from samples_quality

	output:
	file("${name}_fastqc") into fastqc_results

	"""
	mkdir ${name}_fastqc
	fastqc -o ${name}_fastqc -q -t ${params.numCPUs} ${reads}
	"""

}

/*
 Step 1. Split fastq files
*/

process split_Reads{
	input:
	set val(name), file(reads) from samples_analysis

	output:
	set val(name), file('*.gz') into samples_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads1_"
	"""

}

samples_split
	.transpose()
	.set { samples }

/*
 Step 2. Trimming
*/

process trim_reads{
	publishDir "$resDir/qc/trimming", mode: 'copy', pattern: "${name}/*report.txt"

	input:
	set val(name), file(reads) from samples

	output:
	set val(name), file("${name}/*.fq.gz") into trimmed_samples
	file("${name}/*report.txt") into trim_stats

	"""
	trim_galore -q 30 ${reads} -o ${name}
	"""

}

/*
 Step 3. Mapping
*/

process read_mapping{
	publishDir "$resDir/qc/mapping", mode: 'copy', pattern: '*.log'

	input:
	set val(name), file(trimmed_reads) from trimmed_samples

	output:
	set val(name), file('*.bam') into split_mapping
	file('*.log') into mapping_stats

	"""
	file=\$(echo ${trimmed_reads} | sed s/.gz_trimmed.fq.gz//)
	bowtie2 --very-sensitive --end-to-end -p ${params.numCPUs} -x ${params.genomeDirPath} -U ${trimmed_reads} 2> \$file.log | samtools view -b > \$file.bam
	"""

}

/*
 Step 4. Remove mitochondrial and low quality reads
*/

process chrM_MQ_filter{
	publishDir "$resDir/qc/MQ_filter", mode: 'copy', pattern: '*.stats'

	input:
	set val(name), file(mapped_reads) from split_mapping

	output:
	set val(name), file('*.bam') into split_mapping_filtered
	set val(name), file('*.stats') into filtering_stats

	"""
	samtools flagstat -@ ${params.numCPUs} ${mapped_reads} > ${mapped_reads.baseName}_pre-rmQC.stats
	samtools view -q 30 -h ${mapped_reads} -@ ${params.numCPUs} | grep -v chrM | samtools view -bS - > ${mapped_reads.baseName}_filtered.bam
	samtools flagstat -@ ${params.numCPUs} ${mapped_reads.baseName}_filtered.bam > ${mapped_reads.baseName}_post-rmQC.stats
	"""

}

/*
 Step 5. Genotype assignation
*/

process SNPsplit{
	publishDir "$resDir/qc/SNPsplit", mode: 'copy', pattern: "${name}/*.txt"

	input:
	set val(name), file(filtered_reads) from split_mapping_filtered

	output:
	set val(name), file("${name}/*genome1.bam") into split_genome1
	set val(name), file("${name}/*genome2.bam") into split_genome2
	set val(name), file("${name}/*unassigned.bam") into split_unassigned
	file("${name}/*report.txt") into SNPs_report_stats
	file("${name}/*sort.txt") into SNPs_sort_stats

	"""
	SNPsplit --no_sort --snp_file ${varFile} -o ${name} ${filtered_reads}
	"""

}

/*
 Step 6. Merge, sort and filter BAM files per sample
*/

// Organise files

split_genome1.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome1 }

split_genome2.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome2 }

split_unassigned.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { unassigned }

// Process G1 files

process filter_G1_bams{
	publishDir "$resDir/mapping", mode: 'copy', pattern: '*.bam'
	publishDir "$resDir/qc/duplicates", mode: 'copy', pattern: '*.stats'
	label 'filter_bams'

	input:
	set val(name), file(genome1_bam_files) from genome1

	output:
	set val(name), file('*rmdup.bam') into G1_bam_merge, G1_bam
	file('*.stats') into G1_duplicates_stats

	"""
	# Merge
	samtools merge -f ${name}_${params.G1}.bam ${genome1_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G1}.bam O=${name}_${params.G1}_sorted.bam SORT_ORDER=coordinate
	# Remove duplicates
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_${params.G1}_sorted.bam O=${name}_${params.G1}_sorted_rmdup.bam \
		 M=${name}_${params.G1}_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Process G2 files

process filter_G2_bams{
	publishDir "$resDir/mapping", mode: 'copy', pattern: '*.bam'
	publishDir "$resDir/qc/duplicates", mode: 'copy', pattern: '*.stats'
	label 'filter_bams'

	input:
	set val(name), file(genome2_bam_files) from genome2

	output:
	set val(name), file('*rmdup.bam') into G2_bam_merge, G2_bam
	file('*.stats') into G2_duplicates_stats

	"""
	# Merge
	samtools merge -f ${name}_${params.G2}.bam ${genome2_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G2}.bam O=${name}_${params.G2}_sorted.bam SORT_ORDER=coordinate
	# Remove duplicates
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_${params.G2}_sorted.bam O=${name}_${params.G2}_sorted_rmdup.bam \
		 M=${name}_${params.G2}_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Process unassigned files

process filter_unassigned_bams{
	publishDir "$resDir/qc/duplicates", mode: 'copy', pattern: '*.stats'

	input:
	set val(name), file(unassigned_bam_files) from unassigned

	output:
	set val(name), file('*rmdup.bam') into unassigned_bam_merge
	file('*.stats') into unassigned_duplicates_stats

	"""
	# Merge
	samtools merge -f ${name}_unassigned.bam ${unassigned_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_unassigned.bam O=${name}_unassigned_sorted.bam SORT_ORDER=coordinate
	# Remove duplicates
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MarkDuplicates I=${name}_unassigned_sorted.bam O=${name}_unassigned_sorted_rmdup.bam \
		 M=${name}_unassigned_sorted_rmdup.stats REMOVE_DUPLICATES=true ASSUME_SORTED=true
	"""

}

// Merge all reads

process merge_mapped{
	publishDir "$resDir/mapping", mode: 'copy'

	input:
	set val(name), file(genome1_bam), file(genome2_bam), file(unassigned_bam) from G1_bam_merge.join(G2_bam_merge).join(unassigned_bam_merge)

	output:
	set val(name), file('*sorted_rmdup.bam') into all_bam

	"""
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MergeSamFiles MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true \
		I=${genome1_bam} I=${genome2_bam} I=${unassigned_bam} O=${name}_Gall_rmdup.bam
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_Gall_rmdup.bam O=${name}_Gall_sorted_rmdup.bam \
		SORT_ORDER=coordinate
	"""

}

// Set channel with all bam files

G1_bam.mix(G2_bam).mix(all_bam).set { bam_GCnorm }

/*
 Step 8. GC normalization of bam files
*/

process GC_norm{
	publishDir "$resDir/normalizedBams", mode: 'copy', pattern: "*noBL.bam"

	input:
	set val(name), file(bam) from bam_GCnorm

	output:
	set val(name), file("*noBL.bam") into bam_tracks
	set val(name), file("*${params.G1}*noBL.bam") optional true into G1_bam_GCnorm
	set val(name), file("*${params.G2}*noBL.bam") optional true into G2_bam_GCnorm
	set val(name), file('*Gall*noBL.bam') optional true into bam_sizeFactor

	"""
	samtools index ${bam}
	computeGCBias -b ${bam} --effectiveGenomeSize ${params.effectiveSize} -g ${genome} --GCbiasFrequenciesFile ${bam.baseName}_GC-frequencies_noRepeats_noBL.txt \
		--numberOfProcessors ${params.numCPUsDT} --regionSize ${params.binSize} -l 200 --blackListFileName ${repeats} ${blackL}
	correctGCBias -b ${bam}	--effectiveGenomeSize ${params.effectiveSize} -g ${genome} --GCbiasFrequenciesFile ${bam.baseName}_GC-frequencies_noRepeats_noBL.txt \
		--correctedFile ${bam.baseName}_GC_corrected_noRepeats_noBL.bam --binSize ${params.binSize} --numberOfProcessors ${params.numCPUsDT}
	"""

}

/*
 Step 9a. Generate normalized signal tracks
*/

// Calculate size factors

process size_factors{
	input:
	set val(name), file(all_bam_SF) from bam_sizeFactor

	output:
	set val(name), stdout into size_factors

	"""
	samtools view -c ${all_bam_SF} | awk '{print 1000000/\$1}'
	"""

}

// Generate tracks

process signal_tracks{
	label 'tracks'
	publishDir "$resDir/signal", mode: 'copy'

	input:
	set val(name), file(bam_cov), val(sizeFactor) from bam_tracks.combine(size_factors, by: 0)

	output:
	file("*.bw") into bigWig_signal

	"""
	samtools index ${bam_cov}
	bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}.bw --binSize ${params.binSize} --normalizeUsing RPKM \
		--effectiveGenomeSize ${params.effectiveSize} --numberOfProcessors ${params.numCPUs} --outFileFormat bigwig --scaleFactor ${sizeFactor}
	"""

}

/*
 Step 9b. Generate normalized signal tracks
*/

process compare_tracks{
	label 'tracks'
	publishDir "$resDir/signal", mode: 'copy'

	input:
	set val(name), file(bam_G1), file(bam_G2) from G1_bam_GCnorm.combine(G2_bam_GCnorm, by: 0)

	output:
	file("*corrected.bw") into bigWig_log2

	"""
	samtools index ${bam_G1}
	samtools index ${bam_G2}
	bamCompare -b1 ${bam_G1} -b2 ${bam_G2} --outFileName ${name}_log2_${params.G1}_over_${params.G2}_GC_corrected.bw --scaleFactorsMethod readCount \
		--operation log2 --binSize ${params.binSize} --effectiveGenomeSize ${params.effectiveSize} --normalizeUsing None --numberOfProcessors ${params.numCPUs}
	"""

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Main results are saved in $resDir\n" : "There was an error during the execution, check log files." )
}


