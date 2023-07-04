#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2


params.csv = null
params.ref = null
params.bam = null
params.meta = null
params.sample_id = null
params.outdir = "./output"

log.info """\
    S C O M A T I C - N F   P I P E L I N E
    ===================================
    input_csv   : ${params.csv}
    reference   : ${params.ref}
    outdir      : ${params.outdir}
    bam         : ${params.bam}
    metadata    : ${params.meta}
    sample_id    : ${params.sample_id}
    """
    .stripIndent()


process SPLITBAM {
    tag "Sample: $sample_id"

    publishDir params.outdir, mode:'copy'

    input:
        tuple val(sample_id), path(bam), path(bai), path(meta)

    output:
        tuple val(sample_id), path("${sample_id}/Step1_BamCellTypes/*.bam"), path("${sample_id}/Step1_BamCellTypes/*.bai"), path("${sample_id}/Step1_BamCellTypes/")

        script:
        """
        outdir="${sample_id}/Step1_BamCellTypes"

        mkdir -p \${outdir}
        python /scomatic/scripts/SplitBam/SplitBamCellTypes.py \
            --bam $bam \
            --meta $meta \
            --id $sample_id \
            --n_trim 5 \
            --max_nM 5 \
            --max_NH 1 \
            --outdir \${outdir}
        """
}

process BASECOUNTS_SPLIT {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy', overwrite: true
    tag "Sample: $sample_id; BAM: $bam"

    input:
        path ref
        tuple val(sample_id), path(bam), path(bai), path(outdir1)

    output:
        tuple val(sample_id), path ("${sample_id}/Step2_BaseCellCounts/*.tsv"), optional: true

    script:
        """
        outdir="${sample_id}/Step2_BaseCellCounts"
        mkdir -p \${outdir}

        cell_type=\$(basename $bam | awk -F'.' '{print \$(NF-1)}')

        # Temp folder
        temp="${sample_id}/Step2_BaseCellCounts/temp_\${cell_type}"
        mkdir -p \$temp

        # Command line to submit to cluster
        python /scomatic/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
            --ref $ref \
            --chrom all \
            --out_folder "\${outdir}" \
            --min_bq 30 \
            --tmp_dir \$temp \
            --nprocs $task.cpus

        rm -rf \$temp
        """
}


process MERGECOUNTS {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy'

    tag "Sample: $sample_id"

    input:
        tuple val(sample_id), path(tsvs)

    output:
        tuple val(sample_id), path("${sample_id}/Step3_BaseCellCountsMerged/")

    script:
        """
        output_dir3=${sample_id}/Step3_BaseCellCountsMerged
        mkdir -p \$output_dir3

        # Create a tempdir with each sample tsvs
        mkdir -p temp/
        cp ${tsvs} temp/


        python /scomatic/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder temp/ \
        --outfile \${output_dir3}/${sample_id}.BaseCellCounts.AllCellTypes.tsv

        rm -rf temp/
        """
}


process VARIANTCALLING {
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sample_id"

    input:
        path ref
        tuple val(sample_id), path(outdir3)

    output:
        tuple val(sample_id), path("${sample_id}/Step4_VariantCalling/")

    script:
        """
        # Step 4.1
        output_dir4=${sample_id}/Step4_VariantCalling
        mkdir -p \$output_dir4

        python /scomatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
                --infile ${outdir3}/${sample_id}.BaseCellCounts.AllCellTypes.tsv \
                --outfile \${output_dir4}/${sample_id} \
                --ref $ref

        # Step 4.2
        editing=/scomatic/RNAediting/AllEditingSites.hg38.txt
        PON=/scomatic/PoNs/PoN.scRNAseq.hg38.tsv

        python /scomatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
                --infile \${output_dir4}/${sample_id}.calling.step1.tsv \
                --outfile \${output_dir4}/${sample_id} \
                --editing \$editing \
                --pon \$PON

        # check chr prefix before intersecting
        if [ \$(grep -c "^chr" \${output_dir4}/${sample_id}.calling.step2.tsv) -gt 1 ]
        then
            echo "Chr prefix detected"
            repeat_bed="/scomatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed"
        else
            echo "No chr prefix detected"
            repeat_bed="/scomatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker_nochr.bed"
        fi
        echo "Using \${repeat_bed} for hq genomic regions"

        bedtools intersect -header -a \${output_dir4}/${sample_id}.calling.step2.tsv \
            -b \${repeat_bed} | awk '\$1 ~ /^#/ || \$6 == "PASS"' > \${output_dir4}/${sample_id}.calling.step2.pass.tsv
        """
}

process CALLABLESITES {
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sample_id"

    input:
        tuple val(sample_id), path(outdir4)

    output:
        tuple val(sample_id), path("${sample_id}/Step5_CellTypeCallableSites/")

    script:
        """
        # Computing the number of callable sites per cell type
        output_dir5=${sample_id}/Step5_CellTypeCallableSites
        mkdir -p \$output_dir5

        python /scomatic/scripts/GetCallableSites/GetAllCallableSites.py --infile ${outdir4}/${sample_id}.calling.step1.tsv  \
            --outfile \$output_dir5/${sample_id} \
            --max_cov 150 --min_cell_types 2
   """
}


process CALLABLE_PERCT {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sample_id; BAM: $bam"

    input:
        path ref 
        tuple val(sample_id), path(bam), path(bai), path(outdir1), path(outdir4)

    output:
        path("${sample_id}/Step6_UniqueCellCallableSites/*.tsv"), optional: true

    script:
        """
        STEP4_1=$outdir4/${sample_id}.calling.step1.tsv

        output_dir6=${sample_id}/Step6_UniqueCellCallableSites
        mkdir -p \$output_dir6

        cell_type=\$(basename $bam | awk -F'.' '{print \$(NF-1)}')
        echo \$cell_type
        
        temp=\$output_dir6/temp_\${cell_type}
        mkdir -p \$temp

        python /scomatic/scripts/SitesPerCell/SitesPerCell.py --bam $bam    \
            --infile $outdir4/${sample_id}.calling.step1.tsv   \
            --ref $ref \
            --out_folder \$output_dir6 --tmp_dir \$temp --nprocs $task.cpus
        echo
        """
}

process GENOTYPE_CELLS {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sample_id; BAM: $bam"

    input:
        path ref
        tuple val(sample_id), path(bam), path(bai), path(outdir1), path(outdir4), path(meta)

    output:
        path("${sample_id}/Step7_SingleCellAlleles/*.tsv"), optional: true

    script:
        """
        STEP4_2_pass=${outdir4}/${sample_id}.calling.step2.pass.tsv

        output_dir7=${sample_id}/Step7_SingleCellAlleles
        mkdir -p \$output_dir7

        cell_type=\$(basename $bam | awk -F'.' '{print \$(NF-1)}')
        
        temp=\$output_dir7/temp_\${cell_type}
        mkdir -p \$temp

        python /scomatic/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
            --infile \${STEP4_2_pass}   \
            --nprocs $task.cpus   \
            --meta $meta   \
            --outfile \${output_dir7}/\${cell_type}.single_cell_genotype.tsv  \
            --tmp_dir \$temp  \
            --ref $ref

        rm -rf \$temp
        """
}

workflow {

    if (params.csv) {
        input_ch = Channel.fromPath(params.csv, checkIfExists: true).
            splitCsv(header:true).
            map(row -> tuple "${row.sample_id}", "${row.bam}", "${row.bai}", "${row.metadata}")
    } else {        
        input_ch = Channel.from(params.sample_id).
            merge(Channel.fromPath(params.bam, checkIfExists: true)).
            merge(Channel.fromPath(params.bai, checkIfExists: true)).
            merge(Channel.fromPath(params.meta, checkIfExists: true))
    }


    splitbam_outch = SPLITBAM(input_ch)

    basecounts_outch = BASECOUNTS_SPLIT(params.ref, SPLITBAM.out.transpose())

    mergecounts_outch = MERGECOUNTS(BASECOUNTS_SPLIT.out.groupTuple())

    variantcalling_outch = VARIANTCALLING(params.ref, MERGECOUNTS.out)

    callablesites_outch = CALLABLESITES(VARIANTCALLING.out)

    callable_inch = splitbam_outch.join(variantcalling_outch).transpose()

    CALLABLE_PERCT(params.ref, callable_inch)

    gt_inch = splitbam_outch.join(variantcalling_outch).join(input_ch.map {
        sample_id, bam, bai, meta ->
            tuple(sample_id, meta)
    }).transpose()

    GENOTYPE_CELLS(params.ref, gt_inch)
}


workflow.onComplete {
    log.info(workflow.success ? "Success!" : "Something has gone awry...")
}