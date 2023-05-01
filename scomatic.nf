#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2


log.info """\
    S C O M A T I C - N F   P I P E L I N E
    ===================================
    scomatic    : ${params.scomatic}
    reference   : ${params.ref}
    outdir      : ${params.outdir}
    bam         : ${params.bam}
    metadata    : ${params.meta}
    sampleid    : ${params.sampleid}
    """
    .stripIndent()


process SPLITBAM {
    tag "Sample: $sampleid"

    publishDir params.outdir, mode:'copy'

    input:
        tuple val(sampleid), path(meta)
        path bam
        path bai

    output:
        // val sampleid
        // path "${sampleid}/Step1_BamCellTypes"
        // path ("${sampleid}/Step1_BamCellTypes/*.bam")
        // tuple val(sampleid), path("${sampleid}/Step1_BamCellTypes/*.bam"), path("${sampleid}/Step1_BamCellTypes/*.bai")
        tuple val(sampleid), path("${sampleid}/Step1_BamCellTypes/*.bam"), path("${sampleid}/Step1_BamCellTypes/*.bai"), path("${sampleid}/Step1_BamCellTypes/")
        // tuple val(sampleid), path("${sampleid}/Step1_BamCellTypes/")

        script:
        """
        outdir="${sampleid}/Step1_BamCellTypes"

        mkdir -p \${outdir}
        python /scomatic/scripts/SplitBam/SplitBamCellTypes.py \
            --bam $bam \
            --meta $meta \
            --id $sampleid \
            --n_trim 5 \
            --max_nM 5 \
            --max_NH 1 \
            --outdir \${outdir}
        """
}

// TODO: The basecounts step can be split out per cell type rather than in a loop

process BASECOUNTS_SPLIT {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy', overwrite: true
    tag "Sample: $sampleid; BAM: $bam"

    input:
        path ref
        tuple val(sampleid), path(bam), path(bai), path(outdir1)

    output:
        tuple val(sampleid), path ("${sampleid}/Step2_BaseCellCounts/*.tsv"), optional: true
        // path ("${sampleid}/Step2_BaseCellCounts/*.tsv"), optional: true
        // val "${sampleid}/Step2_BaseCellCounts/", emit: step2_dir
        // val sampleid

    script:
        """
        outdir="${sampleid}/Step2_BaseCellCounts"
        mkdir -p \${outdir}

        cell_type=\$(basename $bam | awk -F'.' '{print \$(NF-1)}')

        # Temp folder
        temp="${sampleid}/Step2_BaseCellCounts/temp_\${cell_type}"
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

    tag "Sample: $sampleid"

    input:
        tuple val(sampleid), path(tsvs)

    output:
        tuple val(sampleid), path("${sampleid}/Step3_BaseCellCountsMerged/")

    script:
        """
        output_dir3=${sampleid}/Step3_BaseCellCountsMerged
        mkdir -p \$output_dir3

        # Create a tempdir with each sample tsvs
        mkdir -p temp/
        cp ${tsvs} temp/


        python /scomatic/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder temp/ \
        --outfile \${output_dir3}/${sampleid}.BaseCellCounts.AllCellTypes.tsv

        rm -rf temp/
        """
}


process VARIANTCALLING {
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sampleid"

    input:
        path ref
        tuple val(sampleid), path(outdir3)

    output:
        tuple val(sampleid), path("${sampleid}/Step4_VariantCalling/")

    script:
        """
        # Step 4.1
        output_dir4=${sampleid}/Step4_VariantCalling
        mkdir -p \$output_dir4

        python /scomatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
                --infile ${outdir3}/${sampleid}.BaseCellCounts.AllCellTypes.tsv \
                --outfile \${output_dir4}/${sampleid} \
                --ref $ref

        # Step 4.2
        editing=/scomatic/RNAediting/AllEditingSites.hg38.txt
        PON=/scomatic/PoNs/PoN.scRNAseq.hg38.tsv

        python /scomatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
                --infile \${output_dir4}/${sampleid}.calling.step1.tsv \
                --outfile \${output_dir4}/${sampleid} \
                --editing \$editing \
                --pon \$PON


        bedtools intersect -header -a \${output_dir4}/${sampleid}.calling.step2.tsv \
            -b /scomatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '\$1 ~ /^#/ || \$6 == "PASS"' > \${output_dir4}/${sampleid}.calling.step2.pass.tsv
        """
}

process CALLABLESITES {
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sampleid"

    input:
        tuple val(sampleid), path(outdir4)

    output:
        tuple val(sampleid), path("${sampleid}/Step5_CellTypeCallableSites/")

    script:
        """
        # Computing the number of callable sites per cell type
        output_dir5=${sampleid}/Step5_CellTypeCallableSites
        mkdir -p \$output_dir5

        python /scomatic/scripts/GetCallableSites/GetAllCallableSites.py --infile ${outdir4}/${sampleid}.calling.step1.tsv  \
            --outfile \$output_dir5/${sampleid} \
            --max_cov 150 --min_cell_types 2
   """
}


process CALLABLE_PERCT {
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sampleid; BAM: $bam"

    input:
        path ref 
        tuple val(sampleid), path(bam), path(bai), path(outdir1), path(outdir4)

    output:
        path("${sampleid}/Step6_UniqueCellCallableSites/*.tsv"), optional: true

    script:
        """
        STEP4_1=$outdir4/${sampleid}.calling.step1.tsv

        output_dir6=${sampleid}/Step6_UniqueCellCallableSites
        mkdir -p \$output_dir6

        cell_type=\$(basename $bam | awk -F'.' '{print \$(NF-1)}')
        echo \$cell_type
        
        temp=\$output_dir6/temp_\${cell_type}
        mkdir -p \$temp

        python /scomatic/scripts/SitesPerCell/SitesPerCell.py --bam $bam    \
            --infile $outdir4/${sampleid}.calling.step1.tsv   \
            --ref $ref \
            --out_folder \$output_dir6 --tmp_dir \$temp --nprocs $task.cpus
        echo
        """
}

process GENOTYPE_CELLS {
    label 'multicore'
    publishDir path: params.outdir, mode:'copy'
    tag "Sample: $sampleid; BAM: $bam"

    input:
        path ref
        tuple val(sampleid), path(bam), path(bai), path(outdir1), path(outdir4), path(meta)

    output:
        path("${sampleid}/Step7_SingleCellAlleles/*.tsv"), optional: true

    script:
        """
        STEP4_2_pass=${outdir4}/${sampleid}.calling.step2.pass.tsv

        output_dir7=${sampleid}/Step7_SingleCellAlleles
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
    bam_ch = Channel.fromPath(params.bam, checkIfExists: true)
    bai_ch = Channel.fromPath(params.bai, checkIfExists: true)
    meta_ch = Channel.fromPath(params.meta, checkIfExists: true)
    sampleid_ch = Channel.from(params.sampleid)
    

    // Bind meta channel to sample ids
    samplemeta_ch = sampleid_ch.combine(meta_ch)
    
    splitbam_outch = SPLITBAM(samplemeta_ch, bam_ch, bai_ch)

    basecounts_outch = BASECOUNTS_SPLIT(params.ref, SPLITBAM.out.transpose())

    mergecounts_outch = MERGECOUNTS(BASECOUNTS_SPLIT.out.groupTuple())

    variantcalling_outch = VARIANTCALLING(params.ref, MERGECOUNTS.out)

    callablesites_outch = CALLABLESITES(VARIANTCALLING.out)

    callable_inch = splitbam_outch.join(variantcalling_outch).transpose()

    CALLABLE_PERCT(params.ref, callable_inch)

    gt_inch = splitbam_outch.join(variantcalling_outch).join(samplemeta_ch).transpose()

    GENOTYPE_CELLS(params.ref, gt_inch)
}


workflow.onComplete {
    log.info(workflow.success ? "Success!" : "Something has gone awry...")
}