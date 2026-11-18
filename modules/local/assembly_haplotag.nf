/*
 * Optional (--hap1 / --hap2): haplotag reads by mapping them to the haplotype-resolved assembly.
 * Independent of variant calling and statistical phasing; works for long-read BAMs and for mixed
 * short / long fragments (single-cell DNA).
 */

process ASM_INDEX {
    tag "${params.sample}"
    label 'asm_map'

    input:
    tuple path(hap1), path(hap2)

    output:
    path "diploid.mmi", emit: mmi

    script:
    def preset = params.platform == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    # one diploid reference: contigs of hap1 -> h1_*, hap2 -> h2_*
    { zcat -f ${hap1} | sed 's/^>/>h1_/'; zcat -f ${hap2} | sed 's/^>/>h2_/'; } > diploid.fa
    minimap2 -x ${preset} -t ${task.cpus} -d diploid.mmi diploid.fa
    rm diploid.fa
    """
}

process ASM_MAP {
    tag "${params.sample}"
    label 'asm_map'

    input:
    tuple path(bam), path(bai)
    path mmi

    output:
    path "${params.sample}.asm.byname.bam", emit: bam

    script:
    def preset = params.platform == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    # primary records only, cell barcode / UMI tags carried in the read comment (-T, minimap2 -y)
    samtools fastq -@ 4 -F 0x900 -T ${params.asm_keep_tags} ${bam} \\
      | minimap2 -ax ${preset} -y --secondary=yes -N 10 -p 0.5 -t ${task.cpus} ${params.asm_minimap2_args} ${mmi} - \\
      | samtools view -@ 4 -b -o ${params.sample}.asm.byname.bam -
    """
}

process ASM_ASSIGN {
    tag "${params.sample}"
    label 'haplotag'
    publishDir "${params.outdir}/${params.sample}/05_haplotag/assembly", mode: 'copy',
        saveAs: { it.endsWith('.tsv') || it.endsWith('.tsv.gz') || it.contains('.asm.bam') ? it : null }

    input:
    path bam

    output:
    path "${params.sample}.asm_assignments.tsv.gz", emit: assignments
    path "asm.by_length.tsv", emit: by_length
    path "asm.as_diff.tsv",   emit: as_diff
    tuple path("${params.sample}.asm.bam"), path("${params.sample}.asm.bam.bai"), emit: asm_bam

    script:
    """
    # minimap2 output is grouped by read: decide per read before sorting
    assembly_haplotag.py assign ${bam} ${params.sample}.asm_assignments.tsv.gz asm --min-diff ${params.asm_min_diff}

    # alignments to the assembly, sorted, for IGV / downstream use
    samtools sort -@ ${task.cpus} -o ${params.sample}.asm.bam ${bam}
    samtools index -@ ${task.cpus} ${params.sample}.asm.bam
    """
}

process ASM_TAG {
    tag "${params.sample}"
    label 'haplotag'
    publishDir "${params.outdir}/${params.sample}/05_haplotag/assembly", mode: 'copy'

    input:
    path assignments
    tuple path(bam), path(bai)      // original BAM, or the VCF-haplotagged one (-> agreement table)

    output:
    tuple path("${params.sample}.assembly_haplotagged.bam"), path("${params.sample}.assembly_haplotagged.bam.bai"), emit: bam
    path "asm.bins.tsv",            emit: bins
    path "asm.vs_vcf_haplotag.tsv", emit: agreement
    path "asm.per_cell.tsv.gz",     optional: true, emit: cells

    script:
    """
    assembly_haplotag.py tag ${assignments} ${bam} ${params.sample}.assembly_haplotagged.bam asm \\
        --bin ${((params.haplotag_bin_mb as double) * 1e6) as long}
    samtools index -@ ${task.cpus} ${params.sample}.assembly_haplotagged.bam
    """
}
