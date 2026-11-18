/*
 * Optional (--hap1 / --hap2), UNDER DEVELOPMENT: haplotag reads with the haplotype-resolved assembly.
 *
 * The reads are aligned to hap1 and to hap2 separately; HipHap (Heinz, Marin, Meyerson & Li, bioRxiv 2026)
 * assigns every read to its better-supported haplotype with a confidence score (HapQ). Reads with
 * HapQ >= --asm_min_hapq get HP (+ hq) on the GRCh38 BAM. Independent of variant calling and statistical
 * phasing; works for long-read BAMs and for mixed short / long fragments (single-cell DNA).
 */

process ASM_INDEX {
    tag "${params.sample}"
    label 'asm_map'

    input:
    tuple path(hap1), path(hap2)

    output:
    tuple path("hap1.mmi"), path("hap2.mmi"), emit: mmi

    script:
    def preset = params.platform == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    # HipHap's merged output needs contig names that differ between the haplotypes: prefix h1_ / h2_
    for h in 1 2; do
        f=\$([ \$h = 1 ] && echo ${hap1} || echo ${hap2})
        zcat -f \$f | sed "s/^>/>h\${h}_/" > hap\${h}.fa
        minimap2 -x ${preset} -t ${task.cpus} -d hap\${h}.mmi hap\${h}.fa
        rm hap\${h}.fa
    done
    """
}

process ASM_MAP {
    tag "${params.sample}"
    label 'asm_map'

    input:
    tuple path(bam), path(bai)
    tuple path(mmi1), path(mmi2)

    output:
    tuple path("hap1.bam"), path("hap2.bam"), emit: bams

    script:
    def preset = params.platform == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    # primary records only; cell barcode / UMI tags carried in the read comment (-T, minimap2 -y)
    samtools fastq -@ 4 -F 0x900 -T ${params.asm_keep_tags} ${bam} | bgzip -@ 4 > reads.fq.gz

    # one alignment per haplotype; minimap2 keeps the input read order, which HipHap needs
    for h in 1 2; do
        minimap2 -ax ${preset} -y -t ${task.cpus} ${params.asm_minimap2_args} hap\${h}.mmi reads.fq.gz \\
          | samtools view -@ 4 -b -o hap\${h}.bam -
    done
    rm reads.fq.gz
    """
}

process HIPHAP {
    tag "${params.sample}"
    label 'haplotag'
    publishDir "${params.outdir}/${params.sample}/05_haplotag/assembly", mode: 'copy',
        saveAs: { it.endsWith('.tsv') || it.endsWith('.tsv.gz') || it.contains('.asm.bam') || it.endsWith('.fastq') ? it : null }

    input:
    tuple path(bam1), path(bam2)

    output:
    path "${params.sample}.asm_assignments.tsv.gz", emit: assignments
    path "asm.by_length.tsv", emit: by_length
    path "asm.hapq.tsv",      emit: hapq
    tuple path("${params.sample}.asm.bam"), path("${params.sample}.asm.bam.bai"), emit: asm_bam
    path "*span_chrom.fastq", optional: true, emit: span_chrom

    script:
    """
    # merged alignments, each read on its better-supported haplotype, with HP (1/2) and hq (HapQ)
    hiphap -1 hap1 -2 hap2 -t ${task.cpus} -O bam -o merged.bam ${params.hiphap_args} ${bam1} ${bam2}

    # per-read table (read, HP, HapQ, contig, length) + summaries for the report
    assembly_haplotag.py table merged.bam ${params.sample}.asm_assignments.tsv.gz asm --min-hapq ${params.asm_min_hapq}

    # reads on the assembly, sorted, for IGV / downstream use
    samtools sort -@ ${task.cpus} -o ${params.sample}.asm.bam merged.bam
    samtools index -@ ${task.cpus} ${params.sample}.asm.bam
    rm merged.bam
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
        --min-hapq ${params.asm_min_hapq} --bin ${((params.haplotag_bin_mb as double) * 1e6) as long}
    samtools index -@ ${task.cpus} ${params.sample}.assembly_haplotagged.bam
    """
}
