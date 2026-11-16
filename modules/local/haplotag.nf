/*
 * Haplotag the BAM with the final phased VCF + read counts per haplotype for the report
 */

process HAPLOTAG {
    tag "${combo}"
    label 'haplotag'
    publishDir "${params.outdir}/${params.sample}/05_haplotag", mode: 'copy'

    input:
    tuple val(combo), path(vcf), path(tbi)
    tuple path(bam), path(bai)
    tuple path(ref), path(fai)

    output:
    tuple val(combo), path("${params.sample}.${combo}.haplotagged.bam"), path("${params.sample}.${combo}.haplotagged.bam.bai"), emit: bam
    path "${params.sample}.${combo}.haplotag_list.tsv.gz", optional: true, emit: list
    path "${params.sample}.${combo}.haplotag.log"

    script:
    def out = "${params.sample}.${combo}"
    if (params.haplotag_tool == 'longphase') {
        """
        longphase haplotag -r ${ref} -s ${vcf} -b ${bam} -t ${task.cpus} -o ${out}.haplotagged \\
            > ${out}.haplotag.log 2>&1
        samtools index -@ ${task.cpus} ${out}.haplotagged.bam
        """
    } else {
        """
        whatshap haplotag --ignore-read-groups --skip-missing-contigs \\
            --reference ${ref} --output-threads ${task.cpus} \\
            --output-haplotag-list ${out}.haplotag_list.tsv.gz \\
            -o ${out}.haplotagged.bam ${vcf} ${bam} > ${out}.haplotag.log 2>&1
        samtools index -@ ${task.cpus} ${out}.haplotagged.bam
        """
    }
}

process HAPLOTAG_COUNTS {
    tag "${combo}"
    label 'phase'
    publishDir "${params.outdir}/${params.sample}/05_haplotag", mode: 'copy'

    input:
    tuple val(combo), path(bam), path(bai)

    output:
    path "haplotag_counts.tsv", emit: tsv

    script:
    def bin = (params.haplotag_bin_mb as double) * 1e6
    """
    # primary, mapped reads: chrom, bin, HP (0 = untagged), reads, bases
    samtools view -@ ${task.cpus} -F 0x904 ${bam} \\
      | awk -F'\\t' -v B=${bin as long} '
          { hp = 0
            for (i = 12; i <= NF; i++) if (substr(\$i, 1, 5) == "HP:i:") { hp = substr(\$i, 6); break }
            k = \$3 "\\t" int(\$4 / B) "\\t" hp; n[k]++; b[k] += length(\$10) }
          END { print "chrom\\tbin\\thp\\treads\\tbases"; for (k in n) print k "\\t" n[k] "\\t" b[k] }' \\
      > haplotag_counts.tsv
    """
}
