/*
 * Read-based phasing per chromosome -> locally trustworthy blocks with PS tags
 */

process WHATSHAP_PHASE {
    tag "${chrom}"
    label 'phase'

    input:
    tuple val(chrom), path(het_vcf), path(het_tbi)
    tuple path(bam), path(bai)
    tuple path(ref), path(fai)

    output:
    tuple val('whatshap'), val(chrom), path("reads_whatshap_${chrom}.vcf.gz"), path("reads_whatshap_${chrom}.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -r ${chrom} ${het_vcf} -Oz -o input_${chrom}.vcf.gz
    tabix -p vcf input_${chrom}.vcf.gz

    whatshap phase --reference ${ref} --ignore-read-groups --chromosome ${chrom} \\
        -o reads_whatshap_${chrom}.vcf.gz input_${chrom}.vcf.gz ${bam}
    tabix -f -p vcf reads_whatshap_${chrom}.vcf.gz
    rm -f input_${chrom}.vcf.gz*
    """
}

process LONGPHASE_PHASE {
    tag "${chrom}"
    label 'phase'

    input:
    tuple val(chrom), path(het_vcf), path(het_tbi)
    tuple path(bam), path(bai)
    tuple path(ref), path(fai)

    output:
    tuple val('longphase'), val(chrom), path("reads_longphase_${chrom}.vcf.gz"), path("reads_longphase_${chrom}.vcf.gz.tbi"), emit: vcf

    script:
    def platform = params.platform == 'hifi' ? '--pb' : '--ont'
    """
    bcftools view -r ${chrom} ${het_vcf} -Ov -o input_${chrom}.vcf

    longphase phase -s input_${chrom}.vcf -b ${bam} -r ${ref} -t ${task.cpus} ${platform} \\
        -o reads_longphase_${chrom}

    bcftools sort reads_longphase_${chrom}.vcf -Oz -o reads_longphase_${chrom}.vcf.gz
    tabix -f -p vcf reads_longphase_${chrom}.vcf.gz
    rm -f input_${chrom}.vcf reads_longphase_${chrom}.vcf
    """
}
