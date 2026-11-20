/*
 * Clair3 germline variant calling + germline heterozygous SNP filter
 */

process CLAIR3 {
    tag "${sample}"
    label 'clair3'
    publishDir "${params.outdir}/${params.sample}/01_variant_calling/clair3", mode: 'copy',
        saveAs: { it.endsWith('.vcf.gz') || it.endsWith('.tbi') || it.endsWith('.log') ? file(it).name : null }

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref), path(fai)
    path model_dir

    output:
    tuple val(sample), path("clair3/merge_output.vcf.gz"), path("clair3/merge_output.vcf.gz.tbi"), emit: vcf
    path "clair3/run_clair3.log", optional: true, emit: log

    script:
    def platform = params.platform == 'hifi' ? 'hifi' : 'ont'
    def model    = "${model_dir}/${params.clair3_model}"
    """
    # Clair3 checks for this variable inside the container
    export CONDA_PREFIX=/opt

    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${task.cpus} \\
        --platform=${platform} \\
        --model_path=\$(readlink -f ${model}) \\
        --sample_name=${sample} \\
        --output=\$PWD/clair3 \\
        ${params.clair3_args}
    """
}

process FILTER_HET_SNPS {
    tag "${sample}"
    label 'light'
    publishDir "${params.outdir}/${params.sample}/01_variant_calling", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("germline_het_snps.vcf.gz"), path("germline_het_snps.vcf.gz.tbi"), emit: vcf
    path "qc/het_funnel.tsv",          emit: funnel
    path "qc/het_calls.tsv.gz",        emit: calls
    path "qc/het_kept.tsv.gz",         emit: kept
    path "qc/01_germline_het_snps.stats", emit: stats

    script:
    def vaf = "(FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]))"
    """
    mkdir -p qc
    bcftools view -v snps -m2 -M2 -f PASS ${vcf} -Oz -o pass_snps.vcf.gz
    tabix -p vcf pass_snps.vcf.gz

    bcftools filter -i 'GT="het" && QUAL>=${params.min_qual} && FORMAT/DP>=${params.min_dp} && ${vaf} >= ${params.min_vaf} && ${vaf} <= ${params.max_vaf}' \\
        pass_snps.vcf.gz -Oz -o germline_het_snps.vcf.gz
    tabix -f -p vcf germline_het_snps.vcf.gz
    bcftools stats germline_het_snps.vcf.gz > qc/01_germline_het_snps.stats

    # numbers for the report
    {
      printf "step\\tn\\n"
      printf "all_records\\t%s\\n"            \$(bcftools view -H ${vcf} | wc -l)
      printf "pass\\t%s\\n"                   \$(bcftools view -H -f PASS ${vcf} | wc -l)
      printf "pass_biallelic_snps\\t%s\\n"    \$(bcftools view -H pass_snps.vcf.gz | wc -l)
      printf "het\\t%s\\n"                    \$(bcftools view -H -i 'GT="het"' pass_snps.vcf.gz | wc -l)
      printf "hom_alt\\t%s\\n"                \$(bcftools view -H -i 'GT="AA"' pass_snps.vcf.gz | wc -l)
      printf "het_kept\\t%s\\n"               \$(bcftools view -H germline_het_snps.vcf.gz | wc -l)
    } > qc/het_funnel.tsv

    # every PASS het SNP with the values used by the filter (for the QC plots)
    bcftools query -i 'GT="het"' -f '%CHROM\\t%POS\\t%QUAL\\t[%DP]\\t[%AD{0}]\\t[%AD{1}]\\n' pass_snps.vcf.gz \\
        | awk 'BEGIN{OFS="\\t"; print "chrom","pos","qual","dp","ad_ref","ad_alt"} {print}' | gzip > qc/het_calls.tsv.gz
    bcftools query -f '%CHROM\\t%POS\\n' germline_het_snps.vcf.gz \\
        | awk 'BEGIN{OFS="\\t"; print "chrom","pos"} {print}' | gzip > qc/het_kept.tsv.gz

    rm -f pass_snps.vcf.gz*
    """
}

/*
 * lrsinglecell (IntGenomicsLab) DNA BAMs have duplicates FLAGGED (0x400), not removed.
 * Clair3 counts flagged reads, so drop them before calling / phasing.
 * The original BAM (flags kept) is the one that gets haplotagged.
 */
process DROP_DUPLICATES {
    tag "${params.sample}"
    label 'phase'

    input:
    tuple path(bam), path(bai)

    output:
    tuple path("${params.sample}.nodup.bam"), path("${params.sample}.nodup.bam.bai"), emit: bam
    path "${params.sample}.duplicates.tsv", emit: stats

    script:
    """
    samtools view -@ ${task.cpus} -b -F 0x400 -o ${params.sample}.nodup.bam ${bam}
    samtools index -@ ${task.cpus} ${params.sample}.nodup.bam
    printf "reads\\tn\\nall_primary\\t%s\\nduplicate_flagged\\t%s\\n" \\
        \$(samtools view -@ ${task.cpus} -c -F 0x904 ${bam}) \\
        \$(samtools view -@ ${task.cpus} -c -F 0x904 -f 0x400 ${bam}) > ${params.sample}.duplicates.tsv
    """
}
