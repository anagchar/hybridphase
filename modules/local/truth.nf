/*
 * Optional benchmark against a haplotype-resolved diploid assembly of the same individual.
 *
 *   --hap1 / --hap2          haplotype FASTAs (e.g. hifiasm hap1/hap2 or paternal/maternal)
 *                            -> dipcall -> truth VCF (phased GT) + confident regions BED
 *   --truth_vcf [--truth_bed] an existing truth VCF instead (dipcall is skipped)
 *
 * Every phased VCF of the run (read-only and hybrid) is compared with WhatsHap compare.
 */

process DIPCALL {
    tag "${params.sample}"
    label 'dipcall'
    publishDir "${params.outdir}/${params.sample}/04_qc/truth", mode: 'copy'

    input:
    tuple path(hap1), path(hap2)
    tuple path(ref), path(fai)

    output:
    tuple path("truth.dip.vcf.gz"), path("truth.dip.bed"), emit: truth

    script:
    // -t is per minimap2 process; make -j2 runs the two haplotypes at once
    def t = Math.max(1, (task.cpus / 2) as int)
    """
    run-dipcall -t ${t} ${params.dipcall_args} truth ${ref} ${hap1} ${hap2} > truth.mak
    make -j2 -f truth.mak
    """
}

process TRUTH_CONFIDENT {
    tag "${params.sample}"
    label 'light'
    publishDir "${params.outdir}/${params.sample}/04_qc/truth", mode: 'copy'

    input:
    tuple path(vcf), path(bed)

    output:
    tuple path("truth.confident.vcf.gz"), path("truth.confident.vcf.gz.tbi"), emit: vcf

    script:
    // NO_FILE = no BED given with --truth_vcf
    def regions = bed.name == 'NO_FILE' ? '' : "-T ${bed}"
    """
    # confident calls only: inside the BED and unflagged (dipcall writes "." for good calls)
    bcftools view -f .,PASS ${regions} -Oz -o truth.confident.vcf.gz ${vcf}
    tabix -f -p vcf truth.confident.vcf.gz
    """
}

process COMPARE_TRUTH {
    tag "${id}"
    label 'light'
    publishDir "${params.outdir}/${params.sample}/04_qc/truth", mode: 'copy'

    input:
    tuple val(id), path(vcf), path(tbi), path(truth, stageAs: 'truth/*'), path(truth_tbi, stageAs: 'truth/*')

    output:
    tuple val(id), path("${id}.truth.tsv"), path("${id}.truth_switches.bed"), emit: tsv
    path "${id}.truth.txt"

    script:
    """
    # sample names differ (dipcall: syndip), hence --ignore-sample-name
    whatshap compare --names truth,${id} --ignore-sample-name --only-snvs \\
        --tsv-pairwise ${id}.truth.tsv \\
        --switch-error-bed ${id}.truth_switches.bed \\
        ${truth} ${vcf} > ${id}.truth.txt
    """
}
