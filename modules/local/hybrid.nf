/*
 * Combine read blocks with the population phase, per chromosome; then the genome-wide steps.
 */

// read blocks as soft constraints in the merge, with the same error rate as SHAPEIT4 --use-PS by default
def blockErrorArg() {
    def e = params.read_block_error != null ? params.read_block_error : params.shapeit4_use_ps
    return (e as double) > 0 ? "--block-error ${e}" : ''
}

process HYBRID_MERGE {
    tag "${combo}:${chrom}"
    label 'phase'

    input:
    tuple val(combo), val(pop_phaser), val(chrom), path(reads_vcf), path(reads_tbi), path(pop_vcf, stageAs: "pop/*"), path(pop_tbi, stageAs: "pop/*")
    tuple path(het_vcf), path(het_tbi)
    tuple path(bam), path(bai)
    tuple path(ref), path(fai)

    output:
    tuple val(combo), val(chrom), path("hybrid_${combo}_${chrom}.vcf.gz"), path("hybrid_${combo}_${chrom}.vcf.gz.tbi"), emit: vcf
    tuple val(combo), val(chrom), path("orientation_${combo}_${chrom}.tsv"), emit: orientation

    script:
    def out = "hybrid_${combo}_${chrom}"
    if (pop_phaser == 'none') {
        """
        cp ${reads_vcf} ${out}.vcf.gz
        tabix -f -p vcf ${out}.vcf.gz
        hybrid_merge.py --reads ${reads_vcf} --blocks orientation_${combo}_${chrom}.tsv
        """
    } else if (params.hybrid_merge == 'whatshap') {
        """
        # reads + population VCF as an extra 'pseudo-read';
        # also re-attaches the SNPs that are not in the panel
        bcftools view -r ${chrom} ${het_vcf} -Oz -o input_${chrom}.vcf.gz
        tabix -p vcf input_${chrom}.vcf.gz
        whatshap phase --reference ${ref} --ignore-read-groups --chromosome ${chrom} \\
            -o ${out}.vcf.gz input_${chrom}.vcf.gz ${bam} ${pop_vcf}
        tabix -f -p vcf ${out}.vcf.gz
        hybrid_merge.py --reads ${reads_vcf} --pop ${pop_vcf} --min-snps ${params.orient_min_snps} ${blockErrorArg()} \\
            --blocks orientation_${combo}_${chrom}.tsv
        rm -f input_${chrom}.vcf.gz*
        """
    } else {
        """
        # flip each read block as a whole to agree with the population phase;
        # SNPs outside the panel keep their read phase, panel-only SNPs take the population phase
        hybrid_merge.py --reads ${reads_vcf} --pop ${pop_vcf} --min-snps ${params.orient_min_snps} ${blockErrorArg()} \\
            --out ${out}.vcf.gz --blocks orientation_${combo}_${chrom}.tsv
        tabix -f -p vcf ${out}.vcf.gz
        """
    }
}

process JUNCTION_SPLIT {
    tag "${read_phaser}:${chrom}"
    label 'light'

    input:
    tuple val(read_phaser), val(chrom), path(reads_vcf), path(reads_tbi)

    output:
    tuple val(read_phaser), val(chrom), path("split_${read_phaser}_${chrom}.vcf.gz"), path("split_${read_phaser}_${chrom}.vcf.gz.tbi"), emit: vcf
    tuple val(read_phaser), val(chrom), path("split_${read_phaser}_${chrom}.vcf.gz.splits.tsv"), emit: splits

    script:
    """
    junction_test.py split ${reads_vcf} split_${read_phaser}_${chrom}.vcf.gz --min-snps ${params.junction_min_snps}
    tabix -f -p vcf split_${read_phaser}_${chrom}.vcf.gz
    """
}

process JUNCTION_SCORE {
    tag "${combo}:${chrom}"
    label 'light'

    input:
    tuple val(combo), val(chrom), path(reads_vcf), path(reads_tbi), path(split_vcf), path(split_tbi), path(splits), path(pop_vcf), path(pop_tbi)

    output:
    tuple val(combo), path("junctions_${combo}_${chrom}.tsv"), emit: tsv

    script:
    """
    # re-join the artificially split blocks with the population phase only, then check against the reads
    hybrid_merge.py --reads ${split_vcf} --pop ${pop_vcf} --min-snps ${params.orient_min_snps} ${blockErrorArg()} \\
        --out rejoined.vcf.gz --blocks /dev/null
    junction_test.py score ${splits} ${reads_vcf} rejoined.vcf.gz junctions_${combo}_${chrom}.tsv
    rm -f rejoined.vcf.gz
    """
}

process CONCAT_PHASED {
    tag "${id}"
    label 'light'
    publishDir "${params.outdir}/${params.sample}/${pubdir}", mode: 'copy'

    input:
    tuple val(id), val(pubdir), path(vcfs), path(tbis)

    output:
    tuple val(id), path("${id}.phased.vcf.gz"), path("${id}.phased.vcf.gz.tbi"), emit: vcf
    path "${id}.phases_for_ascat.vcf.gz", emit: ascat
    tuple val(id), path("qc/${id}.phase_stats.tsv"), path("qc/${id}.blocklist.tsv"), emit: stats
    path "qc/${id}.whatshap_stats.*", emit: whatshap_stats

    script:
    // vcfs arrive sorted in chromosome order (see main.nf)
    def ordered = vcfs instanceof List ? vcfs : [vcfs]
    """
    mkdir -p qc
    printf '%s\\n' ${ordered.join(' ')} > files.list
    bcftools concat -f files.list -Oz -o ${id}.phased.vcf.gz
    tabix -f -p vcf ${id}.phased.vcf.gz

    # ASCAT.sc readPhases() uses data.table::fread: no ## header, GT as the only FORMAT field
    bcftools annotate -x INFO,^FORMAT/GT ${id}.phased.vcf.gz | grep -v '^##' | gzip > ${id}.phases_for_ascat.vcf.gz

    whatshap stats --tsv=qc/${id}.whatshap_stats.tsv ${id}.phased.vcf.gz > qc/${id}.whatshap_stats.txt
    phase_stats.py ${id}.phased.vcf.gz qc/${id}.phase_stats.tsv qc/${id}.blocklist.tsv
    """
}

process COMPARE_READS_FINAL {
    tag "${combo}"
    label 'light'
    publishDir "${params.outdir}/${params.sample}/03_hybrid_phasing/${combo}/qc", mode: 'copy'

    input:
    tuple val(combo), path(final_vcf), path(final_tbi), path(reads_vcf), path(reads_tbi)

    output:
    tuple val(combo), path("${combo}.reads_vs_final.tsv"), path("${combo}.reads_vs_final.bed"), emit: tsv
    path "${combo}.reads_vs_final.txt"

    script:
    """
    # within read blocks the final phasing should (almost) never disagree with the reads
    whatshap compare --names reads_only,final --ignore-sample-name --only-snvs \\
        --tsv-pairwise ${combo}.reads_vs_final.tsv \\
        --switch-error-bed ${combo}.reads_vs_final.bed \\
        ${reads_vcf} ${final_vcf} > ${combo}.reads_vs_final.txt
    """
}
