/*
 * Population-based phasing per chromosome.
 *
 *  SHAPEIT4: takes the read-phased VCF and keeps the read blocks intact (--use-PS); only panel SNPs come out.
 *  Beagle  : phases the unphased het SNPs from the panel alone (it cannot use read blocks);
 *            the read blocks are oriented to it afterwards in HYBRID_MERGE.
 *
 *  The panel / map paths are plain strings so that the copies baked into the container
 *  (/opt/ref_data, see containers/hybrid_phasing.def) can be used directly.
 */

def panelPath(chrom) { "${params.panel_dir}/${params.panel_name.replace('{chr}', chrom)}" }
def mapPath(chrom)   { "${params.map_dir}/${params.map_name.replace('{chr}', chrom)}" }

process SHAPEIT4 {
    tag "${read_phaser}:${chrom}:${mode}"
    label 'pop'

    input:
    tuple val(read_phaser), val(chrom), path(vcf), path(tbi), val(mode)   // mode: full | split (junction test)

    output:
    tuple val(read_phaser), val('shapeit4'), val(chrom), val(mode), path("pop_${read_phaser}_shapeit4_${chrom}_${mode}.vcf.gz"), path("pop_${read_phaser}_shapeit4_${chrom}_${mode}.vcf.gz.tbi"), emit: vcf

    script:
    def out = "pop_${read_phaser}_shapeit4_${chrom}_${mode}"
    """
    # SHAPEIT4 needs AC/AN in the target
    bcftools +fill-tags ${vcf} -Ob -o tagged.bcf -- -t AN,AC
    bcftools index -f tagged.bcf

    # keeps the read blocks intact and orients them with the panel (do NOT also pass --scaffold)
    shapeit4 --input tagged.bcf \\
        --reference ${panelPath(chrom)} \\
        --map ${mapPath(chrom)} \\
        --region ${chrom} \\
        --use-PS ${params.shapeit4_use_ps} \\
        --thread ${task.cpus} \\
        ${params.shapeit4_args} \\
        --output ${out}.bcf

    bcftools view ${out}.bcf -Oz -o ${out}.vcf.gz
    tabix -f -p vcf ${out}.vcf.gz
    rm -f tagged.bcf* ${out}.bcf
    """
}

process BEAGLE {
    tag "${chrom}"
    label 'pop'

    input:
    tuple val(chrom), path(het_vcf), path(het_tbi)

    output:
    tuple val(chrom), path("pop_beagle_${chrom}.vcf.gz"), path("pop_beagle_${chrom}.vcf.gz.tbi"), emit: vcf

    script:
    def mem_gb = Math.max(2, (task.memory.toGiga() * 0.85) as int)
    // pre-made Beagle panel (vcf.gz / bref3) if given, else convert the SHAPEIT4 BCF panel
    def beagle_panel = params.beagle_panel_dir ? "${params.beagle_panel_dir}/${params.beagle_panel_name.replace('{chr}', chrom)}" : null
    def ref_file = beagle_panel ?: "ref_${chrom}.vcf.gz"
    def ref_cmd  = beagle_panel ? "" : "bcftools view ${panelPath(chrom)} -Oz -o ref_${chrom}.vcf.gz --threads ${task.cpus}"
    """
    bcftools view -r ${chrom} ${het_vcf} -Oz -o target_${chrom}.vcf.gz
    tabix -p vcf target_${chrom}.vcf.gz

    # reference panel: Beagle reads vcf.gz / bref3, not BCF
    ${ref_cmd}

    # SHAPEIT4 map (pos chr cM, with header) -> PLINK map (chr . cM pos), chromosome named like the VCF
    zcat -f ${mapPath(chrom)} | awk -v c=${chrom} 'NR>1 {print c, ".", \$3, \$1}' OFS='\\t' > ${chrom}.plink.map

    beagle -Xmx${mem_gb}g \\
        gt=target_${chrom}.vcf.gz \\
        ref=${ref_file} \\
        map=${chrom}.plink.map \\
        chrom=${chrom} \\
        impute=false \\
        nthreads=${task.cpus} \\
        out=pop_beagle_${chrom} \\
        ${params.beagle_args}

    tabix -f -p vcf pop_beagle_${chrom}.vcf.gz
    rm -f ref_${chrom}.vcf.gz target_${chrom}.vcf.gz*
    """
}
