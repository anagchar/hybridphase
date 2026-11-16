#!/usr/bin/env nextflow
/*
 * splongget-phasing
 *
 *   BAM ── Clair3 ── germline het SNPs ─┬─ read phasing (WhatsHap | LongPhase)  ── read blocks (PS)
 *                                       └─ population phasing (SHAPEIT4 | Beagle) ─┐
 *                                            read blocks + population phase ── hybrid merge ── chromosome-scale phase
 *                                                                   ── haplotag BAM ── QC + HTML report
 *
 * Every combination of --read_phaser x --pop_phaser is run; the report compares them.
 */

include { CLAIR3; FILTER_HET_SNPS }                      from './modules/local/variant_calling'
include { WHATSHAP_PHASE; LONGPHASE_PHASE }              from './modules/local/read_phasing'
include { SHAPEIT4; BEAGLE }                             from './modules/local/pop_phasing'
include { HYBRID_MERGE; JUNCTION_SPLIT; JUNCTION_SCORE;
          CONCAT_PHASED; COMPARE_READS_FINAL;
          METHOD_CONCORDANCE }                           from './modules/local/hybrid'
include { HAPLOTAG; HAPLOTAG_COUNTS }                    from './modules/local/haplotag'
include { CELLSNP; ASCAT_TABLES }                        from './modules/local/single_cell'
include { REPORT }                                       from './modules/local/report'

// '1-22' / '21,22' / 'chr1-chr5,X' -> [chr1, ...]
def parseChroms(spec, prefix) {
    spec.toString().tokenize(',').collectMany { t ->
        t = t.trim().replaceFirst(/^${prefix}/, '')
        if (t ==~ /\d+-\d+/) {
            def (a, b) = t.tokenize('-').collect { it as int }
            return (a..b).collect { "${prefix}${it}".toString() }
        }
        return ["${prefix}${t}".toString()]
    }
}

def splitList(v) { v.toString().tokenize(',')*.trim().findAll { it } }

def findIndex(String path, List suffixes) {
    def hit = suffixes.collect { file(it(path)) }.find { it.exists() }
    if (!hit) error "No index found for ${path} (tried: ${suffixes.collect { it(path) }.join(', ')})"
    return hit
}

workflow {
    // ---------------------------------------------------------------- parameters
    if (!params.bam) error "Please provide --bam"
    if (!params.ref) error "Please provide --ref"
    if (!params.vcf && !params.clair3_model_dir) error "Please provide --clair3_model_dir (or an existing Clair3 VCF with --vcf)"

    def read_phasers = splitList(params.read_phaser)
    def pop_phasers  = splitList(params.pop_phaser)
    def bad_r = read_phasers - ['whatshap', 'longphase']
    def bad_p = pop_phasers - ['shapeit4', 'beagle', 'none']
    if (bad_r) error "Unknown --read_phaser: ${bad_r}. Use whatshap and/or longphase"
    if (bad_p) error "Unknown --pop_phaser: ${bad_p}. Use shapeit4, beagle and/or none"
    if (!(params.hybrid_merge in ['orient', 'whatshap'])) error "--hybrid_merge must be orient or whatshap"
    if (!(params.haplotag_tool in ['whatshap', 'longphase'])) error "--haplotag_tool must be whatshap or longphase"

    def chroms     = parseChroms(params.chroms, params.chr_prefix)
    def chromIndex = chroms.withIndex().collectEntries { c, i -> [(c): i] }
    def combos     = [read_phasers, pop_phasers].combinations().collect { "${it[0]}_${it[1]}".toString() }
    def ht_combo   = params.haplotag_combo ?: combos[0]
    if (!(ht_combo in combos)) error "--haplotag_combo ${ht_combo} is not one of the combinations run: ${combos}"

    log.info """
    ==============================================================
     splongget-phasing ${workflow.manifest.version}
    ==============================================================
     sample        : ${params.sample}
     bam           : ${params.bam}
     reference     : ${params.ref}
     variant calls : ${params.vcf ?: "Clair3 (${params.clair3_model})"}
     chromosomes   : ${chroms.first()} .. ${chroms.last()} (${chroms.size()})
     combinations  : ${combos.join(', ')}  (merge: ${params.hybrid_merge})
     haplotag      : ${params.haplotag ? "${ht_combo} with ${params.haplotag_tool}" : 'no'}
     outdir        : ${params.outdir}/${params.sample}
    ==============================================================
    """.stripIndent()

    def bam = file(params.bam, checkIfExists: true)
    def bai = findIndex(params.bam, [{ "${it}.bai" }, { it.replaceFirst(/\.bam$/, '.bai') }, { "${it}.csi" }])
    def ref = file(params.ref, checkIfExists: true)
    def fai = file("${params.ref}.fai", checkIfExists: true)
    ch_bam = Channel.value([bam, bai])
    ch_ref = Channel.value([ref, fai])

    // ---------------------------------------------------------------- 1. variant calling
    if (params.vcf) {
        def tbi = findIndex(params.vcf, [{ "${it}.tbi" }, { "${it}.csi" }])
        ch_vcf = Channel.value([params.sample, file(params.vcf, checkIfExists: true), tbi])
    } else {
        CLAIR3(Channel.value([params.sample, bam, bai]), ch_ref, file(params.clair3_model_dir, checkIfExists: true))
        ch_vcf = CLAIR3.out.vcf
    }

    FILTER_HET_SNPS(ch_vcf)
    ch_het       = FILTER_HET_SNPS.out.vcf.map { s, v, t -> [v, t] }.first()
    ch_chrom_het = Channel.fromList(chroms).combine(ch_het)                    // [chrom, vcf, tbi]

    // ---------------------------------------------------------------- 2. read-based phasing
    ch_reads = Channel.empty()                                                 // [read_phaser, chrom, vcf, tbi]
    if ('whatshap' in read_phasers) {
        WHATSHAP_PHASE(ch_chrom_het, ch_bam, ch_ref)
        ch_reads = ch_reads.mix(WHATSHAP_PHASE.out.vcf)
    }
    if ('longphase' in read_phasers) {
        LONGPHASE_PHASE(ch_chrom_het, ch_bam, ch_ref)
        ch_reads = ch_reads.mix(LONGPHASE_PHASE.out.vcf)
    }

    // truth-free junction test: split every large read block in two, hide the link
    def do_junction = params.junction_test && pop_phasers.any { it != 'none' }
    if (do_junction) {
        JUNCTION_SPLIT(ch_reads)
        ch_split = JUNCTION_SPLIT.out.vcf
    } else {
        ch_split = Channel.empty()
    }

    // ---------------------------------------------------------------- 3. population-based phasing
    ch_pop       = Channel.empty()                                             // [read_phaser, pop_phaser, chrom, vcf, tbi]
    ch_pop_split = Channel.empty()
    if ('shapeit4' in pop_phasers) {
        SHAPEIT4(ch_reads.map { rp, c, v, t -> [rp, c, v, t, 'full'] }
                         .mix(ch_split.map { rp, c, v, t -> [rp, c, v, t, 'split'] }))
        ch_pop       = ch_pop.mix(SHAPEIT4.out.vcf.filter { it[3] == 'full' }.map { rp, pp, c, m, v, t -> [rp, pp, c, v, t] })
        ch_pop_split = ch_pop_split.mix(SHAPEIT4.out.vcf.filter { it[3] == 'split' }.map { rp, pp, c, m, v, t -> [rp, pp, c, v, t] })
    }
    if ('beagle' in pop_phasers) {
        BEAGLE(ch_chrom_het)
        // Beagle never sees the reads: one run per chromosome serves every read phaser and the junction test
        ch_beagle    = ch_reads.map { rp, c, v, t -> [c, rp] }.combine(BEAGLE.out.vcf, by: 0)
                               .map { c, rp, v, t -> [rp, 'beagle', c, v, t] }
        ch_pop       = ch_pop.mix(ch_beagle)
        ch_pop_split = ch_pop_split.mix(ch_beagle)
    }
    if ('none' in pop_phasers) {
        ch_pop = ch_pop.mix(ch_reads.map { rp, c, v, t -> [rp, 'none', c, v, t] })
    }

    // ---------------------------------------------------------------- 4. hybrid merge per chromosome
    ch_merge_in = ch_reads.map { rp, c, v, t -> [[rp, c], v, t] }
        .combine(ch_pop.map { rp, pp, c, v, t -> [[rp, c], pp, v, t] }, by: 0)
        .map { k, rv, rt, pp, pv, pt -> ["${k[0]}_${pp}".toString(), pp, k[1], rv, rt, pv, pt] }
    HYBRID_MERGE(ch_merge_in, ch_het, ch_bam, ch_ref)

    // ---------------------------------------------------------------- 5. genome-wide VCFs (read-only and hybrid)
    def byChrom = { id, pubdir, cs, vs, ts ->
        def o = (0..<cs.size()).sort { chromIndex[cs[it]] }
        [id, pubdir, o.collect { vs[it] }, o.collect { ts[it] }]
    }
    ch_concat_in = HYBRID_MERGE.out.vcf
        .groupTuple(by: 0, size: chroms.size())
        .map { combo, cs, vs, ts -> byChrom(combo, "03_hybrid_phasing/${combo}", cs, vs, ts) }
        .mix(ch_reads.map { rp, c, v, t -> ["${rp}_reads".toString(), c, v, t] }
                     .groupTuple(by: 0, size: chroms.size())
                     .map { id, cs, vs, ts -> byChrom(id, "02_read_phasing/${id - ~/_reads$/}", cs, vs, ts) })
    CONCAT_PHASED(ch_concat_in)

    ch_final = CONCAT_PHASED.out.vcf.filter { it[0] in combos }                // [combo, vcf, tbi]
    ch_reads_genome = CONCAT_PHASED.out.vcf.filter { it[0].endsWith('_reads') }
        .map { id, v, t -> [id - ~/_reads$/, v, t] }

    // ---------------------------------------------------------------- 6. truth-free QC
    COMPARE_READS_FINAL(
        ch_final.filter { !it[0].endsWith('_none') }
                .map { id, v, t -> [id.tokenize('_')[0], id, v, t] }
                .combine(ch_reads_genome, by: 0)
                .map { rp, id, v, t, rv, rt -> [id, v, t, rv, rt] })

    if (combos.size() > 1) {
        METHOD_CONCORDANCE(ch_final.toSortedList { a, b -> a[0] <=> b[0] }
                                   .map { l -> [l.collect { it[0] }, l.collect { it[1] }, l.collect { it[2] }] })
        ch_concordance = METHOD_CONCORDANCE.out.tsv
    } else {
        ch_concordance = Channel.empty()
    }

    if (do_junction) {
        ch_junction_in = ch_reads.map { rp, c, v, t -> [[rp, c], v, t] }
            .join(JUNCTION_SPLIT.out.vcf.map { rp, c, v, t -> [[rp, c], v, t] })
            .join(JUNCTION_SPLIT.out.splits.map { rp, c, s -> [[rp, c], s] })
            .combine(ch_pop_split.map { rp, pp, c, v, t -> [[rp, c], pp, v, t] }, by: 0)
            .map { k, rv, rt, sv, st, s, pp, pv, pt -> ["${k[0]}_${pp}".toString(), k[1], rv, rt, sv, st, s, pv, pt] }
        JUNCTION_SCORE(ch_junction_in)
        ch_junctions = JUNCTION_SCORE.out.tsv
            .collectFile(keepHeader: true, storeDir: "${params.outdir}/${params.sample}/04_qc/junction_test") { combo, f -> ["junctions_${combo}.tsv", f] }
    } else {
        ch_junctions = Channel.empty()
    }

    ch_orientation = HYBRID_MERGE.out.orientation
        .collectFile(keepHeader: true, storeDir: "${params.outdir}/${params.sample}/04_qc/block_orientation") { combo, c, f -> ["orientation_${combo}.tsv", f] }

    // ---------------------------------------------------------------- 7. haplotag
    ch_ht_vcf = ch_final.filter { it[0] == ht_combo }
    if (params.haplotag) {
        HAPLOTAG(ch_ht_vcf, ch_bam, ch_ref)
        HAPLOTAG_COUNTS(HAPLOTAG.out.bam)
        ch_ht_counts = HAPLOTAG_COUNTS.out.tsv
    } else {
        ch_ht_counts = Channel.empty()
    }

    // ---------------------------------------------------------------- 8. optional single-cell allele counts
    if (params.barcodes) {
        CELLSNP(ch_ht_vcf, ch_bam, file(params.barcodes, checkIfExists: true))
        ASCAT_TABLES(CELLSNP.out.dir)
        ch_cell_qc = ASCAT_TABLES.out.qc
    } else {
        ch_cell_qc = Channel.empty()
    }

    // ---------------------------------------------------------------- 9. report
    def run_info = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        sample: params.sample, bam: params.bam, ref: params.ref,
        variant_calls: params.vcf ?: "Clair3 ${params.clair3_model}".toString(),
        platform: params.platform, chroms: chroms, combos: combos, hybrid_merge: params.hybrid_merge,
        haplotag_combo: params.haplotag ? ht_combo : null, haplotag_tool: params.haplotag_tool,
        haplotag_bin_mb: params.haplotag_bin_mb,
        filters: [min_qual: params.min_qual, min_dp: params.min_dp, min_vaf: params.min_vaf, max_vaf: params.max_vaf],
        version: workflow.manifest.version, run_name: workflow.runName, started: workflow.start.toString()
    ]))

    ch_report_files = FILTER_HET_SNPS.out.funnel
        .mix(FILTER_HET_SNPS.out.calls, FILTER_HET_SNPS.out.kept)
        .mix(CONCAT_PHASED.out.stats.flatMap { id, s, b -> [s, b] })
        .mix(COMPARE_READS_FINAL.out.tsv.flatMap { id, t, b -> [t, b] })
        .mix(ch_orientation, ch_junctions, ch_concordance, ch_ht_counts, ch_cell_qc)
        .collect()

    REPORT(ch_report_files, fai, run_info)
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Report: ${params.outdir}/${params.sample}/04_qc/phasing_report.html\n"
        : "\nFailed: ${workflow.errorMessage}\n")
}
