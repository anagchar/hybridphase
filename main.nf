#!/usr/bin/env nextflow
/*
 * hybridphase
 *
 *   BAM ── Clair3 ── germline het SNPs ─┬─ read phasing (WhatsHap | LongPhase)  ── read blocks (PS)
 *                                       └─ population phasing (SHAPEIT4 | Beagle) ─┐
 *                                            read blocks + population phase ── hybrid merge ── chromosome-scale phase
 *                                                                   ── haplotag BAM ── QC + HTML report
 *
 * Every combination of --read_phaser x --pop_phaser is run (usually one of each).
 */

include { DROP_DUPLICATES; CLAIR3; FILTER_HET_SNPS }     from './modules/local/variant_calling'
include { WHATSHAP_PHASE; LONGPHASE_PHASE }              from './modules/local/read_phasing'
include { SHAPEIT4; BEAGLE }                             from './modules/local/pop_phasing'
include { HYBRID_MERGE; JUNCTION_SPLIT; JUNCTION_SCORE;
          CONCAT_PHASED; COMPARE_READS_FINAL }                           from './modules/local/hybrid'
include { HAPLOTAG; HAPLOTAG_COUNTS }                    from './modules/local/haplotag'
include { CELLSNP; ASCAT_TABLES }                        from './modules/local/single_cell'
include { DIPCALL; TRUTH_CONFIDENT; COMPARE_TRUTH }       from './modules/local/truth'
include { ASM_INDEX; ASM_MAP; HIPHAP; ASM_TAG }         from './modules/local/assembly_haplotag'
include { REPORT }                                       from './modules/local/report'

// '1-22' / '21,22' / 'chr1-chr5,X' -> [chr1, ...]
def parseChroms(spec, prefix) {
    spec.toString().tokenize(',').collectMany { tok ->
        def t = tok.trim()
        if (prefix) t = t.replace(prefix.toString(), '')
        if (t ==~ '[0-9]+-[0-9]+') {
            def ab = t.tokenize('-').collect { x -> x as int }
            return (ab[0]..ab[1]).collect { n -> "${prefix}${n}".toString() }
        }
        return ["${prefix}${t}".toString()]
    }
}

def splitList(v) { v.toString().tokenize(',')*.trim().findAll { it } }

def findIndex(String path, List candidates) {
    def hit = candidates.collect { c -> file(c) }.find { f -> f.exists() }
    if (!hit) error "No index found for ${path} (tried: ${candidates.join(', ')})"
    return hit
}

// per-chromosome VCFs of one id -> [id, publish dir, vcfs, tbis] in chromosome order
def byChrom(id, pubdir, cs, vs, ts, chromIndex) {
    def o = (0..<cs.size()).sort { i -> chromIndex[cs[i]] }
    return [id, pubdir, o.collect { i -> vs[i] }, o.collect { i -> ts[i] }]
}

def stripSuffix(String s, String suf) { s.endsWith(suf) ? s.substring(0, s.length() - suf.length()) : s }

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
    if ((params.hap1 && !params.hap2) || (params.hap2 && !params.hap1)) error "--hap1 and --hap2 go together (haplotype-resolved assembly)"
    if (params.truth_vcf && params.hap1) error "Give either --hap1/--hap2 or --truth_vcf, not both"
    if (!(params.haplotag_tool in ['whatshap', 'longphase'])) error "--haplotag_tool must be whatshap or longphase"

    def chroms     = parseChroms(params.chroms, params.chr_prefix)
    def chromIndex = chroms.withIndex().collectEntries { c, i -> [(c): i] }
    def combos     = [read_phasers, pop_phasers].combinations().collect { "${it[0]}_${it[1]}".toString() }
    def ht_combo   = params.haplotag_combo ?: combos[0]
    if (!(ht_combo in combos)) error "--haplotag_combo ${ht_combo} is not one of the combinations run: ${combos}"

    log.info """
    ==============================================================
     hybridphase ${workflow.manifest.version}
    ==============================================================
     sample        : ${params.sample}
     bam           : ${params.bam}
     reference     : ${params.ref}
     variant calls : ${params.vcf ?: "Clair3 (${params.clair3_model})"}
     chromosomes   : ${chroms.first()} .. ${chroms.last()} (${chroms.size()})
     combinations  : ${combos.join(', ')}  (merge: ${params.hybrid_merge})
     drop dups     : ${params.drop_duplicates}
     truth         : ${params.hap1 && params.assembly_truth ? 'dipcall on --hap1/--hap2' : params.truth_vcf ?: 'none'}
     assembly tags : ${params.hap1 && params.assembly_haplotag ? 'HipHap on --hap1 / --hap2 (under development)' : 'no'}
     haplotag      : ${params.haplotag ? "${ht_combo} with ${params.haplotag_tool}" : 'no'}
     outdir        : ${params.outdir}/${params.sample}
    ==============================================================
    """.stripIndent()

    def bam = file(params.bam, checkIfExists: true)
    def bai = findIndex(params.bam, ["${params.bam}.bai", "${stripSuffix(params.bam.toString(), '.bam')}.bai", "${params.bam}.csi"])
    def ref = file(params.ref, checkIfExists: true)
    def fai = file("${params.ref}.fai", checkIfExists: true)
    ch_bam_in = Channel.value([bam, bai])           // original BAM: haplotagged and used for single cells
    ch_ref    = Channel.value([ref, fai])

    // ---------------------------------------------------------------- 0. duplicate-flagged reads out
    if (params.drop_duplicates) {
        DROP_DUPLICATES(ch_bam_in)
        ch_bam     = DROP_DUPLICATES.out.bam.first()   // for variant calling and phasing
        ch_dup_tsv = DROP_DUPLICATES.out.stats
    } else {
        ch_bam     = ch_bam_in
        ch_dup_tsv = Channel.empty()
    }

    // ---------------------------------------------------------------- 1. variant calling
    if (params.vcf) {
        def tbi = findIndex(params.vcf, ["${params.vcf}.tbi", "${params.vcf}.csi"])
        ch_vcf = Channel.value([params.sample, file(params.vcf, checkIfExists: true), tbi])
    } else {
        CLAIR3(ch_bam.map { b, i -> [params.sample, b, i] }, ch_ref, file(params.clair3_model_dir, checkIfExists: true))
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
    ch_concat_in = HYBRID_MERGE.out.vcf
        .groupTuple(by: 0, size: chroms.size())
        .map { combo, cs, vs, ts -> byChrom(combo, "03_hybrid_phasing/${combo}", cs, vs, ts, chromIndex) }
        .mix(ch_reads.map { rp, c, v, t -> ["${rp}_reads".toString(), c, v, t] }
                     .groupTuple(by: 0, size: chroms.size())
                     .map { id, cs, vs, ts -> byChrom(id, "02_read_phasing/${stripSuffix(id, '_reads')}", cs, vs, ts, chromIndex) })
    CONCAT_PHASED(ch_concat_in)

    ch_final = CONCAT_PHASED.out.vcf.filter { it[0] in combos }                // [combo, vcf, tbi]
    ch_reads_genome = CONCAT_PHASED.out.vcf.filter { it[0].endsWith('_reads') }
        .map { id, v, t -> [stripSuffix(id, '_reads'), v, t] }

    // ---------------------------------------------------------------- 6. truth-free QC
    COMPARE_READS_FINAL(
        ch_final.filter { !it[0].endsWith('_none') }
                .map { id, v, t -> [id.tokenize('_')[0], id, v, t] }
                .combine(ch_reads_genome, by: 0)
                .map { rp, id, v, t, rv, rt -> [id, v, t, rv, rt] })

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

    // ---------------------------------------------------------------- 6b. optional truth benchmark
    if ((params.hap1 && params.assembly_truth) || params.truth_vcf) {
        if (params.hap1) {
            DIPCALL(Channel.value([file(params.hap1, checkIfExists: true), file(params.hap2, checkIfExists: true)]), ch_ref)
            ch_truth_raw = DIPCALL.out.truth
        } else {
            ch_truth_raw = Channel.value([file(params.truth_vcf, checkIfExists: true),
                                          params.truth_bed ? file(params.truth_bed, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")])
        }
        TRUTH_CONFIDENT(ch_truth_raw)
        // every phased VCF of the run: read-only blocks and all hybrid combinations
        COMPARE_TRUTH(CONCAT_PHASED.out.vcf.combine(TRUTH_CONFIDENT.out.vcf))
        ch_truth_files = COMPARE_TRUTH.out.tsv.flatMap { id, t, b -> [t, b] }
    } else {
        ch_truth_files = Channel.empty()
    }

    // ---------------------------------------------------------------- 7. haplotag
    ch_ht_vcf = ch_final.filter { it[0] == ht_combo }
    if (params.haplotag) {
        HAPLOTAG(ch_ht_vcf, ch_bam_in, ch_ref)
        HAPLOTAG_COUNTS(HAPLOTAG.out.bam)
        ch_ht_counts = HAPLOTAG_COUNTS.out.tsv
    } else {
        ch_ht_counts = Channel.empty()
    }

    // ---------------------------------------------------------------- 7b. optional assembly-based haplotagging
    if (params.hap1 && params.assembly_haplotag) {
        ASM_INDEX(Channel.value([file(params.hap1, checkIfExists: true), file(params.hap2, checkIfExists: true)]))
        ASM_MAP(ch_bam_in, ASM_INDEX.out.mmi)
        HIPHAP(ASM_MAP.out.bams)
        // tag the VCF-haplotagged BAM when there is one: its HP/PS give the agreement table, then get replaced
        ASM_TAG(HIPHAP.out.assignments, params.haplotag ? HAPLOTAG.out.bam.map { c, b, i -> [b, i] } : ch_bam_in)
        ch_asm_files = HIPHAP.out.by_length.mix(HIPHAP.out.hapq, ASM_TAG.out.bins, ASM_TAG.out.agreement,
                                                ASM_TAG.out.cells)
    } else {
        ch_asm_files = Channel.empty()
    }

    // ---------------------------------------------------------------- 8. optional single-cell allele counts
    if (params.barcodes) {
        CELLSNP(ch_ht_vcf, ch_bam_in, file(params.barcodes, checkIfExists: true))
        ASCAT_TABLES(CELLSNP.out.dir)
        ch_cell_qc = ASCAT_TABLES.out.qc
    } else {
        ch_cell_qc = Channel.empty()
    }

    // ---------------------------------------------------------------- 9. report
    def run_info = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        sample: params.sample, bam: params.bam, ref: params.ref,
        variant_calls: params.vcf ?: "Clair3 ${params.clair3_model}".toString(),
        platform: params.platform, drop_duplicates: params.drop_duplicates, chroms: chroms, combos: combos, hybrid_merge: params.hybrid_merge,
        haplotag_combo: params.haplotag ? ht_combo : null, haplotag_tool: params.haplotag_tool,
        haplotag_bin_mb: params.haplotag_bin_mb, asm_min_hapq: params.asm_min_hapq,
        truth: params.hap1 ? "dipcall: ${params.hap1} + ${params.hap2}".toString() : params.truth_vcf,
        filters: [min_qual: params.min_qual, min_dp: params.min_dp, min_vaf: params.min_vaf, max_vaf: params.max_vaf],
        version: workflow.manifest.version, run_name: workflow.runName, started: workflow.start.toString()
    ]))

    ch_report_files = FILTER_HET_SNPS.out.funnel
        .mix(FILTER_HET_SNPS.out.calls, FILTER_HET_SNPS.out.kept)
        .mix(CONCAT_PHASED.out.stats.flatMap { id, s, b -> [s, b] })
        .mix(COMPARE_READS_FINAL.out.tsv.flatMap { id, t, b -> [t, b] })
        .mix(ch_orientation, ch_junctions, ch_ht_counts, ch_cell_qc, ch_dup_tsv, ch_truth_files, ch_asm_files)
        .collect()

    REPORT(ch_report_files, fai, run_info)
}

