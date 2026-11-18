#!/usr/bin/env bash
# Runs the helper scripts (bin/) on synthetic data: no containers, no cluster.
# Needs python3 with pysam, numpy, pandas, matplotlib.
set -euo pipefail
cd "$(dirname "$0")/.."
B=$PWD/bin
T=$PWD/tests/out
rm -rf "$T"
python3 tests/make_synthetic.py "$T"
cd "$T"
idx() { python3 -c "import pysam,sys; pysam.tabix_index(sys.argv[1], preset='vcf', force=True)" "$1"; }

python3 $B/hybrid_merge.py --reads reads.vcf.gz --pop pop.vcf.gz --out final.vcf.gz --blocks data/orientation_whatshap_shapeit4.tsv
idx final.vcf.gz
python3 $B/junction_test.py split reads.vcf.gz split.vcf.gz --min-snps 10
python3 $B/hybrid_merge.py --reads split.vcf.gz --pop pop.vcf.gz --out rejoined.vcf.gz --blocks /dev/null
python3 $B/junction_test.py score split.vcf.gz.splits.tsv reads.vcf.gz rejoined.vcf.gz data/junctions_whatshap_shapeit4.tsv
python3 $B/phase_stats.py reads.vcf.gz data/whatshap_reads.phase_stats.tsv data/whatshap_reads.blocklist.tsv
python3 $B/phase_stats.py final.vcf.gz data/whatshap_shapeit4.phase_stats.tsv data/whatshap_shapeit4.blocklist.tsv

# sanity checks: inside read blocks the merge never changes the relative phase
python3 - <<'PY'
import pysam
r = {(x.chrom, x.pos): (x.samples[0]["GT"], x.samples[0].get("PS")) for x in pysam.VariantFile("reads.vcf.gz") if x.samples[0].phased}
f = {(x.chrom, x.pos): x.samples[0]["GT"] for x in pysam.VariantFile("final.vcf.gz")}
from collections import defaultdict
rel = defaultdict(set)
for k, (gt, ps) in r.items():
    rel[(k[0], ps)].add(gt[0] == f[k][0])
bad = sum(len(v) > 1 for v in rel.values())
assert bad == 0, f"{bad} read blocks were broken by the merge"
print(f"OK: {len(rel)} read blocks kept intact")
PY

python3 - <<'PY'
import json
json.dump(dict(sample="SYNTHETIC", variant_calls="synthetic", chroms=["chr21", "chr22"], combos=["whatshap_shapeit4"],
               hybrid_merge="orient", haplotag_combo="whatshap_shapeit4", haplotag_tool="whatshap",
               run_name="test", started="-", haplotag_bin_mb=1, asm_min_diff=6), open("data/run_info.json", "w"))
PY
# truth benchmark (mocked whatshap compare against the known haplotypes)
python3 ../mock_whatshap_compare.py truth.vcf.gz final.vcf.gz whatshap_shapeit4 data/whatshap_shapeit4.truth.tsv data/whatshap_shapeit4.truth_switches.bed
python3 ../mock_whatshap_compare.py truth.vcf.gz reads.vcf.gz whatshap_reads data/whatshap_reads.truth.tsv data/whatshap_reads.truth_switches.bed
# assembly-based haplotagging on simulated alignments
python3 ../make_synthetic_asm.py asm
python3 $B/assembly_haplotag.py assign asm/asm.byname.bam asm/assign.tsv.gz data/asm --min-diff 6
python3 $B/assembly_haplotag.py tag asm/assign.tsv.gz asm/vcf_tagged.bam asm/tagged.bam data/asm
python3 - <<'PY'
import pysam
n = ok = 0
for r in pysam.AlignmentFile("asm/tagged.bam"):
    if r.has_tag("HP"):
        n += 1
        ok += r.get_tag("HP") == int(r.get_tag("ha")[1])
assert n and ok == n, "HP tag does not match the assigned contig"
print(f"OK: {n} reads tagged from the assembly")
PY
printf "reads\tn\nall_primary\t1000000\nduplicate_flagged\t42000\n" > data/SYNTHETIC.duplicates.tsv
python3 $B/phasing_report.py --data data --fai ref.fa.fai --out report
echo "Report: $T/report/phasing_report.html"
