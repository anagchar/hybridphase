# Changelog

## v0.1.0

First Nextflow version of `legacy/run_splongget.sh`.

- Clair3 variant calling in front (or `--vcf` to reuse calls)
- removal of duplicate-flagged reads (lrsinglecell BAMs flag duplicates instead of removing them)
- read phasing with WhatsHap and/or LongPhase; population phasing with SHAPEIT4 and/or Beagle; all combinations in one run
- block-orientation hybrid merge (`bin/hybrid_merge.py`) next to the original WhatsHap re-phasing step
- truth-free QC: junction masking test, read-block vs population agreement, method agreement, haplotype read fractions
- haplotagging with WhatsHap or LongPhase
- optional cellsnp-lite + ASCAT.sc tables (`--barcodes`)
- self-contained HTML report and tube map (`assets/make_tube_map.py`)
