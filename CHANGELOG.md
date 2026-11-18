# Changelog

## v0.1.0

First release.

- Clair3 variant calling in front (or `--vcf` to reuse calls)
- removal of duplicate-flagged (0x400) reads before calling and phasing
- read phasing with WhatsHap and/or LongPhase; population phasing with SHAPEIT4 and/or Beagle; all combinations in one run
- block-orientation hybrid merge (`bin/hybrid_merge.py`) next to the original WhatsHap re-phasing step
- truth-free QC: junction masking test, read-block vs population agreement, method agreement, haplotype read fractions
- haplotagging with WhatsHap or LongPhase
- optional haplotagging by mapping reads to a haplotype-resolved assembly (`--hap1/--hap2`), with read-length,
  per-cell and VCF-agreement reports
- optional benchmark against a haplotype-resolved diploid assembly (dipcall + WhatsHap compare)
- optional cellsnp-lite + ASCAT.sc tables (`--barcodes`)
- self-contained HTML report and tube map (`assets/make_tube_map.py`)
