#!/usr/bin/env bash
set -euo pipefail

function usage() {
    cat <<EOF
Usage: apptainer run splongget_phasing.sif [OPTIONS]

Required Arguments:
  -b, --bam FILE          Input pseudobulk BAM (with CB:Z tags)
  -r, --ref FILE          Reference genome FASTA
  -v, --vcf FILE          Clair3 output VCF (merge_output.vcf.gz)
  -c, --barcodes FILE     Cell barcodes list (only the first column is used)
  -o, --out-dir DIR       Output directory to save all results & QC

Optional Arguments:
  -t, --threads INT       Number of CPU threads [default: 22]
  -C, --chroms LIST       Autosomes to run, without "chr", e.g. 22 or 21,22 [default: 1-22]
  -h, --help              Show this help message
EOF
    exit 1
}

THREADS=22
CHROMS=$(seq 1 22)
PANEL_DIR="/opt/ref_data/1000G_phase3_hg38"
MAP_DIR="/opt/ref_data/genetic_maps"

BAM=""
REF=""
VCF_INPUT=""
BARCODES=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam) BAM="$2"; shift 2 ;;
        -r|--ref) REF="$2"; shift 2 ;;
        -v|--vcf) VCF_INPUT="$2"; shift 2 ;;
        -c|--barcodes) BARCODES="$2"; shift 2 ;;
        -o|--out-dir) OUT_DIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -C|--chroms) CHROMS="${2//,/ }"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
done

for var in BAM REF VCF_INPUT BARCODES OUT_DIR; do
    if [ -z "${!var:-}" ]; then
        echo "Error: Missing required argument!"
        usage
    fi
done

mkdir -p "${OUT_DIR}/qc" "${OUT_DIR}/phased_bcf" "${OUT_DIR}/cellsnp" "${OUT_DIR}/ascat_allele_counts"

LOG="${OUT_DIR}/pipeline.log"
exec > >(tee -a "${LOG}") 2>&1

echo "=========================================================="
echo " Starting Self-Contained SPLONGGET Phasing Pipeline"
echo " Date: $(date)"
echo " BAM: ${BAM}"
echo " Threads: ${THREADS}"
echo " Chromosomes: ${CHROMS//$'\n'/ }"
echo "=========================================================="

echo "[STEP 1] Filtering Germline Het SNPs..."
bcftools view -v snps -m2 -M2 -f PASS "${VCF_INPUT}" \
  | bcftools filter -i 'GT="het" && QUAL>=20 && FORMAT/DP>=15 && (FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])) >= 0.25 && (FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])) <= 0.75' \
  -Oz -o "${OUT_DIR}/germline_het_snps.vcf.gz"

tabix -f -p vcf "${OUT_DIR}/germline_het_snps.vcf.gz"
bcftools stats "${OUT_DIR}/germline_het_snps.vcf.gz" > "${OUT_DIR}/qc/01_germline_het_snps.stats"

echo "[STEP 2] Phasing: reads first, then SHAPEIT4 connecting the read blocks..."
for CHR in ${CHROMS}; do
    CHR_TAG="chr${CHR}"; P="${OUT_DIR}/phased_bcf"
    echo "--- Processing ${CHR_TAG} ---"

    bcftools view -r "${CHR_TAG}" "${OUT_DIR}/germline_het_snps.vcf.gz" -Ob -o "${P}/input_${CHR_TAG}.bcf"
    bcftools index -f "${P}/input_${CHR_TAG}.bcf"

    # 2a) Reads only -> locally trustworthy blocks with PS tags
    whatshap phase --reference "${REF}" --ignore-read-groups \
        -o "${P}/reads_${CHR_TAG}.vcf.gz" "${P}/input_${CHR_TAG}.bcf" "${BAM}"
    tabix -f -p vcf "${P}/reads_${CHR_TAG}.vcf.gz"

    # 2b) SHAPEIT4 keeps the blocks intact, orients them with the panel
    bcftools +fill-tags "${P}/reads_${CHR_TAG}.vcf.gz" -Ob -o "${P}/reads_${CHR_TAG}.tagged.bcf" -- -t AN,AC
    bcftools index -f "${P}/reads_${CHR_TAG}.tagged.bcf"
    shapeit4 --input "${P}/reads_${CHR_TAG}.tagged.bcf" \
        --reference "${PANEL_DIR}/1kGP_phase3_${CHR_TAG}.bcf" \
        --map "${MAP_DIR}/${CHR_TAG}.b38.gmap.gz" --region "${CHR_TAG}" \
        --sequencing --use-PS 0.001 --thread "${THREADS}" \
        --output "${P}/shapeit4_${CHR_TAG}.phased.bcf"          # do NOT also pass --scaffold
    bcftools view "${P}/shapeit4_${CHR_TAG}.phased.bcf" -Oz -o "${P}/shapeit4_${CHR_TAG}.phased.vcf.gz"
    tabix -f -p vcf "${P}/shapeit4_${CHR_TAG}.phased.vcf.gz"

    # 2c) Needed: SHAPEIT4 only keeps SNPs that are also in the panel, so re-attach the others
    #     (the SHAPEIT4 pseudo-read now agrees with the reads)
    whatshap phase --reference "${REF}" --ignore-read-groups \
        -o "${P}/hybrid_${CHR_TAG}.phased.vcf.gz" \
        "${P}/input_${CHR_TAG}.bcf" "${BAM}" "${P}/shapeit4_${CHR_TAG}.phased.vcf.gz"
    tabix -f -p vcf "${P}/hybrid_${CHR_TAG}.phased.vcf.gz"

    # QC: within read blocks the final phasing should now (almost) never disagree with the reads
    whatshap compare --names reads_only,final \
        --tsv-pairwise "${OUT_DIR}/qc/02_reads_vs_final_${CHR_TAG}.tsv" \
        "${P}/reads_${CHR_TAG}.vcf.gz" "${P}/hybrid_${CHR_TAG}.phased.vcf.gz" \
        > "${OUT_DIR}/qc/02_reads_vs_final_${CHR_TAG}.txt"

    rm -f "${P}/input_${CHR_TAG}.bcf"* "${P}/reads_${CHR_TAG}.tagged.bcf"*
done

echo "[STEP 3] Concatenating Phased Genome..."
for CHR in ${CHROMS}; do
    echo "${OUT_DIR}/phased_bcf/hybrid_chr${CHR}.phased.vcf.gz"
done > "${OUT_DIR}/phased_bcf/autosomes.list"

bcftools concat -f "${OUT_DIR}/phased_bcf/autosomes.list" -Oz -o "${OUT_DIR}/fully_phased_genome.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/fully_phased_genome.vcf.gz"

# Copy for ASCAT.sc (path_to_phases): its readPhases() reads the VCF with data.table::fread,
# which cannot get past the long ## header. Keep only the #CHROM line, and GT as the only FORMAT field.
bcftools annotate -x INFO,^FORMAT/GT "${OUT_DIR}/fully_phased_genome.vcf.gz" \
    | grep -v '^##' | gzip > "${OUT_DIR}/phases_for_ascat.vcf.gz"

whatshap stats --tsv="${OUT_DIR}/qc/03_genome_phasing_stats.tsv" "${OUT_DIR}/fully_phased_genome.vcf.gz" \
    > "${OUT_DIR}/qc/03_genome_phasing_summary.txt"

echo "[STEP 4] Counting Cell-Level Alleles (cellsnp-lite)..."
# flexiplex lists "barcode<TAB>read count"; cellsnp-lite takes the whole line as the
# barcode, so no CB tag matches and every cell gets 0 counts. Keep the barcode only.
cut -f1 "${BARCODES}" > "${OUT_DIR}/barcodes.tsv"

cellsnp-lite \
  -s "${BAM}" \
  -b "${OUT_DIR}/barcodes.tsv" \
  -O "${OUT_DIR}/cellsnp" \
  -R "${OUT_DIR}/fully_phased_genome.vcf.gz" \
  --cellTAG CB \
  --UMItag None \
  --minMAF 0 \
  --minCOUNT 0 \
  -p "${THREADS}" \
  --gzip

echo "[STEP 5] Generating ASCAT.sc Tables & Computing Cell QC..."
Rscript - <<'EOF' "${OUT_DIR}"
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

suppressMessages({
  library(Matrix)
  library(data.table)
})

vcf_path <- file.path(out_dir, "cellsnp", "cellSNP.base.vcf.gz")
barcodes_path <- file.path(out_dir, "cellsnp", "cellSNP.samples.tsv")
ad_path <- file.path(out_dir, "cellsnp", "cellSNP.tag.AD.mtx")
dp_path <- file.path(out_dir, "cellsnp", "cellSNP.tag.DP.mtx")

cmd <- sprintf("gunzip -c %s | grep -v '^##'", vcf_path)
snp_vcf <- fread(cmd = cmd, data.table = FALSE)
barcodes <- fread(barcodes_path, header = FALSE)$V1

# Column-compressed: taking one cell (column) at a time is then fast on the whole genome
ad <- as(readMM(ad_path), "CsparseMatrix")
dp <- as(readMM(dp_path), "CsparseMatrix")

ascat_dir <- file.path(out_dir, "ascat_allele_counts")
qc_rows <- list()

for (i in seq_along(barcodes)) {
  cell_id <- barcodes[i]
  cell_dp <- dp[, i]
  cell_ad <- ad[, i]
  
  covered_snps <- sum(cell_dp > 0)
  total_counts <- sum(cell_dp)
  median_cov   <- if (covered_snps > 0) median(cell_dp[cell_dp > 0]) else 0
  
  qc_rows[[i]] <- data.frame(
    cell_id = cell_id,
    covered_snps = covered_snps,
    total_depth = total_counts,
    median_cov_per_covered_snp = median_cov
  )
  
  ref_counts <- cell_dp - cell_ad
  alt_counts <- cell_ad
  
  df_cell <- data.frame(
    `#CHR`     = snp_vcf$`#CHROM`,
    POS        = snp_vcf$POS,
    Count_A    = 0,
    Count_C    = 0,
    Count_G    = 0,
    Count_T    = 0,
    Good_depth = cell_dp,
    check.names = FALSE
  )
  
  for (base in c("A", "C", "G", "T")) {
    is_ref <- (snp_vcf$REF == base)
    is_alt <- (snp_vcf$ALT == base)
    col_name <- paste0("Count_", base)
    df_cell[[col_name]][is_ref] <- ref_counts[is_ref]
    df_cell[[col_name]][is_alt] <- alt_counts[is_alt]
  }
  
  fwrite(df_cell, file = file.path(ascat_dir, paste0(cell_id, "_alleleCounts.txt")), sep = "\t")
}

fwrite(do.call(rbind, qc_rows), file.path(out_dir, "qc", "04_single_cell_allele_qc.tsv"), sep = "\t")
EOF

echo "=========================================================="
echo " Pipeline Complete!"
echo "=========================================================="