#!/usr/bin/env Rscript
# Per-cell allele count tables for ASCAT.sc + cell QC (run_splongget.sh step 5)
# usage: ascat_tables.R <cellsnp_dir> <out_dir> <qc.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: ascat_tables.R <cellsnp_dir> <out_dir> <qc.tsv>")
cellsnp_dir <- args[1]; ascat_dir <- args[2]; qc_file <- args[3]

suppressMessages({
  library(Matrix)
  library(data.table)
})
dir.create(ascat_dir, showWarnings = FALSE, recursive = TRUE)

snp_vcf  <- fread(cmd = sprintf("gunzip -c %s | grep -v '^##'", file.path(cellsnp_dir, "cellSNP.base.vcf.gz")), data.table = FALSE)
barcodes <- fread(file.path(cellsnp_dir, "cellSNP.samples.tsv"), header = FALSE)$V1
mtx <- function(x) {
  f <- file.path(cellsnp_dir, x)
  if (!file.exists(f)) f <- paste0(f, ".gz")
  as(readMM(f), "CsparseMatrix")   # column-compressed: one cell at a time is fast
}
ad <- mtx("cellSNP.tag.AD.mtx")
dp <- mtx("cellSNP.tag.DP.mtx")

qc_rows <- vector("list", length(barcodes))
for (i in seq_along(barcodes)) {
  cell_dp <- dp[, i]
  cell_ad <- ad[, i]
  covered <- sum(cell_dp > 0)
  qc_rows[[i]] <- data.frame(
    cell_id = barcodes[i],
    covered_snps = covered,
    total_depth = sum(cell_dp),
    median_cov_per_covered_snp = if (covered > 0) median(cell_dp[cell_dp > 0]) else 0
  )
  df <- data.frame(`#CHR` = snp_vcf$`#CHROM`, POS = snp_vcf$POS,
                   Count_A = 0, Count_C = 0, Count_G = 0, Count_T = 0,
                   Good_depth = cell_dp, check.names = FALSE)
  for (base in c("A", "C", "G", "T")) {
    col <- paste0("Count_", base)
    is_ref <- snp_vcf$REF == base
    is_alt <- snp_vcf$ALT == base
    df[[col]][is_ref] <- (cell_dp - cell_ad)[is_ref]
    df[[col]][is_alt] <- cell_ad[is_alt]
  }
  fwrite(df, file.path(ascat_dir, paste0(barcodes[i], "_alleleCounts.txt")), sep = "\t")
}
fwrite(rbindlist(qc_rows), qc_file, sep = "\t")
