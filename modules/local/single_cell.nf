/*
 * Optional (--barcodes): cell-level allele counts on the phased SNPs + ASCAT.sc tables
 */

process CELLSNP {
    tag "${combo}"
    label 'pop'
    publishDir "${params.outdir}/${params.sample}/06_single_cell", mode: 'copy'

    input:
    tuple val(combo), path(vcf), path(tbi)
    tuple path(bam), path(bai)
    path barcodes

    output:
    tuple val(combo), path("cellsnp"), emit: dir
    path "barcodes.tsv"

    script:
    """
    # flexiplex lists "barcode<TAB>read count"; cellsnp-lite would take the whole line as barcode
    cut -f1 ${barcodes} > barcodes.tsv

    cellsnp-lite -s ${bam} -b barcodes.tsv -O cellsnp -R ${vcf} \\
        --cellTAG CB --UMItag None --minMAF 0 --minCOUNT 0 -p ${task.cpus} --gzip
    """
}

process ASCAT_TABLES {
    tag "${combo}"
    label 'pop'
    publishDir "${params.outdir}/${params.sample}/06_single_cell", mode: 'copy'

    input:
    tuple val(combo), path(cellsnp_dir)

    output:
    path "ascat_allele_counts", emit: dir
    path "single_cell_allele_qc.tsv", emit: qc

    script:
    """
    ascat_tables.R ${cellsnp_dir} ascat_allele_counts single_cell_allele_qc.tsv
    """
}
