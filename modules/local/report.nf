process REPORT {
    label 'light'
    publishDir "${params.outdir}/${params.sample}/04_qc", mode: 'copy'

    input:
    path(files, stageAs: 'data/*')
    path fai
    val run_info

    output:
    path "phasing_report.html", emit: html
    path "figures/*",           emit: figures
    path "tables/*",            emit: tables

    script:
    """
    cat > data/run_info.json <<'JSON'
${run_info}
JSON
    phasing_report.py --data data --fai ${fai} --out . \\
        --min-qual ${params.min_qual} --min-dp ${params.min_dp} \\
        --min-vaf ${params.min_vaf} --max-vaf ${params.max_vaf}
    """
}
