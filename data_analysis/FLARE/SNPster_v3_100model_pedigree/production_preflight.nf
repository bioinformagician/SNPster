nextflow.enable.dsl = 2

process AUDIT_REFERENCE {
    tag "chr${chromosome}"
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/references", mode: 'copy'
    input:
    val chromosome
    output:
    path "reference_${chromosome}.json"
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/audit_production_ancestry.py" reference \
        --chromosome ${chromosome} --production-root "${params.production_root}" \
        --output reference_${chromosome}.json
    """
}

process SNAPSHOT_RUNTIME {
    container 'sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/runtime", mode: 'copy'
    output:
    path 'runtime.json'
    script:
    """
    PYTHONPATH=/app python3 "${params.project}/pedigree_benchmark/scripts/audit_production_ancestry.py" \
        runtime --output runtime.json
    """
}

process REPORT_ELIGIBILITY {
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path references
    path runtime
    output:
    path 'preflight.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/audit_production_ancestry.py" report \
        --project "${params.project}" --authority "${params.authority}" \
        --runtime "${runtime}" --references ${references.collect { '"' + it + '"' }.join(' ')} \
        --output preflight.json
    """
}

process REQUIRE_SCIENTIFIC_READINESS {
    cpus 1
    memory '1 GB'
    input:
    path report
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/audit_production_ancestry.py" gate \
        --output "${report}"
    """
}

workflow {
    references = AUDIT_REFERENCE(Channel.fromList((1..22).toList()))
    runtime = SNAPSHOT_RUNTIME()
    report = REPORT_ELIGIBILITY(references.collect(), runtime)
    REQUIRE_SCIENTIFIC_READINESS(report)
}
