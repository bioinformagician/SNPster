nextflow.enable.dsl = 2

process AUDIT_AND_PREPARE {
    cpus 2
    memory '6 GB'
    publishDir "${params.outdir}/${params.mode}/preflight", mode: 'copy'
    input:
    path source_manifest
    output:
    path 'prepared'
    script:
    """
    python3 /workspace/pedigree_benchmark/scripts/flare_nextflow_tasks.py prepare \
      --mode '${params.mode}' --source-manifest '${source_manifest}' --output prepared
    """
}

process CHARACTERIZE_INITIALIZATION {
    cpus 6
    memory '12 GB'
    publishDir "${params.outdir}/${params.mode}/characterization", mode: 'copy'
    input:
    path prepared
    output:
    path 'characterization'
    script:
    """
    python3 /workspace/pedigree_benchmark/scripts/flare_nextflow_tasks.py characterize \
      --prepared '${prepared}' --output characterization
    """
}

process NATIVE_FLARE_CALL {
    tag "${job_id}"
    cpus 6
    memory '12 GB'
    maxForks 2
    publishDir "${params.outdir}/${params.mode}/calls", mode: 'copy'
    input:
    val job_id
    path prepared
    path characterization
    output:
    path "${job_id}"
    script:
    """
    python3 /workspace/pedigree_benchmark/scripts/flare_nextflow_tasks.py infer \
      --prepared '${prepared}' --characterization '${characterization}' \
      --job '${job_id}' --output '${job_id}'
    """
}

process AGGREGATE_AND_VERIFY {
    cpus 2
    memory '6 GB'
    publishDir "${params.outdir}/${params.mode}/summary", mode: 'copy'
    input:
    path prepared
    path calls
    path characterization
    output:
    path 'summary'
    script:
    """
    python3 /workspace/pedigree_benchmark/scripts/flare_nextflow_tasks.py aggregate \
      --prepared '${prepared}' --characterization '${characterization}' --output summary
    """
}

workflow {
    if (!(params.mode in ['smoke', 'full']) || !params.outdir || !params.source_manifest)
        error 'Require --mode smoke|full, --outdir and --source_manifest'
    AUDIT_AND_PREPARE(file(params.source_manifest, checkIfExists: true))
    prepared = AUDIT_AND_PREPARE.out.first()
    CHARACTERIZE_INITIALIZATION(prepared)
    characterization = CHARACTERIZE_INITIALIZATION.out.first()
    jobs = prepared.flatMap { directory ->
        def manifest = new groovy.json.JsonSlurper().parse(directory.resolve('jobs.json').toFile())
        manifest.collect { row -> row.id.toString() }
    }
    NATIVE_FLARE_CALL(jobs, prepared, characterization)
    AGGREGATE_AND_VERIFY(prepared, NATIVE_FLARE_CALL.out.collect(), characterization)
}
