nextflow.enable.dsl=2

process AUDIT_FLARE_INPUTS {
    publishDir "${params.output}/inputs", mode: 'copy'
    input:
    path protocol
    output:
    path 'prepared'
    script:
    """
    python3 '${projectDir}/diagnose_flare.py' prepare --parent '${params.parent}' \
        --protocol '$protocol' --output prepared
    """
}

process TEST_FLARE_MECHANISM {
    tag "$identity"
    publishDir path: { "${params.output}/calls/$identity" }, mode: 'copy'
    input:
    tuple val(identity), path(prepared)
    output:
    path 'result'
    script:
    """
    python3 '${projectDir}/diagnose_flare.py' evaluate --prepared '$prepared' \
        --identity '$identity' --output result
    """
}

process SMOKE_FLARE_MECHANISM {
    publishDir "${params.output}/smoke", mode: 'copy'
    input:
    path prepared
    output:
    path 'smoke'
    script:
    """
    python3 '${projectDir}/diagnose_flare.py' evaluate --prepared '$prepared' \
        --identity 'lonni_HG02954_cap16_baseline_fixed' --output smoke
    """
}

process VERIFY_FLARE_DIAGNOSIS {
    publishDir "${params.output}/summary", mode: 'copy'
    input:
    path prepared
    path results, stageAs: 'results/item????'
    output:
    path 'summary'
    script:
    """
    python3 '${projectDir}/diagnose_flare.py' summarize --prepared '$prepared' \
        --results results/item* --output summary
    """
}

workflow {
    prepared = AUDIT_FLARE_INPUTS(file(params.experiment, checkIfExists: true))
    smoke = SMOKE_FLARE_MECHANISM(prepared)
    jobs = prepared.combine(smoke).flatMap { directory, verified ->
        def jobs = new groovy.json.JsonSlurper().parse(file("$directory/jobs.json").toFile())
        jobs.collect { job -> tuple(job.id, directory) }
    }
    calls = TEST_FLARE_MECHANISM(jobs)
    VERIFY_FLARE_DIAGNOSIS(prepared, calls.collect())
}
