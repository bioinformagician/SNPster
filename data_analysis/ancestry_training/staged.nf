nextflow.enable.dsl=2

params.source_run = null
params.plan_dir = null
params.screen_results = null
params.metadata_stage = 'plan'

process DEVELOPMENT_PLAN {
    publishDir params.output, mode: 'copy'
    input:
    path config
    path samples
    path controller
    path tasks
    path previous
    path traces, stageAs: 'traces/item????'
    output:
    path 'plan'
    script:
    """
    python3 '${projectDir}/staged_comparison.py' plan --config '${config}' \
        --source-run '${params.source_run}' --output plan --durable '${params.output}/plan'
    """
}

process ADVANCE_FLARE {
    publishDir params.output, mode: 'copy'
    input:
    path plan
    path summary
    path audit
    output:
    path 'selection'
    script:
    """
    python3 '${projectDir}/staged_comparison.py' advance --plan '${plan}' \
        --screen-results '${params.screen_results}' --output selection
    """
}

workflow PLAN_ONLY {
    def config = new groovy.json.JsonSlurper().parse(file(params.experiment).toFile())
    def transition = "${params.source_run}/reports/staged-scope-transition"
    DEVELOPMENT_PLAN(
        file(params.experiment, checkIfExists: true),
        file(config.samples, checkIfExists: true),
        file("$transition/controller-after.json", checkIfExists: true),
        file("$transition/owned-tasks-after.json", checkIfExists: true),
        file("${params.source_run}/runtime/reuse-manifest.json", checkIfExists: true),
        Channel.value(file("${params.source_run}/reports/*.smoke.trace.tsv", checkIfExists: true)))
}

workflow ADVANCE_ONLY {
    ADVANCE_FLARE(
        file(params.plan_dir, checkIfExists: true),
        file("${params.screen_results}/summary/summary/summary.json", checkIfExists: true),
        file("${params.screen_results}/independent_verification.json", checkIfExists: true))
}

workflow {
    if (params.metadata_stage == 'plan') {
        PLAN_ONLY()
    } else if (params.metadata_stage == 'advance') {
        ADVANCE_ONLY()
    } else {
        error 'metadata_stage must be plan or advance'
    }
}
