nextflow.enable.dsl=2

params.final_stage = 'plan'
params.root = null

process FINAL_METADATA {
    publishDir params.output, mode: 'copy'
    input:
    path specification
    path audit
    output:
    path 'design'
    script:
    """
    python3 '${projectDir}/final_comparison.py' '${params.final_stage}' \
        --specification '${specification}' --root '${params.root}' --output design
    """
}

process FINAL_REPORT {
    publishDir params.output, mode: 'copy'
    input:
    path audit
    output:
    path 'report'
    script:
    """
    python3 '${projectDir}/final_comparison.py' report --root '${params.root}' --output report
    """
}

workflow {
    if (params.final_stage == 'plan') {
        FINAL_METADATA(file(params.experiment, checkIfExists: true), [])
    } else if (params.final_stage == 'test-plan') {
        FINAL_METADATA(file(params.experiment, checkIfExists: true),
            file("${params.root}/genome/results/independent_verification.json", checkIfExists: true))
    } else if (params.final_stage == 'report') {
        FINAL_REPORT(file("${params.root}/test/results/independent_verification.json", checkIfExists: true))
    } else {
        error 'final_stage must be plan, test-plan or report'
    }
}
