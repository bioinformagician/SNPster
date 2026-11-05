nextflow.enable.dsl=2

params.population_stage = 'prepare'
params.root = null

process POPULATION_DESIGN {
    publishDir params.output, mode: 'copy'
    input:
    path specification
    output:
    path 'design'
    script:
    """
    python3 '${projectDir}/population_followup.py' '${params.population_stage}' \
        --specification '${specification}' --root '${params.root}' --output design
    """
}

process POPULATION_REPORT {
    publishDir params.output, mode: 'copy'
    input:
    path design
    path audit
    output:
    path 'report'
    script:
    """
    python3 '${projectDir}/population_followup.py' report --root '${params.root}' --output report
    """
}

workflow {
    if (params.population_stage in ['prepare', 'plan']) {
        POPULATION_DESIGN(file(params.experiment, checkIfExists: true))
    } else if (params.population_stage == 'report') {
        POPULATION_REPORT(
            file("${params.root}/planning/design", checkIfExists: true),
            file("${params.root}/screen/results/independent_verification.json", checkIfExists: true))
    } else {
        error 'population_stage must be prepare, plan or report'
    }
}
