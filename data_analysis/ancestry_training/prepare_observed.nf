nextflow.enable.dsl=2

process PLAN_OBSERVED_COHORT {
    publishDir "${params.output}/metadata", mode: 'copy'
    input:
    path template
    val source
    output:
    path 'plan'
    script:
    """
    python3 '${projectDir}/prepare_observed.py' plan --source '$source' --template '$template' --output plan
    """
}

process BUILD_OBSERVED_CHROMOSOME {
    tag "chr$chrom"
    publishDir path: { "${params.output}/chromosomes/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(plan)
    output:
    path 'bundle'
    script:
    """
    python3 '${projectDir}/prepare_observed.py' chromosome --plan '$plan' --chromosome $chrom --output bundle
    """
}

process VERIFY_OBSERVED_COHORT {
    publishDir "${params.output}", mode: 'copy'
    input:
    path plan
    path bundles, stageAs: 'bundles/chr??'
    output:
    path 'ready'
    script:
    """
    python3 '${projectDir}/prepare_observed.py' finalize --plan '$plan' --bundles bundles/chr* \
        --destination '${params.output}' --output ready
    """
}

workflow {
    plan = PLAN_OBSERVED_COHORT(file(params.experiment, checkIfExists: true), params.source)
    def config = new groovy.json.JsonSlurper().parse(file(params.experiment).toFile())
    jobs = plan.flatMap { directory ->
        ([22] + config.chromosomes.findAll { it != 22 }).collect { chrom -> tuple(chrom, directory) }
    }
    bundles = BUILD_OBSERVED_CHROMOSOME(jobs)
    VERIFY_OBSERVED_COHORT(plan, bundles.collect())
}
