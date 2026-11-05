nextflow.enable.dsl = 2

process FREEZE_OBSERVED_MODEL_EXPANSION {
    container 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d'
    publishDir "${params.run_dir}/frozen", mode: 'copy'
    output:
    path 'protocol.json'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/observed_dense_model_expansion.py' freeze \
      --source '${params.source}' --plan '${params.plan}' --output protocol.json
    """
}

process PREPARE_OBSERVED_MODEL_INPUTS {
    tag "${mask}/chr${chromosome}"
    cpus 1
    memory '4 GB'
    maxForks 4
    container 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d'
    publishDir path: { "${params.run_dir}/prepared/${mask}_chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(mask), val(chromosome), val(source_task), path(protocol)
    output:
    tuple val(mask), val(chromosome), path('prepared')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/observed_dense_model_expansion.py' prepare \
      --source '${params.source}' --source-task '${source_task}' --protocol '${protocol}' --output prepared
    """
}

process RUN_OBSERVED_MODEL_CANDIDATE {
    tag "${candidate}/${mask}/chr${chromosome}"
    cpus 1
    memory '8 GB'
    maxForks 22
    container 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d'
    publishDir path: { "${params.run_dir}/models/${candidate}/${mask}_chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(mask), val(chromosome), val(candidate), path(prepared)
    output:
    path 'result'
    script:
    def recoveryOption = params.recovery ? "--recovery '${params.recovery}'" : ""
    """
    python3 '${params.project}/pedigree_benchmark/scripts/observed_dense_model_expansion.py' evaluate-candidate \
      --source '${params.source}' --prepared '${prepared}' --protocol '${prepared}/protocol.json' \
      --candidate '${candidate}' ${recoveryOption} --output result
    """
}

process SCORE_OBSERVED_MODEL_EXPANSION {
    cpus 1
    memory '8 GB'
    container 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path(results, stageAs: 'results/result????')
    path protocol
    output:
    path 'summary'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/observed_dense_model_expansion.py' summarize \
      --source '${params.source}' --protocol '${protocol}' --results results/result* \
      --chromosomes ${params.chromosomes.join(' ')} --masks ${params.masks.join(' ')} --output summary
    """
}

workflow {
    def smoke = params.smoke.toString().toLowerCase()
    if (!(smoke in ['true', 'false'])) { error 'smoke must be explicitly true or false' }
    frozen = FREEZE_OBSERVED_MODEL_EXPANSION()
    params.chromosomes = smoke == 'true' ? [22] : (1..22).toList()
    params.masks = smoke == 'true' ? ['lonni'] : ['lonni', 'frederik']
    jobs = Channel.fromList(params.masks).combine(Channel.fromList(params.chromosomes)).combine(frozen).map {
        mask, chromosome, protocol ->
        def sourceTask = "${params.source}/tasks/${mask}_chr${chromosome}"
        tuple(mask, chromosome, sourceTask, protocol)
    }
    prepared = PREPARE_OBSERVED_MODEL_INPUTS(jobs)
    candidates = prepared.flatMap { mask, chromosome, folder ->
        def protocol = new groovy.json.JsonSlurper().parse(file("${folder}/protocol.json"))
        protocol.candidates.keySet().collect { candidate -> tuple(mask, chromosome, candidate, folder) }
    }
    results = RUN_OBSERVED_MODEL_CANDIDATE(candidates)
    SCORE_OBSERVED_MODEL_EXPANSION(results.collect(), frozen)
}
