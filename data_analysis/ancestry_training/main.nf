nextflow.enable.dsl=2

process FREEZE {
    publishDir "${params.output}/frozen", mode: 'copy'
    input:
    path config
    path reuse
    output:
    path 'frozen'
    script:
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' freeze --config '${config}' $reuseArg --output frozen
    """
}

process PREPARE {
    tag "$route/chr$chrom"
    publishDir path: { "${params.output}/panels/$route/chr$chrom" }, mode: 'copy'
    input:
    tuple val(route), val(chrom), path(frozen)
    path reuse
    output:
    tuple val(route), val(chrom), path('panel')
    script:
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' prepare --frozen '${frozen}' --route '$route' --chromosome $chrom $reuseArg --output panel
    """
}

process PHASE {
    tag "$route/chr$chrom"
    publishDir path: { "${params.output}/phased/$route/chr$chrom" }, mode: 'copy'
    input:
    tuple val(route), val(chrom), path(panel)
    path reuse
    output:
    tuple val(route), val(chrom), path('phased')
    script:
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' phase --panel '${panel}' $reuseArg --output phased
    """
}

process BASE {
    tag "$route/chr$chrom/$baseid"
    publishDir path: { "${params.output}/bases/$route/chr$chrom/$baseid" }, mode: 'copy'
    container { family == 'gnomix' ? params.gnomix_image : params.image }
    input:
    tuple val(route), val(chrom), val(baseid), val(family), path(panel)
    path reuse
    output:
    tuple val(route), val(chrom), val(baseid), path('base')
    script:
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' base --panel '${panel}' --base-id '$baseid' $reuseArg --output base
    """
}

process FIT {
    tag "$route/chr$chrom/$model"
    publishDir path: { "${params.output}/models/$route/chr$chrom/$model" }, mode: 'copy'
    container { family == 'gnomix' ? params.gnomix_image : params.image }
    input:
    tuple val(route), val(chrom), val(model), val(family), path(base), path(panel)
    path reuse
    output:
    tuple val(route), val(chrom), val(model), path('fitted')
    script:
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' fit --panel '${panel}' --base '${base}' --model '$model' $reuseArg --output fitted
    """
}

process INFER {
    tag "$route/chr$chrom/$model/queries$first-$end"
    publishDir path: { "${params.output}/predictions/$route/chr$chrom/$model" },
        mode: 'copy', saveAs: { name -> family == 'flare' ? null : name }
    container { family == 'gnomix' ? params.gnomix_image : params.image }
    input:
    tuple val(route), val(chrom), val(model), val(family), val(first), val(end), val(chunkCount), path(fitted), path(base), path(panel), path(phased)
    path reuse
    output:
    tuple val(route), val(chrom), val(model), val(family), val(chunkCount), path('prediction')
    script:
    def range = family == 'flare' ? "--query-start $first --query-end $end" : ''
    def reuseArg = reuse ? "--reuse-manifest '$reuse'" : ''
    """
    python3 '${projectDir}/workflow.py' infer --panel '${panel}' --base '${base}' --fitted '${fitted}' --phased '${phased}' $range $reuseArg --output prediction
    """
}

process COLLECT_INFERENCE {
    tag "$route/chr$chrom/$model"
    publishDir path: { "${params.output}/predictions/$route/chr$chrom/$model" }, mode: 'copy'
    input:
    tuple val(route), val(chrom), val(model), path(chunks, stageAs: 'chunks/item????'), path(panel)
    output:
    tuple val(route), val(model), path('prediction')
    script:
    """
    python3 '${projectDir}/workflow.py' collect-inference --panel '${panel}' --results chunks/item* --output prediction
    """
}

process AGGREGATE {
    tag "$route/$model"
    publishDir path: { "${params.output}/aggregated/$route/$model" }, mode: 'copy'
    input:
    tuple val(route), val(model), path(results, stageAs: 'results/item????')
    output:
    tuple val(route), val(model), path('aggregate')
    script:
    """
    python3 '${projectDir}/workflow.py' aggregate --results results/item* --output aggregate
    """
}

process INTERPRET {
    tag "$route/$model/$interpreter"
    publishDir path: { "${params.output}/interpreters/$route/$model/$interpreter" }, mode: 'copy'
    input:
    tuple val(route), val(model), val(interpreter), path(aggregate)
    output:
    path 'interpreted'
    script:
    """
    python3 '${projectDir}/workflow.py' interpret --aggregate '${aggregate}' --interpreter '$interpreter' --output interpreted
    """
}

process SELECT {
    publishDir "${params.output}/summary", mode: 'copy'
    input:
    path results, stageAs: 'results/item????'
    output:
    path 'summary'
    script:
    """
    python3 '${projectDir}/workflow.py' select --results results/item* --output summary
    """
}

process VERIFY {
    publishDir params.output, mode: 'copy'
    input:
    path summary
    output:
    path 'verification.json'
    script:
    """
    python3 '${projectDir}/workflow.py' verify --aggregate '${summary}' --output verification.json
    """
}

workflow {
    if (!(params.inference_chunk_size.toString() ==~ /[1-9][0-9]*/)) {
        error 'inference_chunk_size must be a positive integer'
    }
    def chunkSize = params.inference_chunk_size.toString().toInteger()
    reuseInputs = Channel.value(params.reuse_manifest ? file(params.reuse_manifest, checkIfExists: true) : [])
    frozen = FREEZE(file(params.experiment, checkIfExists: true), reuseInputs)
    jobs = frozen.flatMap { folder ->
        def config = new groovy.json.JsonSlurper().parse(file("$folder/config.json").toFile())
        config.routes.collectMany { route -> config.chromosomes.collect { chrom -> tuple(route.id, chrom, folder) } }
    }
    panels = PREPARE(jobs, reuseInputs)
    phases = PHASE(panels, reuseInputs)
    // The plan is derived exclusively from the supplied, frozen config, not outcomes.
    def config = new groovy.json.JsonSlurper().parse(file(params.experiment).toFile())
    baseJobs = panels.flatMap { route, chrom, panel ->
        def plan = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
        // Use the same canonical base-policy computation as the task implementation.
        def bases = plan.base_plan.bases
        bases.collect { id, spec -> tuple(route, chrom, id, spec.family, panel) }
    }
    bases = BASE(baseJobs, reuseInputs)
    fittingJobs = bases.combine(panels, by: [0, 1]).flatMap { route, chrom, baseid, base, panel ->
        def plan = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile()).base_plan
        plan.models.findAll { it.base_id == baseid }.collect {
            model -> tuple(route, chrom, model.id, model.family, base, panel)
        }
    }
    fits = FIT(fittingJobs, reuseInputs)
    inferenceJobs = fits.combine(panels, by: [0, 1]).combine(phases, by: [0, 1]).map {
        route, chrom, modelid, fitted, panel, phased ->
        def plan = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile()).base_plan
        def spec = plan.models.find { it.id == modelid }
        tuple(route, chrom, modelid, spec.family, fitted, spec.base_id, panel, phased)
    }
    // Join the exact base identity, not just chromosome or candidate family.
    inferenceWithBase = inferenceJobs.map { route, chrom, model, family, fitted, baseid, panel, phased ->
        tuple(route, chrom, baseid, model, family, fitted, panel, phased)
    }.combine(bases, by: [0, 1, 2]).map {
        route, chrom, baseid, model, family, fitted, panel, phased, base ->
        tuple(route, chrom, model, family, fitted, base, panel, phased)
    }
    chunkJobs = inferenceWithBase.flatMap { route, chrom, model, family, fitted, base, panel, phased ->
        def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
        def size = family == 'flare' ? chunkSize : meta.queries.size()
        def starts = (0..<meta.queries.size()).step(size)
        starts.collect { first ->
            tuple(route, chrom, model, family, first, Math.min(first + size, meta.queries.size()),
                  starts.size(), fitted, base, panel, phased)
        }
    }
    inferred = INFER(chunkJobs, reuseInputs).branch {
        flare: it[3] == 'flare'
        other: true
    }
    completedChunks = inferred.flare.map { route, chrom, model, family, count, folder ->
        tuple(groupKey(tuple(route, chrom, model), count), folder)
    }.groupTuple().map { key, folders ->
        def identity = key.getGroupTarget()
        tuple(identity[0], identity[1], identity[2], folders)
    }.combine(panels, by: [0, 1])
    assembled = COLLECT_INFERENCE(completedChunks)
    predictions = inferred.other.map { route, chrom, model, family, count, folder -> tuple(route, model, folder) }.mix(assembled)
    aggregates = AGGREGATE(predictions.groupTuple(by: [0, 1]))
    interpretationJobs = aggregates.flatMap { route, model, aggregate ->
        config.interpretations.collect { item -> tuple(route, model, item.id, aggregate) }
    }
    interpreted = INTERPRET(interpretationJobs)
    summary = SELECT(interpreted.collect())
    VERIFY(summary)
}
