nextflow.enable.dsl=2

include { FREEZE; PREPARE; BASE; FIT; INFER; COLLECT_INFERENCE; AGGREGATE; INTERPRET; SELECT; VERIFY } from './main.nf'

params.guided_stage = 'run'
params.phase_reuse = null
params.design = null
params.phase_chunk_size = 16
params.java_c1_image = 'sha256:bc523b63202a52dd2aad32ccadb7813c179c37c93fe89168bc3ae2c1c0172547'

process GUIDED_PLAN {
    publishDir params.output, mode: 'copy'
    input:
    path specification
    output:
    path 'design'
    script:
    """
    python3 '${projectDir}/guided_followup.py' plan --config '${specification}' \
        --output design --durable '${params.output}/design'
    """
}

process GUIDED_PHASE {
    tag "$route/chr$chrom/queries$first-$end"
    container params.java_c1_image
    input:
    tuple val(route), val(chrom), val(first), val(end), val(count), path(panel)
    path reuse
    output:
    tuple val(route), val(chrom), val(count), path('phase_chunk')
    script:
    """
    python3 '${projectDir}/guided_followup.py' phase-chunk --panel '${panel}' \
        --reuse '${reuse}' --start $first --end $end --output phase_chunk
    """
}

process COLLECT_GUIDED_PHASE {
    tag "$route/chr$chrom"
    publishDir path: { "${params.output}/phased/$route/chr$chrom" }, mode: 'copy'
    input:
    tuple val(route), val(chrom), path(chunks, stageAs: 'chunks/item????'), path(panel)
    output:
    tuple val(route), val(chrom), path('phased')
    script:
    """
    python3 '${projectDir}/guided_followup.py' collect-phase --panel '${panel}' \
        --folders chunks/item* --output phased
    """
}

process GUIDED_REPORT {
    publishDir params.output, mode: 'copy'
    input:
    path design
    path verification
    output:
    path 'discrimination'
    script:
    """
    python3 '${projectDir}/guided_followup.py' report --results '${params.output}/../screen/results' \
        --design '${design}' --output discrimination
    """
}

workflow {
    if (params.guided_stage == 'plan') {
        GUIDED_PLAN(file(params.experiment, checkIfExists: true))
    } else if (params.guided_stage == 'report') {
        GUIDED_REPORT(file(params.design, checkIfExists: true),
            file("${params.output}/../screen/results/independent_verification.json", checkIfExists: true))
    } else if (params.guided_stage == 'run') {
        def chunkSize = params.inference_chunk_size.toString().toInteger()
        def phaseChunkSize = params.phase_chunk_size.toString().toInteger()
        if (chunkSize < 1 || phaseChunkSize < 1) {
            error 'Inference and phasing chunk sizes must be positive'
        }
        noReuse = Channel.value([])
        nativeReuse = Channel.value(params.reuse_manifest ? file(params.reuse_manifest, checkIfExists: true) : [])
        frozen = FREEZE(file(params.experiment, checkIfExists: true), nativeReuse)
        jobs = frozen.flatMap { folder ->
            def config = new groovy.json.JsonSlurper().parse(file("$folder/config.json").toFile())
            config.routes.collectMany { route -> config.chromosomes.collect { chrom -> tuple(route.id, chrom, folder) } }
        }
        panels = PREPARE(jobs, noReuse)
        phaseJobs = panels.flatMap { route, chrom, panel ->
            def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
            def starts = (0..<meta.queries.size()).step(phaseChunkSize)
            starts.collect { first ->
                tuple(route, chrom, first, Math.min(first + phaseChunkSize, meta.queries.size()), starts.size(), panel)
            }
        }
        phaseChunks = GUIDED_PHASE(phaseJobs, file(params.phase_reuse, checkIfExists: true))
        phaseGroups = phaseChunks.map { route, chrom, count, folder ->
            tuple(groupKey(tuple(route, chrom), count), folder)
        }.groupTuple().map { key, folders ->
            tuple(key.getGroupTarget()[0], key.getGroupTarget()[1], folders)
        }.combine(panels, by: [0, 1])
        phases = COLLECT_GUIDED_PHASE(phaseGroups)
        baseJobs = panels.flatMap { route, chrom, panel ->
            def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
            meta.base_plan.bases.collect { id, spec -> tuple(route, chrom, id, spec.family, panel) }
        }
        bases = BASE(baseJobs, noReuse)
        fitJobs = bases.combine(panels, by: [0, 1]).flatMap { route, chrom, baseid, base, panel ->
            def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
            meta.base_plan.models.findAll { it.base_id == baseid }.collect {
                model -> tuple(route, chrom, model.id, model.family, base, panel)
            }
        }
        fits = FIT(fitJobs, noReuse)
        inferJobs = fits.combine(panels, by: [0, 1]).combine(phases, by: [0, 1]).map {
            route, chrom, model, fitted, panel, phased ->
            def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
            def spec = meta.base_plan.models.find { it.id == model }
            tuple(route, chrom, spec.base_id, model, spec.family, fitted, panel, phased)
        }.combine(bases, by: [0, 1, 2]).flatMap {
            route, chrom, baseid, model, family, fitted, panel, phased, base ->
            def meta = new groovy.json.JsonSlurper().parse(file("$panel/panel.json").toFile())
            def starts = (0..<meta.queries.size()).step(chunkSize)
            starts.collect { first ->
                tuple(route, chrom, model, family, first, Math.min(first + chunkSize, meta.queries.size()),
                      starts.size(), fitted, base, panel, phased)
            }
        }
        inferred = INFER(inferJobs, nativeReuse)
        completed = inferred.map { route, chrom, model, family, count, folder ->
            tuple(groupKey(tuple(route, chrom, model), count), folder)
        }.groupTuple().map { key, folders ->
            def identity = key.getGroupTarget()
            tuple(identity[0], identity[1], identity[2], folders)
        }.combine(panels, by: [0, 1])
        assembled = COLLECT_INFERENCE(completed)
        aggregates = AGGREGATE(assembled.groupTuple(by: [0, 1]))
        interpreted = INTERPRET(aggregates.map { route, model, aggregate ->
            tuple(route, model, 'native', aggregate)
        })
        summary = SELECT(interpreted.collect())
        VERIFY(summary)
    } else {
        error 'guided_stage must be plan, run or report'
    }
}
