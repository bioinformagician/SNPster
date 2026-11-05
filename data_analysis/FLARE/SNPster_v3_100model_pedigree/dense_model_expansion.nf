nextflow.enable.dsl = 2

include { VERIFY_AND_SCORE_COMPLETE_BENCHMARK } from './production_matched_benchmark.nf'

process AWAIT_VERIFIED_UPSTREAM_COMPLETION {
    cpus 1
    memory '512 MB'
    container 'sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'
    containerOptions "${params.container_options} --group-add ${params.socket_gid} -v /var/run/docker.sock:/var/run/docker.sock"
    output:
    path 'upstream.done'
    script:
    """
    echo 'Waiting for ${params.upstream_controller}; existing preprocessing and ten baseline fits are preserved.'
    status=\$(docker wait '${params.upstream_controller}')
    if [[ "\$status" != 0 ]]; then
        echo "Upstream analysis stopped with exit \$status. Expansion cannot bypass its scientific gates." >&2
        docker logs --tail 60 '${params.upstream_controller}'
        exit 1
    fi
    test -s '${params.source_experiment}/full/post_cleanup_verification.json'
    printf '%s\\n' '${params.upstream_controller}' > upstream.done
    """
}

process VERIFY_REUSABLE_PRODUCTION_SCOPE {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/source_verification", mode: 'copy'
    input:
    path ready
    output:
    path 'source_verified.json'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' verify-source \
      --cohort '${params.cohort}' --source-scope '${params.source_scope}' --output source_verified.json
    """
}

process PREPARE_EXPANDED_REFERENCE_PANELS {
    tag "chr${chromosome}"
    memory '8 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/panels/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), path(source_panel), path(gate)
    output:
    tuple val(chromosome), path('panel')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' panels \
      --cohort '${params.cohort}' --folder '${source_panel}' --gate '${gate}' --output panel
    """
}

process FIT_SHARED_GNOMIX_REFERENCE_BASE {
    tag "${window}cM/chr${chromosome}"
    memory '8 GB'
    container 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21'
    beforeScript 'export MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache'
    publishDir path: { "${params.run_dir}/reference_bases/chr${chromosome}/window_${window}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(window), path(panel)
    output:
    tuple val(chromosome), val(window), path('base')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' base \
      --folder '${panel}' --window '${window}' --output base
    """
}

process RUN_EXPANDED_FLARE_OR_REUSED_BASELINE {
    tag "${candidate}/chr${chromosome}"
    memory '8 GB'
    container { candidate.startsWith('gnomix') ? 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21' : 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d' }
    beforeScript 'export MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache'
    publishDir path: { "${params.run_dir}/models/${candidate}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(candidate), path(panel)
    output:
    tuple val(chromosome), val(candidate), path('model')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' evaluate \
      --cohort '${params.cohort}' --folder '${panel}' --candidate '${candidate}' --output model
    """
}

process RUN_EXPANDED_GNOMIX_SMOOTHER {
    tag "${candidate}/chr${chromosome}"
    memory '8 GB'
    container 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21'
    beforeScript 'export MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba-cache'
    publishDir path: { "${params.run_dir}/models/${candidate}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(candidate), path(panel), path(base)
    output:
    tuple val(chromosome), val(candidate), path('model')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' evaluate \
      --cohort '${params.cohort}' --folder '${panel}' --candidate '${candidate}' --base '${base}' --output model
    """
}

process GATE_EXPANDED_WHOLE_GENOME_CAPACITY {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/capacity", mode: 'copy'
    input:
    path certificate
    output:
    path 'source_verified.json', emit: gate
    path 'capacity.json', emit: capacity
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_model_expansion.py' capacity \
      --cohort '${params.cohort}' --source-scope '${params.source_scope}' \
      --folder '${params.run_dir}/../smoke' --output capacity.json
    cp '${certificate}' source_verified.json
    """
}

workflow {
    def mode = params.smoke.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) { error 'smoke must be explicitly true or false' }
    def cohort = new groovy.json.JsonSlurper().parse(file(params.cohort))
    if (cohort.protocol.candidates.size() != 98) { error 'Require the frozen98-configuration expansion' }
    params.chromosomes = mode == 'true' ? [22] : (1..22).toList()
    params.masks = mode == 'true' ? ['lonni'] : cohort.protocol.masks
    def baseline = cohort.candidate_expansion.baseline_candidates.toSet()
    log.info "98 distinct settings:64 FLARE,32 fast Gnomix,2 NNLS; eight unique baselines reused; two old aliases retained in provenance."
    ready = AWAIT_VERIFIED_UPSTREAM_COMPLETION()
    verified = VERIFY_REUSABLE_PRODUCTION_SCOPE(ready)
    if (mode == 'false') {
        gate = GATE_EXPANDED_WHOLE_GENOME_CAPACITY(verified).gate
    } else {
        gate = verified
    }
    inputs = Channel.fromList(params.chromosomes).combine(gate).map { chromosome, certificate ->
        tuple(chromosome, file("${params.source_scope}/panels/chr${chromosome}/panel"), certificate)
    }
    panels = PREPARE_EXPANDED_REFERENCE_PANELS(inputs)
    native_jobs = panels.flatMap { chromosome, folder ->
        cohort.protocol.candidates.findAll { candidate, policy ->
            policy.family != 'gnomix' || baseline.contains(candidate)
        }.collect { candidate, policy -> tuple(chromosome, candidate, folder) }
    }
    base_jobs = panels.flatMap { chromosome, folder ->
        cohort.candidate_expansion.policy.gnomix.window_sizes_cM.collect { window -> tuple(chromosome, window, folder) }
    }
    bases = FIT_SHARED_GNOMIX_REFERENCE_BASE(base_jobs)
    smooth_jobs = bases.combine(panels, by: 0).flatMap { chromosome, window, base, panel ->
        cohort.protocol.candidates.findAll { candidate, policy ->
            policy.family == 'gnomix' && !baseline.contains(candidate) && policy.window_size_cM == window
        }.collect { candidate, policy -> tuple(chromosome, candidate, panel, base) }
    }
    native_models = RUN_EXPANDED_FLARE_OR_REUSED_BASELINE(native_jobs)
    gnomix_models = RUN_EXPANDED_GNOMIX_SMOOTHER(smooth_jobs)
    complete = native_models.mix(gnomix_models).map { chromosome, candidate, folder -> folder }.collect()
    VERIFY_AND_SCORE_COMPLETE_BENCHMARK(complete)
}
