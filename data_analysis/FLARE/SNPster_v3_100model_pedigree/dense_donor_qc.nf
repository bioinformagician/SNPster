nextflow.enable.dsl = 2

process VERIFY_DENSE_QC_INPUTS {
    publishDir "${params.run_dir}/preparation", mode: 'copy'
    input:
    path dense_manifest
    path previous_manifest
    path original_preflight
    path normalization_chain
    output:
    path 'plan.json', emit: plan
    path 'policy', emit: policy
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_donor_qc.py' prepare \
      --project '${params.project}' --dense-run '${params.dense_run}'
    """
}

process MERGE_OBSERVED_IDENTITY_CANDIDATES {
    tag "chr${chromosome}"
    input:
    tuple val(chromosome), path(new_calls, stageAs: 'new/normalized.vcf.gz'),
          path(new_index, stageAs: 'new/normalized.vcf.gz.csi'),
          path(old_calls, stageAs: 'old/normalized.vcf.gz'),
          path(old_index, stageAs: 'old/normalized.vcf.gz.csi'), path(normalization_proof)
    path plan
    output:
    tuple val(chromosome), path("chr${chromosome}.candidates"), path(new_calls), path(new_index),
          path(old_calls), path(old_index), path(normalization_proof)
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_donor_qc.py' merge \
      --plan '${plan}' --chromosome ${chromosome} --new '${new_calls}' --old '${old_calls}' \
      --verification '${normalization_proof}' --output 'chr${chromosome}.candidates'
    """
}

process VERIFY_OBSERVED_CANDIDATE_MERGE {
    tag "chr${chromosome}"
    publishDir "${params.run_dir}/candidate_merges", mode: 'copy'
    input:
    tuple val(chromosome), path(folder), path(new_calls, stageAs: 'new/normalized.vcf.gz'),
          path(new_index, stageAs: 'new/normalized.vcf.gz.csi'),
          path(old_calls, stageAs: 'old/normalized.vcf.gz'),
          path(old_index, stageAs: 'old/normalized.vcf.gz.csi'), path(normalization_proof)
    path plan
    output:
    tuple val(chromosome), path(folder), path("chr${chromosome}.merge.verified.json")
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_donor_qc.py' verify-merge \
      --plan '${plan}' --chromosome ${chromosome} --new '${new_calls}' --old '${old_calls}' \
      --verification '${normalization_proof}' --folder '${folder}' --output 'chr${chromosome}.merge.verified.json'
    """
}

process BUILD_DENSE_IDENTITY_PANEL {
    tag "chr${chromosome}"
    input:
    tuple val(chromosome), path(folder, stageAs: 'candidate_input'), path(verification)
    output:
    tuple val(chromosome), path("chr${chromosome}"), path(folder)
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py' panel \
      --project '${params.project}' --chromosome ${chromosome} \
      --normalized '${folder}/candidates.vcf.gz' --verification '${verification}' \
      --production-root '${params.production_root}' --output 'chr${chromosome}'
    """
}

process VERIFY_EVERY_IDENTITY_GENOTYPE {
    tag "chr${chromosome}"
    publishDir "${params.run_dir}/identity_panels", mode: 'copy'
    input:
    tuple val(chromosome), path(identity_panel), path(merged, stageAs: 'candidate_input')
    path plan
    output:
    tuple val(chromosome), path(identity_panel), path("chr${chromosome}.verified.json")
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_donor_qc.py' verify-panel \
      --plan '${plan}' --chromosome ${chromosome} --folder '${identity_panel}' --merged '${merged}' \
      --production-root '${params.production_root}' \
      --preflight '${params.project}/results_v3_production_preflight_20261021_193442/summary/preflight.json' \
      --output 'chr${chromosome}.verified.json'
    """
}

process SCREEN_ALL_CANDIDATE_KING_PAIRS {
    memory '6 GB'
    publishDir "${params.run_dir}", mode: 'copy', pattern: 'kinship'
    input:
    path panels
    path qc_policy
    output:
    path 'kinship'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py' kinship \
      --project '${qc_policy}' --panels ${panels.collect { '"' + it + '"' }.join(' ')} --output kinship
    """
}

process CLASSIFY_DENSE_IDENTITY_AND_FAMILY_ROLES {
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path plan
    path kinship
    path panels
    path reports
    output:
    path 'genetic_screen.json'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/dense_donor_qc.py' collect \
      --project '${params.project}' --plan '${plan}' --kinship '${kinship}' \
      --panels ${panels.collect { '"' + it + '"' }.join(' ')} \
      --reports ${reports.collect { '"' + it + '"' }.join(' ')} --output genetic_screen.json
    """
}

workflow {
    def mode = params.smoke.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) { error "smoke must be true or false" }
    prepared_inputs = VERIFY_DENSE_QC_INPUTS(
        file("${params.dense_run}/full/published_verification.json", checkIfExists: true),
        file("${params.previous_run}/summary/normalization_complete.json", checkIfExists: true),
        file("${params.project}/results_v3_production_preflight_20261021_193442/summary/preflight.json", checkIfExists: true),
        file('/home/frederik/github_projects/SNPster/data pipeline/standardizer_module/dependencies/hg19ToHg38.over.chain.gz', checkIfExists: true))
    scoped_jobs = prepared_inputs.plan.flatMap { plan_file ->
        def specification = new groovy.json.JsonSlurper().parse(plan_file)
        def chromosomes = mode == 'true' ? [22] : (1..21).toList()
        chromosomes.collect { chromosome ->
            def job = specification.jobs["${chromosome}"]
            def new_path = job['new']
            def old_path = job['old']
            tuple(chromosome, file(new_path, checkIfExists: true), file("${new_path}.csi", checkIfExists: true),
                  file(old_path, checkIfExists: true), file("${old_path}.csi", checkIfExists: true),
                  file(job.verification, checkIfExists: true))
        }
    }
    merged_candidates = MERGE_OBSERVED_IDENTITY_CANDIDATES(scoped_jobs, prepared_inputs.plan)
    verified_candidates = VERIFY_OBSERVED_CANDIDATE_MERGE(merged_candidates, prepared_inputs.plan)
    candidate_identity_panels = BUILD_DENSE_IDENTITY_PANEL(verified_candidates)
    verified_identity_panels = VERIFY_EVERY_IDENTITY_GENOTYPE(candidate_identity_panels, prepared_inputs.plan)
    if (mode == 'false') {
        reused_panel = Channel.value(file("${params.run_dir}/../smoke/identity_panels/chr22", checkIfExists: true))
        all_panels = verified_identity_panels.map { chromosome, directory, proof -> directory }.mix(reused_panel)
        collected_panels = all_panels.collect()
        reused_proof = Channel.value(file("${params.run_dir}/../smoke/identity_panels/chr22.verified.json", checkIfExists: true))
        collected_proofs = verified_identity_panels.map { chromosome, directory, proof -> proof }.mix(reused_proof).collect()
        kinship_outputs = SCREEN_ALL_CANDIDATE_KING_PAIRS(collected_panels, prepared_inputs.policy)
        CLASSIFY_DENSE_IDENTITY_AND_FAMILY_ROLES(prepared_inputs.plan, kinship_outputs, collected_panels, collected_proofs)
    }
}
