nextflow.enable.dsl = 2

process PREPARE_ASSEMBLY_RESOURCES {
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/inputs", mode: 'copy'
    input:
    path input_manifest
    output:
    tuple path('hs37d5.fa'), path('hs37d5.fa.fai'), path('preflight.json')
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/normalize_sgdp.py" prepare \
        --project "${params.project}" --source "${params.source_fasta_zst}" \
        --target "${params.target_fasta}" --chain "${params.chain}"
    """
}

process NORMALIZE_SEQUENCING_CALLS {
    tag "chr${chromosome}"
    cpus 1
    memory '3 GB'
    publishDir "${params.run_dir}/normalized", mode: 'copy', pattern: 'chr*'
    input:
    tuple val(chromosome), path(calls), path(call_index), path(source), path(source_index), path(gate)
    output:
    tuple val(chromosome), path("chr${chromosome}"), path(calls)
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/normalize_sgdp.py" normalize \
        --project "${params.project}" --source "${source}" --target "${params.target_fasta}" \
        --chain "${params.chain}" --calls "${calls}" --chromosome ${chromosome} --output "chr${chromosome}"
    """
}

process VERIFY_NORMALIZED_GENOTYPES {
    tag "chr${chromosome}"
    cpus 1
    memory '3 GB'
    publishDir "${params.run_dir}/verification", mode: 'copy', pattern: '*.verified.json'
    input:
    tuple val(chromosome), path(folder), path(calls), path(source), path(source_index), path(gate)
    output:
    tuple val(chromosome), path(folder), path("chr${chromosome}.verified.json")
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/normalize_sgdp.py" verify \
        --project "${params.project}" --source "${source}" --target "${params.target_fasta}" \
        --chain "${params.chain}" --calls "${calls}" --folder "${folder}" \
        --chromosome ${chromosome} --output "chr${chromosome}.verified.json"
    """
}

process PREPARE_GENETIC_IDENTITY_PANEL {
    tag "chr${chromosome}"
    cpus 1
    memory '3 GB'
    publishDir "${params.run_dir}/identity_panels", mode: 'copy', pattern: 'chr*'
    input:
    tuple val(chromosome), path(folder), path(verification)
    output:
    path "chr${chromosome}"
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py" panel \
        --project "${params.project}" --chromosome ${chromosome} \
        --normalized "${folder}/normalized.vcf.gz" --verification "${verification}" \
        --production-root "${params.production_root}" --output "chr${chromosome}"
    """
}

process AUDIT_CANDIDATE_RELATEDNESS {
    cpus 1
    memory '6 GB'
    publishDir "${params.run_dir}", mode: 'copy', pattern: 'kinship'
    input:
    path panels
    output:
    path 'kinship'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py" kinship \
        --project "${params.project}" --panels ${panels.collect { '"' + it + '"' }.join(' ')} \
        --output kinship
    """
}

process COLLECT_NORMALIZATION_AND_IDENTITY_QC {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path folders, stageAs: 'normalized/*'
    path reports
    path panels, stageAs: 'panels/*'
    path kinship
    output:
    path 'normalization_complete.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py" collect \
        --project "${params.project}" --run-dir "${params.run_dir}" \
        --folders ${folders.collect { '"' + it + '"' }.join(' ')} \
        --reports ${reports.collect { '"' + it + '"' }.join(' ')} \
        --panels ${panels.collect { '"' + it + '"' }.join(' ')} \
        --kinship "${kinship}" --output normalization_complete.json
    """
}

workflow {
    manifest = new groovy.json.JsonSlurper().parse(file(params.input_manifest))
    resources = PREPARE_ASSEMBLY_RESOURCES(file(params.input_manifest, checkIfExists: true))
    chromosomes = params.smoke ? [22] : (1..22).toList()
    jobs = Channel.fromList(chromosomes.collect { chromosome ->
        def root = manifest.data_roots["${chromosome}"]
        tuple(chromosome, file("${root}/basecalls/basecalls.vcf.gz", checkIfExists: true),
              file("${root}/basecalls/basecalls.vcf.gz.csi", checkIfExists: true))
    })
    normalized = NORMALIZE_SEQUENCING_CALLS(jobs.combine(resources))
    verified = VERIFY_NORMALIZED_GENOTYPES(normalized.combine(resources))
    panels = PREPARE_GENETIC_IDENTITY_PANEL(verified)
    if (!params.smoke) {
        relatedness = AUDIT_CANDIDATE_RELATEDNESS(panels.collect())
        COLLECT_NORMALIZATION_AND_IDENTITY_QC(
            normalized.map { chromosome, folder, calls -> folder }.collect(),
            verified.map { chromosome, folder, report -> report }.collect(),
            panels.collect(), relatedness)
    }
}
