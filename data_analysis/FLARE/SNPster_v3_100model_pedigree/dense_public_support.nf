nextflow.enable.dsl = 2

process VERIFY_DENSE_NORMALIZATION_RESOURCES {
    publishDir "${params.run_dir}/resources", mode: 'copy', pattern: 'resource_verification.json'
    input:
    tuple path(source), path(source_index), path(target), path(target_index), path(chain)
    output:
    tuple path(source), path(source_index), path(target), path(target_index), path(chain), emit: assemblies
    path 'resource_verification.json', emit: report
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_dense_public_support.py' resources \
      --project '${params.project}' --source '${source}' --target '${target}' --chain '${chain}' \
      --output resource_verification.json
    """
}

process SELECT_DENSE_PUBLIC_CANDIDATES {
    publishDir "${params.run_dir}/metadata", mode: 'copy'
    input:
    path resource_verification
    output:
    path 'selected.json', emit: people
    path 'recall.json', emit: recall
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_dense_public_support.py' metadata \
      --project '${params.project}' --protocol '${params.project}/config/dense_public_support.json'
    """
}

process DOWNLOAD_DENSE_LIKELIHOODS {
    tag "chr${chromosome}"
    memory '1 GB'
    containerOptions "${params.container_options} --network bridge"
    input:
    tuple val(chromosome), path(selection)
    output:
    tuple val(chromosome), path("chr.sgdp.pub.${chromosome}.bcf"), path('download.json')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py' download \
      --checksums '${params.public_smoke}/metadata/publisher_md5sums.txt' --chromosome ${chromosome} \
      --protocol '${params.project}/config/sgdp_public_import.json'
    """
}

process RECALL_DENSE_OBSERVED_GENOTYPES {
    tag "chr${chromosome}"
    input:
    tuple val(chromosome), path(source), path(download)
    path selection
    output:
    tuple val(chromosome), path('basecalls.vcf.gz'), path('basecalls.json'), path(download), path('likelihood_query.command.json')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py' basecalls \
      --source '${source}' --selected '${selection}' --chromosome ${chromosome} \
      --protocol '${params.project}/config/sgdp_public_import.json'
    """
}

process VERIFY_RAW_SEQUENCING_EVIDENCE {
    tag "chr${chromosome}"
    publishDir path: { "${params.run_dir}/raw_observed/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), path(calls), path(summary), path(download), path(query)
    path selection
    output:
    tuple val(chromosome), path('observed')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_dense_public_support.py' verify-raw \
      --calls '${calls}' --summary '${summary}' --selection '${selection}' \
      --protocol '${params.project}/config/sgdp_public_import.json' --output observed
    cp '${summary}' observed/original_caller.json
    cp '${download}' observed/download.json
    cp '${query}' observed/likelihood_query.command.json
    """
}

process CHECK_ORIGINAL_SEQUENCING_CONTROL {
    tag "chr${chromosome}"
    publishDir path: { "${params.run_dir}/raw_control/chr${chromosome}" }, mode: 'copy', pattern: 'verification.json'
    input:
    tuple val(chromosome), path(observed)
    path selection
    output:
    tuple val(chromosome), path('verification.json')
    script:
    """
    cp '${params.public_smoke}/raw_control/raw_control.vcf.gz' raw_control.vcf.gz
    python3 '${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py' verify \
      --source '${observed}/basecalls.vcf.gz' --control raw_control.vcf.gz --selected '${selection}' \
      --chromosome ${chromosome} --protocol '${params.project}/config/sgdp_public_import.json'
    """
}

process NORMALIZE_DENSE_OBSERVED_GENOTYPES {
    tag "chr${chromosome}"
    publishDir "${params.run_dir}/normalized", mode: 'copy', pattern: 'chr*'
    input:
    tuple val(chromosome), path(observed), path(source), path(source_index), path(target), path(target_index), path(chain)
    output:
    tuple val(chromosome), path("chr${chromosome}"), path(observed), path(source), path(source_index), path(target), path(target_index), path(chain)
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/normalize_sgdp.py' normalize \
      --project '${params.project}' --calls '${observed}/basecalls.vcf.gz' --chromosome ${chromosome} \
      --source '${source}' --target '${target}' --chain '${chain}' --output chr${chromosome}
    """
}

process VERIFY_DENSE_COORDINATES_AND_ALLELES {
    tag "chr${chromosome}"
    publishDir "${params.run_dir}/verification", mode: 'copy', pattern: '*.verified.json'
    input:
    tuple val(chromosome), path(folder), path(observed), path(source), path(source_index), path(target), path(target_index), path(chain)
    output:
    tuple val(chromosome), path(folder), path("chr${chromosome}.verified.json")
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/normalize_sgdp.py' verify \
      --project '${params.project}' --calls '${observed}/basecalls.vcf.gz' --chromosome ${chromosome} \
      --source '${source}' --target '${target}' --chain '${chain}' --folder '${folder}' \
      --output chr${chromosome}.verified.json
    """
}

process AUDIT_REAL_MASK_SUPPORT {
    tag "chr${chromosome}"
    publishDir path: { "${params.run_dir}/density/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), path(folder), path(verification), path(control)
    path selection
    output:
    path 'density'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_dense_public_support.py' density \
      --project '${params.project}' --selection '${selection}' --chromosome ${chromosome} \
      --normalized '${folder}/normalized.vcf.gz' --verification '${verification}' \
      --raw-control '${control}' --output density
    """
}

workflow {
    def mode = params.smoke.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) { error "smoke must be true or false" }
    def chromosomes = mode == 'true' ? [22] : (1..21).toList()
    resource_checks = VERIFY_DENSE_NORMALIZATION_RESOURCES(Channel.value(tuple(
        file(params.source_fasta, checkIfExists: true), file("${params.source_fasta}.fai", checkIfExists: true),
        file(params.target_fasta, checkIfExists: true), file("${params.target_fasta}.fai", checkIfExists: true),
        file(params.chain, checkIfExists: true))))
    metadata = SELECT_DENSE_PUBLIC_CANDIDATES(resource_checks.report)
    people = metadata.people
    recall = metadata.recall
    downloads = DOWNLOAD_DENSE_LIKELIHOODS(Channel.fromList(chromosomes).combine(people))
    called = RECALL_DENSE_OBSERVED_GENOTYPES(downloads, recall)
    verified_observations = VERIFY_RAW_SEQUENCING_EVIDENCE(called, recall)
    control_reports = CHECK_ORIGINAL_SEQUENCING_CONTROL(verified_observations, recall)
    normalized = NORMALIZE_DENSE_OBSERVED_GENOTYPES(verified_observations.combine(resource_checks.assemblies))
    verified = VERIFY_DENSE_COORDINATES_AND_ALLELES(normalized)
    AUDIT_REAL_MASK_SUPPORT(verified.join(control_reports), people)
}
