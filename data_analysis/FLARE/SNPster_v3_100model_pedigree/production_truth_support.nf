nextflow.enable.dsl = 2

process RECOVER_ORIGINAL_ASSAY_CALLS {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/assay", mode: 'copy'
    output:
    path 'assay'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_production_truth_support.py' assay --output assay
    """
}

process RECOVER_OBSERVED_CHROMOSOME_SUPPORT {
    tag "chr${chromosome}"
    memory '8 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/chromosomes/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), path(assay)
    output:
    path 'support'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_production_truth_support.py' chromosome \
      --cohort '${params.cohort}' --chromosome ${chromosome} --assay '${assay}' --run-dir '${params.run_dir}' --output support
    """
}

process FREEZE_ALL22_OBSERVED_SUPPORT {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/frozen", mode: 'copy'
    input:
    path(folders, stageAs: 'chromosomes/support??')
    output:
    path 'frozen'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/prepare_production_truth_support.py' collect \
      --cohort '${params.cohort}' --run-dir '${params.run_dir}' --folders ${folders.join(' ')} --output frozen
    """
}

workflow {
    assay = RECOVER_ORIGINAL_ASSAY_CALLS()
    jobs = assay.flatMap { folder -> (1..22).collect { chromosome -> tuple(chromosome, folder) } }
    support = RECOVER_OBSERVED_CHROMOSOME_SUPPORT(jobs)
    FREEZE_ALL22_OBSERVED_SUPPORT(support.collect())
}
