nextflow.enable.dsl = 2
params.freeze_script = "${params.project}/pedigree_benchmark/scripts/freeze_production_cohort.py"
params.reference_dir = null
params.prepare_references = false

process FREEZE_GENETICALLY_SCREENED_COHORT {
    publishDir "${params.run_dir}/cohort", mode: 'copy'
    cpus 1
    memory '6 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    output:
    path 'frozen'
    script:
    """
    python3 '${params.freeze_script}' --project '${params.project}' --output frozen
    """
}

workflow {
    def mode = params.prepare_references.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) {
        error 'prepare_references must be explicitly true or false'
    }
    if (mode == 'true') {
        PURGE_BENCHMARK_REFERENCES()
    } else {
        FREEZE_GENETICALLY_SCREENED_COHORT()
    }
}

process PREPARE_DONOR_PURGED_IMPUTATION_REFERENCE {
    tag "chr${chromosome}"
    publishDir "${params.reference_dir ?: params.run_dir}/beagle_references", mode: 'copy'
    cpus 1
    memory '6 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    input:
    val chromosome
    output:
    path '*.vcf.gz'
    path '*.vcf.gz.csi'
    path 'chr*.reference.json'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/freeze_dense_split.py' \
      --project '${params.project}' --cohort '${params.cohort}' \
      --reference-chromosome ${chromosome} --output .
    """
}

workflow PURGE_BENCHMARK_REFERENCES {
    def requested = params.chromosomes instanceof Collection ? params.chromosomes : params.chromosomes.toString().tokenize(',')
    def chromosomes = requested.collect { it as int }
    if (!chromosomes || chromosomes.toSet().size() != chromosomes.size() || chromosomes.any { it < 1 || it > 22 }) {
        error 'Require unique autosome numbers in chromosomes'
    }
    PREPARE_DONOR_PURGED_IMPUTATION_REFERENCE(Channel.fromList(chromosomes))
}
