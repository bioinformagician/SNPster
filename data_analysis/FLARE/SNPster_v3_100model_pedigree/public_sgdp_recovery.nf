nextflow.enable.dsl = 2

process VERIFY_PUBLISHED_CHROMOSOME {
    tag "chr${chromosome}"
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/chromosomes", mode: 'copy'
    input:
    tuple val(chromosome), path(calls), path(call_index)
    output:
    path "chr${chromosome}"
    script:
    def root = chromosome == 22 ? params.smoke_dir : "${params.import_dir}/chromosomes/chr${chromosome}"
    def reuse = [1,2,3,4,5,7].contains(chromosome) ?
        "--previous-report \"${params.previous_dir}/chromosomes/chr${chromosome}.json\"" : ''
    """
    python3 "${params.project}/pedigree_benchmark/scripts/recover_sgdp_duplicates.py" recover \
        --root "${root}" --smoke "${params.smoke_dir}" --project "${params.project}" \
        --chromosome ${chromosome} ${reuse} --output "chr${chromosome}"
    """
}

process COLLECT_INDEPENDENT_VERIFICATION {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path folders
    output:
    path 'published_verification.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/recover_sgdp_duplicates.py" collect \
        --import-run "${params.import_dir}" --smoke "${params.smoke_dir}" \
        --previous-run "${params.previous_dir}" --recovery-run "${params.run_dir}" \
        --folders ${folders.collect { '"' + it + '"' }.join(' ')} --output published_verification.json
    """
}

workflow {
    jobs = Channel.fromList((1..22).collect { chromosome ->
        def root = chromosome == 22 ? params.smoke_dir : "${params.import_dir}/chromosomes/chr${chromosome}"
        tuple(chromosome, file("${root}/basecalls/basecalls.vcf.gz", checkIfExists: true),
              file("${root}/basecalls/basecalls.vcf.gz.csi", checkIfExists: true))
    })
    reports = VERIFY_PUBLISHED_CHROMOSOME(jobs)
    COLLECT_INDEPENDENT_VERIFICATION(reports.collect())
}
