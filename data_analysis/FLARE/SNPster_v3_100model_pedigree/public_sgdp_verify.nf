nextflow.enable.dsl = 2

process VERIFY_PUBLISHED_CHROMOSOME {
    tag "chr${chromosome}"
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/chromosomes", mode: 'copy'
    input:
    tuple val(chromosome), path(calls), path(call_index)
    output:
    path "chr${chromosome}.json"
    script:
    def root = chromosome == 22 ? params.smoke_dir : "${params.import_dir}/chromosomes/chr${chromosome}"
    """
    python3 "${params.project}/pedigree_benchmark/scripts/verify_public_sgdp.py" \
        --run "${root}" --metadata-run "${params.smoke_dir}" --chromosome ${chromosome} \
        --output "chr${chromosome}.json"
    """
}

process COLLECT_INDEPENDENT_VERIFICATION {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path reports
    output:
    path 'published_verification.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" collect \
        --import-run "${params.import_dir}" --smoke "${params.smoke_dir}" \
        --verification-run "${params.run_dir}" \
        --reports ${reports.collect { '"' + it + '"' }.join(' ')} \
        --output published_verification.json
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
