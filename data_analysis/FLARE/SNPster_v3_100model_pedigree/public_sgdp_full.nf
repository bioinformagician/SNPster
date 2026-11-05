nextflow.enable.dsl = 2

process REQUIRE_VERIFIED_PUBLIC_SMOKE {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/preflight", mode: 'copy'
    input:
    path verification
    output:
    path 'smoke_gate.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/prepare_sgdp_full_import.py" require \
        --smoke "${params.smoke_dir}" --project "${params.project}"
    """
}

process DOWNLOAD_PUBLIC_LIKELIHOODS {
    tag "chr${chromosome}"
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/chromosomes/chr${chromosome}/source", mode: 'copy', pattern: 'download.json'
    input:
    tuple val(chromosome), path(gate)
    output:
    tuple val(chromosome), path("chr.sgdp.pub.${chromosome}.bcf"), emit: data
    path 'download.json', emit: provenance
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" download \
        --checksums "${params.smoke_dir}/metadata/publisher_md5sums.txt" --chromosome ${chromosome} \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process CALL_RAW_SEQUENCING_GENOTYPES {
    tag "chr${chromosome}"
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/chromosomes/chr${chromosome}/basecalls", mode: 'copy'
    input:
    tuple val(chromosome), path(source)
    output:
    tuple val(chromosome), path('basecalls.vcf.gz'), emit: vcf
    path 'basecalls.vcf.gz.csi'
    path 'basecalls.json'
    path 'likelihood_query.command.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" basecalls \
        --source "${source}" --selected "${params.smoke_dir}/metadata/selected.json" --chromosome ${chromosome} \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process VERIFY_RAW_CALL_CONCORDANCE {
    tag "chr${chromosome}"
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/chromosomes/chr${chromosome}/verification", mode: 'copy'
    input:
    tuple val(chromosome), path(calls)
    output:
    path 'verification.json'
    script:
    """
    cp "${params.smoke_dir}/raw_control/raw_control.vcf.gz" raw_control.vcf.gz
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" verify \
        --source "${calls}" --control raw_control.vcf.gz \
        --selected "${params.smoke_dir}/metadata/selected.json" --chromosome ${chromosome} \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process COLLECT_PUBLIC_DATA_QC {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path reports, stageAs: 'qc??/*'
    output:
    path 'full_import.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/prepare_sgdp_full_import.py" collect \
        --smoke "${params.smoke_dir}" --run "${params.run_dir}" \
        --reports ${reports.collect { '"' + it + '"' }.join(' ')}
    """
}

workflow {
    gate = REQUIRE_VERIFIED_PUBLIC_SMOKE(file("${params.smoke_dir}/published_verification.json", checkIfExists: true))
    jobs = Channel.fromList((1..21).toList()).combine(gate)
    downloads = DOWNLOAD_PUBLIC_LIKELIHOODS(jobs)
    calls = CALL_RAW_SEQUENCING_GENOTYPES(downloads.data)
    reports = VERIFY_RAW_CALL_CONCORDANCE(calls.vcf)
    COLLECT_PUBLIC_DATA_QC(reports.collect())
}
