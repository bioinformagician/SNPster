nextflow.enable.dsl = 2

process FETCH_PUBLIC_METADATA {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/metadata", mode: 'copy'
    output:
    path 'selected.json', emit: people
    path 'publisher_md5sums.txt', emit: checksums
    path 'metadata.original.txt'
    path 'likelihood_readme.txt'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" metadata \
        --project "${params.project}" --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process DOWNLOAD_PUBLIC_LIKELIHOODS {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/source", mode: 'copy', pattern: 'download.json'
    input:
    path checksums
    output:
    path 'chr.sgdp.pub.22.bcf'
    path 'download.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" download \
        --checksums "${checksums}" --chromosome 22 \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process FETCH_ORIGINAL_RAW_CONTROL {
    cpus 1
    memory '1 GB'
    publishDir "${params.run_dir}/raw_control", mode: 'copy'
    output:
    path 'raw_control.vcf.gz', emit: vcf
    path 'raw_control.json'
    path 'raw_control.md5sum'
    path 'public_tar_members.txt'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" raw_control \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process CALL_RAW_SEQUENCING_GENOTYPES {
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/basecalls", mode: 'copy'
    input:
    path source
    path selected
    output:
    path 'basecalls.vcf.gz', emit: vcf
    path 'basecalls.vcf.gz.csi'
    path 'basecalls.json'
    path 'likelihood_query.command.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" basecalls \
        --source "${source}" --selected "${selected}" --chromosome 22 \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

process VERIFY_RAW_CALL_CONCORDANCE {
    cpus 1
    memory '2 GB'
    publishDir "${params.run_dir}/verification", mode: 'copy'
    input:
    path calls
    path control
    path selected
    output:
    path 'verification.json'
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/import_sgdp_public.py" verify \
        --source "${calls}" --control "${control}" --selected "${selected}" --chromosome 22 \
        --protocol "${params.project}/config/sgdp_public_import.json"
    """
}

workflow {
    metadata = FETCH_PUBLIC_METADATA()
    download = DOWNLOAD_PUBLIC_LIKELIHOODS(metadata.checksums)
    control = FETCH_ORIGINAL_RAW_CONTROL()
    calls = CALL_RAW_SEQUENCING_GENOTYPES(download[0], metadata.people)
    VERIFY_RAW_CALL_CONCORDANCE(calls.vcf, control.vcf, metadata.people)
}
