#!/usr/bin/env nextflow

params.dependencies_dir = "/mnt/c/Users/frezz/Desktop/dependencies"
params.genefile_dir     = "/mnt/c/Users/frezz/Desktop/genome_scraping/scraped_genomes/subset"
params.output_dir       = "/mnt/c/Users/frezz/Desktop/docker_testing/nf_output"

process harmonizer {

    errorStrategy 'ignore'

    publishDir params.output_dir, mode: 'copy', overwrite: true

    container 'harmonizer:latest'
    containerOptions "-v ${params.dependencies_dir}:/data -v ${params.genefile_dir}:/input"

    input:
        path microarray_file
    
    output:
        // 1) The parquet mapping file (named channel: parquet)
        path 'harmonization_results/vcf_reference_mapping.parquet', emit: parquet

        // 2) All compressed VCFs from bed_files (named channel: vcfs)
        path 'bed_files/*.vcf.gz', emit: vcfs

    script:
    """
    echo "===== DEBUG START ====="
    echo "PWD:"
    pwd

    echo "Listing here:"
    ls -lah

    echo "Listing /app:"
    ls -lah /app || true

    echo "Listing /input:"
    ls -lah /input || true

    echo "Listing /data:"
    ls -lah /data || true

    echo "===== RUNNING PYTHON ====="

    python /app/main.py \
        --microarray_file /input/"${microarray_file.getName()}" \
        --working_dir "\$PWD"

    echo "===== AFTER PYTHON ====="
    echo "Listing harmonization_results:"
    ls -lah harmonization_results || true

    echo "Listing bed_files:"
    ls -lah bed_files || true
    """
}




process impute {


    // Imputed files will end up under:
    //   /mnt/c/Users/frezz/Desktop/docker_testing/nf_output/impute/...
    publishDir "${params.output_dir}/impute", mode: 'copy', overwrite: true

    container 'imputer:latest'

    // Mount deps and nf_output as /work (like your manual docker run)
    containerOptions "-v ${params.dependencies_dir}:/data -v ${params.output_dir}:/work"

    input:
        // We don't actually need the parquet path in the command,
        // but we use it as a dependency so impute runs *after* harmonizer.
        path vcf_mapping_parquet

    output:
        // Adjust this pattern to match what imputer actually writes.
        // For now we assume it produces gzipped VCFs under /work/impute_results or similar.
        path 'imputed/*.vcf.gz', optional: true

    script:
    """
    echo "===== IMPUTE START ====="
    echo "PWD:"
    pwd

    echo "Listing /data:"
    ls -lah /data || true

    echo "Listing /work:"
    ls -lah /work || true

    echo "vcf_mapping_parquet seen by Nextflow as: ${vcf_mapping_parquet}"

    echo "===== RUNNING IMPUTER ====="

    # This mirrors: python main.py
    # but we use /app/main.py so it doesn't depend on the current working dir.
    python /app/main.py

    echo "===== AFTER IMPUTER ====="
    echo "Listing /work:"
    ls -lah /work || true

    echo "Listing possible impute output dir:"
    ls -lah /work/impute_results || true
    """
}


workflow {

    microarray_ch = Channel.fromPath("${params.genefile_dir}/*")

    harmonizer_out = harmonizer(microarray_ch)

    impute( harmonizer_out.parquet )
}
