#!/usr/bin/env nextflow

params.harmonizer_dependencies = "/mnt/c/Users/frezz/Desktop/harmonizer_dependencies"
params.imputation_dependencies = "/mnt/c/Users/frezz/Desktop/imputer_dependencies"
params.genefile_dir     = "/mnt/c/Users/frezz/Desktop/microarray_data/testing"
params.output_dir       = "/mnt/c/Users/frezz/Desktop/docker_testing/nf_output"

process HARMONIZE {

    //errorStrategy 'ignore'

    publishDir params.output_dir, mode: 'copy', overwrite: true

    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data -v ${params.genefile_dir}:/input"

    input:
        path microarray_file
    
    output:
        // 1) The files needed to be parsed to the imputer (named channel: vcfs)
        path "/${PWD}/bed_files/*.vcf.gz", emit: vcfs

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

    python /app/main.py --microarray_file /input/"${microarray_file.getName()}" --working_dir "/$PWD"

    echo "===== AFTER PYTHON ====="
    echo "Listing harmonization_results:"
    ls -lah harmonization_results || true

    echo "Listing bed_files:"
    ls -lah bed_files || true
    """
}




process IMPUTE {


    publishDir "${params.output_dir}/impute", mode: 'copy', overwrite: true

    container 'imputer:latest'

    // Mount deps and nf_output as /work (like your manual docker run)
    containerOptions "-v ${params.imputation_dependencies}:/data -v ${params.output_dir}:/work"

    input:
        path vcf_files

    output:

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

    HARMONIZE(microarray_ch)

    IMPUTE(HARMONIZE.out.vcfs)
}
