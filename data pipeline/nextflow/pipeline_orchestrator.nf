#!/usr/bin/env nextflow

params.harmonizer_dependencies = "/mnt/c/Users/frezz/Desktop/harmonizer_dependencies"
params.imputation_dependencies = "/mnt/c/Users/frezz/Desktop/imputer_dependencies"
params.genefile_dir     = "/mnt/c/Users/frezz/Desktop/microarray_data/testing"
params.output_dir       = "/mnt/c/Users/frezz/Desktop/docker_testing/nf_output"




process STANDARDIZE {

    errorStrategy 'ignore'
        
    container 'standardizer:latest'
    input:
        path microarray_file
    output:
        path "*.parquet", emit: parquet_file
    script:
    """
    python /app/main.py --microarray_file "${microarray_file}"
    """
}




process HARMONIZE {


    publishDir params.output_dir, mode: 'copy', overwrite: true

    container 'harmonizer-full:latest'
    //containerOptions "-v ${params.harmonizer_dependencies}:/data -v ${params.genefile_dir}:/input"
    //containerOptions "-v ${params.genefile_dir}:/input"

    input:
        path parquet_file
    
    output:
        // 1) The files needed to be parsed to the imputer (named channel: vcfs)
        path "bed_files/*.vcf.gz", emit: vcfs

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

    python /app/main.py --microarray_file "${parquet_file}"

    echo "===== AFTER PYTHON ====="
    echo "Listing harmonization_results:"
    ls -lah harmonization_results || true

    echo "Listing bed_files:"
    ls -lah bed_files || true
    """
}




process IMPUTE {

    debug true  // Print stdout/stderr in real-time

    publishDir "${params.output_dir}/impute", mode: 'copy', overwrite: true

    container 'imputer-full:latest'

    
    //for when mounting huge dependencies as volumes
    //containerOptions "-v ${params.imputation_dependencies}:/data -v ${params.output_dir}:/work --memory=12g"

    // when the dependencies are baked into the image
    //containerOptions "-v ${params.output_dir}:/work --memory=12g"
    containerOptions "--memory=10g"


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

    echo "===== RUNNING IMPUTER ====="

    # VCF files are staged in PWD by Nextflow
    # Pass current directory as the vcf_files location
    python /app/main.py --vcf_files .

    echo "===== AFTER IMPUTER ====="
    echo "Listing /work:"
    ls -lah /work || true

    echo "Listing possible impute output dir:"
    ls -lah /work/impute_results || true
    """
}


workflow {

    microarray_ch = Channel.fromPath("${params.genefile_dir}/*")
    standardized_ch = STANDARDIZE(microarray_ch)
    harmonized_ch = HARMONIZE(standardized_ch.parquet_file)
    IMPUTE(harmonized_ch.vcfs)

    /*STANDARDIZE(microarray_ch)
    HARMONIZE(STANDARDIZE.out.parquet_file)
    IMPUTE(HARMONIZE.out.vcfs)*/
}
