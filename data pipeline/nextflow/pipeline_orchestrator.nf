#!/usr/bin/env nextflow

params.harmonizer_dependencies = "/mnt/c/Users/frezz/Desktop/harmonizer_dependencies"
params.imputation_dependencies = "/mnt/c/Users/frezz/Desktop/imputer_dependencies"
params.genefile_dir     = "/mnt/c/Users/frezz/Desktop/genome_scraping/scraped_genomes/test"
params.output_dir       = "/mnt/c/Users/frezz/Desktop/docker_testing/nf_output"




process STANDARDIZE {

    //errorStrategy 'ignore'
        
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


    //publishDir params.output_dir, mode: 'copy', overwrite: true only for testing

    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"
    //containerOptions "-v ${params.genefile_dir}:/input"

    input:
        path parquet_file
    
    output:
        // 1) The files needed to be parsed to the imputer (named channel: vcfs)
        path "bed_files/*.vcf.gz", emit: vcfs

    script:
    """

    python /app/main.py --microarray_file "${parquet_file}"

    """
}




process IMPUTE {

    publishDir "${params.output_dir}", mode: 'copy'

    container 'imputer:latest'

    
    //for when mounting huge dependencies as volumes
    //containerOptions "-v ${params.imputation_dependencies}:/data -v ${params.output_dir}:/work --memory=12g"

    // when the dependencies are baked into the image
    //containerOptions "-v ${params.output_dir}:/work --memory=12g"
    containerOptions "-v${params.imputation_dependencies}:/data"


    input:
        path vcf_files


    output:

        path 'imputed/*.vcf.gz', emit: imputed_vcfs, optional: true
        path 'output/*.vcf*', emit: output_vcfs, optional: true  
        path 'output/*.csv', emit: samplesheet

    script:
    """

    python /app/main.py --vcf_files .

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
