params.harmonizer_dependencies = "/srv/dependencies/imputation_runner/harmonizer"
params.imputation_dependencies = "/srv/dependencies/imputation_runner/imputer"
params.standardizer_output_dir = "/home/frederik/github_projects/SNPster/data pipeline/temp_test/std"
params.harmonizer_output_dir = "/home/frederik/github_projects/SNPster/data pipeline/temp_test/harm"
params.samplesheet = "/home/frederik/github_projects/SNPster/data pipeline/imputation_runner_module/samplesheet.csv" 


process STANDARDIZE {

    container 'standardizer:latest'


    input:
      tuple val(identifier), val(output_dir), path(microarray_file)

    output:
      tuple val(identifier), val(output_dir), path("*.parquet"), emit: parquet_file

    script:
    """
    python /app/main.py --microarray_file "${microarray_file}" --identifier "${identifier}" --output_dir .
    """
}




process HARMONIZE {
    // harmonize all file in dir, the files are outputted to the dir from the standardizer process, so just use the same output dir for the input of the harmonizer process
    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"

    input:
    tuple val(identifier), val(output_dir), path(parquet_file)

    output:
    tuple val(identifier), val(output_dir), path("output/bed_files/*.vcf.gz"), emit: vcfs

    script:
    """
    mkdir -p output
    python /app/main.py --microarray_file "${parquet_file}" --output_dir output
    """
}


process MAKE_SAMPLESHEETS {
    container 'samplesheet_maker:latest'

    input:
    path "input_vcfs/*"

    output:
    path "*.csv"

    script:
    """
    python /app/create_merging_samplesheets.py --vcf_file_dir input_vcfs --output_dir .
    """
}



process MERGE_VCFS {

    container 'vcf_merger:latest'

    input:
    path samplesheet

    output:
    path "merged_output/*.vcf.gz"

    script:
    """
    mkdir -p merged_output
    python /app/vcf_combiner.py --vcf_samplesheet_path ${samplesheet} --output_dir merged_output
    """
}


process IMPUTE {

    maxForks 11

    container 'imputer:latest'
    containerOptions "-v ${params.imputation_dependencies}:/data"

    input:
    path vcf_file
        
    output:
    path "imputed_output/*.vcf.gz"

    script:
    """
    mkdir -p imputed_output
    python /app/main.py --vcf_file ${vcf_file} --output_dir imputed_output
    """
}


process SPLIT_VCF {

    container 'vcf_splitter:latest'

    input:
    path imputed_vcf

    output:
    path "split_output/*.vcf.gz"

    script:
    """
    mkdir -p split_output
    python /app/vcf_splitter.py --vcf_file ${imputed_vcf} --output_dir split_output
    """
}


process QC_VCFS {

    container 'imputation_qc:latest'

    maxForks 10

    input:
    path split_vcf

    output:
    path "qc_output/*"

    script:
    """
    mkdir -p qc_output
    python /app/main.py --input_file ${split_vcf} --output_dir qc_output
    """
}




workflow {

    // Step 1: Load samples from samplesheet
    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.identifier,
                row.output_dir,
                file(row.file_path)
            )
        }

    // Step 2: Standardize and harmonize each sample
    standardized_ch = STANDARDIZE(samples_ch)
    harmonized_ch = HARMONIZE(standardized_ch.parquet_file)
    
    // Step 3: Collect all harmonized VCF files and create samplesheets
    // Extract just the VCF files from the tuple, flatten, and collect them all
    all_vcfs = harmonized_ch.vcfs
        .map { identifier, output_dir, vcf_files -> vcf_files }
        .flatten()
        .collect()
    
    // Create samplesheets from all collected VCFs
    samplesheets_ch = MAKE_SAMPLESHEETS(all_vcfs)
    
    // Step 4: Flatten samplesheets channel to process each CSV individually
    samplesheets_ch = samplesheets_ch.flatten()
    
    // Step 5: Combine each samplesheet with all VCF files, then merge
    // Use all_vcfs.first() to convert value channel to single-item queue for combine
    merged_input_ch = all_vcfs.combine(samplesheets_ch)
    merged_vcfs_ch = MERGE_VCFS(merged_input_ch)
    
    // Step 6: Flatten merged VCFs and impute each one
    merged_vcfs_flat = merged_vcfs_ch.flatten()
    imputed_vcfs_ch = IMPUTE(merged_vcfs_flat)
    
    // Step 7: Flatten imputed VCFs and split each one
    imputed_vcfs_flat = imputed_vcfs_ch.flatten()
    split_vcfs_ch = SPLIT_VCF(imputed_vcfs_flat)
    
    // Step 8: Flatten split VCFs and QC each one
    split_vcfs_flat = split_vcfs_ch.flatten()
    QC_VCFS(split_vcfs_flat)

}
