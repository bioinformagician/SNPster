params.harmonizer_dependencies = "/srv/dependencies/imputation_runner/harmonizer"
params.imputation_dependencies = "/srv/dependencies/imputation_runner/imputer"
params.standardizer_output_dir = "/path/to/standardizer/output/dir"
params.harmonizer_output_dir = "/path/to/harmonizer/output/dir"
params.imputation_output_dir = "/path/to/output/dir"
params.samplesheet = "/app/imputation_runner_module/input_files.csv" 


process STANDARDIZE {

    container 'standardizer:latest'


    input:
      tuple val(identifier), val(output_dir), path(microarray_file)

    output:
      tuple val(identifier), val(output_dir), path("*.parquet"), emit: parquet_file

    script:
    """
    python /app/main.py --microarray_file "${microarray_file}" --identifier "${identifier}" --output_dir "${params.standardizer_output_dir}"
    """
}




process HARMONIZE {
    // harmonize all file in dir, the files are outputted to the dir from the standardizer process, so just use the same output dir for the input of the harmonizer process
    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"

    input:
    tuple val(identifier), val(output_dir), path(parquet_file)

    output:
    tuple val(identifier), val(output_dir), path("output/*.vcf.gz"), emit: vcfs

    script:
    """
    python /app/main.py --microarray_file "${parquet_file}" --output_dir "${params.harmonizer_output_dir}"

    """
}


process MERGE_VCFS {


    input:
    path file

    script:
    """
    python /app/vcf_combiner.py --vcf_samplesheet_path $file --output_dir $params.output_dir

    mkdir -p merge_output
    cp -a /output/* merge_output/

    """
}


process IMPUTE {

    maxForks 11

    container 'imputer:latest'
    containerOptions "-v ${params.imputation_dependencies}:/data"

    input:
        path file
        
    output:
        path "output/*"

    script:
    """
    python /app/main.py --vcf_file $file --output_dir $params.imputation_output_dir

    mkdir -p imputation_output
    cp -a /output/* imputation_output/
        

    """
}


process SPLIT_VCF {

  input:
    path file

  script:
  """
  python /app/vcf_splitter.py --vcf_samplesheet_path $file --output_dir $params.output_dir
  """
}


process QC_VCFS {

    container 'imputation_qc:latest'

    maxForks = 10

    input:
        path file


    script:
    """
    python /app/main.py --input_file $file --output_dir $params.output_dir
    """
}




workflow {


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

    standardized_ch = STANDARDIZE(samples_ch)
    harmonized_ch = HARMONIZE(standardized_ch.parquet_file)
    
    // Collect all VCF files from all harmonizer processes into a single list
    all_vcfs = harmonized_ch.vcfs
        .map { identifier, output_dir, vcf_files -> vcf_files }
        .flatten()
        .collect()
        .map { vcf_list -> tuple("combined", "./output", vcf_list) }
    
    IMPUTE(all_vcfs)

}
