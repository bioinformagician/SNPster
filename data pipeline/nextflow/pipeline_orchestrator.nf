params.harmonizer_dependencies = "/srv/dependencies/imputation_runner/harmonizer"
params.imputation_dependencies = "/srv/dependencies/imputation_runner/imputer"
params.samplesheet = "/app/imputation_runner_module/input_files.csv" 



process STANDARDIZE {

    container 'standardizer:latest'
    

    input:
      tuple val(identifier), val(output_dir), path(microarray_file)

    output:
      tuple val(identifier), val(output_dir), path("*.parquet"), emit: parquet_file

    script:
    """
    python /app/main.py --microarray_file "${microarray_file}" --identifier "${identifier}"
    """
}




process HARMONIZE {

    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"

    input:
      tuple val(identifier), val(output_dir), path(parquet_file)

    output:
      tuple val(identifier), val(output_dir), path("bed_files/*.vcf.gz"), emit: vcfs

    script:
    """
    python /app/main.py --microarray_file "${parquet_file}"
    """
}




process IMPUTE {

  publishDir path: { output_dir }, mode: 'copy'

  container 'imputer:latest'
  containerOptions "-v ${params.imputation_dependencies}:/data -e HEAP_GB=8"

  input:
    tuple val(identifier), val(output_dir), path(vcf_files)

  output:
    path "output/*"

  script:
  """
  python /app/main.py --vcf_files .

  mkdir -p output
  cp -a /output/* output/
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
    IMPUTE(harmonized_ch.vcfs)

}
