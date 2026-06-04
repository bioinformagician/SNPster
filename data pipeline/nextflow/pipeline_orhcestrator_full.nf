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
        path("*.parquet"), emit: parquets

    script:
    """
    python /app/main.py --microarray_file "${microarray_file}" --imputation_id "${identifier}" --output_dir .
    """
}




process HARMONIZE {
    // harmonize all file in dir, the files are outputted to the dir from the standardizer process, so just use the same output dir for the input of the harmonizer process
    container 'harmonizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"

    input:
    path(parquet_file)

    output:
    path("output/*.vcf.gz"), emit: vcfs

    script:
    """
    mkdir -p output
    python /app/main.py --microarray_file "${parquet_file}" --output_dir output
    """
}


process MERGE_VCFS {

    container 'vcf_merger:latest'

    input:
    tuple val(chr), path(vcf_files)

    output:
    path "merged_output/*.vcf.gz"

    
    script:
    """
        echo "full_vcf_path,chrom" > samplesheet.csv

        for f in ${vcf_files.join(' ')}; do
            echo "\$PWD/\$f,${chr}" >> samplesheet.csv
        done

        python /app/vcf_combiner.py \
            --vcf_samplesheet_path samplesheet.csv
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
    harmonized_ch = HARMONIZE(standardized_ch.parquets.flatten())

    all_vcfs = harmonized_ch.vcfs

    vcf_bundles_ch = all_vcfs
        .flatten()
        .map { vcf_path ->
            def m = vcf_path.name =~ /IMPID\d+\.chr([^\.]+)\..*\.vcf\.gz$/
            if (!m) {
                throw new IllegalArgumentException("Could not parse chromosome from: ${vcf_path.name}")
            }

            def chr = "chr${m[0][1]}"
            tuple(chr, vcf_path)
        }
        .groupTuple()

    merged_ch = MERGE_VCFS(vcf_bundles_ch)

    imputed_ch = IMPUTE(merged_ch.flatten())


}
