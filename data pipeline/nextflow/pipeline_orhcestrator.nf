params.harmonizer_dependencies = "/srv/dependencies/imputation_runner/harmonizer"
params.imputation_dependencies = "/srv/dependencies/imputation_runner/imputer"
params.ancestry_dependencies = "/srv/dependencies/imputation_runner/imputer/beagle_references"
params.standardizer_output_dir = "/home/frederik/github_projects/SNPster/data pipeline/temp_test/std"
params.harmonizer_output_dir = "/home/frederik/github_projects/SNPster/data pipeline/temp_test/harm"
params.samplesheet = "/home/frederik/github_projects/SNPster/data pipeline/imputation_runner_module/samplesheet.csv" 


process STANDARDIZE {

    errorStrategy 'ignore'

    container 'standardizer:latest'
    containerOptions "-v ${params.harmonizer_dependencies}:/data"


    input:
            tuple val(identifier), val(file_id), val(output_dir), path(microarray_file)

    output:
                tuple val(identifier), val(output_dir), path("*.parquet"), emit: standardized

    script:
    """

    python /app/main.py --microarray_file "${microarray_file}" --imputation_id "${identifier}" --file_id "${file_id}" --output_dir .

    """
}


process COMBINE_MICROARRAY {

    maxRetries 2
    errorStrategy 'ignore'

    container 'file_combiner:latest'

    input:
    tuple val(identifier), val(output_dir), path(parquet_files)

    output:
    path("combined_output/*.parquet"), emit: parquets

    script:
    """

    mkdir -p combined_output
    python /app/main.py \
        --parquet_files ${parquet_files.collect { "\"${it}\"" }.join(' ')} \
        --imputation_id "${identifier}" \
        --output_dir combined_output

    """
}




process HARMONIZE {

    maxRetries 5
    errorStrategy 'retry'


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

    maxRetries 2
    errorStrategy 'retry'

    container 'vcf_merger:latest'

    input:
    tuple val(chr), path(vcf_files)

    output:
    path "merged_output/*.vcf.gz"

    
    script:
    """

        mkdir -p merged_output
        echo "full_vcf_path,chrom" > samplesheet.csv

        for f in ${vcf_files.join(' ')}; do
            echo "\$PWD/\$f,${chr}" >> samplesheet.csv
        done

        python /app/vcf_combiner.py \
            --vcf_samplesheet_path samplesheet.csv \
            --output_dir merged_output

    """
}



process CALCULATE_ANCESTRY {

    maxRetries 2
    errorStrategy 'retry'
    
    cpus 4

    container 'ancestry:latest'
    containerOptions "-v ${params.ancestry_dependencies}:/data --network docker_default --cpus=4 -e DB_HOST=postgres -e DB_PORT=5432 -e DB_NAME=snpster_db -e DB_USER=postgres -e DB_PASSWORD=zod50902"

    input:
    path merged_chr22_vcf

    output:
    path "ancestry_output/ancestry_success.txt", emit: success_marker


    script:
    """
    mkdir -p ancestry_output

    python /app/main.py \
        --vcf_file ${merged_chr22_vcf}

    # Only created if ancestry completed successfully.
    echo "ok" > ancestry_output/ancestry_success.txt
    """
}




process IMPUTE {

    maxRetries 2
    errorStrategy 'retry'

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

    maxRetries 2
    errorStrategy 'retry'

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


process QC_VCF {

    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    time '15m'
    container 'imputation_qc:latest'


    publishDir path: { output_dir }, mode: 'copy', overwrite: true

    input:
    tuple val(identifier), val(output_dir), path(split_vcf)

    output:
    path "qc_output/*"

    script:
    """

    # Try default engine first for speed; fall back to pandas on native crashes.
    python /app/main.py --input_file ${split_vcf} --output_dir qc_output || \
    ENGINE=pandas python /app/main.py --input_file ${split_vcf} --output_dir qc_output

    """
}




workflow {

    samples_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.identifier,
                row.file_id,
                row.output_dir,
                file(row.file_path)
            )
        }

    standardized_ch = STANDARDIZE(samples_ch)

    standardized_records_ch = standardized_ch.standardized
        .map { identifier, output_dir, parquet_file ->
            tuple(identifier.toString(), output_dir.toString(), parquet_file)
        }

    combiner_input_ch = standardized_records_ch
        .groupTuple()
        .map { grouped ->
            def identifier = grouped[0]
            def output_dirs = grouped[1]
            def parquet_files = grouped[2]
            def unique_output_dirs = output_dirs.unique()
            if (unique_output_dirs.size() != 1) {
                throw new IllegalArgumentException(
                    "Identifier ${identifier} has multiple output directories: ${unique_output_dirs}"
                )
            }
            tuple(identifier, unique_output_dirs[0], parquet_files)
        }

    id_output_dir_ch = combiner_input_ch
        .map { identifier, output_dir, parquet_files ->
            tuple(identifier.toString(), output_dir)
        }

    combined_ch = COMBINE_MICROARRAY(combiner_input_ch)
    harmonized_ch = HARMONIZE(combined_ch.parquets.flatten())

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

    merged_flat_ch = merged_ch.flatten()

    chr22_merged_ch = merged_flat_ch
        .filter { vcf_file -> vcf_file.name == "chr22.merged.vcf.gz" }

    CALCULATE_ANCESTRY(chr22_merged_ch)

    imputed_ch = IMPUTE(merged_flat_ch)

    splitter_channel = SPLIT_VCF(imputed_ch)

    split_with_id_ch = splitter_channel
        .flatten()
        .map { split_vcf ->
            def m = split_vcf.name =~ /(?i)^IMPID(\d+)(?:_chr[^\.]+|\.chr[^\.]+)\.split\.vcf\.gz$/
            if (!m) {
                throw new IllegalArgumentException("Could not parse imputation ID from split VCF: ${split_vcf.name}")
            }

            tuple(m[0][1], split_vcf)
        }

    qc_input_ch = split_with_id_ch
        .combine(id_output_dir_ch, by: 0)
        .map { identifier, split_vcf, output_dir ->
            tuple(identifier, output_dir, split_vcf)
        }

    QC_VCF(qc_input_ch)


}
