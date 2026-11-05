nextflow.enable.dsl = 2
params.harmonizer_options = ''
params.preprocessing_script = "${params.project}/pedigree_benchmark/scripts/production_matched_preprocessing.py"
params.imputation_root = '/srv/dependencies/imputation_runner/imputer'

process PREPARE_REAL_RSID_TRUTH_PANEL {
    cpus 1
    memory '10 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/preparation", mode: 'copy'
    input:
    val chromosomes
    output:
    path 'rsids'
    script:
    """
    python3 '${params.preprocessing_script}' rsids \
      --cohort '${params.cohort}' --chromosomes ${chromosomes.join(' ')} --output rsids
    """
}

process PREPARE_OBSERVED_FOUNDERS {
    tag "chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/founders/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), path(markers)
    output:
    tuple val(chromosome), path('founders')
    script:
    """
    python3 '${params.preprocessing_script}' founders \
      --cohort '${params.cohort}' --chromosome ${chromosome} --markers '${markers}' --output founders
    """
}

process PHASE_SIMULATION_FOUNDERS_WITH_PRODUCTION {
    tag "chr${chromosome}"
    memory '8 GB'
    container 'sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'
    containerOptions "${params.container_options} -v ${params.imputation_root}:/data:ro"
    input:
    tuple val(chromosome), path(founders)
    output:
    tuple val(chromosome), path(founders), path('phased')
    script:
    """
    mkdir phased
    python3 /app/main.py --vcf_file '${founders}/founders_chr${chromosome}.vcf.gz' --output_dir phased
    """
}

process INDEX_PHASED_FOUNDERS {
    tag "chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/founder_phase/chr${chromosome}" }, mode: 'copy', pattern: 'phased'
    input:
    tuple val(chromosome), path(founders), path(phased, stageAs: 'phase_input')
    output:
    tuple val(chromosome), path(founders), path('phased')
    script:
    """
    mkdir phased
    cp -a phase_input/. phased/
    bcftools index phased/imputed_chr${chromosome}.risk.vcf.gz
    """
}

process PREPARE_MENDELIAN_ARRAY_CASES {
    tag "chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/truth", mode: 'copy'
    input:
    tuple val(chromosome), path(founders), path(phased)
    output:
    tuple val(chromosome), path("chr${chromosome}")
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/production_matched_preprocessing.py' cases \
      --cohort '${params.cohort}' --chromosome ${chromosome} --folder '${founders}' \
      --phased '${phased}/imputed_chr${chromosome}.risk.vcf.gz' --masks ${params.masks.join(' ')} --output chr${chromosome}
    """
}

process PREPARE_SYNTHETIC_PEDIGREE_ARRAY_CASES {
    tag "chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/truth", mode: 'copy'
    input:
    tuple val(chromosome), path(founders)
    output:
    tuple val(chromosome), path("chr${chromosome}")
    script:
    """
    python3 '${params.preprocessing_script}' cases \
      --cohort '${params.cohort}' --chromosome ${chromosome} --folder '${founders}' \
      --masks ${params.masks.join(' ')} --output chr${chromosome}
    """
}

process ASSEMBLE_RAW_MICROARRAY_UPLOADS {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/raw", mode: 'copy'
    input:
    path folders
    output:
    path 'raw'
    script:
    """
    python3 '${params.preprocessing_script}' assemble \
      --cohort '${params.cohort}' --folders ${folders.join(' ')} --masks ${params.masks.join(' ')} --output raw
    """
}

process PRODUCTION_STANDARDIZE {
    tag "${mask}/${identity}"
    memory '10 GB'
    container 'sha256:85e06b985dbaa682af024e5ad5af3cbbacc6eba24d5721fda47b355f6fa70544'
    containerOptions "${params.container_options} -e PVAR_REF_FILE=/srv/dependencies/imputation_runner/harmonizer/all_hg38.pvar.zst -e HG19_TO_HG38_CHAIN_FILE=/app/dependencies/hg19ToHg38.over.chain.gz -e HG18_TO_HG38_CHAIN_FILE=/app/dependencies/hg18ToHg38.over.chain.gz"
    publishDir path: { "${params.run_dir}/standardized/${mask}/${identity}" }, mode: 'copy'
    input:
    tuple val(mask), val(identity), path(raw)
    output:
    tuple val(mask), val(identity), path('std/*.parquet')
    script:
    """
    mkdir std
    python3 /app/main.py --microarray_file '${raw}' --imputation_id '${identity}' --file_id 1 --output_dir std
    """
}

process PRODUCTION_COMBINE_FILES {
    tag "${mask}/${identity}"
    container 'sha256:8b10363aadff5187ef09eff2b5b5d09b8b850173392dede0eadebf2885575ad4'
    publishDir path: { "${params.run_dir}/combined/${mask}/${identity}" }, mode: 'copy'
    input:
    tuple val(mask), val(identity), path(parquet)
    output:
    tuple val(mask), val(identity), path('combined/*.parquet')
    script:
    """
    mkdir combined
    python3 /app/main.py --parquet_files '${parquet}' --imputation_id '${identity}' --output_dir combined
    """
}

process PRODUCTION_HARMONIZE {
    tag "${mask}/${identity}/chr${chromosome}"
    container 'sha256:cd7907c83ad8d41017ce1e235172ace6c5c00d04565a11f75c0b185260c7a74f'
    containerOptions "${params.container_options} ${params.harmonizer_options} -e PLINK_REFERENCE_FASTA=/srv/dependencies/imputation_runner/harmonizer/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    publishDir path: { "${params.run_dir}/harmonized/${mask}/chr${chromosome}/${identity}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), val(identity), path(parquet)
    output:
    tuple val(chromosome), val(mask), val(identity), path('harm/*.vcf.gz')
    script:
    """
    mkdir harm
    python3 /app/main.py --microarray_file '${parquet}' --output_dir harm
    """
}

process PRODUCTION_MERGE_CHROMOSOMES {
    tag "${mask}/chr${chromosome}"
    container 'sha256:f7f4bd94ac09a56125dcbccaa67e1854b97fbdd3164cc5b103ff642319b7483d'
    publishDir path: { "${params.run_dir}/merged/${mask}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), path(vcfs)
    output:
    tuple val(chromosome), val(mask), path('merged/*.vcf.gz')
    script:
    def records = vcfs.collect { "'${it},chr${chromosome}'" }.join(' ')
    """
    mkdir merged
    printf '%s\\n' 'full_vcf_path,chrom' ${records} > samplesheet.csv
    python3 /app/vcf_combiner.py --vcf_samplesheet_path samplesheet.csv --output_dir merged
    """
}

process PRODUCTION_BEAGLE_IMPUTE {
    tag "${mask}/chr${chromosome}"
    memory '8 GB'
    container 'sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'
    containerOptions "${params.container_options} -v ${params.imputation_root}:/data:ro"
    input:
    tuple val(chromosome), val(mask), path(merged)
    output:
    tuple val(chromosome), val(mask), path(merged), path("imputed_${mask}_chr${chromosome}")
    script:
    """
    mkdir imputed_${mask}_chr${chromosome}
    python3 /app/main.py --vcf_file '${merged}' --output_dir imputed_${mask}_chr${chromosome}
    """
}

process INDEX_PRODUCTION_IMPUTATION {
    tag "${mask}/chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/imputed/${mask}/chr${chromosome}" }, mode: 'copy', pattern: 'imputed_*'
    input:
    tuple val(chromosome), val(mask), path(merged), path(imputed, stageAs: 'imputed_input')
    output:
    tuple val(chromosome), val(mask), path(merged), path("imputed_${mask}_chr${chromosome}")
    script:
    """
    mkdir imputed_${mask}_chr${chromosome}
    cp -a imputed_input/. imputed_${mask}_chr${chromosome}/
    bcftools index imputed_${mask}_chr${chromosome}/imputed_chr${chromosome}.risk.vcf.gz
    """
}

process PRODUCTION_SPLIT_SAMPLES {
    tag "${mask}/chr${chromosome}"
    container 'sha256:3db53a67ae2ce73ddd3fadb485a2d19c7e18dbc58ad6c86cc77ef88f330d1ccd'
    publishDir path: { "${params.run_dir}/split/${mask}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), path(merged), path(imputed)
    output:
    tuple val(chromosome), val(mask), path('split/*.vcf.gz')
    script:
    """
    cp '${imputed}/imputed_chr${chromosome}.risk.vcf.gz' split_source.vcf.gz
    mkdir split
    python3 /app/vcf_splitter.py --vcf_file split_source.vcf.gz --output_dir split
    """
}

process PRODUCTION_QC_SAMPLES {
    tag "${mask}/${identity}/chr${chromosome}"
    container 'sha256:5d3d0fda0b7fd310efe5c112314d3923cdee2bd49d7abe95f5b6de3f895421e8'
    publishDir path: { "${params.run_dir}/qc/${mask}/chr${chromosome}/${identity}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), val(identity), path(vcf)
    output:
    tuple val(chromosome), val(mask), val(identity), path('qc/*.vcf.gz')
    script:
    """
    python3 /app/main.py --input_file '${vcf}' --output_dir qc || \
      ENGINE=pandas python3 /app/main.py --input_file '${vcf}' --output_dir qc
    """
}

process INDEX_QC_OUTPUTS {
    tag "${mask}/${identity}/chr${chromosome}"
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/indexed_qc/${mask}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), val(identity), path(vcf)
    output:
    tuple val(chromosome), val(mask), path("IMPID${identity}.vcf.gz"), path("IMPID${identity}.vcf.gz.csi")
    script:
    """
    cp '${vcf}' IMPID${identity}.vcf.gz
    bcftools index IMPID${identity}.vcf.gz
    """
}

process VERIFY_PRODUCTION_GENOTYPES {
    tag "${mask}/chr${chromosome}"
    memory '6 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/preprocessing/${mask}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(mask), path(merged), path(imputed), path(qc), path(indices), path(truth)
    output:
    tuple val(chromosome), val(mask), path('verified/preprocessing_verified.json')
    script:
    """
    python3 '${params.preprocessing_script}' verify \
      --cohort '${params.cohort}' --chromosome ${chromosome} --mask '${mask}' --folder '${truth}' \
      --merged '${merged}' --imputed '${imputed}/imputed_chr${chromosome}.risk.vcf.gz' --qc ${qc.join(' ')} --output verified
    """
}

process GATE_ALL_PREPROCESSING {
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir "${params.run_dir}/preprocessing_gate", mode: 'copy'
    input:
    path(reports, stageAs: 'reports/report??.json')
    output:
    path 'preprocessing_passed.json'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/verify_production_matched.py' gate \
      --cohort '${params.cohort}' --chromosomes ${params.chromosomes.join(' ')} --masks ${params.masks.join(' ')} \
      --reports ${reports.join(' ')} --output preprocessing_passed.json
    """
}

process PREPARE_ANCESTRY_MARKER_PANELS {
    tag "chr${chromosome}"
    memory '8 GB'
    container 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
    publishDir path: { "${params.run_dir}/panels/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(masks), path(imputed), path(merged, stageAs: 'merged??/*'), path(gate), path(truth)
    output:
    tuple val(chromosome), path('panel')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/production_matched_models.py' panels \
      --cohort '${params.cohort}' --chromosome ${chromosome} --masks ${masks.join(' ')} \
      --folders ${imputed.join(' ')} --merged ${merged.join(' ')} --truth '${truth}' --output panel
    """
}

process RUN_NATIVE_ANCESTRY_MODEL {
    tag "${candidate}/chr${chromosome}"
    memory '8 GB'
    container { candidate.startsWith('gnomix') ? 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21' : 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d' }
    publishDir path: { "${params.run_dir}/models/${candidate}/chr${chromosome}" }, mode: 'copy'
    input:
    tuple val(chromosome), val(candidate), path(panel)
    output:
    tuple val(chromosome), val(candidate), path('model')
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/production_matched_models.py' evaluate \
      --cohort '${params.cohort}' --folder '${panel}' --candidate '${candidate}' --output model
    """
}

process VERIFY_AND_SCORE_COMPLETE_BENCHMARK {
    memory '8 GB'
    container 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21'
    publishDir "${params.run_dir}/summary", mode: 'copy'
    input:
    path(models, stageAs: 'models/model???')
    output:
    path 'summary'
    script:
    """
    python3 '${params.project}/pedigree_benchmark/scripts/verify_production_matched.py' score \
      --cohort '${params.cohort}' --run-dir '${params.run_dir}' --chromosomes ${params.chromosomes.join(' ')} \
      --masks ${params.masks.join(' ')} --models ${models.join(' ')} --output summary
    """
}

workflow {
    def mode = params.smoke.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) { error 'smoke must be explicitly true or false' }
    def cohort = new groovy.json.JsonSlurper().parse(file(params.cohort))
    params.chromosomes = mode == 'true' ? [22] : (1..22).toList()
    params.masks = mode == 'true' ? ['lonni'] : cohort.protocol.masks
    log.info "Production-matched scope: chromosomes=${params.chromosomes}; masks=${params.masks}; configurations=10"
    rsids = PREPARE_REAL_RSID_TRUTH_PANEL(params.chromosomes)
    founders = PREPARE_OBSERVED_FOUNDERS(rsids.flatMap { folder ->
        params.chromosomes.collect { chromosome -> tuple(chromosome, file("${folder}/chr${chromosome}.tsv")) }
    })
    if (cohort.truth_design == 'synthetic_reference_haplotype_development') {
        cases = PREPARE_SYNTHETIC_PEDIGREE_ARRAY_CASES(founders)
    } else {
        phased = INDEX_PHASED_FOUNDERS(PHASE_SIMULATION_FOUNDERS_WITH_PRODUCTION(founders))
        cases = PREPARE_MENDELIAN_ARRAY_CASES(phased)
    }
    raw = ASSEMBLE_RAW_MICROARRAY_UPLOADS(cases.map { chromosome, folder -> folder }.collect())
    standardized = PRODUCTION_STANDARDIZE(raw.flatMap { folder ->
        def jobs = new groovy.json.JsonSlurper().parse(file("${folder}/jobs.json"))
        jobs.collect { row -> tuple(row.mask, row.production_id, file("${folder}/${row.raw}")) }
    })
    combined = PRODUCTION_COMBINE_FILES(standardized)
    harmonization = combined.flatMap { mask, identity, files ->
        (files instanceof List ? files : [files]).collect { path ->
            def match = path.name =~ /\.chr([0-9]+)\./
            if (!match) { error "Invalid production combined filename: ${path}" }
            tuple(match[0][1] as int, mask, identity, path)
        }
    }
    harmonized = PRODUCTION_HARMONIZE(harmonization)
    batches = harmonized.groupTuple(by: [0, 1]).map { chromosome, mask, identities, files ->
        if (identities.toSet().size() != cohort.people.size()) { error 'Incomplete/duplicate production merge batch' }
        tuple(chromosome, mask, files.sort { a, b -> a.name <=> b.name })
    }
    merged_batches = PRODUCTION_MERGE_CHROMOSOMES(batches)
    imputed = INDEX_PRODUCTION_IMPUTATION(PRODUCTION_BEAGLE_IMPUTE(merged_batches))
    split = PRODUCTION_SPLIT_SAMPLES(imputed)
    qc_jobs = split.flatMap { chromosome, mask, files ->
        (files instanceof List ? files : [files]).collect { path ->
            def match = path.name =~ /IMPID([0-9]+)\./
            if (!match) { error "Invalid production split filename: ${path}" }
            tuple(chromosome, mask, match[0][1], path)
        }
    }
    qc = PRODUCTION_QC_SAMPLES(qc_jobs)
    indexed = INDEX_QC_OUTPUTS(qc).groupTuple(by: [0, 1])
    checks = imputed.join(indexed, by: [0, 1]).combine(cases, by: 0)
    verified = VERIFY_PRODUCTION_GENOTYPES(checks)
    gate = GATE_ALL_PREPROCESSING(verified.map { chromosome, mask, report -> report }.collect())
    panel_inputs = imputed.map { chromosome, mask, merged, folder -> tuple(chromosome, mask, folder, merged) }
        .groupTuple(by: 0).combine(gate).combine(cases, by: 0)
    panels = PREPARE_ANCESTRY_MARKER_PANELS(panel_inputs)
    jobs = panels.flatMap { chromosome, folder ->
        cohort.protocol.candidates.keySet().collect { candidate -> tuple(chromosome, candidate, folder) }
    }
    models = RUN_NATIVE_ANCESTRY_MODEL(jobs)
    VERIFY_AND_SCORE_COMPLETE_BENCHMARK(models.map { chromosome, candidate, folder -> folder }.collect())
}
