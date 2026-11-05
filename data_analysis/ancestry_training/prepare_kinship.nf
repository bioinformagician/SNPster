nextflow.enable.dsl=2

process CATALOG_PEDIGREES {
    publishDir "${params.output}/catalog", mode: 'copy'
    input:
    val parent
    path pedigree
    path authority
    output:
    path 'catalog'
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' catalog --parent '$parent' \
        --pedigree '$pedigree' --authority '$authority' --output catalog
    """
}

process SMOKE_REAL_CHROMOSOME {
    tag 'chr22 complete genotype + PLINK/LD smoke'
    publishDir "${params.output}/genotypes/all/chr22", mode: 'copy'
    input:
    path catalog
    val dense
    output:
    tuple val(22), path('master')
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' chromosome --catalog '$catalog' \
        --chromosome 22 --dense-root '$dense' --output master
    """
}

process ASSEMBLE_OBSERVED_CHROMOSOME {
    tag "chr$chrom"
    publishDir path: { "${params.output}/genotypes/all/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(catalog)
    val dense
    output:
    tuple val(chrom), path('master')
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' chromosome --catalog '$catalog' \
        --chromosome $chrom --dense-root '$dense' --output master
    """
}

process GENOME_WIDE_KINSHIP {
    publishDir "${params.output}/kinship", mode: 'copy'
    input:
    path bundles, stageAs: 'bundles/chr??'
    output:
    path 'screen'
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' kinship --bundles bundles/chr* --output screen
    """
}

process ASSIGN_FAMILY_COMPONENTS {
    publishDir "${params.output}/metadata", mode: 'copy'
    input:
    path catalog
    path screen
    output:
    path 'split'
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' split --catalog '$catalog' --screen '$screen' --output split
    """
}

process RELEASE_PIPELINE_CHROMOSOME {
    tag "chr$chrom"
    publishDir path: { "${params.output}/genotypes/pipeline/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(master)
    path split
    output:
    tuple val(chrom), path('bundle')
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' release --master '$master' \
        --split '$split' --chromosome $chrom --output bundle
    """
}

process AUDIT_RELEASE_CHROMOSOME {
    tag "chr$chrom exact GT, phase, sample and mask axes"
    publishDir "${params.output}/verification/chromosomes", mode: 'copy'
    input:
    tuple val(chrom), path(master), path(bundle)
    path split
    output:
    path "chr${chrom}.json"
    script:
    """
    python3 '${projectDir}/audit_kinship.py' chromosome --master '$master' --bundle '$bundle' \
        --samples '$split/samples.tsv' --chromosome $chrom --output chr${chrom}.json
    """
}

process FINALIZE_REUSABLE_DATASET {
    publishDir "${params.output}", mode: 'copy'
    input:
    path catalog
    path split
    path bundles, stageAs: 'bundles/chr??'
    path audits
    output:
    path 'configs'
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' finalize --catalog '$catalog' --split '$split' \
        --bundles bundles/chr* --destination '${params.output}' --output configs
    """
}

workflow BUILD_DATASET {
    catalog = CATALOG_PEDIGREES(params.parent, file(params.pedigree, checkIfExists: true),
                               file(params.authority, checkIfExists: true))
    smoke = SMOKE_REAL_CHROMOSOME(catalog, params.dense)
    // Gate all other autosomes on a complete chr22 real-container assembly and PLINK run.
    jobs = catalog.combine(smoke.map { chrom, directory -> true }).flatMap { directory, passed ->
        (1..21).collect { chrom -> tuple(chrom, directory) }
    }

    remaining = ASSEMBLE_OBSERVED_CHROMOSOME(jobs, params.dense)
    masters = smoke.mix(remaining)
    screen = GENOME_WIDE_KINSHIP(masters.map { chrom, directory -> directory }.collect())
    split = ASSIGN_FAMILY_COMPONENTS(catalog, screen)
    released = RELEASE_PIPELINE_CHROMOSOME(masters, split)
    audits = AUDIT_RELEASE_CHROMOSOME(masters.join(released), split)
    FINALIZE_REUSABLE_DATASET(catalog, split,
        released.map { chrom, directory -> directory }.collect(), audits.collect())
}

workflow RECOVER_RELEASE {
    catalog = Channel.value(file("${params.output}/catalog/catalog", checkIfExists: true))
    split = Channel.value(file("${params.output}/metadata/split", checkIfExists: true))
    bundles = Channel.value((1..22).collect { chrom ->
        file("${params.output}/genotypes/pipeline/chr$chrom/bundle", checkIfExists: true)
    })
    audits = Channel.value((1..22).collect { chrom ->
        file("${params.output}/verification/chromosomes/chr${chrom}.json", checkIfExists: true)
    })
    FINALIZE_REUSABLE_DATASET(catalog, split, bundles, audits)
}

workflow {
    if (params.recover_release) {
        RECOVER_RELEASE()
    } else {
        BUILD_DATASET()
    }
}
