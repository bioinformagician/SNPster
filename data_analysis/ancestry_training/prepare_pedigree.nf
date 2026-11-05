nextflow.enable.dsl=2

process PLAN_DISJOINT_PEDIGREES {
    publishDir "${params.output}/plan", mode: 'copy'
    input:
    val parent
    path design
    output:
    path 'plan'
    script:
    """
    python3 '${projectDir}/prepare_pedigree.py' plan --parent '$parent' --design '$design' --output plan
    """
}

process SIMULATE_PEDIGREE_CHROMOSOME {
    tag "chr$chrom: founder haplotypes + three generations"
    publishDir path: { "${params.output}/genotypes/all/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(plan)
    output:
    tuple val(chrom), path('master')
    script:
    """
    python3 '${projectDir}/prepare_pedigree.py' simulate --plan '$plan' --chromosome $chrom --output master
    """
}

process EXPORT_PEDIGREE_CHROMOSOME {
    tag "chr$chrom: model-ready phased training + harmonized masks"
    publishDir path: { "${params.output}/genotypes/pipeline/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(master)
    path plan
    output:
    tuple val(chrom), path('bundle')
    script:
    """
    python3 '${projectDir}/prepare_kinship.py' release --master '$master' --split '$plan' \
        --chromosome $chrom --output bundle
    """
}

process AUDIT_PEDIGREE_CHROMOSOME {
    tag "chr$chrom: independent exact founder inheritance, truth + model exports"
    publishDir "${params.output}/verification/chromosomes", mode: 'copy'
    input:
    tuple val(chrom), path(master), path(bundle)
    path plan
    output:
    path "chr${chrom}.json"
    script:
    """
    python3 '${projectDir}/audit_pedigree.py' chromosome --master '$master' --bundle '$bundle' \
        --plan '$plan' --chromosome $chrom --output chr${chrom}.json
    """
}

process FINALIZE_PEDIGREE_DATASET {
    publishDir "${params.output}", mode: 'copy'
    input:
    path plan
    path masters, stageAs: 'masters/chr??'
    path audits
    output:
    path 'metadata'
    path 'configs'
    script:
    """
    python3 '${projectDir}/prepare_pedigree.py' finalize --plan '$plan' --bundles masters/chr* \
        --destination '${params.output}' --output ready
    mv ready/metadata metadata
    mv ready/configs configs
    """
}

workflow {
    smoke_only = params.smoke_only.toString().toBoolean()
    plan = PLAN_DISJOINT_PEDIGREES(params.parent, file(params.design, checkIfExists: true))
    chromosomes = smoke_only ? [22] : (1..22).toList()
    // A separate immutable chr22 smoke run must pass the complete inheritance/export audit first.
    jobs = plan.flatMap { directory -> chromosomes.collect { chrom -> tuple(chrom, directory) } }
    masters = SIMULATE_PEDIGREE_CHROMOSOME(jobs)
    exports = EXPORT_PEDIGREE_CHROMOSOME(masters, plan)
    audits = AUDIT_PEDIGREE_CHROMOSOME(masters.join(exports), plan)
    if (!smoke_only) {
        FINALIZE_PEDIGREE_DATASET(plan, masters.map { chrom, directory -> directory }.collect(), audits.collect())
    }
}
