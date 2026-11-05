process BASE_PANEL {
    tag "$panel/chr$chrom: matched complete-case panel + disjoint phase-only"
    publishDir path: { "${params.output}/panels/$panel/input/chr$chrom" }, mode: 'copy'
    input:
    tuple val(panel), val(chrom), path(projection)
    path cohorts
    path inventory
    output:
    tuple val(panel), val(chrom), path('base')
    script:
    """
    python3 '${projectDir}/master_resources.py' base --inventory '$inventory' --cohort '$cohorts/$panel' \
        --projection '$projection' --chromosome $chrom --output base
    """
}

process SIMULATE_PANEL {
    tag "$panel/chr$chrom: preserve G0 and simulate G1-G3"
    publishDir path: { "${params.output}/panels/$panel/genotypes/all/chr$chrom" }, mode: 'copy'
    input:
    tuple val(panel), val(chrom), path(base)
    path cohorts
    output:
    tuple val(panel), val(chrom), path('master')
    script:
    """
    python3 '${projectDir}/master_resources.py' simulate --cohort '$cohorts/$panel' --base '$base' \
        --chromosome $chrom --output master
    """
}

process EXPORT_PANEL {
    tag "$panel/chr$chrom: training/reference and harmonized query interfaces"
    publishDir path: { "${params.output}/panels/$panel/genotypes/pipeline/chr$chrom" }, mode: 'copy'
    input:
    tuple val(panel), val(chrom), path(master)
    path cohorts
    output:
    tuple val(panel), val(chrom), path('bundle')
    script:
    """
    python3 '${projectDir}/master_resources.py' export --cohort '$cohorts/$panel' --master '$master' \
        --chromosome $chrom --output bundle
    """
}

process AUDIT_PANEL {
    tag "$panel/chr$chrom: original GT, exact founder homologs, ancestry and exports"
    publishDir path: { "${params.output}/panels/$panel/verification/chromosomes" }, mode: 'copy'
    input:
    tuple val(panel), val(chrom), path(projection), path(base), path(master), path(bundle)
    path cohorts
    output:
    tuple val(panel), val(chrom), path("chr${chrom}.json")
    path "chr${chrom}.source.json"
    script:
    """
    python3 '${projectDir}/audit_resources.py' panel --projection '$projection' --base '$base' \
        --master '$master' --bundle '$bundle' --cohort '$cohorts/$panel' --chromosome $chrom --output chr${chrom}.json
    """
}
