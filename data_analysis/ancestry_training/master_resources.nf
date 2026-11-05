nextflow.enable.dsl=2

include {
    BASE_PANEL as SMOKE_BASE
    SIMULATE_PANEL as SMOKE_SIMULATE
    EXPORT_PANEL as SMOKE_EXPORT
    AUDIT_PANEL as SMOKE_AUDIT
    BASE_PANEL
    SIMULATE_PANEL
    EXPORT_PANEL
    AUDIT_PANEL
} from './master_panels'

process INVENTORY_ALL_RESOURCES {
    publishDir "${params.output}", mode: 'copy'
    input:
    path specification
    output:
    path 'inventory'
    script:
    """
    python3 '${projectDir}/master_resources.py' inventory --specification '$specification' --output inventory
    """
}

process PROJECT_RESOURCE_CHROMOSOME {
    tag "chr$chrom: all people; original GT and missingness"
    publishDir path: { "${params.output}/archive/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(inventory)
    output:
    tuple val(chrom), path('projection')
    script:
    """
    python3 '${projectDir}/master_resources.py' project --inventory '$inventory' --chromosome $chrom --output projection
    """
}

process AUDIT_RESOURCE_CHROMOSOME {
    tag "chr$chrom: every archived observation against original source"
    publishDir "${params.output}/verification/resources", mode: 'copy'
    input:
    tuple val(chrom), path(projection)
    path inventory
    output:
    path "chr${chrom}.json"
    script:
    """
    python3 '${projectDir}/audit_resources.py' project --inventory '$inventory' --projection '$projection' --output chr${chrom}.json
    """
}

process SCREEN_ALL_RESOURCE_RELATIONSHIPS {
    publishDir "${params.output}", mode: 'copy'
    input:
    path bundles, stageAs: 'bundles/chr??'
    path inventory
    path audits
    output:
    path 'kinship'
    script:
    """
    python3 '${projectDir}/master_resources.py' kinship --inventory '$inventory' --bundles bundles/chr* --output kinship
    """
}

process IMPORT_VERIFIED_RESOURCE_PREPARATION {
    publishDir "${params.output}", mode: 'copy'
    input:
    val prepared
    path inventory
    output:
    path 'archive', emit: archive
    path 'kinship', emit: kinship
    path 'verification', emit: audits
    script:
    """
    python3 '${projectDir}/master_resources.py' import-preparation --prepared '$prepared' \
        --inventory '$inventory' --output imported
    mv imported/archive archive
    mv imported/kinship kinship
    mv imported/verification verification
    """
}

process FREEZE_BIOLOGICAL_FAMILY_ROLES {
    publishDir "${params.output}", mode: 'copy'
    input:
    path inventory
    path kinship
    output:
    path 'cohorts'
    script:
    """
    python3 '${projectDir}/master_resources.py' split --inventory '$inventory' --screen '$kinship' \
        --destination '${params.output}' --output cohorts
    """
}

process KEEP_UNSCORED_CASES {
    tag "chr$chrom: unknown affinity + related-reference/replicate stress"
    publishDir path: { "${params.output}/unscored/chr$chrom" }, mode: 'copy'
    input:
    tuple val(chrom), path(projection)
    path cohorts
    path inventory
    output:
    path 'cases'
    script:
    """
    python3 '${projectDir}/master_resources.py' cases --inventory '$inventory' --cohorts '$cohorts' \
        --projection '$projection' --chromosome $chrom --output cases
    """
}

process FINALIZE_PLATFORM_COHORT {
    tag "$panel: all 22 chromosomes independently checked"
    publishDir path: { "${params.output}/panels/$panel" }, mode: 'copy'
    input:
    tuple val(panel), val(chromosomes), path(masters, stageAs: 'masters/chr??'), path(audits)
    path cohorts
    output:
    path 'metadata', emit: metadata
    path 'configs', emit: configs
    script:
    """
    python3 '${projectDir}/master_resources.py' finalize-cohort --cohort '$cohorts/$panel' \
        --bundles masters/chr* --destination '${params.output}/panels/$panel' --output ready
    mv ready/metadata metadata
    mv ready/configs configs
    """
}

process FINALIZE_RESOURCE_MASTER {
    publishDir "${params.output}", mode: 'copy'
    input:
    path inventory
    path cohorts
    path bundles, stageAs: 'bundles/chr??'
    path releases, stageAs: 'releases/panel??'
    path config_releases, stageAs: 'configs/panel??'
    path cases, stageAs: 'cases/chr??'
    output:
    path 'master'
    script:
    """
    python3 '${projectDir}/master_resources.py' finalize-master --inventory '$inventory' --cohorts '$cohorts' \
        --bundles bundles/chr* --releases releases/panel* --config-releases configs/panel* \
        --destination '${params.output}' --output master
    """
}

workflow {
    smoke_only = params.smoke_only.toString().toBoolean()
    inventory = INVENTORY_ALL_RESOURCES(file(params.specification, checkIfExists: true))
    if (params.prepared_resources) {
        imported = IMPORT_VERIFIED_RESOURCE_PREPARATION(params.prepared_resources, inventory)
        projections = imported.archive.flatMap { directory ->
            (1..22).collect { chrom -> tuple(chrom, directory.resolve("chr$chrom/projection")) }
        }
        source_audits = imported.audits
        kinship = imported.kinship
    } else {
        projections = PROJECT_RESOURCE_CHROMOSOME(inventory.flatMap { directory ->
            (smoke_only ? [22, 4] : (1..22).toList()).collect { chrom -> tuple(chrom, directory) }
        })
        source_audits = AUDIT_RESOURCE_CHROMOSOME(projections, inventory)
        if (!smoke_only) {
            kinship = SCREEN_ALL_RESOURCE_RELATIONSHIPS(
                projections.map { chrom, directory -> directory }.collect(), inventory, source_audits.collect())
        }
    }
    if (!smoke_only) {
        cohorts = FREEZE_BIOLOGICAL_FAMILY_ROLES(inventory, kinship)
        panels = ['observed', 'gse77508', 'gse73996', 'sgdp']
        smoke_jobs = projections.filter { chrom, directory -> chrom == 22 }.flatMap { chrom, directory ->
            panels.collect { panel -> tuple(panel, chrom, directory) }
        }
        smoke_bases = SMOKE_BASE(smoke_jobs, cohorts, inventory)
        smoke_masters = SMOKE_SIMULATE(smoke_bases, cohorts)
        smoke_exports = SMOKE_EXPORT(smoke_masters, cohorts)
        smoke_audits = SMOKE_AUDIT(
            smoke_jobs.join(smoke_bases, by: [0, 1]).join(smoke_masters, by: [0, 1]).join(smoke_exports, by: [0, 1]), cohorts)
        smoke_gate = smoke_audits[0].map { panel, chrom, audit -> audit }.collect().map { checked ->
            if (checked.size() != 4) { throw new IllegalStateException('All four real chr22 platform audits must pass') }
            true
        }
        full_jobs = projections.filter { chrom, directory -> chrom != 22 }.combine(smoke_gate)
            .flatMap { chrom, directory, gate -> panels.collect { panel -> tuple(panel, chrom, directory) } }
        bases = BASE_PANEL(full_jobs, cohorts, inventory)
        masters = SIMULATE_PANEL(bases, cohorts)
        exports = EXPORT_PANEL(masters, cohorts)
        audits = AUDIT_PANEL(
            full_jobs.join(bases, by: [0, 1]).join(masters, by: [0, 1]).join(exports, by: [0, 1]), cohorts)
        all_masters = smoke_masters.mix(masters)
        all_audits = smoke_audits[0].mix(audits[0]).map { panel, chrom, audit -> tuple(panel, audit) }.groupTuple()
        cohort_releases = FINALIZE_PLATFORM_COHORT(all_masters.groupTuple().join(all_audits), cohorts)
        cases = KEEP_UNSCORED_CASES(projections, cohorts, inventory)
        FINALIZE_RESOURCE_MASTER(inventory, cohorts, projections.map { chrom, directory -> directory }.collect(),
                                 cohort_releases.metadata.collect(), cohort_releases.configs.collect(), cases.collect())
    }
}
