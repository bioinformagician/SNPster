nextflow.enable.dsl = 2

include {
    PREPARE_ASSEMBLY_RESOURCES
    NORMALIZE_SEQUENCING_CALLS
    VERIFY_NORMALIZED_GENOTYPES
    AUDIT_CANDIDATE_RELATEDNESS
    COLLECT_NORMALIZATION_AND_IDENTITY_QC
} from './public_sgdp_normalization'
include { PREPARE_GENETIC_IDENTITY_PANEL_V2 } from './public_sgdp_normalization_v2'

workflow {
    def mode = params.smoke.toString().toLowerCase()
    if (!(mode in ['true', 'false'])) {
        error "smoke must be explicitly true or false, not '${params.smoke}'"
    }
    def smokeOnly = mode == 'true'
    def chromosomes = smokeOnly ? [22] : (1..22).toList()
    log.info "Validated normalization scope: smoke=${smokeOnly}; chromosomes=${chromosomes.join(',')}"
    manifest = new groovy.json.JsonSlurper().parse(file(params.input_manifest))
    resources = PREPARE_ASSEMBLY_RESOURCES(file(params.input_manifest, checkIfExists: true))
    jobs = Channel.fromList(chromosomes.collect { chromosome ->
        def root = manifest.data_roots["${chromosome}"]
        tuple(chromosome, file("${root}/basecalls/basecalls.vcf.gz", checkIfExists: true),
              file("${root}/basecalls/basecalls.vcf.gz.csi", checkIfExists: true))
    })
    normalized = NORMALIZE_SEQUENCING_CALLS(jobs.combine(resources))
    verified = VERIFY_NORMALIZED_GENOTYPES(normalized.combine(resources))
    panels = PREPARE_GENETIC_IDENTITY_PANEL_V2(verified)
    if (!smokeOnly) {
        relatedness = AUDIT_CANDIDATE_RELATEDNESS(panels.collect())
        COLLECT_NORMALIZATION_AND_IDENTITY_QC(
            normalized.map { chromosome, folder, calls -> folder }.collect(),
            verified.map { chromosome, folder, report -> report }.collect(),
            panels.collect(), relatedness)
    }
}
