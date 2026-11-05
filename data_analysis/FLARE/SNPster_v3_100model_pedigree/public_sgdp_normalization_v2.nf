nextflow.enable.dsl = 2

include {
    PREPARE_ASSEMBLY_RESOURCES
    NORMALIZE_SEQUENCING_CALLS
    VERIFY_NORMALIZED_GENOTYPES
    AUDIT_CANDIDATE_RELATEDNESS
    COLLECT_NORMALIZATION_AND_IDENTITY_QC
} from './public_sgdp_normalization'

process PREPARE_GENETIC_IDENTITY_PANEL_V2 {
    tag "chr${chromosome}"
    cpus 1
    memory '3 GB'
    publishDir "${params.run_dir}/identity_panels", mode: 'copy', pattern: 'chr*'
    input:
    tuple val(chromosome), path(folder, stageAs: 'normalized_input'), path(verification)
    output:
    path "chr${chromosome}"
    script:
    """
    python3 "${params.project}/pedigree_benchmark/scripts/sgdp_genetic_qc.py" panel \
        --project "${params.project}" --chromosome ${chromosome} \
        --normalized "${folder}/normalized.vcf.gz" --verification "${verification}" \
        --production-root "${params.production_root}" --output "chr${chromosome}"
    """
}

workflow {
    manifest = new groovy.json.JsonSlurper().parse(file(params.input_manifest))
    resources = PREPARE_ASSEMBLY_RESOURCES(file(params.input_manifest, checkIfExists: true))
    chromosomes = params.smoke ? [22] : (1..22).toList()
    jobs = Channel.fromList(chromosomes.collect { chromosome ->
        def root = manifest.data_roots["${chromosome}"]
        tuple(chromosome, file("${root}/basecalls/basecalls.vcf.gz", checkIfExists: true),
              file("${root}/basecalls/basecalls.vcf.gz.csi", checkIfExists: true))
    })
    normalized = NORMALIZE_SEQUENCING_CALLS(jobs.combine(resources))
    verified = VERIFY_NORMALIZED_GENOTYPES(normalized.combine(resources))
    panels = PREPARE_GENETIC_IDENTITY_PANEL_V2(verified)
    if (!params.smoke) {
        relatedness = AUDIT_CANDIDATE_RELATEDNESS(panels.collect())
        COLLECT_NORMALIZATION_AND_IDENTITY_QC(
            normalized.map { chromosome, folder, calls -> folder }.collect(),
            verified.map { chromosome, folder, report -> report }.collect(),
            panels.collect(), relatedness)
    }
}
