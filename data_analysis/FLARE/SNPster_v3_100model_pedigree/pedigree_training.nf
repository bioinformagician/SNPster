nextflow.enable.dsl=2

params.evaluation_plan = params.evaluation_plan ?: '/workspace/results_v3_pedigree_pilot_10/pedigree/plans/plans/expected_ancestry.tsv'
params.pedigree_chromosomes = params.pedigree_chromosomes ?: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
params.pedigree_plink_map = params.pedigree_plink_map ?: params.genetic_map
params.pedigree_families = params.pedigree_families ?: 108
params.pedigree_batch_size = params.pedigree_batch_size ?: 10
params.pedigree_seed = params.pedigree_seed ?: 20261014
params.pedigree_reuse_dir = params.pedigree_reuse_dir ?: null
params.frozen_candidates_manifest = params.frozen_candidates_manifest ?: null

include { VALIDATE_FROZEN_INPUTS } from './frozen_inputs.nf'
include {
    PDB_SIMULATE; PDB_AGGREGATE_TRUTH; PDB_MERGE_PANEL; PDB_FLARE_INFER
} from './main.nf'

process PREPARE_TRAINING_PLANS {
    label 'pdb_light'
    publishDir "${params.outdir}/pedigree/plans", mode:'copy'
    input:
      path training_meta
      path validation_meta
      path labels
      path founder_vcf
      path reference_vcf
      path reference_panel
      path evaluation_plan
      path validated
    output:
      path 'plans', emit: plans
    script:
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/prepare_training_pedigrees.py' \\
      --training-meta '${training_meta}' --validation-meta '${validation_meta}' \\
      --label-meta '${labels}' --founder-vcf '${founder_vcf}' \\
      --reference-vcf '${reference_vcf}' --reference-panel '${reference_panel}' \\
      --evaluation-plan '${evaluation_plan}' --outdir plans \\
      --families ${params.pedigree_families} --batch-size ${params.pedigree_batch_size} \\
      --seed ${params.pedigree_seed}
    """
}

process SPLIT_TRAINING_FOUNDER_VCF {
    tag "training founders chr${chr}"
    label 'pdb_light'
    input:
      tuple val(chr), path(vcf), path(index)
      path plans
    output:
      tuple val(chr), path("training_founders.chr${chr}.vcf.gz"), emit: chromosomes
    script:
    """
    set -euo pipefail
    bcftools view -r '${chr}' -m2 -M2 -v snps \\
      -Oz -o 'training_founders.chr${chr}.vcf.gz' '${vcf}'
    """
}

workflow {
    if (params.pedigree_design != 'full' || params.pedigree_chromosomes != '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22')
        error 'Training screen requires explicitly full design and all 22 autosomes'
    if ((params.pedigree_families as int) < 9 || (params.pedigree_families as int) % 9 != 0)
        error 'Training cohort requires balanced complete sets of nine classes'
    if (!params.pedigree_reuse_dir || !params.frozen_candidates_manifest)
        error 'Training cohort requires frozen v1 and validated reusable static inputs'
    VALIDATE_FROZEN_INPUTS(
        file(params.frozen_candidates_manifest, checkIfExists:true),
        file(params.pedigree_reuse_dir, checkIfExists:true),
        file(params.pedigree_plink_map, checkIfExists:true),
        file(params.flare_ref, checkIfExists:true)
    )
    rows = VALIDATE_FROZEN_INPUTS.out.manifest.splitCsv(header:true, sep:'\t')
    models = rows.map { r ->
        if (r.representation != 'v1_frozen_k18' || r.candidate_id != 'v1_frozen_k18')
            error 'Only frozen v1 is allowed in the training screen'
        tuple(r.representation.toString(), file(r.model_path, checkIfExists:true),
              r.min_maf as double, r.min_mac as int, r.flare_seed as int)
    }.unique()
    PREPARE_TRAINING_PLANS(
        file(params.pedigree_training_meta, checkIfExists:true),
        file(params.pedigree_validation_meta, checkIfExists:true),
        file(params.pedigree_label_meta, checkIfExists:true),
        file(params.source_train_vcf, checkIfExists:true),
        file(params.flare_ref, checkIfExists:true),
        file(params.flare_panel, checkIfExists:true),
        file(params.evaluation_plan, checkIfExists:true),
        VALIDATE_FROZEN_INPUTS.out.manifest
    )
    batches = PREPARE_TRAINING_PLANS.out.plans.flatMap { dir ->
        dir.toFile().listFiles().findAll { it.isDirectory() && it.name.startsWith('batch_') }.sort { it.name }
            .collect { folder -> tuple(folder.name, file(folder)) }
    }
    founder_jobs = Channel.fromList((1..22).collect { chr ->
        tuple(chr, file(params.source_train_vcf, checkIfExists:true),
              file("${params.source_train_vcf}.tbi", checkIfExists:true))
    })
    SPLIT_TRAINING_FOUNDER_VCF(founder_jobs, PREPARE_TRAINING_PLANS.out.plans)
    maps = Channel.value(file("${params.pedigree_reuse_dir}/maps", checkIfExists:true))
    sites = Channel.value(file("${params.pedigree_reuse_dir}/panel_sites.tsv", checkIfExists:true))
    jobs = batches.combine(SPLIT_TRAINING_FOUNDER_VCF.out.chromosomes).combine(maps).combine(sites)
        .map { batch,plan,chr,vcf,mapdir,panel -> tuple(batch,chr,plan,vcf,mapdir,panel) }
    PDB_SIMULATE(jobs)
    data = PDB_SIMULATE.out.panel_truth
    PDB_AGGREGATE_TRUTH(data.map { batch,chr,vcf,truth -> truth }.collect())
    PDB_MERGE_PANEL(data.map { batch,chr,vcf,truth -> tuple(batch,chr,vcf) }.groupTuple())
    inference = models.combine(PDB_MERGE_PANEL.out.genomes)
        .map { rep,model,maf,mac,seed,batch,genome -> tuple(rep,model,maf,mac,seed,batch,genome) }
    PDB_FLARE_INFER(inference)
}
