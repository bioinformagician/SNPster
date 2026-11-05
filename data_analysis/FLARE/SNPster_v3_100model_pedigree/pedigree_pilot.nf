nextflow.enable.dsl=2

params.frozen_candidates_manifest = params.frozen_candidates_manifest ?: null
params.pedigree_reuse_dir = params.pedigree_reuse_dir ?: null
params.pedigree_families = params.pedigree_families ?: 27
params.pedigree_batch_size = params.pedigree_batch_size ?: 10
params.pedigree_seed = params.pedigree_seed ?: 20261008
params.pedigree_chromosomes = params.pedigree_chromosomes ?: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
params.pedigree_allow_reuse = params.pedigree_allow_reuse == null ? true : params.pedigree_allow_reuse
params.pedigree_plink_map = params.pedigree_plink_map ?: params.genetic_map

include {
    PDB_PREPARE_FOUNDERS; PDB_BUILD_MAPS; PDB_PANEL_SITES; PDB_MAKE_PLANS
    PDB_SIMULATE; PDB_AGGREGATE_TRUTH; PDB_MERGE_PANEL; PDB_FLARE_INFER
    PDB_EVALUATE_FINALIST; PDB_SUMMARIZE; PDB_SELECT_UNIFIED_FINALISTS
    AGGREGATE_FINAL
} from './main.nf'

include { VALIDATE_FROZEN_INPUTS } from './frozen_inputs.nf'

workflow {
    if (!params.frozen_candidates_manifest)
        error 'Frozen pilot requires --frozen_candidates_manifest'
    def chromosomes = params.pedigree_chromosomes.toString().tokenize(',').collect { it.trim() as int }
    if (!chromosomes || chromosomes.toSet().size() != chromosomes.size() || !chromosomes.every { it >= 1 && it <= 22 })
        error 'Pilot chromosomes must be unique autosomes 1..22'
    if ((params.pedigree_families as int) < 9 || (params.pedigree_batch_size as int) < 1)
        error 'Pilot requires at least nine families and a positive batch size'
    if (!(params.pedigree_design.toString() in ['full','reduced']))
        error 'Pilot design must be full or reduced'

    VALIDATE_FROZEN_INPUTS(
        file(params.frozen_candidates_manifest, checkIfExists:true),
        params.pedigree_reuse_dir ? file(params.pedigree_reuse_dir, checkIfExists:true) : [],
        file(params.pedigree_plink_map, checkIfExists:true),
        file(params.flare_ref, checkIfExists:true)
    )
    rows = VALIDATE_FROZEN_INPUTS.out.manifest.splitCsv(header:true, sep:'\t')
    candidates = rows.map { r ->
        tuple(r.representation.toString(), r.candidate_id.toString(), file(r.candidate_dir, checkIfExists:true))
    }
    models = rows.map { r ->
        tuple(r.representation.toString(), file(r.model_path, checkIfExists:true),
              r.min_maf as double, r.min_mac as int, r.flare_seed as int)
    }.unique()
    candidate_dirs = candidates.map { rep,cid,dir -> dir }

    PDB_PREPARE_FOUNDERS(
        file(params.pedigree_validation_meta, checkIfExists:true),
        file(params.pedigree_training_meta, checkIfExists:true),
        params.pedigree_additional_meta ? file(params.pedigree_additional_meta, checkIfExists:true) : [],
        file(params.pedigree_label_meta, checkIfExists:true),
        file(params.pedigree_mosaic_meta, checkIfExists:true),
        file(params.flare_ref, checkIfExists:true),
        file(params.flare_panel, checkIfExists:true),
        VALIDATE_FROZEN_INPUTS.out.manifest
    )
    PDB_MAKE_PLANS(PDB_PREPARE_FOUNDERS.out.donors,
                   file(params.pedigree_training_meta, checkIfExists:true),
                   VALIDATE_FROZEN_INPUTS.out.manifest)
    if (params.pedigree_reuse_dir) {
        maps = Channel.value(file("${params.pedigree_reuse_dir}/maps", checkIfExists:true))
        sites = Channel.value(file("${params.pedigree_reuse_dir}/panel_sites.tsv", checkIfExists:true))
    } else {
        PDB_BUILD_MAPS(file(params.pedigree_plink_map, checkIfExists:true), VALIDATE_FROZEN_INPUTS.out.manifest)
        PDB_PANEL_SITES(file(params.flare_ref, checkIfExists:true), VALIDATE_FROZEN_INPUTS.out.manifest)
        maps = PDB_BUILD_MAPS.out.maps
        sites = PDB_PANEL_SITES.out.sites
    }
    batches = PDB_MAKE_PLANS.out.plans.flatMap { dir ->
        dir.toFile().listFiles().findAll { it.isDirectory() }.sort { it.name }
            .collect { folder -> tuple(folder.name, file(folder)) }
    }
    references = Channel.fromList(chromosomes.collect { chr ->
        tuple(chr, file("${params.pedigree_vcf_dir}/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz",
                        checkIfExists:true))
    })
    jobs = batches.combine(references).combine(maps).combine(sites)
        .map { batch,plan,chr,vcf,mapdir,panel -> tuple(batch,chr,plan,vcf,mapdir,panel) }
    PDB_SIMULATE(jobs)
    data = PDB_SIMULATE.out.panel_truth
    PDB_AGGREGATE_TRUTH(data.map { batch,chr,vcf,truth -> truth }.collect())
    PDB_MERGE_PANEL(data.map { batch,chr,vcf,truth -> tuple(batch,chr,vcf) }.groupTuple())
    inference = models.combine(PDB_MERGE_PANEL.out.genomes)
        .map { rep,model,maf,mac,seed,batch,genome -> tuple(rep,model,maf,mac,seed,batch,genome) }
    PDB_FLARE_INFER(inference)
    features = PDB_FLARE_INFER.out.features.map { rep,batch,f -> tuple(rep,f) }.groupTuple()
    evaluations = candidates.combine(features, by:0).combine(PDB_AGGREGATE_TRUTH.out.truth)
        .map { rep,cid,dir,files,truth -> tuple(rep,cid,dir,files,truth) }
    PDB_EVALUATE_FINALIST(evaluations, PDB_PREPARE_FOUNDERS.out.summary,
        PDB_MAKE_PLANS.out.plans.map { dir -> file("${dir}/manifest.json", checkIfExists:true) })
    PDB_SUMMARIZE(PDB_EVALUATE_FINALIST.out.metrics.collect())
    AGGREGATE_FINAL(candidate_dirs.collect())
    PDB_SELECT_UNIFIED_FINALISTS(candidate_dirs.collect(), PDB_EVALUATE_FINALIST.out.predictions.collect())
}
