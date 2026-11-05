nextflow.enable.dsl=2

params.search_space           = params.search_space           ?: '/workspace/config/search_space.json'
params.generated_panel_configs= params.generated_panel_configs?: '/workspace/config/generated/panel_configs.tsv'
params.generated_model_configs= params.generated_model_configs?: '/workspace/config/generated/model_configs/*.json'
params.fixed_models           = params.fixed_models           ?: '/workspace/config/fixed_models.tsv'
params.dev_datasets           = params.dev_datasets           ?: '/workspace/config/dev_datasets.tsv'
params.diagnostic_datasets    = params.diagnostic_datasets    ?: '/workspace/config/diagnostic_datasets.tsv'
params.full_datasets          = params.full_datasets          ?: '/workspace/config/full_datasets.tsv'
params.real_benchmarks        = params.real_benchmarks        ?: '/workspace/config/real_benchmarks.tsv'
params.stage1_criteria        = params.stage1_criteria        ?: '/workspace/config/stage1_criteria.json'
params.stage2_criteria        = params.stage2_criteria        ?: '/workspace/config/stage2_criteria.json'
params.stage3_criteria        = params.stage3_criteria        ?: '/workspace/config/stage3_criteria.json'
params.outdir                 = params.outdir                 ?: '/workspace/results_v3_sweep_100_disjoint_reduced'
params.flare_jar              = params.flare_jar              ?: '/opt/flare/flare.jar'
params.flare_ref              = params.flare_ref              ?: '/data/flare2_test/source_model/source_reference.autosomes.vcf.gz'
params.flare_panel            = params.flare_panel            ?: '/data/flare2_test/source_model/source_reference.ref.panel'
params.source_train_vcf       = params.source_train_vcf       ?: '/data/flare2_test/source_model/source_train.autosomes.vcf.gz'
params.genetic_map            = params.genetic_map            ?: '/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map'
params.max_stage2_candidates  = params.max_stage2_candidates  ?: 0
params.config_fingerprint      = params.config_fingerprint      ?: 'unversioned'
params.resume_stage2_manifest = params.resume_stage2_manifest ?: null
params.resume_source_outdir = params.resume_source_outdir ?: params.outdir
// Pedigree settings explicitly initialized: prior zip left these undefined.
params.pedigree_enabled = params.pedigree_enabled == null ? true : params.pedigree_enabled
params.pedigree_families = params.pedigree_families ?: 90
params.pedigree_batch_size = params.pedigree_batch_size ?: 10
params.pedigree_seed = params.pedigree_seed ?: 20261008
params.pedigree_chromosomes = params.pedigree_chromosomes ?: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
params.pedigree_allow_reuse = params.pedigree_allow_reuse == null ? true : params.pedigree_allow_reuse
params.pedigree_keep_high_density = params.pedigree_keep_high_density == null ? false : params.pedigree_keep_high_density
params.pedigree_vcf_dir = params.pedigree_vcf_dir ?: '/data/hgdp_1kg_references'
params.pedigree_training_meta = params.pedigree_training_meta ?: '/data/flare2_test/source_model/source_train_metadata.tsv'
params.pedigree_validation_meta = params.pedigree_validation_meta ?: '/data/flare2_test/source_model/source_validation_metadata.tsv'
params.pedigree_additional_meta = params.pedigree_additional_meta ?: null
params.pedigree_label_meta = params.pedigree_label_meta ?: '/data/flare2_test/source_model/ancestry_role_audit.tsv'
params.pedigree_mosaic_meta = params.pedigree_mosaic_meta ?: '/data/flare2_test/source_model/training_mosaics8/training_mosaics8_metadata.tsv'
params.pedigree_design = params.pedigree_design ?: 'reduced'
params.pedigree_plink_map = params.pedigree_plink_map ?: '/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map'
params.pedigree_pedsim = params.pedigree_pedsim ?: '/opt/ped-sim/ped-sim'
params.pedigree_min_per_source = params.pedigree_min_per_source ?: 2
params.pedigree_interference = params.pedigree_interference ?: ''
params.pedigree_fingerprint = params.pedigree_fingerprint ?: 'unversioned'
params.unified_fingerprint = params.unified_fingerprint ?: 'unversioned'
params.unified_selection_config = params.unified_selection_config ?: '/workspace/config/unified_selection.json'


process PANEL_PROBS {
    tag "${panel_id}"
    label 'flare_train'
    publishDir "${params.outdir}/panel_probs", mode: 'copy', overwrite: true

    input:
    tuple val(panel_id), val(panel_cm), val(min_maf), val(min_mac)

    output:
    tuple val(panel_id), path("${panel_id}.panels"), emit: panels
    path "${panel_id}.log", emit: logs

    script:
    """
    set -euo pipefail
    java -Xmx48g -jar "${params.flare_jar}" \\
      gt="${params.source_train_vcf}" \\
      ref="${params.flare_ref}" \\
      ref-panel="${params.flare_panel}" \\
      map="${params.genetic_map}" \\
      panel-probs=true \\
      panel-cm=${panel_cm} \\
      min-maf=${min_maf} \\
      min-mac=${min_mac} \\
      array=false \\
      nthreads=${task.cpus} \\
      seed=20261002 \\
      out="${panel_id}"
    """
}

process BUILD_INITIAL_MODEL {
    tag "${rep_id}"
    label 'model_build'
    publishDir "${params.outdir}/models", mode: 'copy', overwrite: true

    input:
    tuple val(panel_id), val(rep_id), path(model_config), val(min_maf), val(min_mac), val(flare_seed), val(refinement_mode), val(refine_em), path(panels_file)

    output:
    tuple val(rep_id), path("${rep_id}.initial.model"), path("${rep_id}.config.json"), val(min_maf), val(min_mac), val(flare_seed), val(refinement_mode), val(refine_em), emit: initial
    path "${rep_id}.initial.log"
    path "${rep_id}.initial.json"

    script:
    """
    set -euo pipefail
    OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} \\
    python3 "${projectDir}/scripts/create_model_file_strategy.py" \\
      --config "${model_config}" \\
      --panels "${panels_file}" \\
      --output-prefix "${rep_id}.initial"

    cp "${model_config}" "${rep_id}.config.json"
    """
}

process REFINE_MODEL {
    tag "${rep_id}"
    label 'flare_train'
    publishDir "${params.outdir}/models", mode: 'copy', overwrite: true

    input:
    tuple val(rep_id), path(initial_model), path(config_json, stageAs: 'refinement_input.config.json'), val(min_maf), val(min_mac), val(flare_seed), val(refinement_mode), val(refine_em)

    output:
    tuple val(rep_id), path("${rep_id}.final.model"), val(min_maf), val(min_mac), val(flare_seed), path("${rep_id}.config.json"), emit: models
    path "${rep_id}.train.log"

    script:
    """
    set -euo pipefail
    if [[ "${refinement_mode}" == "none" ]]; then
      cp "${initial_model}" "${rep_id}.final.model"
      echo "FLARE refinement skipped by configuration" > "${rep_id}.train.log"
    elif [[ "${refinement_mode}" == "update_p" ]]; then
      java -Xmx64g -jar "${params.flare_jar}" \\
        gt="${params.source_train_vcf}" \\
        ref="${params.flare_ref}" \\
        ref-panel="${params.flare_panel}" \\
        map="${params.genetic_map}" \\
        model="${initial_model}" \\
        update-p=true \\
        em=${refine_em} \\
        probs=false \\
        array=false \\
        min-maf=${min_maf} \\
        min-mac=${min_mac} \\
        nthreads=${task.cpus} \\
        seed=${flare_seed} \\
        out="${rep_id}.train"

      cp "${rep_id}.train.model" "${rep_id}.final.model"
      rm -f "${rep_id}.train.anc.vcf.gz" "${rep_id}.train.global.anc.gz"
    else
      echo "Unsupported refinement_mode=${refinement_mode}" >&2
      exit 2
    fi
    cp refinement_input.config.json "${rep_id}.config.json"
    """
}
process IMPORT_FIXED_MODEL {
    tag "${rep_id}"
    label 'light'
    publishDir "${params.outdir}/models", mode: 'copy', overwrite: true

    input:
    tuple val(rep_id), val(model_path), val(min_maf), val(min_mac), val(flare_seed), val(description)

    output:
    tuple val(rep_id), path("${rep_id}.final.model"), val(min_maf), val(min_mac), val(flare_seed), path("${rep_id}.config.json"), emit: models

    script:
    """
    set -euo pipefail
    cp "${model_path}" "${rep_id}.final.model"
    python3 - <<'PY'
import json
from pathlib import Path
Path("${rep_id}.config.json").write_text(json.dumps({
  "rep_id":"${rep_id}", "source":"fixed", "original_model":"${model_path}",
  "min_maf":${min_maf}, "min_mac":${min_mac}, "flare_seed":${flare_seed},
  "description":"${description}"
}, indent=2, sort_keys=True))
PY
    """
}

process RUN_REPRESENTATION_DEV {
    tag "${rep_id} :: ${dataset_id}"
    label 'flare_infer'
    publishDir "${params.outdir}/features/dev", mode: 'copy', overwrite: true

    input:
    tuple val(dataset_id), val(role), val(input_kind), val(input_path), val(truth_path), val(truth_mode), val(weight_tag), val(platform), val(array_flag),
          val(rep_id), path(model_file), val(min_maf), val(min_mac), val(flare_seed), path(model_meta)

    output:
    tuple val(rep_id), val(dataset_id), path("${rep_id}__${dataset_id}.features.tsv"), emit: features
    path "${rep_id}__${dataset_id}.log", optional: true

    script:
    def arr = array_flag.toString().toLowerCase()
    """
    set -euo pipefail
    if [[ "${input_kind}" == "features" ]]; then
      cp "${input_path}" "${rep_id}__${dataset_id}.features.tsv"
      exit 0
    fi
    if [[ "${input_kind}" == "global_anc" ]]; then
      cp "${input_path}" global.anc.gz
    elif [[ "${input_kind}" == "vcf" ]]; then
      java -Xmx32g -jar "${params.flare_jar}" \\
        gt="${input_path}" ref="${params.flare_ref}" ref-panel="${params.flare_panel}" map="${params.genetic_map}" \\
        model="${model_file}" update-p=false em=false probs=false array=${arr} \\
        min-maf=${min_maf} min-mac=${min_mac} nthreads=${task.cpus} seed=${flare_seed} out=run
      cp run.global.anc.gz global.anc.gz
      cp run.log "${rep_id}__${dataset_id}.log"
    else
      echo "Unsupported input_kind=${input_kind}" >&2; exit 2
    fi
    python3 "${projectDir}/scripts/extract_global_anc.py" --input global.anc.gz --output "${rep_id}__${dataset_id}.features.tsv"
    rm -f run.anc.vcf.gz run.anc.vcf.gz.tbi run.model run.global.anc.gz global.anc.gz
    """
}

process SCREEN_REPRESENTATION {
    tag "${rep_id}"
    label 'decoder'
    publishDir "${params.outdir}/stage1_screen", mode: 'copy', overwrite: true

    input:
    tuple val(rep_id), val(dataset_ids), path(feature_files)

    output:
    tuple val(rep_id), path("screen_${rep_id}"), emit: result_dirs

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p "screen_${rep_id}"
    python3 "${projectDir}/scripts/make_screen_candidate.py" --rep-id "${rep_id}" --search-space "${params.search_space}" --output screen.json
    SKLEARN_N_JOBS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} \\
    python3 "${projectDir}/scripts/fit_evaluate_candidate.py" \\
      --candidate screen.json --datasets "${params.dev_datasets}" --features-dir . \\
      --real-benchmarks "${params.real_benchmarks}" --success-criteria "${params.stage1_criteria}" \\
      --outdir "screen_${rep_id}"
    """
}

process SELECT_REPRESENTATIONS {
    tag 'select representations'
    label 'light'
    publishDir "${params.outdir}/selection", mode: 'copy', overwrite: true

    input:
    path screen_dirs

    output:
    path 'selected_representations.tsv', emit: selected

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    python3 "${projectDir}/scripts/select_representations.py" --root . --search-space "${params.search_space}" \\
      --fixed-models "${params.fixed_models}" --model-config-dir "${projectDir}/config/generated/model_configs" --output selected_representations.tsv
    """
}

process GENERATE_STAGE2_CANDIDATES {
    tag 'generate stage2 candidates'
    label 'light'
    publishDir "${params.outdir}/selection", mode: 'copy', overwrite: true

    input:
    path selected_reps

    output:
    path 'stage2_candidates/*.json', emit: candidates

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p stage2_candidates
    python3 "${projectDir}/scripts/generate_decoder_candidates.py" --selected-reps "${selected_reps}" \\
      --search-space "${params.search_space}" --model-config-dir "${projectDir}/config/generated/model_configs" --outdir stage2_candidates --max-candidates ${params.max_stage2_candidates}
    """
}

process FIT_STAGE2_CANDIDATE {
    tag { candidate_config.baseName }
    label 'decoder'
    publishDir "${params.outdir}/stage2_candidates", mode: 'copy', overwrite: true

    input:
    path candidate_config
    path feature_files

    output:
    path "candidate_*", emit: result_dirs

    script:
    def cid=candidate_config.baseName
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p "candidate_${cid}"
    SKLEARN_N_JOBS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} \\
    python3 "${projectDir}/scripts/fit_evaluate_candidate.py" \\
      --candidate "${candidate_config}" --datasets "${params.dev_datasets}" --features-dir . \\
      --real-benchmarks "${params.real_benchmarks}" --success-criteria "${params.stage2_criteria}" \\
      --outdir "candidate_${cid}"
    """
}

process SELECT_STAGE2_CANDIDATES {
    tag 'select stage2 candidates'
    label 'light'
    publishDir "${params.outdir}/selection", mode: 'copy', overwrite: true

    input:
    path results_manifest

    output:
    path 'selected_stage2/*.json', emit: candidates
    path 'selected_stage2.tsv', emit: table

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p selected_stage2
    if [[ '${params.pedigree_enabled}' == 'true' ]]; then
      python3 "${projectDir}/pedigree_benchmark/scripts/select_stage2_diverse.py" --manifest '${results_manifest}' --search-space "${params.search_space}" \\
        --outdir selected_stage2 --table selected_stage2.tsv
    else
      python3 "${projectDir}/scripts/select_candidates.py" --manifest '${results_manifest}' --search-space "${params.search_space}" \\
        --outdir selected_stage2 --table selected_stage2.tsv
    fi
    """
}

process RUN_REPRESENTATION_DIAGNOSTIC {
    tag "${rep_id} :: ${dataset_id}"
    label 'flare_infer'
    publishDir "${params.outdir}/features/diagnostic", mode: 'copy', overwrite: true

    input:
    tuple val(dataset_id), val(role), val(input_kind), val(input_path), val(truth_path), val(truth_mode), val(weight_tag), val(platform), val(array_flag),
          val(rep_id), path(model_file), val(min_maf), val(min_mac), val(flare_seed), path(model_meta)

    output:
    tuple val(rep_id), val(dataset_id), path("${rep_id}__${dataset_id}.features.tsv"), emit: features

    script:
    def arr=array_flag.toString().toLowerCase()
    """
    set -euo pipefail
    if [[ "${input_kind}" == "vcf" ]]; then
      java -Xmx32g -jar "${params.flare_jar}" \\
        gt="${input_path}" ref="${params.flare_ref}" ref-panel="${params.flare_panel}" map="${params.genetic_map}" \\
        model="${model_file}" update-p=false em=false probs=false array=${arr} \\
        min-maf=${min_maf} min-mac=${min_mac} nthreads=${task.cpus} seed=${flare_seed} out=run
      python3 "${projectDir}/scripts/extract_global_anc.py" --input run.global.anc.gz --output "${rep_id}__${dataset_id}.features.tsv"
    elif [[ "${input_kind}" == "global_anc" ]]; then
      python3 "${projectDir}/scripts/extract_global_anc.py" --input "${input_path}" --output "${rep_id}__${dataset_id}.features.tsv"
    elif [[ "${input_kind}" == "features" ]]; then
      cp "${input_path}" "${rep_id}__${dataset_id}.features.tsv"
    else
      echo "Unsupported input_kind=${input_kind}" >&2; exit 2
    fi
    rm -f run.anc.vcf.gz run.anc.vcf.gz.tbi run.model run.global.anc.gz
    """
}

process FULL_EVALUATE_CANDIDATE {
    tag { candidate_config.baseName }
    label 'decoder'
    publishDir "${params.outdir}/stage3_full", mode: 'copy', overwrite: true

    input:
    path candidate_config
    path feature_files

    output:
    path "candidate_*", emit: result_dirs

    script:
    def cid=candidate_config.baseName
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p "candidate_${cid}"
    SKLEARN_N_JOBS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} \\
    python3 "${projectDir}/scripts/fit_evaluate_candidate.py" \\
      --candidate "${candidate_config}" --datasets "${params.full_datasets}" --features-dir . \\
      --real-benchmarks "${params.real_benchmarks}" --success-criteria "${params.stage3_criteria}" \\
      --outdir "candidate_${cid}" --save-model
    """
}

process AGGREGATE_FINAL {
    tag 'final leaderboard'
    label 'light'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path result_dirs

    output:
    path 'leaderboard.tsv'
    path 'leaderboard.html'
    path 'all_metrics.json'
    path 'strategy_summary.tsv'
    path 'model_family_summary.tsv'

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    python3 "${projectDir}/scripts/aggregate_results_v3.py" --root . --outdir aggregate --title 'SNPster ancestry strategy benchmark v3'
    cp aggregate/leaderboard.tsv .
    cp aggregate/leaderboard.html .
    cp aggregate/all_metrics.json .
    cp aggregate/strategy_summary.tsv .
    cp aggregate/model_family_summary.tsv .
    """
}

process SELECT_FINALISTS {
    tag 'select finalists'
    label 'light'
    publishDir "${params.outdir}/final", mode: 'copy', overwrite: true

    input:
    path result_dirs

    output:
    path 'finalists.tsv', emit: table
    path 'finalists/*.json', emit: finalists

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    python3 - <<'PY'
import json, shutil
from pathlib import Path
import pandas as pd
search=json.loads(Path("${params.search_space}").read_text()); n=int(search['stage3']['retain_finalists'])
rows=[]
for p in Path('.').glob('candidate_*/metrics.json'):
 m=json.loads(p.read_text()); rows.append({'candidate_id':m['candidate_id'],'pass':bool(m.get('pass',False)),'real_hard_pass':bool(m.get('real_hard_pass',False)),'score':float(m.get('selection_score',999)),'dir':str(p.parent)})
df=pd.DataFrame(rows).sort_values(['pass','real_hard_pass','score'],ascending=[False,False,True]).head(n)
Path('finalists').mkdir(exist_ok=True)
for _,r in df.iterrows(): shutil.copy2(Path(r.dir)/'candidate.json', Path('finalists')/(r.candidate_id+'.json'))
df.drop(columns=['dir']).to_csv('finalists.tsv',sep='\t',index=False)
print(df.drop(columns=['dir']).to_string(index=False))
PY
    """
}


/*
 * PEDIGREE UNIFIED SELECTION (enabled by default)
 * Runs after stage-2 shortlist, before finalist selection. Evaluates all
 * saved stage-3 decoders using the SAME ancestry-vector TV objective as
 * reference-source validation, with scenario/combination/family balancing.
 * Old stage-3 leaderboard remains available for diagnostic comparison.
 * Input is clean, high-density phased founder genotypes, projected to FLARE
 * panel sites for inference. This is NOT consumer-array/Beagle validation.
 */
process PDB_PREPARE_FOUNDERS {
    tag 'heldout high-density founders'
    publishDir "${params.outdir}/pedigree/inputs", mode:'copy'
    label 'pdb_light'
    input:
      path validation_meta
      path training_meta
      path additional_meta
      path label_meta
      path mosaic_meta
      path reference_vcf
      path reference_panel
      path finalized
    output:
      path 'eligible_heldout_founders.tsv', emit: donors
      path 'eligible_heldout_founders.tsv.summary.json', emit: summary
    script:
    def additionalOpt = additional_meta ? "--additional-metadata '${additional_meta}'" : ''
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/prepare_founders.py' \\
      --validation '${validation_meta}' --training '${training_meta}' \\
      --vcf '${params.pedigree_vcf_dir}/hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz' \\
      --min-per-source ${params.pedigree_min_per_source} \\
      --reference-vcf '${reference_vcf}' --reference-panel '${reference_panel}' \\
      --label-metadata '${label_meta}' --mosaic-metadata '${mosaic_meta}' \\
      --benchmark-outdir '${params.outdir}' --design '${params.pedigree_design}' \\
      --output eligible_heldout_founders.tsv ${additionalOpt}
    """
}

process PDB_BUILD_MAPS {
    tag 'GRCh38 Ped-sim genetic maps'
    publishDir "${params.outdir}/pedigree/maps", mode:'copy'
    label 'pdb_light'
    input:
      path plink_map
      path finalized
    output:
      path 'maps', emit: maps
    script:
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/convert_plink_map.py' \\
      --plink-map '${plink_map}' --outdir maps
    """
}

process PDB_MAKE_PLANS {
    tag "${params.pedigree_families} pedigrees"
    publishDir "${params.outdir}/pedigree/plans", mode:'copy'
    label 'pdb_light'
    input:
      path donors
      path training_meta
      path finalized  // dependency on stage-2 shortlisting, not final selection
    output:
      path 'plans', emit: plans
    script:
    def reuse = params.pedigree_allow_reuse.toString().toBoolean() ? '--allow-reuse-across-families' : ''
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/make_plans.py' \\
      --metadata '${donors}' --exclude '${training_meta}' \\
      --families ${params.pedigree_families} --batch-size ${params.pedigree_batch_size} \\
      --seed ${params.pedigree_seed} --design '${params.pedigree_design}' --outdir plans ${reuse}
    """
}

process PDB_PANEL_SITES {
    tag 'FLARE panel sites for clean validation'
    label 'pdb_light'
    input:
      path reference_vcf
      path finalized
    output:
      path 'panel_sites.tsv', emit: sites
    script:
    """
    set -euo pipefail
    bcftools query -f '%CHROM\\t%POS\\t%POS\\n' '${reference_vcf}' > panel_sites.tsv
    test -s panel_sites.tsv
    """
}

process PDB_SIMULATE {
    tag "${batch} chr${chr}"
    label 'pdb_sim'
    publishDir "${params.outdir}/pedigree/chromosomes", mode:'copy', pattern: 'pdb_*'
    input:
      tuple val(batch), val(chr), path(plan), path(vcf), path(maps), path(sites)
    output:
      tuple val(batch), val(chr), path("pdb_${batch}_chr${chr}.panel.vcf.gz"), path("pdb_${batch}_chr${chr}.truth.tsv"), emit: panel_truth
      path("pdb_${batch}_chr${chr}.bp"), emit: breakpoints
      path("pdb_${batch}_chr${chr}.full.vcf.gz"), optional: true, emit: full_genome
    script:
    def interferenceOpt = params.pedigree_interference ? "--intf '${params.pedigree_interference}'" : '--pois'
    def simSeed = ((params.pedigree_seed as long) + (chr as long) * 15427L +
                   ((batch.replace('batch_','') as long) * 1000003L)) % 2000000000L
    def keepFull = params.pedigree_keep_high_density.toString().toBoolean()
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    bcftools view -S '${plan}/founder_samples.txt' -m2 -M2 -v snps '${vcf}' -Ou \\
      | bcftools view -i 'F_MISSING=0' -Oz -o founders.vcf.gz
    python3 '${projectDir}/pedigree_benchmark/scripts/check_founder_vcf.py' \\
      --vcf founders.vcf.gz --samples '${plan}/founder_samples.txt' \\
      --chromosome ${chr} --map '${maps}/chr${chr}.simmap'
    '${params.pedigree_pedsim}' -d '${plan}/pedigrees.def' -m '${maps}/chr${chr}.simmap' \\
       -i founders.vcf.gz -o sim_${batch}_${chr} \\
       --set_founders '${plan}/set_founders.tsv' --founder_ids --bp_all --keep_phase \\
       --nogz --err_rate 0 --miss_rate 0 --seed ${simSeed} ${interferenceOpt}
    python3 '${projectDir}/pedigree_benchmark/scripts/verify_founder_assignments.py' \\
      --expected '${plan}/set_founders.tsv' --actual sim_${batch}_${chr}.ids
    # FLARE-clean variant projection: preserve all descendant genotypes at reference markers.
    bcftools view -S '${plan}/descendant_samples.txt' -T '${sites}' \\
      -Oz -o 'pdb_${batch}_chr${chr}.panel.vcf.gz' sim_${batch}_${chr}.vcf
    test "\$(bcftools index -n 'pdb_${batch}_chr${chr}.panel.vcf.gz' 2>/dev/null || bcftools view -H 'pdb_${batch}_chr${chr}.panel.vcf.gz' | wc -l)" -gt 0
    mv sim_${batch}_${chr}.bp pdb_${batch}_chr${chr}.bp
    python3 '${projectDir}/pedigree_benchmark/scripts/realized_ancestry.py' \\
      --plan '${plan}/pedigrees.json' --bp 'pdb_${batch}_chr${chr}.bp' --chr ${chr} \\
      --output 'pdb_${batch}_chr${chr}.truth.tsv'
    if [[ '${keepFull}' == 'true' ]]; then
      bcftools view -S '${plan}/descendant_samples.txt' -Oz \\
        -o 'pdb_${batch}_chr${chr}.full.vcf.gz' sim_${batch}_${chr}.vcf
    fi
    """
}

process PDB_MERGE_PANEL {
    tag "FLARE-clean autosomes ${batch}"
    label 'pdb_light'
    publishDir "${params.outdir}/pedigree/panel_genomes", mode:'copy'
    input:
      tuple val(batch), val(chrs), path(vcfs)
    output:
      tuple val(batch), path("${batch}.autosomes.vcf.gz"), emit: genomes
    script:
    def ordered = (0..<chrs.size()).toList().sort { a,b -> (chrs[a] as int) <=> (chrs[b] as int) }
                           .collect { i -> "'${vcfs[i]}'" }.join(' ')
    """
    set -euo pipefail
    bcftools concat -Oz -o '${batch}.autosomes.vcf.gz' ${ordered}
    tabix -p vcf '${batch}.autosomes.vcf.gz'
    """
}

process PDB_AGGREGATE_TRUTH {
    tag 'realized multigeneration ancestry'
    label 'pdb_light'
    publishDir "${params.outdir}/pedigree/truth", mode:'copy'
    input:
      path truths
    output:
      path 'all_samples_ancestry.tsv', emit: truth
    script:
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/aggregate_truth.py' \\
      --inputs *.truth.tsv --chromosomes '${params.pedigree_chromosomes}' \\
      --output all_samples_ancestry.tsv
    """
}

process PDB_FLARE_INFER {
    tag "${rep_id} : ${batch} [clean]"
    label 'pdb_flare'
    publishDir "${params.outdir}/pedigree/features", mode:'copy'
    input:
      tuple val(rep_id), path(model_file), val(min_maf), val(min_mac), val(flare_seed),
            val(batch), path(genome)
    output:
      tuple val(rep_id), val(batch), path("${rep_id}__${batch}.features.tsv"), emit: features
    script:
    """
    set -euo pipefail
    java -Xmx32g -jar '${params.flare_jar}' \\
      gt='${genome}' ref='${params.flare_ref}' ref-panel='${params.flare_panel}' \\
      map='${params.genetic_map}' model='${model_file}' \\
      update-p=false em=false probs=false array=false \\
      min-maf=${min_maf} min-mac=${min_mac} nthreads=${task.cpus} seed=${flare_seed} out=run
    python3 '${projectDir}/scripts/extract_global_anc.py' \\
      --input run.global.anc.gz --output '${rep_id}__${batch}.features.tsv'
    """
}

process PDB_EVALUATE_FINALIST {
    tag "${candidate_id} [unified pedigree evaluation]"
    label 'pdb_decoder'
    publishDir "${params.outdir}/pedigree/evaluations", mode:'copy'
    input:
      tuple val(rep_id), val(candidate_id), path(candidate_dir), path(feature_files), path(truth)
      path founder_audit
      path plan_manifest
    output:
      path("${candidate_id}.pedigree_predictions.tsv"), emit: predictions
      path("${candidate_id}.pedigree_metrics.json"), emit: metrics
    script:
    """
    set -euo pipefail
    echo '${params.pedigree_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/evaluate_finalist.py' \\
      --model '${candidate_dir}/model.joblib' --truth '${truth}' \\
      --candidate '${candidate_id}' --representation '${rep_id}' \\
      --feature-files *.features.tsv --chromosomes '${params.pedigree_chromosomes}' \\
      --founder-audit '${founder_audit}' --plan-manifest '${plan_manifest}' \\
      --predictions '${candidate_id}.pedigree_predictions.tsv' \\
      --metrics '${candidate_id}.pedigree_metrics.json'
    """
}

process PDB_SELECT_UNIFIED_FINALISTS {
    tag 'unified ancestry-proportion finalist selection'
    label 'pdb_light'
    publishDir "${params.outdir}/final", mode:'copy', overwrite:true
    input:
      path candidate_dirs
      path pedigree_prediction_files
    output:
      path 'finalists.tsv', emit: table
      path 'finalists/*.json', emit: finalists
      path 'unified_leaderboard.tsv', emit: leaderboard
      path 'unified_diagnostics.json', emit: diagnostics
    script:
    """
    set -euo pipefail
    echo '${params.unified_fingerprint}' >/dev/null
    python3 '${projectDir}/pedigree_benchmark/scripts/select_unified_ancestry.py' \\
      --candidate-dirs candidate_*/ \\
      --ped-predictions *.pedigree_predictions.tsv \\
      --config '${params.unified_selection_config}' \\
      --outdir .
    """
}

process PDB_SUMMARIZE {
    tag 'pedigree validation leaderboard'
    label 'pdb_light'
    publishDir "${params.outdir}/pedigree/final", mode:'copy'
    input:
      path metric_files
    output:
      path 'pedigree_leaderboard.tsv', emit: leaderboard
    script:
    """
    set -euo pipefail
    python3 '${projectDir}/pedigree_benchmark/scripts/summarize_metrics.py' \\
      --metrics *.pedigree_metrics.json --output pedigree_leaderboard.tsv
    """
}

workflow {
    if (params.resume_stage2_manifest) {
        // Recovery requires completed results and published successful model/feature artifacts.
        // It DOES NOT run any initial model build, refinement, stage1 screen or stage2 fit.
        recovered = file(params.resume_stage2_manifest, checkIfExists:true)
        selected_stage2_file = file("${params.resume_source_outdir}/selection/selected_representations.tsv", checkIfExists:true)
        selected_id_ch = Channel.fromPath(selected_stage2_file)
            .splitCsv(header:true,sep:'\t')
            .map { r -> tuple(r.rep_id.toString(),true) }

        all_models = Channel.fromPath("${params.resume_source_outdir}/models/*.final.model",checkIfExists:true)
            .map { model ->
                def rep=model.name.replaceFirst(/\.final\.model$/, '')
                def meta=file("${params.resume_source_outdir}/models/${rep}.config.json",checkIfExists:true)
                def cfg=new groovy.json.JsonSlurper().parse(meta)
                if (cfg.rep_id.toString()!=rep) error "Recovery metadata mismatch for ${rep}"
                tuple(rep,model,(cfg.min_maf as double),(cfg.min_mac as int),(cfg.flare_seed as int),meta)
            }
        // Limit feature staging to successful representations in the published selection.
        recovery_dev_features = Channel.fromPath("${params.resume_source_outdir}/features/dev/*.features.tsv",checkIfExists:true)
            .map { f -> tuple(f.name.tokenize('__')[0],f) }
            .join(selected_id_ch,by:0)
            .map { rep,f,keep -> f }
        all_dev_features = recovery_dev_features.collect()
        SELECT_STAGE2_CANDIDATES(Channel.value(recovered))
    } else {
        panel_cfg_ch = Channel.fromPath(params.generated_panel_configs)
            .splitCsv(header:true,sep:'\t')
            .map { r -> tuple(r.panel_id,r.panel_cm as double,r.min_maf as double,r.min_mac as int) }
        PANEL_PROBS(panel_cfg_ch)
        model_cfg_ch = Channel.fromPath(params.generated_model_configs)
            .map { f ->
                def r=new groovy.json.JsonSlurper().parse(f)
                tuple(r.panel_id.toString(),r.rep_id.toString(),f,
                      r.min_maf as double,r.min_mac as int,r.flare_seed as int,
                      r.refinement_mode.toString(),r.refine_em as boolean)
            }
        BUILD_INITIAL_MODEL(model_cfg_ch.combine(PANEL_PROBS.out.panels,by:0))
        REFINE_MODEL(BUILD_INITIAL_MODEL.out.initial)
        fixed_model_ch = Channel.fromPath(params.fixed_models)
            .splitCsv(header:true,sep:'\t')
            .filter { r -> r.enabled.toString().toLowerCase()!='false' }
            .map { r -> tuple(r.rep_id,r.model_path,r.min_maf as double,r.min_mac as int,r.flare_seed as int,r.description) }
        IMPORT_FIXED_MODEL(fixed_model_ch)
        all_models = REFINE_MODEL.out.models.mix(IMPORT_FIXED_MODEL.out.models)
        dev_dataset_ch = Channel.fromPath(params.dev_datasets)
            .splitCsv(header:true,sep:'\t')
            .filter { r -> r.enabled.toString().toLowerCase()!='false' }
            .map { r -> tuple(r.dataset_id,r.role,r.input_kind,r.input_path,r.truth_path,r.truth_mode,
                              r.weight_tag,r.platform,r.array) }
        RUN_REPRESENTATION_DEV(dev_dataset_ch.combine(all_models))
        rep_feature_groups = RUN_REPRESENTATION_DEV.out.features.groupTuple(by:0)
        SCREEN_REPRESENTATION(rep_feature_groups)
        SELECT_REPRESENTATIONS(SCREEN_REPRESENTATION.out.result_dirs.map { rep_id,d -> d }.collect())
        GENERATE_STAGE2_CANDIDATES(SELECT_REPRESENTATIONS.out.selected)
        stage2_candidate_ch = GENERATE_STAGE2_CANDIDATES.out.candidates.flatten()
        all_dev_features = RUN_REPRESENTATION_DEV.out.features.map { rep_id,dataset_id,f -> f }.collect()
        FIT_STAGE2_CANDIDATE(stage2_candidate_ch,all_dev_features)
        // Collect paths AS TEXT. The selector only stages one small manifest,
        // never tens of thousands of symlinked candidate directories.
        stage2_manifest = FIT_STAGE2_CANDIDATE.out.result_dirs
            .map { d -> d.toString() }
            .collectFile(name:'stage2_result_dirs.txt',newLine:true,sort:true)
        SELECT_STAGE2_CANDIDATES(stage2_manifest)
        selected_id_ch = SELECT_REPRESENTATIONS.out.selected
            .splitCsv(header:true,sep:'\t')
            .map { r -> tuple(r.rep_id.toString(),true) }
    }
    selected_models = all_models.join(selected_id_ch,by:0)
        .map { rep_id,model_file,min_maf,min_mac,flare_seed,model_meta,keep ->
            tuple(rep_id,model_file,min_maf,min_mac,flare_seed,model_meta) }
    diagnostic_dataset_ch = Channel
      .fromPath(params.diagnostic_datasets)
      .splitCsv(header:true, sep:'\t')
      .filter { r -> r.enabled.toString().toLowerCase() != 'false' }
      .map { r -> tuple(r.dataset_id, r.role, r.input_kind, r.input_path, r.truth_path, r.truth_mode,
                        r.weight_tag, r.platform, r.array) }

    diagnostic_jobs = diagnostic_dataset_ch.combine(selected_models)
    RUN_REPRESENTATION_DIAGNOSTIC(diagnostic_jobs)

    diagnostic_features = RUN_REPRESENTATION_DIAGNOSTIC.out.features.map { rep_id, dataset_id, f -> f }
    if (params.resume_stage2_manifest) {
        stage3_features = all_dev_features.combine(diagnostic_features.collect())
            .map { dev,diag -> dev + diag }
    } else {
        stage3_features = RUN_REPRESENTATION_DEV.out.features.map { rep_id,dataset_id,f -> f }
            .mix(diagnostic_features).collect()
    }

    selected_stage2_candidate_ch = SELECT_STAGE2_CANDIDATES.out.candidates.flatten()
    FULL_EVALUATE_CANDIDATE(selected_stage2_candidate_ch, stage3_features)

    final_dirs_for_aggregate = FULL_EVALUATE_CANDIDATE.out.result_dirs.collect()
    AGGREGATE_FINAL(final_dirs_for_aggregate)

    // Retain the original v3 leaderboard for diagnostic and reproducibility only.
    // No model is promoted solely on the reference-only selection score when
    // unified pedigree selection is enabled.
    if (params.pedigree_enabled.toString().toBoolean()) {
        def pdbChrs = params.pedigree_chromosomes.toString().tokenize(',').collect { it.trim() as int }
        if (!pdbChrs || pdbChrs.unique().size() != pdbChrs.size() || !pdbChrs.every { it >= 1 && it <= 22 })
            error 'Invalid --pedigree_chromosomes (expected unique autosomes 1..22)'
        if ((params.pedigree_families as int) < 9 || (params.pedigree_batch_size as int) < 1)
            error 'Unified selection requires --pedigree_families >= 9 (all nine simulation classes) and positive batch size'
        if (!(params.pedigree_design.toString() in ['full', 'reduced']))
            error 'Invalid --pedigree_design (expected full or reduced)'

        // Pedigree simulation depends on successful stage-2 shortlisting,
        // not the previously reference-only SELECT_FINALISTS.
        PDB_PREPARE_FOUNDERS(
            file(params.pedigree_validation_meta, checkIfExists:true),
            file(params.pedigree_training_meta, checkIfExists:true),
            params.pedigree_additional_meta ? file(params.pedigree_additional_meta, checkIfExists:true) : [],
            file(params.pedigree_label_meta, checkIfExists:true),
            file(params.pedigree_mosaic_meta, checkIfExists:true),
            file(params.flare_ref, checkIfExists:true),
            file(params.flare_panel, checkIfExists:true),
            SELECT_STAGE2_CANDIDATES.out.table
        )
        PDB_BUILD_MAPS(file(params.pedigree_plink_map, checkIfExists:true),
                       SELECT_STAGE2_CANDIDATES.out.table)
        PDB_PANEL_SITES(file(params.flare_ref, checkIfExists:true),
                        SELECT_STAGE2_CANDIDATES.out.table)
        PDB_MAKE_PLANS(PDB_PREPARE_FOUNDERS.out.donors,
                       file(params.pedigree_training_meta, checkIfExists:true),
                       SELECT_STAGE2_CANDIDATES.out.table)

        pdb_batches = PDB_MAKE_PLANS.out.plans.flatMap { dir ->
            dir.toFile().listFiles().findAll { it.isDirectory() }.sort { it.name }
                .collect { folder -> tuple(folder.name, file(folder)) }
        }
        pdb_refs = Channel.fromList(pdbChrs.collect { chr ->
            def vcfPath = "${params.pedigree_vcf_dir}/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"
            tuple(chr, file(vcfPath, checkIfExists:true))
        })
        pdb_jobs = pdb_batches.combine(pdb_refs)
            .combine(PDB_BUILD_MAPS.out.maps)
            .combine(PDB_PANEL_SITES.out.sites)
            .map { batch, plan, chr, vcf, maps, sites -> tuple(batch, chr, plan, vcf, maps, sites) }
        PDB_SIMULATE(pdb_jobs)
        pdb_chrom_data = PDB_SIMULATE.out.panel_truth
        pdb_truth_ch = pdb_chrom_data.map { batch,chr,vcf,truth -> truth }.collect()
        PDB_AGGREGATE_TRUTH(pdb_truth_ch)
        pdb_for_merge = pdb_chrom_data.map { batch,chr,vcf,truth -> tuple(batch, chr, vcf) }
            .groupTuple()
        PDB_MERGE_PANEL(pdb_for_merge)

        // Every successful full-stage candidate participates in the final
        // ancestry-proportion selection, not just legacy finalists.
        pdb_stage3_configs = selected_stage2_candidate_ch.map { f ->
            def config = new groovy.json.JsonSlurper().parse(f)
            tuple(config.candidate_id.toString(), config.representation.toString())
        }
        pdb_model_dirs = FULL_EVALUATE_CANDIDATE.out.result_dirs.map { folder ->
            tuple(folder.baseName.toString().replaceFirst('^candidate_', ''), folder)
        }
        pdb_candidate_dirs = pdb_stage3_configs.join(pdb_model_dirs, by:0)
            .map { cid, rep, folder -> tuple(rep, cid, folder) }
        pdb_unique_rep_ids = pdb_stage3_configs.map { cid,rep -> tuple(rep, true) }.unique()
        pdb_rep_models = pdb_unique_rep_ids.join(all_models, by:0)
            .map { rep, keep, model_file, min_maf, min_mac, flare_seed, model_meta ->
                tuple(rep, model_file, min_maf, min_mac, flare_seed)
            }
        pdb_inference_jobs = pdb_rep_models.combine(PDB_MERGE_PANEL.out.genomes)
            .map { rep, model_file, min_maf, min_mac, flare_seed, batch, genome ->
                tuple(rep, model_file, min_maf, min_mac, flare_seed, batch, genome)
            }
        PDB_FLARE_INFER(pdb_inference_jobs)
        pdb_rep_features = PDB_FLARE_INFER.out.features
            .map { rep,batch,f -> tuple(rep,f) }.groupTuple()
        pdb_evaluation_jobs = pdb_candidate_dirs.combine(pdb_rep_features, by:0)
            .combine(PDB_AGGREGATE_TRUTH.out.truth)
            .map { rep, cid, folder, features, truth -> tuple(rep,cid,folder,features,truth) }
        PDB_EVALUATE_FINALIST(pdb_evaluation_jobs, PDB_PREPARE_FOUNDERS.out.summary,
                             PDB_MAKE_PLANS.out.plans.map { dir -> file("${dir}/manifest.json", checkIfExists:true) })
        PDB_SUMMARIZE(PDB_EVALUATE_FINALIST.out.metrics.collect())

        PDB_SELECT_UNIFIED_FINALISTS(
            FULL_EVALUATE_CANDIDATE.out.result_dirs.collect(),
            PDB_EVALUATE_FINALIST.out.predictions.collect()
        )
    } else {
        // Legacy v3 selection available explicitly for comparison only.
        SELECT_FINALISTS(FULL_EVALUATE_CANDIDATE.out.result_dirs.collect())
    }
}
