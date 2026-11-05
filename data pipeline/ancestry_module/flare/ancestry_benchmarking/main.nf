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
params.outdir                 = params.outdir                 ?: '/workspace/results_v3'
params.flare_jar              = params.flare_jar              ?: '/opt/flare/flare.jar'
params.flare_ref              = params.flare_ref              ?: '/data/flare2_test/source_model/source_reference.autosomes.vcf.gz'
params.flare_panel            = params.flare_panel            ?: '/data/flare2_test/source_model/source_reference.ref.panel'
params.source_train_vcf       = params.source_train_vcf       ?: '/data/flare2_test/source_model/source_train.autosomes.vcf.gz'
params.genetic_map            = params.genetic_map            ?: '/data/flare2_test/model_sweep_allchr/plink.autosomes.GRCh38.map'
params.max_stage2_candidates  = params.max_stage2_candidates  ?: 0
params.config_fingerprint      = params.config_fingerprint      ?: 'unversioned'
params.resume_stage2_manifest  = params.resume_stage2_manifest  ?: null

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
    cp "${config_json}" "${rep_id}.config.json"
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
    path result_manifest

    output:
    path 'selected_stage2/*.json', emit: candidates
    path 'selected_stage2.tsv', emit: table

    script:
    """
    set -euo pipefail
    echo "config_fingerprint=${params.config_fingerprint}" >/dev/null
    mkdir -p selected_stage2
    python3 "${projectDir}/scripts/select_candidates.py" --manifest "${result_manifest}" --search-space "${params.search_space}" \\
      --outdir selected_stage2 --table selected_stage2.tsv
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
    path 'finalists.tsv'
    path 'finalists/*.json'

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

workflow {
    if (params.resume_stage2_manifest) {
        log.info "Resuming at SELECT_STAGE2_CANDIDATES using ${params.resume_stage2_manifest}"
        stage2_result_manifest = Channel.value(file(params.resume_stage2_manifest, checkIfExists: true))
        selected_reps = Channel.value(file("${params.outdir}/selection/selected_representations.tsv", checkIfExists: true))
        recovered_rep_ids = selected_reps
          .splitCsv(header:true, sep:'\t')
          .map { r -> r.rep_id }
        all_models = recovered_rep_ids.map { rep_id ->
            def model_file = file("${params.outdir}/models/${rep_id}.final.model", checkIfExists: true)
            def model_meta = file("${params.outdir}/models/${rep_id}.config.json", checkIfExists: true)
            def r = new groovy.json.JsonSlurper().parse(model_meta)
            if (r.rep_id.toString() != rep_id) {
                error "Recovered model metadata does not match representation ${rep_id}"
            }
            tuple(rep_id, model_file, r.min_maf as double, r.min_mac as int, r.flare_seed as int, model_meta)
        }
        recovered_dev_datasets = Channel
          .fromPath(params.dev_datasets, checkIfExists: true)
          .splitCsv(header:true, sep:'\t')
          .filter { r -> r.enabled.toString().toLowerCase() != 'false' }
          .map { r -> r.dataset_id }
        dev_features = recovered_rep_ids.combine(recovered_dev_datasets)
          .map { rep_id, dataset_id ->
              file("${params.outdir}/features/dev/${rep_id}__${dataset_id}.features.tsv", checkIfExists: true)
          }
    } else {
        panel_cfg_ch = Channel
          .fromPath(params.generated_panel_configs)
          .splitCsv(header:true, sep:'\t')
          .map { r -> tuple(r.panel_id, r.panel_cm as double, r.min_maf as double, r.min_mac as int) }

        PANEL_PROBS(panel_cfg_ch)

        model_cfg_ch = Channel
          .fromPath(params.generated_model_configs)
          .map { f ->
              def r = new groovy.json.JsonSlurper().parse(f)
              tuple(r.panel_id.toString(), r.rep_id.toString(), f,
                    r.min_maf as double, r.min_mac as int, r.flare_seed as int,
                    r.refinement_mode.toString(), r.refine_em as boolean)
          }

        build_jobs = model_cfg_ch.combine(PANEL_PROBS.out.panels, by:0)
        BUILD_INITIAL_MODEL(build_jobs)
        REFINE_MODEL(BUILD_INITIAL_MODEL.out.initial)

        fixed_model_ch = Channel
          .fromPath(params.fixed_models)
          .splitCsv(header:true, sep:'\t')
          .filter { r -> r.enabled.toString().toLowerCase() != 'false' }
          .map { r -> tuple(r.rep_id, r.model_path, r.min_maf as double, r.min_mac as int, r.flare_seed as int, r.description) }
        IMPORT_FIXED_MODEL(fixed_model_ch)

        all_models = REFINE_MODEL.out.models.mix(IMPORT_FIXED_MODEL.out.models)

        dev_dataset_ch = Channel
          .fromPath(params.dev_datasets)
          .splitCsv(header:true, sep:'\t')
          .filter { r -> r.enabled.toString().toLowerCase() != 'false' }
          .map { r -> tuple(r.dataset_id, r.role, r.input_kind, r.input_path, r.truth_path, r.truth_mode,
                            r.weight_tag, r.platform, r.array) }

        dev_jobs = dev_dataset_ch.combine(all_models)
        RUN_REPRESENTATION_DEV(dev_jobs)

        rep_feature_groups = RUN_REPRESENTATION_DEV.out.features.groupTuple(by:0)
        SCREEN_REPRESENTATION(rep_feature_groups)

        screen_dirs = SCREEN_REPRESENTATION.out.result_dirs.map { rep_id, d -> d }.collect()
        SELECT_REPRESENTATIONS(screen_dirs)

        GENERATE_STAGE2_CANDIDATES(SELECT_REPRESENTATIONS.out.selected)
        stage2_candidate_ch = GENERATE_STAGE2_CANDIDATES.out.candidates.flatten()
        all_dev_features = RUN_REPRESENTATION_DEV.out.features.map { rep_id, dataset_id, f -> f }.collect()
        FIT_STAGE2_CANDIDATE(stage2_candidate_ch, all_dev_features)

        stage2_result_manifest = FIT_STAGE2_CANDIDATE.out.result_dirs
          .map { result_dir -> result_dir.toString() }
          .collectFile(name: 'stage2_result_dirs.txt', newLine: true, sort: false)
        selected_reps = SELECT_REPRESENTATIONS.out.selected
        dev_features = RUN_REPRESENTATION_DEV.out.features.map { rep_id, dataset_id, f -> f }
    }

    SELECT_STAGE2_CANDIDATES(stage2_result_manifest)

    selected_id_ch = selected_reps
      .splitCsv(header:true, sep:'\t')
      .map { r -> tuple(r.rep_id, true) }

    selected_models = all_models.join(selected_id_ch, by:0)
      .map { rep_id, model_file, min_maf, min_mac, flare_seed, model_meta, keep ->
          tuple(rep_id, model_file, min_maf, min_mac, flare_seed, model_meta) }

    diagnostic_dataset_ch = Channel
      .fromPath(params.diagnostic_datasets)
      .splitCsv(header:true, sep:'\t')
      .filter { r -> r.enabled.toString().toLowerCase() != 'false' }
      .map { r -> tuple(r.dataset_id, r.role, r.input_kind, r.input_path, r.truth_path, r.truth_mode,
                        r.weight_tag, r.platform, r.array) }

    diagnostic_jobs = diagnostic_dataset_ch.combine(selected_models)
    RUN_REPRESENTATION_DIAGNOSTIC(diagnostic_jobs)

    diagnostic_features = RUN_REPRESENTATION_DIAGNOSTIC.out.features.map { rep_id, dataset_id, f -> f }
    stage3_features = dev_features.mix(diagnostic_features)
      .collect()

    selected_stage2_candidate_ch = SELECT_STAGE2_CANDIDATES.out.candidates.flatten()
    FULL_EVALUATE_CANDIDATE(selected_stage2_candidate_ch, stage3_features)

    final_dirs_for_aggregate = FULL_EVALUATE_CANDIDATE.out.result_dirs.collect()
    AGGREGATE_FINAL(final_dirs_for_aggregate)

    final_dirs_for_selection = FULL_EVALUATE_CANDIDATE.out.result_dirs.collect()
    SELECT_FINALISTS(final_dirs_for_selection)
}
