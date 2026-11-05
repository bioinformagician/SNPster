#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 || "$1" != --run ]]; then
    echo "Usage: bash run_dense_public_support.sh --run" >&2; exit 1
fi
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(cd "$project/../../.." && pwd)"
active="$(docker ps --filter name=snpster-ancestry-v2- --format '{{.Names}}')"
if [[ -n "$active" ]]; then
    printf 'Refusing a concurrent ancestry controller:\n%s\n' "$active" >&2; exit 1
fi
suffix="$(date +%Y%m%d_%H%M%S)"
run="$project/results_v3_dense_public_support_$suffix"
name="snpster-ancestry-v2-dense-support-$suffix"
mkdir "$run"
mkdir -p "$run/smoke/launch" "$run/full/launch" "$run/controller_home"
sha256sum "$project/dense_public_support.nf" "$project/dense_public_support.config" \
    "$project/production_matched.config" "$project/run_dense_public_support.sh" \
    "$project/config/dense_public_support.json" "$project/config/sgdp_public_import.json" \
    "$project/config/sgdp_normalization.json" \
    "$project/pedigree_benchmark/scripts/prepare_dense_public_support.py" \
    "$project/pedigree_benchmark/scripts/import_sgdp_public.py" \
    "$project/pedigree_benchmark/scripts/normalize_sgdp.py" \
    "$project/pedigree_benchmark/scripts/recover_sgdp_duplicates.py" \
    "$project/pedigree_benchmark/scripts/verify_public_sgdp.py" \
    "$project/pedigree_benchmark/scripts/verify_sgdp_full_manifest.py" \
    "$project/pedigree_benchmark/scripts/sgdp_genetic_qc.py" \
    "$project/pedigree_benchmark/scripts/audit_production_ancestry.py" \
    "$project/results_v3_sgdp_normalization_full_v2_20261022_095828/inputs/hs37d5.fa" \
    "$project/results_v3_sgdp_normalization_full_v2_20261022_095828/inputs/hs37d5.fa.fai" \
    "$project/results_v3_sgdp_normalization_full_v2_20261022_095828/summary/normalization_complete.json" \
    "$project/results_v3_public_sgdp_smoke_20261021_203540/metadata/metadata.original.txt" \
    "$project/results_v3_public_sgdp_smoke_20261021_203540/metadata/publisher_md5sums.txt" \
    "$project/results_v3_public_sgdp_smoke_20261021_203540/raw_control/raw_control.vcf.gz" \
    "$project/results_v3_production_preflight_20261021_193442/summary/preflight.json" \
    "$project/results_v3_production_cohort_20261022_110441/cohort/frozen/cohort.json" \
    "$project/results_v3_model_family_screen_05/inputs/reference_labels.tsv" \
    "$project/results_v3_model_family_screen_05/inputs/pure_train_labels.tsv" \
    "$project/results_v3_phasefree_array_validation_01/lonni.observed_keys.tsv" \
    "$project/results_v3_phasefree_array_validation_01/frederik.observed_keys.tsv" \
    "$repository/data pipeline/standardizer_module/dependencies/hg19ToHg38.over.chain.gz" \
    /srv/dependencies/imputation_runner/harmonizer/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    /srv/dependencies/imputation_runner/harmonizer/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
    "$HOME/.nextflow/framework/26.04.4/nextflow-26.04.4-one.jar" > "$run/source_sha256.txt"
printf '%s\n' "$name" > "$run/controller.name"
docker run --detach -t --name "$name" --init --network none --cpus=4 --memory=4g \
    --user "$(id -u):$(id -g)" --group-add "$(stat -c '%g' /var/run/docker.sock)" \
    -e NXF_HOME="$run/controller_home" -e NXF_VER=26.04.4 -e NXF_OFFLINE=true \
    -e NXF_DISABLE_CHECK_LATEST=true -e PYTHONDONTWRITEBYTECODE=1 \
    -e PROJECT="$project" -e REPOSITORY="$repository" -e EXPERIMENT="$run" \
    -v "$repository:$repository:ro" -v "$run:$run:rw" -v /srv/dependencies:/srv/dependencies:ro \
    -v "$HOME/.nextflow/framework/26.04.4:$run/controller_home/framework/26.04.4:ro" \
    -v /var/run/docker.sock:/var/run/docker.sock -w "$run" \
    --entrypoint /bin/bash sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1 \
    -euo pipefail -c '
        stty cols 180 rows 60
        options="--user $(id -u):$(id -g) -v /srv/dependencies:/srv/dependencies:ro"
        for folder in pedigree_benchmark/scripts config \
            results_v3_sgdp_normalization_full_v2_20261022_095828 \
            results_v3_public_sgdp_smoke_20261021_203540 \
            results_v3_production_preflight_20261021_193442 \
            results_v3_production_cohort_20261022_110441 \
            results_v3_model_family_screen_05 \
            results_v3_phasefree_array_validation_01; do
            options="$options -v $PROJECT/$folder:$PROJECT/$folder:ro"
        done
        for mode in smoke full; do
            run="$EXPERIMENT/$mode"
            if [[ "$mode" == smoke ]]; then smoke=true; else smoke=false; fi
            cd "$run/launch"
            echo "Recovering dense observed support ($mode), never release-imputed GT/DS/GP; production strategy unchanged."
            nextflow -log "$run/nextflow.log" -C "$PROJECT/dense_public_support.config" \
                run "$PROJECT/dense_public_support.nf" -ansi-log true -work-dir "$run/work" \
                -with-trace "$run/trace.tsv" -with-report "$run/report.html" \
                -with-timeline "$run/timeline.html" -with-dag "$run/pipeline.html" \
                --project "$PROJECT" --run_dir "$run" --smoke "$smoke" \
                --container_options "$options"
            python3 "$PROJECT/pedigree_benchmark/scripts/prepare_dense_public_support.py" finalize --run "$run"
            session="$(awk -F "\t" "END {print \$6}" .nextflow/history)"
            [[ "$session" =~ ^[a-f0-9-]{36}$ ]]
            nextflow clean -n "$session" > "$run/cache_cleanup.plan"
            if ! grep -q "^Would remove " "$run/cache_cleanup.plan" ||
                grep "^Would remove " "$run/cache_cleanup.plan" | grep -v -F "Would remove $run/work/" >/dev/null; then
                echo "Missing or unsafe exact-session cleanup plan" >&2; exit 1
            fi
            nextflow clean -f "$session" > "$run/cache_cleanup.log"
            find "$run/work" -depth -mindepth 1 -type d -empty -delete
            python3 "$PROJECT/pedigree_benchmark/scripts/prepare_dense_public_support.py" finalize --run "$run" --post-cleanup
        done
        echo "Dense observed support completed; all-autosome genetic role/family screening and new cohort freeze remain required before ancestry testing."
    ' > "$run/controller.id"
echo "Controller: $name"
echo "Monitor: docker logs -f $name"
echo "Output: $run"
docker inspect --format '{{.Name}} {{.State.Status}}' "$name"
