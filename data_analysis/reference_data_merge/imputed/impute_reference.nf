nextflow.enable.dsl = 2

/*
 * Phase + impute three external reference datasets
 * against HGDP+1kGP GRCh38.
 *
 * Datasets:
 *   GSE77508    196 samples
 *   GSE73996    381 samples
 *   Kazakh      218 samples
 *
 * 66 tasks total.
 * Maximum 11 concurrent tasks.
 * 2 CPUs and 14 GB RAM per task.
 */

params.gse77508_dir =
    "/srv/dependencies/ancestry/gse77508"

params.gse73996_dir =
    "/srv/dependencies/ancestry/harmonized/gse73996/selected/phased"

params.kazakh_dir =
    "/srv/dependencies/ancestry/harmonized/kazakh/vcf/phased"

params.ref_dir =
    "/srv/dependencies/ancestry/hgdp_1kg_references"

params.map_dir =
    "/srv/dependencies/imputation_runner/imputer/plink.GRCh38.map"

params.out_dir =
    "${projectDir}/output"

// Test chr1 on all three datasets before running everything.
params.test_chr1 = false


process IMPUTE_REFERENCE {

    tag "${dataset} chr${chr}"

    executor 'local'

    cpus 2
    memory '14 GB'
    maxForks 11

    errorStrategy 'retry'
    maxRetries 2

    publishDir { "${params.out_dir}/${dataset}" },
        mode: 'copy',
        overwrite: true

    input:

    tuple val(dataset),
          val(chr),
          path(gt),
          path(gt_index),
          path(ref),
          path(ref_index),
          path(genetic_map)

    output:

    tuple val(dataset),
          val(chr),
          path("${dataset}.chr${chr}.imputed.vcf.gz"),
          path("${dataset}.chr${chr}.imputed.vcf.gz.csi"),
          path("${dataset}.chr${chr}.imputed.log"),
          path("${dataset}.chr${chr}.markers.txt"),
          emit: imputed

    script:

    """
    set -euo pipefail

    # Limit the Java heap to 12 GB, leaving ~2 GB
    # for other memory usage within the 14 GB allocation.
    export _JAVA_OPTIONS="-Xmx12g"

    beagle \\
        gt=${gt} \\
        ref=${ref} \\
        map=${genetic_map} \\
        impute=true \\
        nthreads=${task.cpus} \\
        out=${dataset}.chr${chr}.imputed

    # Index the imputed VCF.
    bcftools index -f \\
        ${dataset}.chr${chr}.imputed.vcf.gz

    # Record the number of resulting markers.
    bcftools index -n \\
        ${dataset}.chr${chr}.imputed.vcf.gz \\
        > ${dataset}.chr${chr}.markers.txt

    echo "Completed ${dataset} chr${chr}"
    """
}


workflow {

    datasets = [
        [
            name: 'gse77508',
            dir: params.gse77508_dir,
            prefix: 'GSE77508'
        ],
        [
            name: 'gse73996',
            dir: params.gse73996_dir,
            prefix: 'gse73996'
        ],
        [
            name: 'kazakh',
            dir: params.kazakh_dir,
            prefix: 'kazakh'
        ]
    ]

    chromosomes = params.test_chr1
        ? [1]
        : (1..22).toList()

    inputs = Channel
        .fromList(datasets)
        .flatMap { dataset ->

            chromosomes.collect { chr ->

                def gt = file(
                    "${dataset.dir}/${dataset.prefix}.chr${chr}.phased.vcf.gz",
                    checkIfExists: true
                )

                // Accept either CSI or TBI input indexes.
                def gt_index = file(
                    "${gt}.csi"
                )

                if (!gt_index.toFile().exists()) {
                    gt_index = file(
                        "${gt}.tbi",
                        checkIfExists: true
                    )
                }

                def ref = file(
                    "${params.ref_dir}/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz",
                    checkIfExists: true
                )

                def ref_index = file(
                    "${ref}.tbi",
                    checkIfExists: true
                )

                def genetic_map = file(
                    "${params.map_dir}/plink.chr${chr}.GRCh38.map",
                    checkIfExists: true
                )

                tuple(
                    dataset.name,
                    chr,
                    gt,
                    gt_index,
                    ref,
                    ref_index,
                    genetic_map
                )
            }
        }

    IMPUTE_REFERENCE(inputs)

    QC_REFERENCE(IMPUTE_REFERENCE.out.imputed)
}


/*
 * Post-imputation QC.
 *
 * Keep:
 *   - observed SNPs
 *   - imputed SNPs with DR2 >= 0.8
 *
 * Remove:
 *   - imputed SNPs with DR2 < 0.8 or missing DR2
 *   - indels
 *   - multiallelic sites
 */

process QC_REFERENCE {

    tag "${dataset} chr${chr}"

    executor 'local'

    cpus 2
    memory '3 GB'
    maxForks 11

    errorStrategy 'retry'
    maxRetries 2

    publishDir { "${params.out_dir}/${dataset}/qc" },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(dataset),
          val(chr),
          path(imputed_vcf),
          path(imputed_index),
          path(beagle_log),
          path(marker_count)

    output:
    tuple val(dataset),
          val(chr),
          path("${dataset}.chr${chr}.qc.vcf.gz"),
          path("${dataset}.chr${chr}.qc.vcf.gz.csi"),
          path("${dataset}.chr${chr}.qc.tsv"),
          emit: filtered

    script:
    """
    set -euo pipefail

    PREFIX=${dataset}.chr${chr}

    echo "QC: ${dataset} chr${chr}"

    # Check the expected Beagle annotations.
    bcftools view -h ${imputed_vcf} |
        grep -q '^##INFO=<ID=DR2,'

    bcftools view -h ${imputed_vcf} |
        grep -q '^##INFO=<ID=IMP,'

    # Keep observed genotypes OR high-quality imputed variants.
    # Restrict output to biallelic SNPs.
    bcftools view \\
        -m2 -M2 \\
        -v snps \\
        -i 'INFO/IMP=0 || (INFO/IMP=1 && INFO/DR2>=0.8)' \\
        --threads ${task.cpus} \\
        -Oz \\
        -o \${PREFIX}.qc.vcf.gz \\
        ${imputed_vcf}

    bcftools index -f -c \${PREFIX}.qc.vcf.gz

    RAW=\$(bcftools index -n ${imputed_vcf})
    KEPT=\$(bcftools index -n \${PREFIX}.qc.vcf.gz)

    OBSERVED=\$(bcftools view \\
        -H -i 'INFO/IMP=0' \\
        \${PREFIX}.qc.vcf.gz | wc -l)

    IMPUTED=\$(bcftools view \\
        -H -i 'INFO/IMP=1' \\
        \${PREFIX}.qc.vcf.gz | wc -l)

    REMOVED=\$((RAW - KEPT))

    SAMPLES=\$(bcftools query -l \\
        \${PREFIX}.qc.vcf.gz | wc -l)

    if [ "\$KEPT" -eq 0 ]; then
        echo "ERROR: QC removed all variants."
        exit 1
    fi

    printf "dataset\\tchromosome\\traw_variants\\tqc_variants\\tobserved_kept\\timputed_kept\\tremoved\\tsamples\\n" \\
        > \${PREFIX}.qc.tsv

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${dataset}" "${chr}" \\
        "\$RAW" "\$KEPT" \\
        "\$OBSERVED" "\$IMPUTED" \\
        "\$REMOVED" "\$SAMPLES" \\
        >> \${PREFIX}.qc.tsv

    echo "Raw variants: \$RAW"
    echo "QC variants:  \$KEPT"
    echo "Observed:     \$OBSERVED"
    echo "Imputed:      \$IMPUTED"
    echo "Removed:      \$REMOVED"
    """
}
