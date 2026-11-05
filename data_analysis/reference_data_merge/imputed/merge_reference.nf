nextflow.enable.dsl = 2

/*
 * Build imputed FLARE reference v2.
 *
 * 4 reference datasets:
 *   HGDP+1kGP
 *   GSE77508 (QC)
 *   GSE73996 (QC)
 *   Kazakh   (QC)
 *
 * 22 chromosome tasks.
 * Maximum 11 concurrent tasks.
 * 2 CPUs / 14 GB per task.
 */

params.hgdp_dir =
    "/srv/dependencies/ancestry/hgdp_1kg_references"

params.imputed_dir =
    "${projectDir}/output"

params.previous_ref =
    "/srv/dependencies/ancestry/flare_reference_v1/output"

params.out_dir =
    "${projectDir}/merged_reference/output"

// SNP selection parameters
params.min_maf = 0.01
params.min_spacing_bp = 10000


process MERGE_REFERENCE {

    tag "chr${chr}"

    executor 'local'

    cpus 2
    memory '14 GB'
    maxForks 11

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chr),
          path(hgdp),
          path(hgdp_idx),
          path(gse77508),
          path(gse77508_idx),
          path(gse73996),
          path(gse73996_idx),
          path(kazakh),
          path(kazakh_idx),
          path(original_names),
          path(sanitized_names),
          path(flare_panel)

    output:
    tuple val(chr),
          path("flare_reference.chr${chr}.vcf.gz"),
          path("flare_reference.chr${chr}.vcf.gz.csi"),
          path("backbone.chr${chr}.vcf.gz"),
          path("backbone.chr${chr}.vcf.gz.csi"),
          path("validation.chr${chr}.tsv"),
          emit: reference

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo "Building imputed reference chr${chr}"
    echo "========================================"

    # -------------------------------------------------------
    # 1. Four-way exact allele intersection
    #
    # Calculate MAF from original HGDP genotypes.
    # Remove variants with MAF below threshold.
    # Thin to one SNP per minimum physical spacing.
    #
    # Output is a sites-only VCF.
    # -------------------------------------------------------

    bcftools isec \\
        -n=4 \\
        -c none \\
        -w1 \\
        ${hgdp} \\
        ${gse77508} \\
        ${gse73996} \\
        ${kazakh} \\
        -Ou \\
    | bcftools +fill-tags -Ou -- -t MAF \\
    | bcftools view \\
        -i 'INFO/MAF>=${params.min_maf}' \\
        -Ou \\
    | bcftools view -G -Ov \\
    | awk -v gap=${params.min_spacing_bp} '
        /^#/ {print; next}
        {
            if (\$1 != chrom || \$2 - last >= gap) {
                print
                chrom = \$1
                last = \$2
            }
        }
      ' \\
    | bcftools view \\
        -Oz \\
        -o backbone.chr${chr}.vcf.gz

    bcftools index -f -c backbone.chr${chr}.vcf.gz

    BACKBONE=\$(bcftools index -n backbone.chr${chr}.vcf.gz)

    if [ "\$BACKBONE" -eq 0 ]; then
        echo "ERROR: No backbone SNPs on chr${chr}"
        exit 1
    fi

    echo "Backbone SNPs: \$BACKBONE"

    # -------------------------------------------------------
    # 2. Extract exact variants from each reference dataset
    #
    # isec -c none ensures REF/ALT match the backbone.
    # -------------------------------------------------------

    bcftools isec \\
        -n=2 -c none -w1 \\
        ${hgdp} backbone.chr${chr}.vcf.gz \\
        -Oz -o hgdp.subset.vcf.gz

    bcftools isec \\
        -n=2 -c none -w1 \\
        ${gse77508} backbone.chr${chr}.vcf.gz \\
        -Oz -o gse77508.subset.vcf.gz

    bcftools isec \\
        -n=2 -c none -w1 \\
        ${gse73996} backbone.chr${chr}.vcf.gz \\
        -Oz -o gse73996.subset.vcf.gz

    bcftools isec \\
        -n=2 -c none -w1 \\
        ${kazakh} backbone.chr${chr}.vcf.gz \\
        -Oz -o kazakh.subset.vcf.gz

    for file in \\
        hgdp.subset.vcf.gz \\
        gse77508.subset.vcf.gz \\
        gse73996.subset.vcf.gz \\
        kazakh.subset.vcf.gz
    do
        bcftools index -f -c "\$file"

        N=\$(bcftools index -n "\$file")

        if [ "\$N" -ne "\$BACKBONE" ]; then
            echo "ERROR: \$file has \$N markers, expected \$BACKBONE"
            exit 1
        fi
    done

    # -------------------------------------------------------
    # 3. Merge samples across all four datasets
    # -------------------------------------------------------

    bcftools merge \\
        --threads ${task.cpus} \\
        -m none \\
        hgdp.subset.vcf.gz \\
        gse77508.subset.vcf.gz \\
        gse73996.subset.vcf.gz \\
        kazakh.subset.vcf.gz \\
        -Oz \\
        -o merged.raw.vcf.gz

    # -------------------------------------------------------
    # 4. Verify original sample ordering
    #
    # We can safely reuse the previous sanitized names
    # only if sample identity and order have not changed.
    # -------------------------------------------------------

    bcftools query -l merged.raw.vcf.gz > current_names.txt

    if ! cmp -s current_names.txt ${original_names}; then
        echo "ERROR: Sample names/order differ from previous reference."
        diff -u ${original_names} current_names.txt || true
        exit 1
    fi

    # -------------------------------------------------------
    # 5. Apply sanitized sample IDs
    # -------------------------------------------------------

    bcftools reheader \\
        -s ${sanitized_names} \\
        -o flare_reference.chr${chr}.vcf.gz \\
        merged.raw.vcf.gz

    bcftools index -f -c \\
        flare_reference.chr${chr}.vcf.gz

    # -------------------------------------------------------
    # 6. Validation
    # -------------------------------------------------------

    FINAL_MARKERS=\$(bcftools index -n \\
        flare_reference.chr${chr}.vcf.gz)

    SAMPLES=\$(bcftools query -l \\
        flare_reference.chr${chr}.vcf.gz | wc -l)

    test "\$FINAL_MARKERS" -eq "\$BACKBONE"
    test "\$SAMPLES" -eq 4886

    # Validate every sample appears exactly once in panel.
    bcftools query -l flare_reference.chr${chr}.vcf.gz \\
        | LC_ALL=C sort > merged.samples.sorted

    awk '{print \$1}' ${flare_panel} \\
        | LC_ALL=C sort > panel.samples.sorted

    cmp merged.samples.sorted panel.samples.sorted

    # Check for missing or unphased genotypes.
    # Exit immediately if any invalid genotype is found.
    bcftools query \\
        -f '[%GT\\t]\\n' \\
        flare_reference.chr${chr}.vcf.gz \\
    | awk '
        index(\$0, "/") || index(\$0, ".") {
            print "ERROR: missing/unphased genotype" > "/dev/stderr"
            exit 1
        }
    '

    printf "chromosome\\tmarkers\\tsamples\\tstatus\\n" \\
        > validation.chr${chr}.tsv

    printf "chr${chr}\\t%s\\t%s\\tOK\\n" \\
        "\$FINAL_MARKERS" "\$SAMPLES" \\
        >> validation.chr${chr}.tsv

    echo "chr${chr}: \$FINAL_MARKERS markers, \$SAMPLES samples, OK"
    """
}


process SUMMARY {

    cpus 1
    memory '1 GB'

    publishDir "${params.out_dir}",
        mode: 'copy',
        overwrite: true

    input:
    path reports

    output:
    path "validation_summary.tsv"

    script:
    """
    printf "chromosome\\tmarkers\\tsamples\\tstatus\\n" \\
        > validation_summary.tsv

    for file in validation.chr*.tsv; do
        tail -n +2 "\$file"
    done | sort -V -k1,1 >> validation_summary.tsv
    """
}


workflow {

    chromosome_inputs = Channel
        .fromList((1..22).toList())
        .map { chr ->

            def hgdp = file(
                "${params.hgdp_dir}/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz",
                checkIfExists: true
            )

            def a = file(
                "${params.imputed_dir}/gse77508/qc/gse77508.chr${chr}.qc.vcf.gz",
                checkIfExists: true
            )

            def b = file(
                "${params.imputed_dir}/gse73996/qc/gse73996.chr${chr}.qc.vcf.gz",
                checkIfExists: true
            )

            def c = file(
                "${params.imputed_dir}/kazakh/qc/kazakh.chr${chr}.qc.vcf.gz",
                checkIfExists: true
            )

            tuple(
                chr,

                hgdp,
                file("${hgdp}.tbi", checkIfExists: true),

                a,
                file("${a}.csi", checkIfExists: true),

                b,
                file("${b}.csi", checkIfExists: true),

                c,
                file("${c}.csi", checkIfExists: true),

                file(
                    "${params.previous_ref}/original_sample_names.txt",
                    checkIfExists: true
                ),

                file(
                    "${params.previous_ref}/sanitized/sample_names.txt",
                    checkIfExists: true
                ),

                file(
                    "${params.previous_ref}/sanitized/flare_reference.ref.panel",
                    checkIfExists: true
                )
            )
        }

    MERGE_REFERENCE(chromosome_inputs)

    reports = MERGE_REFERENCE.out.reference
        .map { chr, vcf, idx, backbone, backbone_idx, report ->
            report
        }
        .collect()

    SUMMARY(reports)
}
