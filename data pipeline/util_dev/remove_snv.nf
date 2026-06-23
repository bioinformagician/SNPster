nextflow.enable.dsl = 2

params.input_dir  = "/srv/dependencies/imputation_runner/imputer/beagle_references"
params.output_dir = "${params.input_dir}/snv_indel_biallelic"

process FILTER_REFERENCE {

    tag { vcf_file.simpleName }

    publishDir params.output_dir, mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path "*.SNV_INDEL_biallelic.vcf.gz"
    path "*.SNV_INDEL_biallelic.vcf.gz.tbi"

    script:
    """
    set -euo pipefail

    base=\$(basename "${vcf_file}" .vcf.gz)

    # Avoid duplicating suffix if rerun on already-filtered files
    base=\${base/.SNV_INDEL_biallelic/}

    output="\${base}.SNV_INDEL_biallelic.vcf.gz"

    bcftools view \\
      -m2 -M2 \\
      -v snps,indels \\
      "${vcf_file}" \\
      -Oz -o "\$output"

    tabix -f -p vcf "\$output"
    """
}

workflow {

    vcf_ch = Channel
        .fromPath("${params.input_dir}/*.vcf.gz")
        .filter { it.name.contains("numericCHR") }
        .filter { !it.name.contains("SNV_INDEL_biallelic") }

    FILTER_REFERENCE(vcf_ch)
}