nextflow.enable.dsl = 2

params.input_dir = "/srv/dependencies/imputation_runner/imputer/beagle_references"
params.output_dir = "${params.input_dir}/numeric_chroms"

process RENAME_CHROMS {

    tag { vcf_file.simpleName }

    publishDir params.output_dir, mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path "*.numericCHR.vcf.gz"
    path "*.numericCHR.vcf.gz.tbi"

    script:
    """
    set -euo pipefail

    for chr in {1..22}; do
      echo -e "chr\${chr}\\t\${chr}"
    done > rename_from_chr.txt

    base=\$(basename "${vcf_file}" .vcf.gz)
    output="\${base}.numericCHR.vcf.gz"

    bcftools annotate \\
      --rename-chrs rename_from_chr.txt \\
      -Oz \\
      -o "\$output" \\
      "${vcf_file}"

    tabix -f -p vcf "\$output"
    """
}

workflow {

    vcf_ch = Channel
        .fromPath("${params.input_dir}/*.vcf.gz")
        .filter { !it.name.endsWith(".tbi") }

    RENAME_CHROMS(vcf_ch)
}