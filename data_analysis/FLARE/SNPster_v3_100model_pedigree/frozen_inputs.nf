nextflow.enable.dsl=2

process VALIDATE_FROZEN_INPUTS {
    label 'pdb_light'
    publishDir "${params.outdir}/pilot_inputs", mode:'copy'
    input:
      path manifest
      path reusable
      path genetic_map
      path reference_vcf
    output:
      path 'validated_candidates.tsv', emit: manifest
      path 'pilot_input_audit.json', emit: audit
    script:
    def reuseOpt = reusable ? "--reusable-inputs '${reusable}'" : ''
    """
    set -euo pipefail
    python3 '${projectDir}/pedigree_benchmark/scripts/validate_frozen_candidates.py' \\
      --manifest '${manifest}' --output validated_candidates.tsv --report pilot_input_audit.json \\
      --genetic-map '${genetic_map}' --reference-vcf '${reference_vcf}' ${reuseOpt}
    """
}
