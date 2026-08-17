nextflow.enable.dsl = 2

params.input_dir = "/srv/dependencies/imputation_runner/imputer/beagle_references_hgdp"
params.output_dir = "${params.input_dir}/prepared_vcf"
params.panel_output_dir = "/srv/dependencies/imputation_runner"
params.download_raw = true
params.reference_url = "https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
params.metadata_url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"

process DOWNLOAD_HGDP_1KGP_REFERENCE {

    tag "chr${chromosome}"

    publishDir params.input_dir, mode: 'copy', overwrite: false

    input:
    val chromosome

    output:
    path "hgdp1kgp_chr${chromosome}.filtered.SNV_INDEL.phased.shapeit5.bcf", emit: bcf
    path "hgdp1kgp_chr${chromosome}.filtered.SNV_INDEL.phased.shapeit5.bcf.csi"

    script:
    """
    set -euo pipefail

    filename="hgdp1kgp_chr${chromosome}.filtered.SNV_INDEL.phased.shapeit5.bcf"

    curl \\
      --fail \\
      --location \\
      --retry 5 \\
      --retry-all-errors \\
      --output "\${filename}" \\
      "${params.reference_url}/\${filename}"

    curl \\
      --fail \\
      --location \\
      --retry 5 \\
      --retry-all-errors \\
      --output "\${filename}.csi" \\
      "${params.reference_url}/\${filename}.csi"
    """
}

process DOWNLOAD_METADATA {

    publishDir params.input_dir, mode: 'copy', overwrite: false

    output:
    path "gnomad_meta_updated.tsv", emit: metadata

    script:
    """
    set -euo pipefail

    curl \\
      --fail \\
      --location \\
      --retry 5 \\
      --retry-all-errors \\
      --output "gnomad_meta_updated.tsv" \\
      "${params.metadata_url}"
    """
}

process CONVERT_METADATA_TO_PANEL {

    publishDir params.panel_output_dir, mode: 'copy', overwrite: true

    input:
    path metadata_file

    output:
    path "hgdp_1kg_panel.txt", emit: panel

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    # Mapping from gnomAD population codes to standard superpopulations
    POPULATION_INFERENCE_MAPPING = {
        'afr': 'AFR', 'amr': 'AMR', 'eas': 'EAS', 'sas': 'SAS',
        'nfe': 'EUR', 'fin': 'EUR', 'mid': 'EUR', 'oth': 'OTH',
    }
    
    print("Reading metadata from ${metadata_file}...")
    df = pd.read_csv("${metadata_file}", sep='\\t', low_memory=False)
    
    # Extract columns
    panel_df = pd.DataFrame({
        'sample': df['s'],
        'pop': df['population'],
        'super_pop_raw': df['population_inference.pop'],
        'gender': df['project_meta.sex']
    })
    
    # Map to standard superpopulations
    panel_df['super_pop'] = panel_df['super_pop_raw'].map(POPULATION_INFERENCE_MAPPING)
    
    # Handle unmapped values
    unmapped = panel_df['super_pop'].isna()
    if unmapped.any():
        print(f"Warning: {unmapped.sum()} samples with unmapped superpopulation")
        panel_df.loc[unmapped, 'super_pop'] = 'OTH'
    
    # Filter out 'OTH' superpopulation
    before = len(panel_df)
    panel_df = panel_df[panel_df['super_pop'] != 'OTH']
    print(f"Filtered out {before - len(panel_df)} samples with OTH superpopulation")
    
    # Standardize gender
    panel_df['gender'] = panel_df['gender'].astype(str).str.lower()
    panel_df['gender'] = panel_df['gender'].replace({'nan': 'unknown', 'na': 'unknown'})
    
    # Select final columns
    panel_df = panel_df[['sample', 'pop', 'super_pop', 'gender']]
    
    print(f"Total samples in panel: {len(panel_df)}")
    print(f"Superpopulation distribution:")
    print(panel_df['super_pop'].value_counts().to_string())
    
    # Save to file
    panel_df.to_csv('hgdp_1kg_panel.txt', sep='\\t', index=False)
    print("Panel file created successfully!")
    """
}

process PREPARE_HGDP_1KGP_REFERENCE {

    tag { bcf_file.simpleName }

    publishDir params.output_dir, mode: 'copy', overwrite: false

    input:
    path bcf_file

    output:
    path "*.SNV_biallelic.numericCHR.vcf.gz"
    path "*.SNV_biallelic.numericCHR.vcf.gz.tbi"

    script:
    """
    set -euo pipefail

    for chromosome in {1..22}; do
        printf 'chr%s\\t%s\\n' "\${chromosome}" "\${chromosome}"
    done > rename_from_chr.txt

    base=\$(basename "${bcf_file}" .bcf)
    output="\${base}.SNV_biallelic.numericCHR.vcf.gz"

    # Keep only biallelic SNPs; this intentionally removes indels.
    bcftools annotate \\
      --rename-chrs rename_from_chr.txt \\
      -Ou \\
      "${bcf_file}" | \\
    bcftools view \\
      --min-alleles 2 \\
      --max-alleles 2 \\
      --types snps \\
      -Oz \\
      -o "\${output}"

    tabix -f -p vcf "\${output}"
    """
}

workflow {
  if (params.download_raw) {
    bcf_ch = DOWNLOAD_HGDP_1KGP_REFERENCE(Channel.fromList((1..22))).bcf
    metadata_ch = DOWNLOAD_METADATA()
  } else {
    bcf_ch = Channel.fromPath("${params.input_dir}/*.bcf", checkIfExists: true)
    metadata_ch = Channel.fromPath("${params.input_dir}/gnomad_meta_updated.tsv", checkIfExists: true)
  }

    PREPARE_HGDP_1KGP_REFERENCE(bcf_ch)
    CONVERT_METADATA_TO_PANEL(metadata_ch)
}