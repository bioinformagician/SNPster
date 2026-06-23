params.vcf_samplesheet_dir = "/path/to/samplesheet/dir"
params.output_dir = "/path/to/output/dir"

process MERGE_VCFS {
  
  maxRetries 2
  errorStrategy 'retry'

  cpus 2
  memory '20 GB'
  maxForks 4

  input:
    path file

  script:
  """
  python /app/vcf_combiner.py --vcf_samplesheet_path $file --output_dir $params.output_dir
  """
}

workflow {
  samplesheet_ch = Channel.fromPath("${params.vcf_samplesheet_dir}/vcf_merge_sheet_chr*.csv")
  MERGE_VCFS(samplesheet_ch)
}