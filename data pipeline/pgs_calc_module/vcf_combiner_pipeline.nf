params.samplsheet_dir = "/path/to/samplesheet/dir"
params.output_dir = "/path/to/output/dir"

process MERGE_VCFS {

  input:
    path file

  script:
  """
  python /app/vcf_combiner.py --vcf_samplesheet_path $file --output_dir $params.output_dir
  """
}


workflow {

    samplesheet_ch = Channel.fromPath("$params.samplsheet_dir/*.csv")
    MERGE_VCFS(samplesheet_ch)

}