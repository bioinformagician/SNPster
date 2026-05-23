params.vcf_file_dir = "/path/to/samplesheet/dir"
params.output_dir = "/path/to/output/dir"

process IMPUTE_VCFS {

  maxForks 11

  input:
    path file


  script:
  """
  python /app/main.py --vcf_file $file --output_dir $params.output_dir
  """
}


workflow {

  vcf_file_ch = Channel.fromPath("$params.vcf_file_dir/*.vcf.gz")
    IMPUTE_VCFS(vcf_file_ch)

}