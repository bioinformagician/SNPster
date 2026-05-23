params.vcf_file_dir = "/path/to/samplesheet/dir"
params.output_dir = "/path/to/output/dir"

process QC_VCFS {

  maxForks = 10

  input:
    path file


  script:
  """
  python /app/main.py --input_file $file --output_dir $params.output_dir
  """
}


workflow {

    vcf_file_ch = Channel.fromPath("$params.vcf_file_dir/split*.vcf.gz")
    QC_VCFS(vcf_file_ch)

}