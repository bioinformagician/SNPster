params.vcf_dir = "/path/to/vcf/dir"
params.output_dir = "/path/to/output/dir"

process SPLIT_VCF {

  errorStrategy 'ignore'

  input:
    path file

  script:
  """
  python /app/vcf_splitter.py --vcf_file $file --output_dir $params.output_dir
  """
}


workflow {

    vcf_ch = Channel.fromPath("$params.vcf_dir/*.vcf.gz")
    SPLIT_VCF(vcf_ch)

}