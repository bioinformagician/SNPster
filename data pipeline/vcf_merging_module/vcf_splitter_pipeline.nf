params.samplsheet_dir = "/path/to/samplesheet/dir"
params.output_dir = "/path/to/output/dir"

process SPLIT_VCF {

  input:
    path file

  script:
  """
  python /app/vcf_splitter.py --vcf_samplesheet_path $file --output_dir $params.output_dir
  """
}


workflow {

    samplesheet_ch = Channel.fromPath("$params.samplsheet_dir/vcf_splitting_chr*.csv")
    SPLIT_VCF(samplesheet_ch)

}