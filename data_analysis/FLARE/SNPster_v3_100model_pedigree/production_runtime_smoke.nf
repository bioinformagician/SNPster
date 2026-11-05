nextflow.enable.dsl = 2
params.harmonizer_options = ''

process PROBE_ORIGINAL_PRODUCTION_RUNTIME {
    tag "${name}"
    container { image }
    containerOptions { "${params.container_options} ${name == 'harmonizer' ? params.harmonizer_options : ''}" }
    publishDir "${params.run_dir}/probes", mode: 'copy'
    input:
    tuple val(name), val(image)
    output:
    path "${name}.txt"
    script:
    def checks = [
        standardizer: 'plink2 --version',
        combiner: 'python3 /app/main.py --help',
        harmonizer: '''test -x /app/dependencies/plink && test -x /app/dependencies/plink2 &&
          /app/dependencies/plink --version && /app/dependencies/plink2 --version &&
          sha256sum /app/dependencies/plink /app/dependencies/plink2 &&
          python3 -c 'import pysam; f = pysam.FastaFile("/srv/dependencies/imputation_runner/harmonizer/Homo_sapiens.GRCh38.dna.primary_assembly.fa"); print("FASTA chr22:", f.get_reference_length("22")); f.close()' ''',
        merger: 'bcftools --version',
        imputer: 'java -version && test -r "$BEAGLE_JAR" && sha256sum "$BEAGLE_JAR"',
        splitter: 'bcftools --version && bcftools plugin -l | grep -x split',
        qc: 'bgzip --version && tabix --version',
        gnomix: 'python3 -c "import numpy, scipy, sklearn; print(numpy.__version__, scipy.__version__, sklearn.__version__)"',
        flare: 'java -version && bcftools --version',
        preparation: 'bcftools --version && bgzip --version && tabix --version && plink2 --version'
    ]
    assert checks.containsKey(name)
    """
    set -o pipefail
    ps -e -o pid= -o ppid= > '${name}.txt'
    test -s '${name}.txt'
    {
      ${checks[name]}
    } > tools.txt 2>&1 || { cat tools.txt >&2; exit 1; }
    cat tools.txt >> '${name}.txt'
    python3 -c 'print("Original subprocess tools and real procps monitoring passed: ${name}")'
    """
}

workflow {
    Channel.of(
        tuple('standardizer', 'sha256:85e06b985dbaa682af024e5ad5af3cbbacc6eba24d5721fda47b355f6fa70544'),
        tuple('combiner', 'sha256:8b10363aadff5187ef09eff2b5b5d09b8b850173392dede0eadebf2885575ad4'),
        tuple('harmonizer', 'sha256:cd7907c83ad8d41017ce1e235172ace6c5c00d04565a11f75c0b185260c7a74f'),
        tuple('merger', 'sha256:f7f4bd94ac09a56125dcbccaa67e1854b97fbdd3164cc5b103ff642319b7483d'),
        tuple('imputer', 'sha256:eb8107e12137619252a4300f672e1e5d7ebbd493a4ab3a76809e00ad162194b1'),
        tuple('splitter', 'sha256:3db53a67ae2ce73ddd3fadb485a2d19c7e18dbc58ad6c86cc77ef88f330d1ccd'),
        tuple('qc', 'sha256:5d3d0fda0b7fd310efe5c112314d3923cdee2bd49d7abe95f5b6de3f895421e8'),
        tuple('gnomix', 'sha256:c4fd14a99e9d84839245b9242f36cf84306c32e67fa93b65fa5dbe5408d30e21'),
        tuple('flare', 'sha256:08a7c550ecd61aa3709379b3cf610b4d650cb655b59550c435d36e774c83868d'),
        tuple('preparation', 'sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b')
    ).set { images }
    PROBE_ORIGINAL_PRODUCTION_RUNTIME(images)
}
