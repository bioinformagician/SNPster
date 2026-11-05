#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKFLOW="${1:-$ROOT/main.nf}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

for chr in 1 2 10; do
  {
    printf '##fileformat=VCFv4.2\n'
    for contig in 1 2 10; do
      printf '##contig=<ID=%s,length=1000>\n' "$contig"
    done
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild_a\tchild_b\n'
    printf '%s\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\n' "$chr"
    printf '%s\t20\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|0\n' "$chr"
  } | bcftools view -Oz -o "$TMP/chr${chr}.vcf.gz"
  test ! -e "$TMP/chr${chr}.vcf.gz.tbi"
  test ! -e "$TMP/chr${chr}.vcf.gz.csi"
done

# Exercise the real process, isolated from the expensive benchmark workflow.
awk '
  /^process PDB_MERGE_PANEL \{/ { copying = 1 }
  copying { print }
  copying && /^}/ { exit }
' "$WORKFLOW" > "$TMP/main.nf"
grep -q '^process PDB_MERGE_PANEL {' "$TMP/main.nf"
cat >> "$TMP/main.nf" <<'NF'

workflow {
    inputs = Channel.of(tuple('test_batch', [10, 2, 1],
        [file("${params.inputs}/chr10.vcf.gz"),
         file("${params.inputs}/chr2.vcf.gz"),
         file("${params.inputs}/chr1.vcf.gz")]))
    PDB_MERGE_PANEL(inputs)
}
NF

cd "$TMP"
nextflow -C /dev/null run main.nf -work-dir "$TMP/work" \
  --inputs "$TMP" --outdir "$TMP/results"

test -s "$TMP/results/pedigree/panel_genomes/test_batch.autosomes.vcf.gz"
GENOME="$(find "$TMP/work" -type f -name test_batch.autosomes.vcf.gz)"
test -s "$GENOME"
test -s "$GENOME.tbi"
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$GENOME" > actual.tsv
{
  for chr in 1 2 10; do
    printf '%s\t10\t0|1\t1|0\n' "$chr"
    printf '%s\t20\t1|1\t0|0\n' "$chr"
  done
} > expected.tsv
diff -u expected.tsv actual.tsv
printf 'child_a\nchild_b\n' > expected_samples.txt
bcftools query -l "$GENOME" > actual_samples.txt
diff -u expected_samples.txt actual_samples.txt
test "$(bcftools index -n "$GENOME")" -eq 6
test "$(bcftools view -H -r 10:10-10 "$GENOME" | wc -l)" -eq 1
echo 'PASS: unindexed chromosome inputs merge in numeric order, retaining samples and phased genotypes; output is indexed.'
