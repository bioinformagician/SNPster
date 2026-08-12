import argparse
import subprocess
import time
from pathlib import Path


def get_number_of_variants_in_vcf(vcf_file: str) -> int:
	"""Count variant records by streaming `bcftools view -H` output."""
	vcf_path = Path(vcf_file)
	if not vcf_path.exists():
		raise FileNotFoundError(f"VCF file does not exist: {vcf_file}")

	try:
		proc = subprocess.Popen(
			["bcftools", "view", "-H", str(vcf_path)],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			text=True,
		)
	except FileNotFoundError as exc:
		raise RuntimeError("bcftools is not installed or not available in PATH") from exc

	count = 0
	assert proc.stdout is not None
	for _ in proc.stdout:
		count += 1

	stderr_text = proc.stderr.read().strip() if proc.stderr is not None else ""
	return_code = proc.wait()
	if return_code != 0:
		raise RuntimeError(
			f"bcftools failed while counting variants for {vcf_file}: {stderr_text or f'exit code {return_code}'}"
		)

	return count


def main() -> None:
	parser = argparse.ArgumentParser(description="Benchmark bcftools variant counting on VCF files")
	parser.add_argument(
		"vcf_dir",
		nargs="?",
		default="/srv/imputed/frederik_myheritage/306/qc_output",
		help="Directory containing *.vcf.gz files",
	)
	args = parser.parse_args()

	vcf_dir = Path(args.vcf_dir)
	if not vcf_dir.is_dir():
		raise NotADirectoryError(f"Directory does not exist: {vcf_dir}")

	vcf_files = sorted(vcf_dir.glob("*.vcf.gz"))
	if not vcf_files:
		raise FileNotFoundError(f"No .vcf.gz files found in {vcf_dir}")

	print(f"Found {len(vcf_files)} files in {vcf_dir}")

	total_variants = 0
	total_seconds = 0.0
	start_all = time.perf_counter()

	for vcf_file in vcf_files:
		start = time.perf_counter()
		n_variants = get_number_of_variants_in_vcf(str(vcf_file))
		elapsed = time.perf_counter() - start

		total_variants += n_variants
		total_seconds += elapsed

		print(f"{vcf_file.name}: {n_variants} variants ({elapsed:.4f}s)")

	elapsed_all = time.perf_counter() - start_all
	print("\nSummary")
	print(f"Total variants: {total_variants}")
	print(f"Sum of per-file times: {total_seconds:.4f}s")
	print(f"Wall time (entire run): {elapsed_all:.4f}s")
	print(f"Avg time per file: {total_seconds / len(vcf_files):.4f}s")


if __name__ == "__main__":
	main()
