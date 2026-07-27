import pandas as pd
from pathlib import Path
from typing import List
import pyarrow as pa
import pyarrow.parquet as pq
from config import MAXIMUM_OVERLAP_PERCENTAGE





class Combiner:
    def __init__(self, file_paths: List[Path],
                 combined_data: pd.DataFrame | None = None,
                 output_dir: str | Path | None = None,
                 imputation_id: str | int | None = None):
        self.file_paths = [Path(file_path) for file_path in file_paths]
        self.combined_data = combined_data if combined_data is not None else pd.DataFrame()
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.imputation_id = imputation_id
        self.parquet_metadata: dict[bytes, bytes] = {}

    def combine_microarray_data(self) -> pd.DataFrame:
        
        for file in self.file_paths:
            if not file.exists():
                raise FileNotFoundError(f"Missing parquet file: {file}")
        
        dataframes: list[pd.DataFrame] = []
        expected_file_imputation_id: str | None = None

        for idx, file in enumerate(self.file_paths):
            table = pq.read_table(file)
            dataframes.append(table.to_pandas())

            file_metadata = table.schema.metadata or {}

            file_imputation_id_bytes = file_metadata.get(b'imputation_id')
            if file_imputation_id_bytes is None:
                raise ValueError(
                    f"Missing parquet metadata key b'imputation_id' in file {file}."
                )

            try:
                file_imputation_id = file_imputation_id_bytes.decode('utf-8')
            except UnicodeDecodeError as exc:
                raise ValueError(
                    f"Invalid UTF-8 imputation_id metadata in file {file}."
                ) from exc

            if expected_file_imputation_id is None:
                expected_file_imputation_id = file_imputation_id
            elif file_imputation_id != expected_file_imputation_id:
                raise ValueError(
                    "Conflicting imputation_id across parquet files: "
                    f"expected {expected_file_imputation_id!r}, found {file_imputation_id!r} in file {file}."
                )

            if self.imputation_id is not None and file_imputation_id != str(self.imputation_id):
                raise ValueError(
                    "Provided imputation_id does not match file metadata: "
                    f"expected {self.imputation_id!r}, found {file_imputation_id!r} in file {file}."
                )

            if idx == 0:
                self.parquet_metadata = dict(file_metadata)
                continue

            # Keep metadata stable across files (except vendor and file_id).
            for key, value in file_metadata.items():
                if key in {b'vendor', b'file_id'}:
                    continue
                if key not in self.parquet_metadata:
                    self.parquet_metadata[key] = value
                elif self.parquet_metadata[key] != value:
                    raise ValueError(
                        f"Conflicting parquet metadata for key {key!r} in file {file}."
                    )

        self.combined_data = pd.concat(dataframes, ignore_index=True)
        return self.combined_data
    
    def validate_combined_data(self) -> pd.DataFrame:
        if len(self.file_paths) == 1:
            print("Single input parquet detected; skipping cross-file genotype validation.")
            return self.combined_data

        # Normalize genotype order: TA -> AT, GA -> AG, etc.
        self.combined_data["genotype"] = (
            self.combined_data["genotype"]
            .astype(str)
            .apply(lambda x: "".join(sorted(x)))
        )

        # Get rows where the RSID occurs more than once
        duplicate_rsids = self.combined_data[
            self.combined_data.duplicated(subset="# rsid", keep=False)
        ].copy()

        # For each duplicated RSID, determine whether genotypes agree
        genotype_check = (
            duplicate_rsids
            .groupby("# rsid")["genotype"]
            .nunique()
            .reset_index(name="n_genotypes")
        )

        genotype_check["genotype_status"] = genotype_check["n_genotypes"].apply(
            lambda n: "same" if n == 1 else "different"
        )

        n_intersecting_rsids = len(genotype_check)
        n_disagreeing_rsids = (genotype_check["genotype_status"] == "different").sum()
        
        percent_disagreeing = (n_disagreeing_rsids / n_intersecting_rsids * 100) if n_intersecting_rsids > 0 else 0.0
        
        if percent_disagreeing > MAXIMUM_OVERLAP_PERCENTAGE:
            raise ValueError(
                f"More than {MAXIMUM_OVERLAP_PERCENTAGE}% of intersecting RSIDs have conflicting genotypes: "
                f"{percent_disagreeing:.2f}% ({n_disagreeing_rsids}/{n_intersecting_rsids})"
                f"Data likely from different individuals. Please check input."
            )
        
        else:
            
            #remove disagreeing RSIDs from combined_data
            disagreeing_rsids = genotype_check.loc[genotype_check["genotype_status"] == "different", "# rsid"]
            self.combined_data = self.combined_data[~self.combined_data["# rsid"].isin(disagreeing_rsids)].reset_index(drop=True)
            
            print(f"Total variants after removing disagreeing RSIDs: {len(self.combined_data)}")
                
            

    def write_parquet_output(self) -> None:
        """Naming of file should follow {IMPIDx}.chr{CHR}.{STAGE}.fileextension"""
        if self.combined_data.empty:
            raise ValueError("No combined data available. Run load_microarray_data() first.")

        if self.output_dir is None:
            raise ValueError("output_dir is not set for Combiner.")

        output_dir = Path(self.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        metadata = dict(self.parquet_metadata)
        metadata[b'vendor'] = b'combined'
        metadata.pop(b'file_id', None)

        resolved_imputation_id = self.imputation_id
        if resolved_imputation_id is None and b'imputation_id' in metadata:
            resolved_imputation_id = metadata[b'imputation_id'].decode('utf-8')
        if resolved_imputation_id is None:
            raise ValueError("imputation_id is required (constructor argument or parquet metadata).")

        metadata[b'imputation_id'] = str(resolved_imputation_id).encode('utf-8')

        if "chromosome" in self.combined_data.columns:
            grouped_data = self.combined_data.groupby("chromosome", sort=False)
        else:
            grouped_data = [("ALL", self.combined_data)]

        for chrom, df_chrom in grouped_data:
            output_path = output_dir / f"IMPID{resolved_imputation_id}.chr{chrom}.standardizedMicroarray.parquet"

            table = pa.Table.from_pandas(df_chrom)
            table = table.replace_schema_metadata(metadata)

            pq.write_table(table, output_path)
            print(f"Standardized data written to {output_path}")
            print(
                "Metadata: "
                f"vendor=combined, "
                f"genome_build={metadata.get(b'genome_build', b'').decode('utf-8')}, "
                f"is_forward_strand={metadata.get(b'is_forward_strand', b'').decode('utf-8')}, "
                f"imputation_id={metadata.get(b'imputation_id', b'').decode('utf-8')}"
            )
                