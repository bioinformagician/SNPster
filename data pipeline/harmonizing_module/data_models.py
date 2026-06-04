from typing import List, Optional
from dataclasses import dataclass
import pandas as pd
import re


@dataclass
class DataContainer:
    imputation_id: Optional[str] = None
    vendor: Optional[str] = None
    genome_build: Optional[str] = None
    is_forward_strand: Optional[bool] = True
    microarray_data: Optional[pd.DataFrame] = None
    reference_data: Optional[pd.DataFrame] = None
    harmonized_data: Optional[pd.DataFrame] = None
    harmonization_stats: Optional[pd.DataFrame] = None
    

    
    def _fix_one_allele(self, allele, ref, alts):
        """
        Return (new_allele, flipped) where:
        - new_allele is the original allele or its complement if that matches
        - flipped is True iff we changed to the complement
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        if pd.isna(allele):
            return allele, False

        allele = str(allele).upper()

        if allele == ref or allele in alts:
            return allele, False

        comp = complement.get(allele)
        if comp == ref or comp in alts:
            return comp, True

        return allele, False


    def flip_alleles_to_match_ref_alt(self) -> None:
        """
        For each row, ensure allele1/allele2 match REF or ALT (multi-allelic supported).
        If not, try complement; if complement matches, replace allele and mark *_flipped=True.
        """
        out = self.harmonized_data

        # Normalize strings (only for columns that exist)
        for col in ("REF", "ALT", "allele1", "allele2"):
            if col in out.columns:
                out[col] = out[col].astype("string").str.upper()
                out[col] = out[col].replace(["NAN", "NONE", "NULL", ""], pd.NA)

        def upper_or_empty(value):
            return "" if pd.isna(value) else str(value).upper()

        def allele_ok(allele, ref, alts):
            return (not pd.isna(allele)) and (allele == ref or allele in alts)

        def process_row(r):
            ref = upper_or_empty(r.get("REF"))
            alt_str = upper_or_empty(r.get("ALT"))
            alts = {x for x in alt_str.split(",") if x and x != "."}

            a1_new, a1_flip = self._fix_one_allele(r.get("allele1"), ref, alts)
            a2_new, a2_flip = self._fix_one_allele(r.get("allele2"), ref, alts)

            a1_ok = allele_ok(a1_new, ref, alts)
            a2_ok = allele_ok(a2_new, ref, alts)

            return pd.Series({
                "allele1": a1_new,
                "allele2": a2_new,
                "allele1_flipped": a1_flip,
                "allele2_flipped": a2_flip,
                "allele1_ok": a1_ok,
                "allele2_ok": a2_ok,
                "both_ok": a1_ok and a2_ok,
            })

        cols = [
            "allele1", "allele2",
            "allele1_flipped", "allele2_flipped",
            "allele1_ok", "allele2_ok", "both_ok"
        ]
        out[cols] = out.apply(process_row, axis=1)

        self.harmonized_data = out
    
    def set_harmonization_stats(self) -> None:
        total = len(self.harmonized_data)

        # Booleans sum to counts; cast to int just in case columns are 0/1 ints
        flipped_allele1 = int(self.harmonized_data["allele1_flipped"].sum())
        flipped_allele2 = int(self.harmonized_data["allele2_flipped"].sum())
        ok_allele1      = int(self.harmonized_data["allele1_ok"].sum())
        ok_allele2      = int(self.harmonized_data["allele2_ok"].sum())
        ok_both         = int(self.harmonized_data["both_ok"].sum())

        pct = (lambda x: (x / total * 100.0) if total else 0.0)

        data = pd.DataFrame()


        data = pd.DataFrame([{
            "total_snps": total,
            "flipped_allele1": flipped_allele1,
            "flipped_allele1_pct": round(pct(flipped_allele1), 2),
            "flipped_allele2": flipped_allele2,
            "flipped_allele2_pct": round(pct(flipped_allele2), 2),
            "ok_allele1": ok_allele1,
            "ok_allele1_pct": round(pct(ok_allele1), 2),
            "ok_allele2": ok_allele2,
            "ok_allele2_pct": round(pct(ok_allele2), 2),
            "ok_both": ok_both,
            "ok_both_pct": round(pct(ok_both), 2),
        }])

        self.harmonization_stats = data
    
    def harmonize_data(self) -> None:
        
        

        self.microarray_data[['allele1', 'allele2']] = self.microarray_data['genotype'].str.extract(r'([ACGT])([ACGT])')

    
        #left join reference_data onto user_data on rsid
        self.harmonized_data = pd.merge(self.microarray_data, self.reference_data, left_on="# rsid", right_on="ID", how="left")


        self.flip_alleles_to_match_ref_alt()
        self.set_harmonization_stats()


        self.harmonized_data['genotype'] = self.harmonized_data['allele1'] + self.harmonized_data['allele2']
        self.harmonized_data = self.harmonized_data[self.harmonized_data["both_ok"] == True].drop(columns=["ID", "allele1_flipped", "allele2_flipped", "allele1_ok", "allele2_ok", "both_ok", "REF", "ALT", "allele1", "allele2"])
        
        # Use the reference genome positions (POS, CHROM) instead of user's original positions
        # Drop old position/chromosome if they exist, then rename reference columns
        if 'position' in self.harmonized_data.columns:
            self.harmonized_data = self.harmonized_data.drop(columns=['position'])
        if 'chromosome' in self.harmonized_data.columns:
            self.harmonized_data = self.harmonized_data.drop(columns=['chromosome'])
        
        # Rename reference columns to standard names
        self.harmonized_data = self.harmonized_data.rename(columns={'POS': 'position', 'CHROM': 'chromosome'})
        
        # Filter out alternate contigs, patches, and unplaced sequences
        # Keep only standard chromosomes: 1-22, X, Y, MT
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']
        valid_chromosomes_with_prefix = valid_chromosomes + [f'chr{c}' for c in valid_chromosomes]
        
        initial_count = len(self.harmonized_data)
        self.harmonized_data['chrom_clean'] = self.harmonized_data['chromosome'].astype(str).str.replace('chr', '')
        self.harmonized_data = self.harmonized_data[self.harmonized_data['chrom_clean'].isin(valid_chromosomes)]
        self.harmonized_data = self.harmonized_data.drop(columns=['chrom_clean'])
        filtered_count = initial_count - len(self.harmonized_data)
        
        if filtered_count > 0:
            print(f"Filtered out {filtered_count} variants on non-standard chromosomes (alt contigs, patches, etc.)")
        
        # Keep only the 4 columns needed for PLINK --23file format
        self.harmonized_data = self.harmonized_data[['# rsid', 'chromosome', 'position', 'genotype']]




