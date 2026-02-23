import subprocess

class EnvironmentHandler:
    def __init__(self, samplesheet_path: str, 
                 output_dir: str, 
                 reference_data_path: str,
                 pgs_id_file: str,
                 low_memory : bool, #split analysis into multiple runs of pgs calcs, calling NF pipeline multiple times (cache will be used)
                 ):
        
        self.samplesheet_path = samplesheet_path
        self.output_dir = output_dir
        self.reference_data_path = reference_data_path
        self.pgs_id_file = pgs_id_file
        self.low_memory = low_memory



class PGSCalculator_Config:
    def __init__(self, 
                 environment_handler: EnvironmentHandler,
                 pgs_id_str: str,
                 target_build: str = "GRCh38"):
        self.environment_handler = environment_handler
        self.pgs_id_str = pgs_id_str
        self.target_build = target_build
        self.set_pgs_id_str()
        
    
    def set_pgs_id_str(self):
        
        with open(self.environment_handler.pgs_id_file) as f:
            first_line = f.readline().strip()
        
        #check if pgs substring contained in first line
        if 'PGS' not in first_line:
            raise ValueError("PGS ID string not found in the first line of the file.")
        
        else:
            self.pgs_id_str = first_line
        
        
        
        


class PGSCalculator:
    def __init__(self, 
                 environment_handler: EnvironmentHandler,
                 pgscalculator_config: PGSCalculator_Config):
        self.environment_handler = environment_handler
        self.pgscalculator_config = pgscalculator_config
    
    def run_pgs_calculation(self):
        
        """nextflow run pgscatalog/pgsc_calc \
            -profile <docker/singularity/conda> \
            --input samplesheet.csv --target_build GRCh37 \
            --pgs_id PGS001229 \
            --run_ancestry pgsc_HGDP+1kGP_v1.tar.zst"""
        
        command = [
            "nextflow", "run", "pgscatalog/pgsc_calc",
            "-profile", "singularity",
            "--input", self.environment_handler.samplesheet_path,
            "--target_build", self.pgscalculator_config.target_build,
            "--pgs_id", self.pgscalculator_config.pgs_id_str,
            "--run_ancestry", self.environment_handler.reference_data_path,
            "--outdir", self.environment_handler.output_dir,
            "--min_overlap", "0.5"
        ]


        if self.environment_handler.low_memory == "true":

            pgs_id_list = self.pgscalculator_config.pgs_id_str.split(",")

            for pgs_id in pgs_id_list:

                command = [
                    "nextflow", "run", "pgscatalog/pgsc_calc",
                    "-profile", "singularity",
                    "--input", self.environment_handler.samplesheet_path,
                    "--target_build", self.pgscalculator_config.target_build,
                    "--pgs_id", pgs_id,
                    "--run_ancestry", self.environment_handler.reference_data_path,
                    "--outdir", self.environment_handler.output_dir,
                    "--min_overlap", "0.5"
                ]

                print(f"Running command {command}")

                try:
                    result = subprocess.run(command, check=True, capture_output=True, text=True)
                    print("PGS calculation completed successfully.")
                    print("Output:", result.stdout)

                except subprocess.CalledProcessError as e:
                    print("Error during PGS calculation:")
                    print("Return code:", e.returncode)
                    print("Output:", e.output)
                    print("Error message:", e.stderr)
        


        else:
            
            command = [
            "nextflow", "run", "pgscatalog/pgsc_calc",
            "-profile", "singularity",
            "--input", self.environment_handler.samplesheet_path,
            "--target_build", self.pgscalculator_config.target_build,
            "--pgs_id", self.pgscalculator_config.pgs_id_str,
            "--run_ancestry", self.environment_handler.reference_data_path,
            "--outdir", self.environment_handler.output_dir,
            "--min_overlap", "0.5"
        ]
            
            print(f"Running command {command}")
        
            try:
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                print("PGS calculation completed successfully.")
                print("Output:", result.stdout)

            except subprocess.CalledProcessError as e:
                print("Error during PGS calculation:")
                print("Return code:", e.returncode)
                print("Output:", e.output)
                print("Error message:", e.stderr)
        
    