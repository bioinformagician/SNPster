from pgs_classes import EnvironmentHandler, PGSCalculator_Config, PGSCalculator
from config import SAMPLESHEET_PATH, OUTPUT_DIR, REFERENCE_DATA_PATH, PGS_ID_FILE


environment_handler = EnvironmentHandler(
    samplesheet_path=SAMPLESHEET_PATH,
    output_dir=OUTPUT_DIR,
    reference_data_path=REFERENCE_DATA_PATH,
    pgs_id_file=PGS_ID_FILE
)

pgs_calculator_config = PGSCalculator_Config(
    environment_handler=environment_handler,
    pgs_id_str=None)

pgs_calculator = PGSCalculator(
    environment_handler=environment_handler,
    pgscalculator_config=pgs_calculator_config
)

pgs_calculator.run_pgs_calculation()

