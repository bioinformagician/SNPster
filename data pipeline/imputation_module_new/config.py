import os

BEAGLE_JAR = os.getenv("BEAGLE_JAR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\dependencies\beagle.27Feb25.75f.jar")
JAVA_EXE = os.getenv("JAVA_EXE", r"C:\Program Files\Java\jre1.8.0_471\bin\java.exe")
OUTPUT_DIR = os.getenv("OUTPUT_DIR", r"C:\Users\frezz\Desktop\imputer_output\testing")
PLINK_MAP_DIR = os.getenv("PLINK_MAP_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\big_dependencies\plink.GRCh38.map")
BEAGLE_REFERENCE_DIR = os.getenv("BEAGLE_REFERENCE_DIR", r"C:\Users\frezz\Downloads\snpster\data pipeline\imputation_module\big_dependencies\beagle_references")
HEAP_GB = int(os.getenv("HEAP_GB", "8"))
THREADS = int(os.getenv("THREADS", "1"))
