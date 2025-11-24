import os

BEAGLE_JAR = os.getenv("BEAGLE_JAR", "/data/beagle.27Feb25.75f.jar")

JAVA_EXE = os.getenv("JAVA_EXE", "/usr/bin/java")



HEAP_GB = 8
THREADS = 2
GP_MIN = 0.90
DS_TOL = 0.20
SNPS_ONLY = True
BIALLELIC_ONLY = True