import os

#use docker environment variables or default values

PLINK_PREFIX = os.getenv("PLINK_PREFIX", "plink_temp")


PLINK_REFERENCE_FASTA = os.getenv(
    "PLINK_REFERENCE_FASTA",
    "/data/hs37d5.fa.zst",
)

TEST_FILE = os.getenv(
    "TEST_FILE",
    "/app/dependencies/genome_Frederik_FangelTolberg_v5_Full_20241117223640.txt",
)
PVAR_REF_FILE = os.getenv(
    "PVAR_REF_FILE",
    "/data/all_phase3.pvar.zst",
)

PLINK_1_9_PATH = os.getenv(
    "PLINK_1_9_PATH",
    "/app/dependencies/plink",
)

PLINK_1_9_MEMORY_MB = os.getenv(
    "PLINK_1_9_MEMORY_MB",
    "1000",
)

PLINK_1_9_THREADS = os.getenv(
    "PLINK_1_9_THREADS",
    "1",
)


PLINK_2_0_PATH = os.getenv(
    "PLINK_2_0_PATH",
    "/app/dependencies/plink2",
)


ACCEPTED_VENDORS_DICT = { #setup is key = vendor and value = identifier substring
    '23andme': '23andme',
    'AncestryDNA': 'ancestry',
    'ftdna': 'ftdna',
    'MyHeritage': 'myheritage',
    'livingdna': 'livingdna',
    'SelfDecode': 'selfdecode',
    'Genes for Good': 'genesforgood',
}

GENOME_BUILD_DICT = {
    'build 38' : 'GRCh38',
    'human reference build 38' : 'GRCh38',
    'human assembly build 38' : 'GRCh38',
    'GRCh38' : 'GRCh38',
    'human assembly build 37' : 'GRCh37',
    'GRCh37' : 'GRCh37',
    'human reference build 37' : 'GRCh37',
    'Reference Build 37' :  'GRCh37',
    'build 37' : 'GRCh37',
    'build 36' : 'GRCh36',
    'human reference build 36' : 'GRCh36',
    'human assembly build 36' : 'GRCh36',
    'GRCh36' : 'GRCh36',
}