import os


TEST_FILE = os.getenv("TEST_FILE", r"C:\Users\frezz\Downloads\snpster\data pipeline\standardizer_module\dependencies\test_file\23andMe_Full_20140502164209.zip")
HG18_TO_HG38_CHAIN_FILE = os.getenv("HG18_TO_HG38_CHAIN_FILE", r"C:\Users\frezz\Downloads\snpster\data pipeline\standardizer_module\dependencies\hg18ToHg38.over.chain.gz")
HG19_TO_HG38_CHAIN_FILE = os.getenv("HG19_TO_HG38_CHAIN_FILE", r"C:\Users\frezz\Downloads\snpster\data pipeline\standardizer_module\dependencies\hg19ToHg38.over.chain.gz")


ACCEPTED_VENDORS_DICT = { #setup is key = vendor and value = identifier substring in user upload files
    '23andme': '23andme',
    'AncestryDNA': 'ancestry',
    'ftdna': 'ftdna',
    'MyHeritage': 'myheritage',
    'livingdna': 'livingdna',
    'SelfDecode': 'selfdecode',
    'Genes for Good': 'genesforgood',
}

FORWARD_STRAND_VENDORS = ['23andme', 'AncestryDNA', 'ftdna', 'MyHeritage', 'livingdna', 'SelfDecode', 'Genes for Good']

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

GRCH_TO_HG_IDENTIFIER_DICT = {
    'GRCh36': 'hg18',
    'GRCh37': 'hg19',
    'GRCh38': 'hg38',
}

CHAIN_FILE_DICT={
                    'GRCh36': HG18_TO_HG38_CHAIN_FILE,
                    'GRCh37': HG19_TO_HG38_CHAIN_FILE,
                    'GRCh38': None
                    }



PVAR_REF_FILE = os.getenv(
    "PVAR_REF_FILE",
    r"C:\Users\frezz\Desktop\harmonizer_dependencies\all_hg38.pvar.zst",
)
