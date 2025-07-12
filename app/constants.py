"""
Константы для аннотации TAS2R38 и работы с Ensembl API.
"""

BITTER_GENE = 'TAS2R38'
BITTER_POSITIONS = [49, 262, 296]
ENSEMBL_SERVER = 'https://rest.ensembl.org'
GENE_LOOKUP_URL = (
    ENSEMBL_SERVER +
    '/lookup/symbol/homo_sapiens/{gene_name}?expand=1;content-type=application/json'
)
VARIATION_REGION_URL = (
    ENSEMBL_SERVER +
    '/overlap/region/human/{region}?feature=variation;content-type=application/json'
)
VEP_ID_URL = ENSEMBL_SERVER + '/vep/human/id'
ENSEMBL_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}
