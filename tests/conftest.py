import pytest
import src.webservice
import pandas as pd

@pytest.fixture()
def uniprot_array():
    id_list = ['V9S9Y9', 'Q97UB2']
    columns = ['organism', 'mass', 'length']
    uniprot_dict = webservice.getUniProtDict(id_list, columns)
    uniprot_array = pd.DataFrame.from_dict(uniprot_dict, orient='index')
    uniprot_array['organism_split'] = uniprot_array.apply(lambda row: " ".join(row.organism.split(" ")[0:2]), axis=1)
    uniprot_array.index.names = ['id']

    return uniprot_array

