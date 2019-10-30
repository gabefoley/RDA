import pytest

def test_uniprot_array(uniprot_array):
    for k,v in uniprot_array.items():
        print (k, v)
        print()
        assert uniprot_array.columns.all(['organism', 'mass', 'length', 'organism_split'])
