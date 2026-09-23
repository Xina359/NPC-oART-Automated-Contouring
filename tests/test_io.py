import pandas as pd
import pytest
from oart_analysis.io import require_numeric, require_unique, read_source


def test_numeric_text_and_boolean_are_rejected():
    for value in ('0.7',True,None,float('nan')):
        with pytest.raises(ValueError): require_numeric(pd.DataFrame({'x':[value]}),['x'],'synthetic')


def test_missing_and_duplicate_keys():
    for values in (['A','A'],['A',''],['A',None]):
        with pytest.raises(ValueError): require_unique(pd.DataFrame({'id':values}),['id'],'synthetic')


def test_missing_workbook_explains_data_availability(tmp_path):
    with pytest.raises(FileNotFoundError,match='explicitly supplied'):
        read_source(tmp_path/'missing.xlsx')
