import copy
import json
import shutil
import zipfile
from pathlib import Path
import pytest
from oart_analysis import cli, plot_data, plotting


def test_bundled_coverage_and_real_case():
    data, metadata = plot_data.load_bundled()
    assert {key:len(rows) for key,rows in data.items()} == {
        'summary_s4':90,'summary_s5':315,'performance':2430,'timing':135,'gtvn':108,'case':21,
    }
    assert metadata['input_mode']=='bundled'
    assert metadata['full_patient_level_analysis'] is False
    assert {r['study_id'] for r in data['case']}=={'MC_R023'}


def test_bundled_validation_rejects_missing_pairs_and_bad_timing():
    data,_=plot_data.load_bundled()
    incomplete=copy.deepcopy(data)
    incomplete['performance'].pop()
    with pytest.raises(ValueError,match='135 complete'): plot_data.validate(incomplete)
    wrong=copy.deepcopy(data)
    wrong['timing'][0]['combined_s']+=1
    with pytest.raises(ValueError,match='timing sum'): plot_data.validate(wrong)


def test_manifest_detects_modified_csv(tmp_path,monkeypatch):
    original=Path(plot_data.files('oart_analysis').joinpath('data'))
    shutil.copytree(original,tmp_path/'data')
    path=tmp_path/'data/prospective_timing.csv'
    path.write_bytes(path.read_bytes()+b'\n')
    monkeypatch.setattr(plot_data,'files',lambda package:tmp_path)
    with pytest.raises(ValueError,match='checksum mismatch'): plot_data.load_bundled()


def test_default_ignores_local_workbook(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path/'data').mkdir()
    (tmp_path/'data/Source_data.xlsx').write_text('invalid workbook',encoding='utf-8')
    (tmp_path/'Source_data.xlsx').write_text('also invalid',encoding='utf-8')
    def unexpected_read(path): raise AssertionError('Default mode must not read a workbook')
    def lightweight_render(data,figure,output):
        folder=output/f'Figure{figure}';folder.mkdir(parents=True)
        return folder
    monkeypatch.setattr(cli,'read_source',unexpected_read)
    monkeypatch.setattr(plotting,'render',lightweight_render)
    cli.main(['plot','--figure','S4','--output',str(tmp_path/'out')])
    metadata=json.loads((tmp_path/'out/data_source.json').read_text())
    assert metadata['input_mode']=='bundled'
    assert not metadata['full_patient_level_analysis']


@pytest.mark.parametrize('file_exists',[False,True])
def test_explicit_invalid_workbook_never_falls_back(tmp_path,monkeypatch,file_exists):
    path=tmp_path/'explicit.xlsx'
    if file_exists:path.write_text('invalid workbook',encoding='utf-8')
    def unexpected_bundle(): raise AssertionError('Explicit workbook mode must not fall back')
    monkeypatch.setattr(plot_data,'load_bundled',unexpected_bundle)
    with pytest.raises((FileNotFoundError,zipfile.BadZipFile)):
        cli.main(['plot','--figure','S4','--source-data',str(path),'--output',str(tmp_path/'out')])


def test_full_analysis_requires_explicit_workbook():
    with pytest.raises(SystemExit) as error:cli.main(['analyse'])
    assert error.value.code==2


def test_explicit_workbook_routes_to_analysis(tmp_path,monkeypatch):
    data,_=plot_data.load_bundled()
    calls=[]
    monkeypatch.setattr(cli,'read_source',lambda path:calls.append('read') or {'example':True})
    monkeypatch.setattr(cli,'analyse',lambda source,output:calls.append('analyse') or {'table_cells_verified':1,'composite_values_verified':1})
    monkeypatch.setattr(plot_data,'from_workbook',lambda source:(data,dict(input_mode='workbook',full_patient_level_analysis=True)))
    cli.main(['analyse','--source-data',str(tmp_path/'selected.xlsx'),'--output',str(tmp_path/'out')])
    assert calls==['read','analyse']
    assert json.loads((tmp_path/'out/data_source.json').read_text())['input_mode']=='workbook'
