import numpy as np
import pandas as pd
import pytest
from oart_analysis.statistics import describe, paired_test, benjamini_hochberg, orient_effect, precision_sample_size
from oart_analysis.analysis import compute_composites, paired_values
from oart_analysis.constants import SMALL


def test_paired_t_and_ci():
    result=paired_test(np.arange(5),np.arange(5)+np.arange(1,6))
    assert result['difference']['mean']==3
    assert result['difference']['sd']==pytest.approx(np.sqrt(2.5))
    assert result['t']==pytest.approx(4.242640687119285)
    assert result['p']==pytest.approx(0.013235599563682695)
    assert result['difference']['ci_lower']==pytest.approx(3-2.7764451051977987*np.sqrt(.5))


def test_bh_restores_order_and_uses_family_size():
    p=np.array([.2,.01,.6,.03,.02])
    assert benjamini_hochberg(p)==pytest.approx([.25,.05,.6,.05,.05])
    assert benjamini_hochberg([.01,.01,.8])==pytest.approx([.015,.015,.8])


def test_distance_ci_sign_and_order():
    assert orient_effect(-2,-3,-1,'HD95')==(2,1,3)
    assert orient_effect(.2,.1,.3,'DSC')==(.2,.1,.3)


def test_degenerate_pairs_and_invalid_inputs():
    assert paired_test([1,2,3],[1,2,3])['p']==1
    assert paired_test([1,2,3],[2,3,4])['p']==0
    with pytest.raises(ValueError): describe([1,np.nan,2])
    with pytest.raises(ValueError): benjamini_hochberg([.1,-1])
    with pytest.raises(ValueError): paired_test([1,2],[1,2,3])
    assert describe([-1,1])['cv'] is None


def synthetic_oars():
    return pd.DataFrame([{'Centre':'Example','Model state':'Example','Study ID':'SYNTHETIC_1',
        'ROI group':'Small OARs','ROI':roi,'Metric':'DSC','Value':i/20,
        'source_row':i+2} for i,roi in enumerate(SMALL)])


def test_composite_complete_equal_weight():
    frame=synthetic_oars()
    result=compute_composites(frame)
    assert result.Value.iloc[0]==pytest.approx(.225)
    assert result['Number of OARs averaged'].iloc[0]==10
    with pytest.raises(ValueError): compute_composites(frame.iloc[:-1])
    with pytest.raises(ValueError): compute_composites(pd.concat([frame,frame.iloc[:1]]))


def test_pairing_is_by_id_not_row_order():
    rows=[]
    for i in range(12):
        for state,value in [('Conventional U-Net comparator',i),('Locked source model',i+2)]:
            rows.append({'Centre':'Main centre','Model state':state,'Study ID':f'SYNTHETIC_{i:02d}',
                         'ROI':'GTVn','Metric':'DSC','Value':value/20})
    frame=pd.DataFrame(rows).sample(frac=1,random_state=7)
    ids,before,after=paired_values(frame,'Main centre','GTVn','DSC')
    assert after-before==pytest.approx(np.full(12,.1))
    with pytest.raises(ValueError): paired_values(frame.iloc[:-1],'Main centre','GTVn','DSC')


def test_sample_size_rounding():
    result=precision_sample_size()
    assert result['evaluable_n']==114
    assert result['recruitment_n']==135
