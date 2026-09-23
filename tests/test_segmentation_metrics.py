import numpy as np
import pytest
from scipy.spatial.distance import cdist
from oart_analysis.segmentation_metrics import segmentation_metrics, union_masks


def test_identical_mask():
    x=np.zeros((9,10,11),bool);x[2:7,2:8,2:9]=True
    result=segmentation_metrics(x,x,(2,1,.5))
    assert result['DSC']==1
    assert result['HD95_mm']==0
    assert result['ASD_mm']==0
    assert result['reference_volume_cm3']==pytest.approx(5*6*7/1000)


def test_single_points_use_anisotropic_spacing():
    x=np.zeros((8,8,8),bool);y=x.copy()
    x[1,1,1]=True;y[2,3,4]=True
    distance=np.sqrt(2**2+6**2+12**2)
    result=segmentation_metrics(x,y,(2,3,4))
    assert result['DSC']==0
    assert result['HD95_mm']==pytest.approx(distance)
    assert result['ASD_mm']==pytest.approx(distance)


def test_pooled_weighting_and_percentile():
    x=np.zeros((12,12,12),bool);y=x.copy()
    a=np.array([[1,1,1],[7,1,1]])
    b=np.array([[1,2,1],[4,4,1],[8,3,1],[10,10,10]])
    x[tuple(a.T)]=True;y[tuple(b.T)]=True
    distance=cdist(a,b)
    both=np.r_[distance.min(axis=0),distance.min(axis=1)]
    result=segmentation_metrics(x,y,(1,1,1))
    assert result['ASD_mm']==pytest.approx(both.mean())
    assert result['HD95_mm']==pytest.approx(np.percentile(both,95,method='linear'))
    assert not np.isclose(both.mean(),(distance.min(axis=0).mean()+distance.min(axis=1).mean())/2)


def test_union_and_empty_masks():
    a=np.zeros((8,8,8),bool);b=a.copy()
    a[1:3,1:3,1:3]=True;b[5:7,5:7,5:7]=True
    whole=union_masks([a,b])
    assert whole.sum()==16
    assert segmentation_metrics(whole,whole,(1,1,1))['DSC']==1
    with pytest.raises(ValueError): segmentation_metrics(a,np.zeros_like(a),(1,1,1))
    with pytest.raises(ValueError): segmentation_metrics(a,a,(1,0,1))
    with pytest.raises(ValueError): segmentation_metrics(a,a[:2],(1,1,1))
    with pytest.raises(ValueError): segmentation_metrics(a.astype(float)*.5,a,(1,1,1))
