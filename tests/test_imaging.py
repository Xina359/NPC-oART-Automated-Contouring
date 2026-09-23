import numpy as np
import pytest
sitk=pytest.importorskip('SimpleITK')
from oart_analysis.imaging import labelmap_pair, require_same_grid, rtstruct_pair
from oart_analysis.segmentation_metrics import segmentation_metrics


def test_labelmaps_geometry_axis_order_and_bilateral_union(tmp_path):
    mask=np.zeros((8,9,10),np.uint8)
    mask[2:4,2:4,2:4]=1;mask[5:7,5:7,5:7]=2
    image=sitk.GetImageFromArray(mask);image.SetSpacing((.5,1,2))
    a=tmp_path/'ref.nii.gz';b=tmp_path/'pred.nii.gz'
    sitk.WriteImage(image,str(a));sitk.WriteImage(image,str(b))
    ref,pred,spacing=labelmap_pair(a,b,[1,2],[1,2])
    assert spacing==(2,1,.5)
    assert ref.sum()==16
    assert segmentation_metrics(ref,pred,spacing)['DSC']==1
    with pytest.raises(ValueError): labelmap_pair(a,b,[1,3],[1,2])
    shifted=sitk.Image(image);shifted.SetOrigin((1,0,0))
    with pytest.raises(ValueError): require_same_grid(image,shifted)


def test_rtstruct_roundtrip(tmp_path):
    pydicom=pytest.importorskip('pydicom')
    rt_utils=pytest.importorskip('rt_utils')
    from pydicom.dataset import FileDataset, FileMetaDataset
    from pydicom.uid import generate_uid, CTImageStorage, ExplicitVRLittleEndian
    ct=tmp_path/'ct';ct.mkdir()
    study,series,frame=generate_uid(),generate_uid(),generate_uid()
    for k in range(5):
        sop=generate_uid();meta=FileMetaDataset()
        meta.MediaStorageSOPClassUID=CTImageStorage;meta.MediaStorageSOPInstanceUID=sop
        meta.TransferSyntaxUID=ExplicitVRLittleEndian
        ds=FileDataset(str(ct/f'{k}.dcm'),{},file_meta=meta,preamble=b'\0'*128)
        ds.SOPClassUID=CTImageStorage;ds.SOPInstanceUID=sop;ds.Modality='CT'
        ds.StudyInstanceUID=study;ds.SeriesInstanceUID=series;ds.FrameOfReferenceUID=frame
        ds.PatientName='SYNTHETIC';ds.PatientID='SYNTHETIC';ds.PatientBirthDate='';ds.PatientSex=''
        ds.StudyDate='20000101';ds.StudyTime='120000';ds.SeriesDate='20000101';ds.SeriesTime='120000'
        ds.StudyID='1';ds.SeriesNumber=1;ds.InstanceNumber=k+1
        ds.StudyDescription='Synthetic test';ds.SeriesDescription='Synthetic CT'
        ds.Rows=16;ds.Columns=16;ds.PixelSpacing=[1,1.5];ds.SliceThickness=2
        ds.ImageOrientationPatient=[1,0,0,0,1,0];ds.ImagePositionPatient=[0,0,k*2]
        ds.BitsAllocated=16;ds.BitsStored=16;ds.HighBit=15;ds.PixelRepresentation=1
        ds.SamplesPerPixel=1;ds.PhotometricInterpretation='MONOCHROME2'
        ds.RescaleIntercept=0;ds.RescaleSlope=1
        ds.PixelData=np.zeros((16,16),np.int16).tobytes()
        ds.save_as(str(ct/f'{k}.dcm'),enforce_file_format=True)
    mask=np.zeros((16,16,5),bool);mask[3:9,4:10,1:4]=True
    rt=rt_utils.RTStructBuilder.create_new(dicom_series_path=str(ct))
    rt.add_roi(mask=mask,name='SYNTHETIC_ROI',approximate_contours=False)
    path=tmp_path/'structure.dcm';rt.save(str(path))
    ref,pred,spacing=rtstruct_pair(ct,path,path,['SYNTHETIC_ROI'],['SYNTHETIC_ROI'])
    assert spacing==(1,1.5,2)
    assert np.array_equal(ref,mask)
    assert segmentation_metrics(ref,pred,spacing)['DSC']==1
