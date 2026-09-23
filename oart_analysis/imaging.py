"""Optional local imaging adapter. Install the [imaging] extra to use this module."""
from pathlib import Path
import json
import numpy as np
from .segmentation_metrics import segmentation_metrics, union_masks
from .io import write_csv, write_json


def require_same_grid(reference, prediction):
    if reference.GetDimension()!=3 or prediction.GetDimension()!=3:
        raise ValueError('Only scalar three-dimensional labelmaps are supported')
    if reference.GetNumberOfComponentsPerPixel()!=1 or prediction.GetNumberOfComponentsPerPixel()!=1:
        raise ValueError('Labelmaps must be scalar images')
    if reference.GetSize()!=prediction.GetSize():
        raise ValueError('Mask sizes differ; align the images explicitly before evaluation')
    for attr in ('GetSpacing','GetOrigin','GetDirection'):
        if not np.allclose(getattr(reference,attr)(),getattr(prediction,attr)(),rtol=0,atol=1e-6):
            raise ValueError(f'Mask geometry differs: {attr}; no automatic registration or resampling is performed')
    for image in (reference,prediction):
        direction=np.array(image.GetDirection()).reshape(3,3)
        if not np.allclose(direction.T @ direction,np.eye(3),atol=1e-6):
            raise ValueError('Nonorthogonal image axes are unsupported')


def labelmap_pair(reference_path, prediction_path, reference_labels, prediction_labels):
    import SimpleITK as sitk
    reference, prediction=sitk.ReadImage(str(reference_path)),sitk.ReadImage(str(prediction_path))
    require_same_grid(reference,prediction)
    arrays=[]
    for image,labels in ((reference,reference_labels),(prediction,prediction_labels)):
        if not labels or any(not isinstance(i,int) or isinstance(i,bool) or i<=0 for i in labels):
            raise ValueError('Select one or more positive integer labels; use both label IDs for bilateral structures')
        array=sitk.GetArrayFromImage(image)
        if not np.isfinite(array).all() or not np.equal(array,np.round(array)).all() or (array<0).any():
            raise ValueError('Expected a nonnegative integer labelmap, not unthresholded probabilities')
        missing=set(labels)-set(np.unique(array))
        if missing: raise ValueError(f'Requested label IDs are absent: {sorted(missing)}')
        arrays.append(np.isin(array,labels))
    # SimpleITK arrays are z,y,x; GetSpacing() is x,y,z.
    return *arrays,tuple(reversed(reference.GetSpacing()))


def _validate_dicom_series(series):
    if len(series)<2:
        raise ValueError('At least two regularly spaced CT slices are required')
    first=series[0]
    orientation=np.array(first.ImageOrientationPatient,float)
    row,column=orientation[:3],orientation[3:]
    normal=np.cross(row,column)
    axes=np.stack((row,column,normal))
    if not np.allclose(axes @ axes.T,np.eye(3),atol=1e-5):
        raise ValueError('Invalid DICOM orientation cosines')
    spacings=np.array(first.PixelSpacing,float)
    if np.any(spacings<=0): raise ValueError('Invalid DICOM PixelSpacing')
    origin=np.array(first.ImagePositionPatient,float)
    positions=[]
    for ds in series:
        if ds.Modality!='CT' or ds.SeriesInstanceUID!=first.SeriesInstanceUID or ds.FrameOfReferenceUID!=first.FrameOfReferenceUID:
            raise ValueError('CT directory must contain one series in one frame of reference')
        if (ds.Rows,ds.Columns)!=(first.Rows,first.Columns) or not np.allclose(ds.PixelSpacing,spacings,atol=1e-6):
            raise ValueError('Inconsistent CT dimensions or spacing')
        if not np.allclose(ds.ImageOrientationPatient,orientation,atol=1e-6):
            raise ValueError('Inconsistent CT orientations')
        displacement=np.array(ds.ImagePositionPatient,float)-origin
        if not np.allclose(displacement-np.dot(displacement,normal)*normal,0,atol=1e-3):
            raise ValueError('Tilted/sheared CT grids need explicit preprocessing')
        positions.append(np.dot(displacement,normal))
    steps=np.diff(positions)
    if not (np.all(steps>0) or np.all(steps<0)) or not np.allclose(steps,steps[0],rtol=1e-4,atol=1e-3):
        raise ValueError('Duplicate, irregular or unsorted CT slice positions')
    return float(spacings[0]),float(spacings[1]),float(abs(steps[0]))


def rtstruct_pair(ct_directory, reference_path, prediction_path, reference_rois, prediction_rois):
    import pydicom
    from rt_utils import RTStructBuilder
    # Use a directory containing only the selected CT series. RTSTRUCT files
    # should be outside that directory. These are filesystem placeholders in
    # config/imaging.example.json; replace them with your own local paths.
    ref=RTStructBuilder.create_from(dicom_series_path=str(ct_directory),rt_struct_path=str(reference_path))
    pred=RTStructBuilder.create_from(dicom_series_path=str(ct_directory),rt_struct_path=str(prediction_path))
    spacing=_validate_dicom_series(ref.series_data)
    _validate_dicom_series(pred.series_data)
    uids=[ds.SOPInstanceUID for ds in ref.series_data]
    if uids!=[ds.SOPInstanceUID for ds in pred.series_data]:
        raise ValueError('The two RTSTRUCTs were not rasterized on identical slices')
    output=[]
    for rt,path,names in ((ref,reference_path,reference_rois),(pred,prediction_path,prediction_rois)):
        if not names or len(names)!=len(set(names)): raise ValueError('Provide unique ROI names')
        dataset=pydicom.dcmread(str(path),stop_before_pixels=True)
        frames={s.FrameOfReferenceUID for s in dataset.ReferencedFrameOfReferenceSequence}
        if frames!={ref.series_data[0].FrameOfReferenceUID}:
            raise ValueError('RTSTRUCT frame of reference does not match CT')
        available=rt.get_roi_names()
        if len(available)!=len(set(available)):
            raise ValueError('RTSTRUCT contains ambiguous duplicate ROI names')
        masks=[]
        for name in names:
            if name not in available: raise ValueError(f'ROI is absent from RTSTRUCT: {name}')
            mask=rt.get_roi_mask_by_name(name)
            if not mask.any(): raise ValueError(f'Empty RTSTRUCT component: {name}')
            masks.append(mask)
        # rt-utils uses row,column,slice array order.
        output.append(union_masks(masks))
    return *output,spacing


def evaluate_config(config_path, output):
    config_path=Path(config_path)
    config=json.loads(config_path.read_text(encoding='utf-8'))
    cases=config.get('cases',[])
    if not cases: raise ValueError('Add cases and replace the local input paths in the imaging configuration')
    def local_path(value):
        if not value or 'REPLACE_' in str(value): raise ValueError('Replace every REPLACE_ path in the imaging configuration')
        p=Path(value)
        p=p if p.is_absolute() else config_path.parent/p
        if not p.exists(): raise FileNotFoundError(p)
        return p
    results=[]
    keys=set()
    for case in cases:
        identifier=case['study_id']
        for roi in case['rois']:
            key=(identifier,case.get('model_state',''),roi['name'])
            if key in keys: raise ValueError(f'Duplicate imaging evaluation key: {key}')
            keys.add(key)
            if case['format']=='labelmap':
                ref,pred,spacing=labelmap_pair(local_path(case['reference_path']),local_path(case['prediction_path']),roi['reference_labels'],roi['prediction_labels'])
            elif case['format']=='rtstruct':
                ref,pred,spacing=rtstruct_pair(local_path(case['ct_directory']),local_path(case['reference_path']),local_path(case['prediction_path']),roi['reference_rois'],roi['prediction_rois'])
            else: raise ValueError('format must be labelmap or rtstruct')
            results.append(dict(study_id=identifier,model_state=case.get('model_state',''),roi=roi['name'],**segmentation_metrics(ref,pred,spacing)))
    if not results: raise ValueError('No ROIs configured')
    output=Path(output)
    write_csv(results,output/'segmentation_metrics.csv')
    write_json({'cases':len(cases),'evaluated_rois':len(results),'convention':'6-neighbour voxel boundary centres; pooled distances; linear percentile; ASD weighted by surface-point counts; no automatic empty-mask handling'},output/'metric_convention.json')
    return len(results)
