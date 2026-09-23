import pandas as pd
from oart_analysis.plotting import summary_rows


def test_summary_mapping_uses_column_names_not_column_order():
    row={'Centre':'Synthetic centre','ROI group':'Synthetic group','ROI':'Synthetic ROI',
         'Metric':'DSC','Unit':'unitless','Paired n':12,'Before model':'before','After model':'after',
         'Mean difference (after-before)':.1,'SD of paired differences':.02,
         '95% CI lower':.08,'95% CI upper':.12,
         'Two-sided paired t-test P':'<0.0001','BH-adjusted q':.001,'source_row':2}
    for state in ('Before','After'):
        for label in ('mean','SD','minimum','maximum','Q1','median','Q3'):
            row[f'{state} {label}']=.4
    frame=pd.DataFrame([row])
    result=summary_rows(frame[sorted(frame.columns,reverse=True)])[0]
    assert result['centre']=='Synthetic centre'
    assert result['mean_difference']==.1
    assert result['before_mean']==.4
    assert result['p_numeric']==.0001
    assert result['q_numeric']==.001
