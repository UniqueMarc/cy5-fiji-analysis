"""Portable exact-pixel summaries and Prism-ready columns. Python >= 3.9, stdlib only.

Punctum-level Student tests provide exploratory comparisons. Puncta
are nested in fields/cells; these P values are NOT biological-replicate tests.
"""
import argparse
import csv
import json
import math
import statistics as st
from pathlib import Path
from student_distribution import student_two_sided_p, stars

METRICS = ['Area_um2', 'Integrated_density_raw']

def write_csv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0] if rows else {})
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def student(x, y):
    df = len(x) + len(y) - 2
    if len(x) < 2 or len(y) < 2:
        return None, None, None
    variance = ((len(x)-1)*st.variance(x)+(len(y)-1)*st.variance(y))/df
    delta = st.mean(x)-st.mean(y)
    if variance == 0:
        return (0.0, df, 1.0) if delta == 0 else (None, df, 0.0)
    t = delta / math.sqrt(variance*(1/len(x)+1/len(y)))
    return t, df, student_two_sided_p(t, df)

def holm(rows, key):
    valid = [r for r in rows if r['P_raw'] is not None]
    last = 0.0
    for i, row in enumerate(sorted(valid, key=lambda r:r['P_raw'])):
        last = max(last, min(1.0, (len(valid)-i)*row['P_raw']))
        row[key] = last
        row[key+'_stars'] = stars(last)

def calculate(data, group_order):
    rows, images = data['regions'], data['images']
    assert len({r['Unique_ID'] for r in rows}) == len(rows), 'Duplicate punctum ID'
    by_image = {im['Image_Index']:im for im in images}
    assert len(by_image) == len(images), 'Duplicate image ID'
    for r in rows:
        assert r['Group'] == by_image[r['Image_Index']]['Group'], 'Mismatched group'
        assert all(math.isfinite(r[k]) and r[k] > 0 for k in METRICS)
    arrays = {g:{k:[r[k] for r in rows if r['Group']==g] for k in METRICS} for g in group_order}
    summary=[]
    for g in group_order:
        ims=[im for im in images if im['Group']==g]
        regions=[r for r in rows if r['Group']==g]
        assert sum(im['Regions'] for im in ims)==len(regions)
        rec=dict(Group=g,Images=len(ims),Automated_nuclear_regions=sum(im['Nuclei'] for im in ims),Cy5_regions=len(regions))
        for k in METRICS:
            a=arrays[g][k]
            rec.update({k+'_mean':st.mean(a) if a else None,k+'_median':st.median(a) if a else None,k+'_SD':st.stdev(a) if len(a)>1 else None})
        summary.append(rec)
    comparisons=[]
    pairs=[('mRNA',g) for g in group_order if g!='mRNA'] + [(g,g.replace('_mRNA','_chole_mRNA')) for g in ['DOPC_mRNA','DOPG_mRNA','66DOPC33DOPG_mRNA'] if g in group_order and g.replace('_mRNA','_chole_mRNA') in group_order]
    for k in METRICS:
        local=[]
        for a,b in pairs:
            if a not in arrays or b not in arrays:continue
            x,y=arrays[a][k],arrays[b][k]
            if len(x)<2 or len(y)<2:continue
            t,df,p=student(x,y)
            local.append(dict(Metric=k,Group_A=a,Group_B=b,n_A=len(x),n_B=len(y),Mean_A=st.mean(x),Mean_B=st.mean(y),t=t,df=df,P_raw=p,Stars_raw=stars(p),Observation_unit='pooled Cy5 region; exploratory'))
        holm(local,'P_Holm_within_metric')
        comparisons.extend(local)
    holm(comparisons,'P_Holm_all_metrics')
    return summary,comparisons,arrays

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_json',type=Path,help='Geometry_QC_results.json or frozen_exact_measurements.json')
    parser.add_argument('output_folder',type=Path,help='New folder; will not overwrite')
    parser.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'parameters/batch_20260930.json')
    args=parser.parse_args()
    data=json.loads(args.input_json.read_text(encoding='utf-8'))
    groups=json.loads(args.config.read_text(encoding='utf-8'))['group_order']
    summary,comparisons,arrays=calculate(data,groups)
    args.output_folder.mkdir(parents=True,exist_ok=False)
    write_csv(args.output_folder/'Exact_regions_long.csv',data['regions'])
    write_csv(args.output_folder/'Images.csv',data['images'])
    write_csv(args.output_folder/'Group_summary.csv',summary)
    write_csv(args.output_folder/'Exploratory_Student_tests.csv',comparisons)
    for k in METRICS:
        count=max((len(arrays[g][k]) for g in groups),default=0)
        write_csv(args.output_folder/(k+'_wide.csv'),[{g:arrays[g][k][i] if i<len(arrays[g][k]) else '' for g in groups} for i in range(count)],groups)
    # Preserve field-level summaries; do not invent biological replicate IDs.
    field=[]
    for im in data['images']:
        vals=[r for r in data['regions'] if r['Image_Index']==im['Image_Index']]
        rec={k:im[k] for k in ['Image_Index','Group','Regions','Nuclei']}
        for k in METRICS:rec[k+'_mean']=st.mean(r[k] for r in vals) if vals else None
        field.append(rec)
    write_csv(args.output_folder/'Per_image_summary.csv',field)
    print(f"Exported {len(data['images'])} fields, {len(data['regions'])} Cy5 regions -> {args.output_folder}")

if __name__=='__main__':main()
