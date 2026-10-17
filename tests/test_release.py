"""Regression checks against the independently saved per-cell figure data."""
from pathlib import Path
import copy
import csv
import json
import math
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from summarize import calculate, export, METRICS


class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'results/per_cell/per_cell_measurements.json').read_text())
        cls.groups=json.loads((ROOT/'parameters/batch_20260930.json').read_text())['group_order']
        cls.expected=json.loads((ROOT/'tests/expected_per_cell.json').read_text())
        cls.cells,cls.puncta,cls.summary,cls.arrays=calculate(cls.data,cls.groups)

    def test_all_cells_match_independent_figure_data(self):
        actual={(r['Image_Index'],r['Cell_ID']):r for r in self.cells}
        self.assertEqual(len(actual),537)
        for e in self.expected:
            a=actual[e['Image_Index'],e['Cell_ID']]
            for k,v in e.items():
                if isinstance(v,float): self.assertTrue(math.isclose(a[k],v,rel_tol=1e-12,abs_tol=1e-12),(k,a,e))
                else: self.assertEqual(a[k],v,(k,a,e))

    def test_population_and_manifest(self):
        manifest=json.loads((ROOT/'data_manifest.json').read_text())
        by_id={m['Image_Index']:m for m in manifest}
        self.assertEqual(len(self.data['images']),115)
        self.assertEqual(len(self.puncta),7152)
        self.assertEqual(sum(p['Cell_ID']>0 for p in self.puncta),6507)
        self.assertEqual(sum(p['Cell_ID']==0 for p in self.puncta),645)
        self.assertEqual([r['Nonempty_regions'] for r in self.summary],[73,56,56,51,71,67,56,75])
        for im in self.data['images']:
            for key in ['Group','SHA256','Raw_relative_path']:
                self.assertEqual(im[key],by_id[im['Image_Index']][key])

    def test_exact_components_unchanged(self):
        exact=json.loads((ROOT/'results/final_20260930/frozen_exact_measurements.json').read_text())
        by_id={p['Unique_ID']:p for p in exact['regions']}
        self.assertEqual(len(by_id),len(self.puncta))
        for p in self.puncta:
            e=by_id[p['Unique_ID']]
            for key in ['Image_Index','Group','Region_ID','Area_pixels','Integrated_density_raw']:
                self.assertEqual(p[key],e[key])
            self.assertAlmostEqual(p['Area_um2'],e['Area_um2'],places=11)

    def test_zero_puncta_not_zero_means(self):
        empty=[c for c in self.cells if c['Puncta_count']==0]
        self.assertEqual(len(empty),32)
        for c in empty:
            for key in METRICS: self.assertIsNone(c[key])
            self.assertEqual(c['Sum_punctum_area_um2'],0)
            self.assertEqual(c['Sum_punctum_integrated_intensity_raw'],0)
        self.assertEqual(sum(len(self.arrays[g][METRICS[0]]) for g in self.groups),505)

    def test_integer_raw_sums_fractional_means(self):
        self.assertTrue(all(p['Integrated_density_raw']==int(p['Integrated_density_raw']) for p in self.puncta))
        self.assertTrue(any(c[METRICS[1]] is not None and c[METRICS[1]]!=int(c[METRICS[1]]) for c in self.cells))
        self.assertEqual(sum(c['Sum_punctum_integrated_intensity_raw'] for c in self.cells),
                         sum(p['Integrated_density_raw'] for p in self.puncta if p['Cell_ID']))

    def test_duplicate_and_unknown_owner_rejected(self):
        duplicate=copy.deepcopy(self.data);duplicate['puncta'].append(duplicate['puncta'][0])
        with self.assertRaisesRegex(ValueError,'Duplicate punctum'): calculate(duplicate,self.groups)
        unknown=copy.deepcopy(self.data);unknown['puncta'][0]['Cell_ID']=99999
        with self.assertRaisesRegex(ValueError,'missing cell'): calculate(unknown,self.groups)

    def test_invalid_values_and_population_rejected(self):
        for key,value in [('Area_um2',float('nan')),('Area_pixels',0),('Integrated_density_raw',1.5)]:
            bad=copy.deepcopy(self.data);bad['puncta'][0][key]=value
            with self.assertRaises(ValueError): calculate(bad,self.groups)
        bad=copy.deepcopy(self.data);bad['images'][0]['Assigned_puncta']+=1
        with self.assertRaisesRegex(ValueError,'counts mismatch'): calculate(bad,self.groups)

    def test_export_roundtrip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)/'tables';export(self.data,self.groups,folder)
            for path in folder.glob('*.csv'):
                self.assertEqual(path.read_bytes(),(ROOT/'results/per_cell/tables'/path.name).read_bytes())
            for key in METRICS:
                with (folder/(key+'_wide.csv')).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
                for g in self.groups: self.assertEqual([float(r[g]) for r in rows if r[g]],self.arrays[g][key])
            with self.assertRaises(FileExistsError): export(self.data,self.groups,folder)

    def test_legacy_punctum_input_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'per_cell_measurements'):
            calculate({'regions':[],'images':[]},self.groups)


if __name__=='__main__': unittest.main(verbosity=2)
