"""Read-only checks of frozen measurements and portable statistics."""
from pathlib import Path
import csv, json, math, sys, unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from summarize import calculate, student
from student_distribution import student_two_sided_p

class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'results/final_20260930/frozen_exact_measurements.json').read_text(encoding='utf-8'))
        cls.expected=json.loads((ROOT/'tests/expected_final.json').read_text(encoding='utf-8'))
        cls.manifest=json.loads((ROOT/'data_manifest.json').read_text(encoding='utf-8'))
        cls.groups=json.loads((ROOT/'parameters/batch_20260930.json').read_text())['group_order']
        cls.summary,cls.tests,cls.arrays=calculate(cls.data,cls.groups)

    def test_population_and_manifest(self):
        self.assertEqual(len(self.manifest),115)
        self.assertEqual(len(self.data['regions']),7152)
        self.assertEqual(len(self.data['images']),115)
        self.assertEqual(sum(i['Nuclei'] for i in self.data['images']),537)
        self.assertTrue(all(i['Included']==1 for i in self.manifest))
        self.assertEqual({i['Image_Index'] for i in self.data['images']},{i['Image_Index'] for i in self.manifest if i['Included']})
        self.assertTrue(all(len(i['SHA256'])==64 for i in self.manifest))

    def test_final_means_counts_and_tests(self):
        for s,e in zip(self.summary,self.expected['Summary']):
            self.assertEqual(s['Group'],e['Group'])
            self.assertEqual(s['Cy5_regions'],e['Regions_after'])
            self.assertEqual(s['Automated_nuclear_regions'],e['Nuclei_after'])
            for k in ['Area_um2','Integrated_density_raw']:
                for measure in ['mean','median','SD']:
                    self.assertTrue(math.isclose(s[k+'_'+measure],e[k+'_'+measure],rel_tol=1e-12))
        self.assertEqual(len(self.tests),20)
        for t,e in zip(self.tests,self.expected['Comparisons']):
            self.assertEqual((t['Metric'],t['Group_A'],t['Group_B']),(e['Metric'],e['Group_A'],e['Group_B']))
            self.assertTrue(math.isclose(t['P_raw'],e['P_raw'],rel_tol=1e-9,abs_tol=1e-14),(t,e))
            self.assertEqual(t['Stars_raw'],e['Stars_raw'])
            for k,old in [('P_Holm_within_metric','P_Holm_10'),('P_Holm_all_metrics','P_Holm_20')]:
                self.assertTrue(math.isclose(t[k],e[old],rel_tol=1e-9,abs_tol=1e-14))

    def test_geometry_and_intensity_definitions(self):
        ims={i['Image_Index']:i for i in self.data['images']}
        for r in self.data['regions']:
            self.assertAlmostEqual(r['Area_um2'],r['Area_pixels']*ims[r['Image_Index']]['Pixel_area_um2'],places=11)
            self.assertAlmostEqual(r['Integrated_density_raw'],r['Mean_raw_intensity']*r['Area_pixels'],places=7)
            self.assertLessEqual(r['Area_pixels'],r['Original_pixels'])

    def test_distribution_and_zero_variance(self):
        for t in [0.01,1,10]:
            self.assertAlmostEqual(student_two_sided_p(t,1),1-2*math.atan(t)/math.pi,places=12)
        self.assertEqual(student([1,1,1],[1,1,1]),(0.0,4,1.0))
        self.assertEqual(student([1,1,1],[2,2,2]),(None,4,0.0))

    def test_wide_tables(self):
        for k in ['Area_um2','Integrated_density_raw']:
            path=ROOT/'results/final_20260930/tables'/(k+'_wide.csv')
            with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
            for g in self.groups:
                self.assertEqual([float(r[g]) for r in rows if r[g]!=''],self.arrays[g][k])

if __name__=='__main__':unittest.main(verbosity=2)
