"""Export per-cell means and punctum assignments. Python >= 3.9, stdlib only.

One nonempty nucleus-associated territory contributes one observation to each
mean-punctum distribution. No significance tests are run.
"""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics as st

METRICS = ('Mean_punctum_area_um2', 'Mean_punctum_integrated_intensity_raw')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value, name, minimum=0):
    require(isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value >= minimum, f'Invalid {name}: {value!r}')


def integer(value, name, minimum=0):
    finite(value, name, minimum)
    require(int(value) == value, f'{name} must be an integer')


def calculate(data, group_order):
    require('cells' in data and 'puncta' in data,
            'Expected per_cell_measurements.json. Run all three Fiji stages first.')
    require(len(set(group_order)) == len(group_order), 'Duplicate group names')
    images = {i['Image_Index']: i for i in data['images']}
    require(len(images) == len(data['images']), 'Duplicate image ID')
    for i, im in images.items():
        integer(i, 'Image_Index', 1)
        require(im['Group'] in group_order, f'Unknown group: {im["Group"]}')
    cells = {}
    for cell in data['cells']:
        key = cell['Image_Index'], cell['Cell_ID']
        integer(key[1], 'Cell_ID', 1)
        require(key[0] in images and key not in cells, 'Missing image or duplicate cell ID')
        cells[key] = cell
    by_cell, by_image = defaultdict(list), defaultdict(list)
    puncta, ids, region_ids = [], set(), set()
    for p in data['puncta']:
        image_id, cell_id = p['Image_Index'], p['Cell_ID']
        integer(cell_id, 'Cell_ID')
        integer(p['Region_ID'], 'Region_ID', 1)
        require(image_id in images, 'Punctum refers to missing image')
        require(cell_id == 0 or (image_id, cell_id) in cells, 'Punctum refers to missing cell')
        require(p['Unique_ID'] not in ids and (image_id, p['Region_ID']) not in region_ids,
                'Duplicate punctum assignment')
        ids.add(p['Unique_ID']); region_ids.add((image_id, p['Region_ID']))
        integer(p['Area_pixels'], 'Area_pixels', 1)
        integer(p['Integrated_density_raw'], 'Integrated_density_raw')
        finite(p['Area_um2'], 'Area_um2')
        im = images[image_id]
        expected_area = p['Area_pixels'] * im['Pixel_width_um'] * im['Pixel_height_um']
        require(math.isclose(p['Area_um2'], expected_area, rel_tol=1e-10), 'Pixel area mismatch')
        by_image[image_id].append(p)
        if cell_id:
            by_cell[image_id, cell_id].append(p)
        puncta.append(dict(p, Group=im['Group'], Cell_key=f'{image_id:04d}_C{cell_id:03d}' if cell_id else '',
                           Assignment='assigned' if cell_id else 'unassigned'))
    results = []
    for (image_id, cell_id), cell in sorted(cells.items(), key=lambda x: (group_order.index(images[x[0][0]]['Group']), *x[0])):
        points = by_cell[image_id, cell_id]
        n = len(points)
        area = sum(p['Area_um2'] for p in points)
        intensity = sum(int(p['Integrated_density_raw']) for p in points)
        pixels = sum(p['Area_pixels'] for p in points)
        results.append(dict(cell, Group=images[image_id]['Group'], Cell_key=f'{image_id:04d}_C{cell_id:03d}',
            Puncta_count=n, Mean_punctum_area_um2=area/n if n else None,
            Mean_punctum_integrated_intensity_raw=intensity/n if n else None,
            Sum_punctum_area_um2=area, Sum_punctum_integrated_intensity_raw=intensity,
            Mean_intensity_over_puncta_pixels_raw=intensity/pixels if pixels else None))
    for image_id, im in images.items():
        image_cells = [c for c in results if c['Image_Index'] == image_id]
        image_puncta = by_image[image_id]
        assigned = sum(p['Cell_ID'] > 0 for p in image_puncta)
        require(len(image_cells) == im['Nuclear_regions'], 'Nuclear population mismatch')
        require(len(image_puncta) == im['All_puncta'], 'Punctum population mismatch')
        require(assigned == im['Assigned_puncta'] and len(image_puncta)-assigned == im['Unassigned_puncta'],
                'Assigned/unassigned counts mismatch')
    require(sum(c['Sum_punctum_integrated_intensity_raw'] for c in results)
            == sum(int(p['Integrated_density_raw']) for p in puncta if p['Cell_ID']), 'Intensity conservation failed')
    arrays, summary = {}, []
    for group in group_order:
        cs = [c for c in results if c['Group'] == group]
        ims = [im for im in images.values() if im['Group'] == group]
        arrays[group] = {k: [c[k] for c in cs if c[k] is not None] for k in METRICS}
        row = dict(Group=group, Images=len(ims), Nuclear_regions=len(cs),
                   Nonempty_regions=sum(c['Puncta_count'] > 0 for c in cs),
                   Zero_puncta_regions=sum(c['Puncta_count'] == 0 for c in cs),
                   Assigned_puncta=sum(c['Puncta_count'] for c in cs),
                   Unassigned_puncta=sum(im['Unassigned_puncta'] for im in ims))
        for key in METRICS:
            values = arrays[group][key]
            row.update({key+'_mean': st.mean(values) if values else None,
                        key+'_median': st.median(values) if values else None,
                        key+'_SD': st.stdev(values) if len(values) > 1 else None})
        summary.append(row)
    return results, puncta, summary, arrays


def write_csv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def export(data, groups, folder):
    cells, puncta, summary, arrays = calculate(data, groups)
    folder.mkdir(parents=True, exist_ok=False)
    write_csv(folder/'Per_cell_results.csv', cells)
    write_csv(folder/'Puncta_assignments.csv', puncta)
    write_csv(folder/'Group_summary.csv', summary)
    write_csv(folder/'Images.csv', data['images'])
    for key in METRICS:
        length = max((len(arrays[g][key]) for g in groups), default=0)
        write_csv(folder/(key+'_wide.csv'),
                  [{g: arrays[g][key][i] if i < len(arrays[g][key]) else '' for g in groups}
                   for i in range(length)], groups)
    return cells, puncta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_json', type=Path)
    parser.add_argument('output_folder', type=Path, help='New folder; will not overwrite')
    parser.add_argument('--config', type=Path, default=Path(__file__).resolve().parents[1]/'parameters/batch_20260930.json')
    args = parser.parse_args()
    data = json.loads(args.input_json.read_text(encoding='utf-8-sig'))
    groups = json.loads(args.config.read_text(encoding='utf-8-sig'))['group_order']
    cells, puncta = export(data, groups, args.output_folder)
    print(f'Exported {len(cells)} nuclear territories and {len(puncta)} puncta to {args.output_folder}')


if __name__ == '__main__':
    main()
