# Per-cell measurements

`per_cell_measurements.json` contains 115 images, 537 nucleus-associated territories and all 7,152 puncta with stable image, region and cell IDs. The `tables/` directory contains CSV outputs from `scripts/summarize.py`.

Of 537 territories, 505 contain puncta and contribute a mean to the area and fluorescence distributions. The remaining 32 have zero counts and sums but missing punctum means. The 645 unassigned puncta are preserved in the detail table and do not enter per-cell means.

Run `python scripts/summarize.py results/per_cell/per_cell_measurements.json NEW_OUTPUT` from the repository root to reproduce the tables. These nuclear territories are analysis regions, not measured membrane boundaries.
