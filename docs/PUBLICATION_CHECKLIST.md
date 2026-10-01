# Before public release

- [x] Portable source, actual batch configuration and small reproducible example.
- [x] Final labels and valid-dataset manifest.
- [x] Raw-file hashes and historical source hashes.
- [x] Final wide-format measurements and explicit statistical scope.
- [ ] Confirm all author names/order, institution and contact.
- [ ] Approve code license and data license with the laboratory/rightsholders.
- [ ] Add a valid `CITATION.cff` with confirmed authors, repository URL and
  release date. Do not use guessed author lists, ORCIDs or placeholder DOIs.
- [ ] Verify/complete biological-replicate and timepoint metadata from records.
- [ ] Deposit the 115 valid full-resolution ND2 stacks used for analysis.
- [ ] Deposit full derived masks/label TIFFs if required for the paper; only
  a small example is bundled here. Re-running the code regenerates them.
- [ ] Add verified dataset URL/DOI and GitHub URL in README and manuscript.
- [ ] Review source/data release permission and remove private local-only files.
- [ ] Tag the final release and archive its exact code version with the dataset.

Suggested raw archive layout:

```text
raw_data/
  <Group>/<current ND2 filename>
  data_manifest.csv
```

`Raw_relative_path` is relative to `raw_data`; use that directory as the
Fiji input parent. Large ND2 stacks are not part of this code package.
The public CSV/JSON tables contain no machine-specific absolute paths.

Nothing has been uploaded or published by preparation of this local folder.
