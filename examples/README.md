# Synthetic 48-sample TPM example

`synthetic_tpm_48.tsv` is fabricated test data, not a biological dataset. Species A, B and C appear together in the first 12 samples at 333333.333333 TPM each. Species D appears alone at 1000000 TPM in the remaining 36 samples. All 48 samples are retained in the analysis denominator.

```bash
python "Computing co-occurrence_Single.py" examples/synthetic_tpm_48.tsv --input-kind tpm --output results/demo
```

With the default criteria, the selected communities are AB, AC, BC and ABC: pairwise observed frequency is 12/48, expected count is 3 and fold enrichment is 4; the triple's expected count is 0.75 and fold enrichment is 16. The four-member group has zero observed occurrences and is not selected. This construction checks known arithmetic, not statistical calibration or biological accuracy on real data.

A relative-abundance threshold of 0.001 corresponds to **strictly greater than 1000 TPM**; 0.0001 corresponds to strictly greater than 100 TPM. The default minimum co-occurrence of 10 observations requires appearance in at least 10 of 48 samples (approximately 20.8%). Report any changes to these settings.

No valid SBML models are bundled. Tests create temporary XML filename placeholders only to check the runner's input contract; actual simulations require valid models and a configured solver.
