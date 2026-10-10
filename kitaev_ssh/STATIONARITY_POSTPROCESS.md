# Post-run thermalization cut

The sampler can save every measured block by running with `--warmup 0`. The
whole requested sweep budget is sampled, so no early measurements are hidden
from later analysis. Keep the usual `--measure-every` and `--block-samples`
values, and make `--sweeps` an integer multiple of their product.

For a full momentum-grid run, for example:

```sh
./kitaev_ssh_mc ... --warmup 0 --sweeps 60000 \
  --measure-every 1 --block-samples 25 --vbs 1 --output run
```

After that seed finishes, analyze its own block file:

```sh
python3 kitaev_ssh/detect_equilibration.py run_blocks.csv \
  --window-blocks 24 --min-tail-blocks 48
```

The default diagnostics track `X_x/y/z`, matter energy, each `Szz(q)`, and a
per-block all-q maximum. The plot shows block and cumulative all-q-max curves,
the individual momentum channels, and the recommended first retained block.
The JSON records the detector settings, sweep range, and retained phase-ratio
summary. The `_equilibrated_blocks.csv` contains the retained original complex
phase and numerator sums, so it can be reblocked without losing reweighting
information. Run one file per independent seed; do not concatenate chains
before selecting their cuts.

The cut requires each selected block series' suffix mean to agree with the
final comparison window within `z` standard errors and its fitted linear trend
to be below the same `z` threshold. Defaults are `z=2`, a comparison window of
roughly one eighth of the blocks (bounded to 4–32), and a tail of at least two
windows. Pass `--observable NAME` one or more times to select the series used
for the cut. If no suffix passes, the diagnostic JSON and plot are still
written, but no retained-block file is emitted.

This screen helps choose a reproducible cut; it does not prove equilibrium.
Block means can remain correlated, the max over q is nonlinear, and different
observables can relax at different rates. Inspect the plot and repeat across
independent seeds. Treat the generated single-chain jackknife error as a
diagnostic; final uncertainties should preserve chain boundaries and use
blocking appropriate to the retained data. If a previous run used a positive
`--warmup`, its discarded configurations were never saved, so they cannot be
recovered by post-processing.
