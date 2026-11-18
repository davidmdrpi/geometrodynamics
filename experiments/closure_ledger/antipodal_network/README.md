# Antipodal clock and two-port network

[Derivations, results and limitations](../../../docs/antipodal_clock_network.md) · [Archived runs](../runs/20261117_antipodal_network/)

Run from the repository root with Python 3.12, NumPy, SciPy, and Matplotlib. The original runs used NumPy 2.3.5 and SciPy 1.17.0. The full-domain study refines to 196,608 cover tetrahedra. The analytic compact-clock experiment and the real-time two-port experiment use different spacetimes.

```bash
export OPENBLAS_NUM_THREADS=1
module=experiments.closure_ledger.antipodal_network
output=$(mktemp -d)
for stage in clock gravity network interference refinement transport round_exterior; do
    python -m "$module.$stage" > "$output/$stage.json" || exit 1
done
python -m "$module.receiver" "$output/refinement.json" > "$output/receiver.json"
python -m "$module.figures" "$output" "$output/refinement.png"
python -m unittest discover -s tests -p test_antipodal_clock_network.py
```

The receiver uses the supplied run’s calibration. Each time-domain step balances work independently at both rims. The frequency-domain solve eliminates the interior of each region in local optical time and applies the clock phase at its attachment. The pulse and driven receiver use different boundary preparations. The sampled near-zero responses and exact-null eigenvalues depend on floating-point libraries.

| Calculation | Observable |
|---|---|
| clock, gravity | Allowed clock offsets, Maxwell waist parity, Pell pairs, continuously scalable action, linear TT and its quadratic obstruction |
| network, transport | Free antipodal tube exchange, nonorientability, separate rim work, signed histories under refinement |
| interference, refinement | Tube-only versus full-domain refinement of a calibrated dark port, including the unresolved exterior response |
| round_exterior | Round geodesic rims, consistent-mass Dirichlet poles, full uniform-port power, and a fixed-clock fine-grid test |
| receiver | Passive absorbed power, its split between the tube and the exterior, and the independently calculable product-tube response |
