# Parent Throat Action — Fail-Fast Classifier

## Purpose

The postprocessor gives raw diagnostics. The next practical layer is a machine
verdict:

```text
PASS / FAIL / INCOMPLETE
```

with explicit reasons.

That is what
[cfd_runtime_failfast.py](/home/trevnorris/Downloads/em_projected/cfd_runtime_failfast.py)
does.

---

## 1. Inputs

The classifier consumes the JSON summary produced by
[cfd_runtime_monitor_postprocess.py](/home/trevnorris/Downloads/em_projected/cfd_runtime_monitor_postprocess.py).

The key metrics are:

- `rms_S_rho`
- `rms_R_cont`
- `rms_R_pois_exact`
- `Q_r_tail_cv`
- `mu_eff2_tail_median`
- `alpha_fit_tail_mean`
- `alpha_fit_tail_std`

---

## 2. Default fast-screen thresholds

The default early-kill thresholds are:

\[
\frac{\|R_{\rm cont}\|_{\rm rms}}{\|S_\rho\|_{\rm rms}} \le 0.05,
\qquad
\frac{\|R_{\rm Pois}^{\rm exact}\|_{\rm rms}}{\|S_\rho\|_{\rm rms}} \le 0.05,
\]
\[
{\rm CV}(Q_r^{\rm tail}) \le 0.05,
\qquad
|\mu_{\rm eff,tail}^2| \le 0.25,
\]
\[
|\alpha_{\rm fit,tail}-2| \le 0.1,
\qquad
{\rm std}(\alpha_{\rm fit,tail}) \le 0.1.
\]

These are deliberately fail-fast rather than publication-grade tolerances.

For the \(Q_r\) plateau check, the CLI also supports a resolution-derived
threshold:

```bash
python cfd_runtime_failfast.py summary.json \
  --q-tail-cv-max-from-resolution
```

This flag is mutually exclusive with an explicit `--q-tail-cv-max` override;
the CLI raises instead of silently choosing one threshold source.

That mode reads `tail_n_points` from the monitor summary and sets

\[
{\rm CV}(Q_r^{\rm tail})_{\max}
=\sqrt{\frac{2}{\max(N_{\rm tail},4)}}.
\]

The heuristic treats the exterior tail shells as near-stationary independent
samples. A mean or CV-type tail statistic has sampling noise that scales like
\(1/\sqrt{N_{\rm tail}}\); applying a factor-of-2 safety margin gives the
\(\sqrt{2/N_{\rm tail}}\) factor. The \(\max(N_{\rm tail},4)\) clamp avoids
over-tightening or over-trusting very short tails. This is still a fail-fast
rule, but it ties at least one tolerance to the number of shell samples in the
exterior readout instead of using only a fixed magic number.

If the source scale is numerically zero, the classifier skips the normalized
projection test and warns that the projection channel was not assessed on that
snapshot. By default that makes the verdict `INCOMPLETE`; a synthetic exterior
test may opt out by marking the snapshot as source-free exterior-only.

If the optics channel is missing entirely, the verdict is `INCOMPLETE` unless a
different failure already triggered.

---

## 3. Synthetic self-test

[step_35_parent_throat_action_failfast_classifier_sympy.py](/home/trevnorris/Downloads/em_projected/step_35_parent_throat_action_failfast_classifier_sympy.py)
drives the classifier on the synthetic cases below plus one CLI-precedence
guard:

1. **Newton-like exterior**: should `PASS`.
2. **Yukawa exterior**: should `FAIL`.
3. **Bad optics exterior**: should `FAIL`.
4. **Projection-broken snapshot**: should `FAIL`.
5. **Missing optics snapshot**: should `INCOMPLETE`.
6. **Near-zero source snapshot**: should `INCOMPLETE`.
7. **Default source-free exterior**: should `INCOMPLETE`, proving that the
   source-free opt-out must be explicit.
8. **Q_r threshold-boundary snapshot**: should `PASS`.
9. **Q_r just-outside-threshold snapshot**: should `FAIL`.
10. **mu_eff2 threshold-boundary snapshot**: should `PASS`.
11. **mu_eff2 just-outside-threshold snapshot**: should `FAIL`.
12. **Resolution-derived Q_r threshold-boundary snapshot** at
    `tail_n_points = 200`: should `PASS`.
13. **Resolution-derived Q_r just-outside snapshot**: should `FAIL`.
14. **CLI Q_r threshold conflict guard**: passing both `--q-tail-cv-max` and
    `--q-tail-cv-max-from-resolution` should raise.

So the classifier is not just checking that “some numbers exist.” It is
explicitly discriminating the failure classes we care about.

---

## 4. CLI

```bash
python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json
python cfd_runtime_failfast.py summary.json --output-json verdict.json
python cfd_runtime_failfast.py summary.json \
  --q-tail-cv-max-from-resolution \
  --output-json verdict_resolution.json
```

That is now the shortest operational path from a real CFD dump to a direct
gravity-analog falsification verdict.
