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

If the source scale is numerically zero, the classifier skips the normalized
projection test and warns that the projection channel was not assessed on that
snapshot. By default that makes the verdict `INCOMPLETE`; a synthetic exterior
test may opt out by marking the snapshot as source-free exterior-only.

If the optics channel is missing entirely, the verdict is `INCOMPLETE` unless a
different failure already triggered.

---

## 3. Synthetic self-test

[step_35_parent_throat_action_failfast_classifier_sympy.py](/home/trevnorris/Downloads/em_projected/step_35_parent_throat_action_failfast_classifier_sympy.py)
drives the classifier on six synthetic cases:

1. **Newton-like exterior**: should `PASS`.
2. **Yukawa exterior**: should `FAIL`.
3. **Bad optics exterior**: should `FAIL`.
4. **Projection-broken snapshot**: should `FAIL`.
5. **Missing optics snapshot**: should `INCOMPLETE`.
6. **Near-zero source snapshot**: should `INCOMPLETE`.

So the classifier is not just checking that “some numbers exist.” It is
explicitly discriminating the failure classes we care about.

---

## 4. CLI

```bash
python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json
python cfd_runtime_failfast.py summary.json --output-json verdict.json
```

That is now the shortest operational path from a real CFD dump to a direct
gravity-analog falsification verdict.
