---
unit_id: 178
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 178

Light confirmation of a pass-2 CLEAN stage (0 findings, no Codex edits; only the
orchestrator's reliability re-run refreshed the committed `.txt` outputs).

## Per-finding outcomes

No findings were raised by the auditor (`findings_count: 0`, verdict `clean`), so
there is nothing to re-resolve. Per the light-confirmation protocol I re-confirmed
(a) the clean verdict and `.wl` independence, (b) deliverable residuals, (c) banners
and exit codes, and (d) an empty script diff.

### (a) Audit verdict / `.wl` independence — confirmed

The auditor flagged no `mathematica_transliteration` because the `.wl` carries a
genuinely independent route the `.py` lacks. Confirmed by direct read:

- **`.wl` (wl:105-109):** `nuFromData = FullSimplify[Coefficient[Normal[Series[Log[pA^2/dA^2], {eps,0,1}]], eps*lam], …]; expectZero["nu via log-data vs slippage", nuFromData - nuExpected]` — `nu_r` extracted as the `eps*lam` coefficient of the **log-series** of the raw ratio `P_A^2/Δ_A^2`, built straight from the component drifts `pA/dA`, bypassing `pExpected/dExpected`.
- **`.py` reaches `nu` differently:** a ratio-series `nu_from_series = (NAr_lin/N0r - 1)/(eps*lam)` (py:50-52) and `nu_direct` via the closed-form weights `2*p_expected - 2*d_expected` (py:128); there is no `Log[]`-of-ratio coefficient extraction anywhere in the `.py` (grep confirms no `log` route).

Logarithmic-derivative vs ratio-series is a distinct algebraic path → the `.wl` is
NOT a port; second-engine policy satisfied. The other three `.wl` checks are close
parallels, but the independent fourth route is sufficient.

### (b) Deliverable values — all verify

Every deliverable residual is 0 / PASS in the refreshed committed outputs:
`nu_r - 2(p_r-d_r)=0`, `Xi_1 - (nubar_N - kappa_1)=0`, `p_r formula=0`, `d_r formula=0`,
`alpha+beta-1=0`, `chi-zeta-1=0`, `nu_r - [kappa1 + sigma_r]=0`, and the `.wl`-only
`nu via log-data vs slippage=0`. Mathematica prints PASS for all eight; SymPy emits all
zero residuals and exits 0. No FAIL/ERROR/traceback in either committed `.txt`.

### (c) Banners / engines / no FAIL

Both committed outputs and both exec logs open with banner `STAGE 178 —
OUTGOING-PORT CO-LOADING THEOREM`. SymPy exit_code 0, Mathematica exit_code 0
(closing line `Stage 178 Mathematica audit passed.`). No FAIL anywhere.

### (d) Script diff empty

`redteam/pass2/exec_logs/stage_178_diff.patch` is 0 bytes; `git diff HEAD` on both
script files is empty → scripts unchanged, as expected for a clean stage.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `nu_r - 2(p_r-d_r) = 0`,
`Xi_1 - (nubar_N - kappa_1) = 0`, `p_r formula = 0` / `d_r formula = 0`,
`alpha+beta-1 = 0` / `chi-zeta-1 = 0`, `nu_r - [kappa1 + sigma_r] = 0`.

**Mathematica:** exit=0. Notable lines: `PASS: nu_r - 2(p_r-d_r)`,
`PASS: Xi_1 - (nubar_N - kappa_1)`, `PASS: nu via log-data vs slippage`,
`Stage 178 Mathematica audit passed.`

**Output freshness:** confirmed. `.txt` mtimes 2026-06-09 00:20 are newer than the
`.py`/`.wl` mtimes 2026-05-30 01:15 → outputs were regenerated post-run.

## Material-change assessment

`material_change`: false. No edits to any script; deliverable closed forms are
unchanged from the audited state. No downstream unit is affected by this re-run.

## Side observations (non-blocking)

None.

## Verdict justification

This is a clean, unmodified stage. The `.wl` retains an independent log-series route
for `nu_r` (wl:105-109) distinct from the `.py`'s ratio-series / closed-form path, so
the no-transliteration verdict holds. All eight deliverable residuals are 0 / PASS in
the freshly regenerated outputs (txt mtimes > script mtimes), both engines carry the
STAGE 178 banner and exit 0 with no FAIL, and the script diff is empty. Verdict:
**verified**.
