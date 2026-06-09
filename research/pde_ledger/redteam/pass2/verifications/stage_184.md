---
unit_id: 184
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 184

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Output-only finding; no source authoring. Both committed transcripts were regenerated
(`scripts/output/...sympy_audit.txt` and `mathematica/output/...mathematica_audit.txt`,
mtime Jun 09 00:24, newer than both scripts at Jun 08 22:43-22:44). Line 3 of each now reads
`STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES`.

**Assessment:**
Correct. The stale pre-renumber `STAGE 167` banner is gone from both committed `.txt`; all
PASS lines present. Fresh-run trigger satisfied.

### F2 — tautological_check (both engines)

**Classification:** resolved

**What changed:**
- SymPy (py:59-77): `Rtarget = Rtarget0*exp(small*R1)` — now an INDEPENDENT perturbed object
  with its own free drift `R1`; the complement is routed through a separately-named identity
  object `selected_branch_E = 1 - eps_eta_var`. The product-drift residual
  `dln_Rtarget_T2 - dln_selected_branch_E` is printed BEFORE applying the branch law.
- Mathematica (wl:46-69): `rTarget` enters `branchVars`/`branchVelocities` as its own variable
  with velocity `r1*rTarget`; `selectedBranchClosedForm = 1 - epsEta` is separately named;
  `productDriftResidual = logDrift[rTarget*tShape] - logDrift[selectedBranchClosedForm]`.

**Assessment:**
Non-tautological in BOTH engines, confirmed from the printed pre-law residuals (NOT a literal
`0` from identical-term cancellation):
- sympy .txt:16 → `(-Sigma_eta*eps_eta + (R_1 + Xi_1)*(eps_eta - 1))/(eps_eta - 1)`
- math .txt:16 → `r1 - (epsEta*sigmaEta)/(-1 + epsEta) + xi1`
Both are genuine symbolic functions of R1/Ξ₁/Σ_η/ε_η that can fail for arbitrary `R1`; the
residual vanishes only after substituting the branch law `R1 = -Xi1 - eps_eta/(1-eps_eta)·Σ_η`.
`Rtarget` is no longer `Lam0*(1-eps_eta_var)/T2`, so the F2 self-reference is removed in both.
Complement check (A5) preserved and still PASS.

### F3 — mathematica_transliteration (user-authorized re-author)

**Classification:** resolved

**What changed:**
The `.wl` was re-authored (wl:46-91). `SeriesCoefficient[Log[ratio],{small,0,1}]` is entirely
gone (grep: 0 matches). Drifts are now obtained by an independent first-variation construction:
`branchVars={rTr,tShape,epsEta,rTarget}`, `branchVelocities={theta1*rTr, xi1*tShape,
sigmaEta*epsEta, r1*rTarget}`, `firstVariation[expr]=(D[expr,#]&/@branchVars).branchVelocities`
(wl:49-52), `logDrift[expr]=firstVariation[expr]/expr` (wl:54-57). D2 `logDrift[tTr]` (wl:74),
D3 `logDrift[nTr]` (wl:80), D4 `logDrift[epsEta]` (wl:85) all route through this directional-
derivative / chain-rule field.

**Assessment:**
This is a genuinely distinct route — a Lie-derivative-along-velocity-field linearization —
versus py's unchanged exponential-ansatz + `series(log(ratio)).coeff(small,1)` (py:66,82,89,94).
Not a relabel or a port: the engines now disagree in method while agreeing on every residual.
All `expectZero` checks PASS, banner `STAGE 184`, `math -script` exit 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `product-drift residual before branch law = (-Sigma_eta*eps_eta + (R_1 + Xi_1)*(eps_eta - 1))/(eps_eta - 1)` (symbolic, falsifiable)
- `R_target T^2 drift - delta ln(1 - eps_eta) = 0`; `delta ln T_* - Sigma_tr = 0`; `delta ln N_* - Sigma_nt = 0`; `delta ln eps_eta - Sigma_eta = 0`; `selected-branch complement identity = 0`

**Mathematica:** exit=0. Notable lines:
- `product-drift residual before branch law = r1 - (epsEta*sigmaEta)/(-1 + epsEta) + xi1` (symbolic, falsifiable)
- `PASS: R_target T^2 drift - delta ln(1 - eps_eta)`; `PASS: delta ln T_* - Sigma_tr`; `PASS: delta ln N_* - Sigma_nt`; `PASS: delta ln eps_eta - Sigma_eta`; `PASS: selected-branch complement identity`
- `Stage 184 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` mtimes (Jun 09 00:24:31) are newer than their
scripts (`.py` Jun 08 22:43:58, `.wl` Jun 08 22:44:47). Line 3 of each reads `STAGE 184 …`.

## Material-change assessment

`material_change`: false. Every deliverable value is preserved (B_*, C_*, the three branch-
coordinate laws, the complement drift) — the diff changed only the verification *method* (F2
de-toothing + F3 re-author) and the `.py` change is confined to the F2 product-identity block
(py:56-77 plus the `Ecomp` rename at py:98). No derived result that downstream units depend on
changed. The 9 deliverable laws/constants reconcile identically to the pre-fix outputs.

## Side observations (non-blocking)

- The Mathematica complement uses the closed form `selectedBranchClosedForm = 1 - epsEta`
  (background-level), while the product-drift check varies `rTarget` independently; both are
  consistent with the notes' two-object routing the directive asked for. No issue.
- `R1_branch_law` (py:73) / `targetDriftLaw` (wl:66) are hand-supplied to *close* the check;
  the falsifiability is demonstrated by the printed pre-substitution residual, which is correct.

## Verdict justification

All three findings are resolved. F2 is genuinely non-tautological in both engines: `R_target`
enters as an independent perturbed symbol and the printed pre-law product-drift residual is a
live symbolic function of R1/Ξ₁/Σ_η that vanishes only under the branch law, with the
complement routed through a separately-named `(1-ε_η)` object. F3's `.wl` is re-authored to a
first-variation (directional-derivative) construction with `SeriesCoefficient[Log[…]]` fully
removed, so the two engines are now method-independent (173-lesson satisfied). F1's transcripts
are fresh with the `STAGE 184` banner. Both engines exit 0, all checks PASS, deliverable values
preserved (method-only change). Verdict: verified.
