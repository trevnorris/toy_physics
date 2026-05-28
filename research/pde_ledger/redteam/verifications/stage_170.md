---
unit_id: 170
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: true
---

# Verification — unit 170

Three findings were tracked through resolution: F1 (paper_misalignment, user-resolved to
direction (a): add the missing Sec. 5 check to both engines, no paper edit), F2
(mathematica_transliteration), and F3 (cosmetic banner relabel — promoted from the original
report's verdict-justification note into a numbered directive item). All three were applied by
the orchestrator per the directive's `## Applied:` blocks. I checked the current file state, the
git diff vs HEAD, and the post-fix exec logs.

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
A new "Weak-axisymmetric signature (1, 1/2, -1) and scalar amplitudes" block was added to BOTH
engines:
- SymPy `scripts/...sympy_audit.py:131-166` (Section 5).
- Mathematica `mathematica/...mathematica_audit.wl:129-159` (Section 5).
The block defines the verified linear outlet maps as helpers (`kappa_map`/`gamma_map` in py L144-148;
`kappaMap`/`gammaMap` in wl L138-139), the standalone amplitude literals
`kappa1 = 3(1-σ)(D2_1 + D0_1/9)/(σ D0)` and `gamma1 = -(1-σ)(N0_1 - P0 D0_1)/(9 σ N0)`
(py L150-151 / wl L140-141), feeds lane-scaled bundle defects `δD_(A,n)=eps·λ_A·D_n^(1)`,
`δN_(A,0)=eps·λ_A·N_0^(1)` with `λ=(1,1/2,-1)` through those maps, and asserts both the closed-form
amplitudes (`δκ_W^(2A)=eps·λ_A·κ1`, `δγ_W^(2A)=eps·λ_A·γ1`) and the bare signature ratios
(`21 = (1/2)·20`, `22 = -20` for both κ and γ). +10 checks per engine, matching the `## Applied: F1`
summary.

**Assessment:**
Correct and addresses the finding. The check is non-tautological. The amplitude literals
`kappa1`/`gamma1` are written independently of the maps (py L150-151 are bare expressions, not
calls to `kappa_map`). The maps themselves are anchored upstream: `kappa_map`'s literal form
`3(1-σ)(dD2+dD0/9)/(σ D0)` is the exact target verified in Section 2 (py L73-76 against the
derived `dkappa_from_du2`, which comes from the genuine eps-linearization of `u2=-D2/D0`); likewise
`gamma_map` is anchored at py L77-80 / wl L81-84. So a wrong amplitude literal — e.g. `D0_1/8`
instead of `D0_1/9` in `kappa1` while the map kept `/9` — would leave a non-vanishing residual and
trip `expect_zero`/`expectZero`. The signature-ratio checks (`dkW[22] + dkW[20]`, etc.) directly
exercise the lane-factor propagation. The exec logs show all 10 new lines per engine at `= 0` /
`PASS`. The only nit (cosmetic, non-blocking): the SymPy f-string label prints the literal word
"lambda" for every lane, so log lines 44-49 read identically for 20/21/22 — but each is a distinct
assertion with the correct `lam` substituted, and the `.wl` labels them distinctly. No effect on
correctness.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
In `mathematica/...mathematica_audit.wl`:
- L52-56: the `Series`+`Coefficient`+`Normal` linearization was replaced by the derivative route
  `du2/du4/dP0 = FullSimplify[(D[*Full, eps] /. eps -> 0), ...]` — a distinct mechanism from the
  SymPy `series().coeff(eps,1)` it previously mirrored.
- L66-74: the SymPy-only `du2sym`/`dP0sym` dummy-symbol Solve indirection was removed; inversion is
  now direct `Solve[du2 == du2Hyb, dkappa]` / `Solve[dP0/P0 == dP0OverP0Hyb, dgamma]`.

**Assessment:**
Matches the directive's required change exactly. Confirmed by grep: no `Series`/`Coefficient`
remain in the `.wl`; the only `du2sym`/`dP0sym` occurrence is in an explanatory comment (L66), not
live code. All `expectZero` targets are unchanged (still the paper's boxed formulas). The exec log
shows every check PASS and the trace/anomaly residuals still match the SymPy ones (e.g.
`a_kappa from map = -1/3*((aD0 + 9*aD2)*(-1 + sigma))/(D0*sigma)`). The second engine now reaches
the same boxed targets by an independent linearization-and-inversion path, satisfying the
second-engine-independence policy.

### F3 — cosmetic stage-label correction (low)

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:30`: `banner("STAGE 153 — ...")` → `banner("STAGE 170 — ...")`.
- `mathematica/...mathematica_audit.wl:27`: `banner["STAGE 153 — ...")` → `banner["STAGE 170 — ...")`.

**Assessment:**
Correct in both engines; both transcripts now print `STAGE 170 — LINEAR GROUPED-P2 DIRECT OUTLET
MAP` (output line 3 in each log). No effect on any assertion.

## Exec log assessment

**SymPy:** exit=0 (inferred — every residual prints `= 0`, the run reaches the final
carry-forward block, and `expect_zero` raises `AssertionError` on any nonzero, so a clean
end-to-end run implies exit 0; grep for error/traceback/assertion returned nothing). Notable lines:
- `delta u2 + (dD2 + dD0/9)/D0 = 0`
- `delta kappa_W - 3(1-sigma)(dD2 + dD0/9)/(sigma D0) = 0`
- `kappa signature: 21 = (1/2) 20 = 0` / `kappa signature: 22 = -20 = 0`
- `gamma signature: 21 = (1/2) 20 = 0` / `gamma signature: 22 = -20 = 0`

**Mathematica:** exit=0 (the `.wl` ends with `Exit[0]`; `expectZero` calls `fail[...]`→`Exit[1]`
on any nonzero; grep for FAIL returned nothing, so reaching the end with all PASS implies exit 0).
Notable lines:
- `PASS: delta kappa_W - 3(1-sigma)(dD2 + dD0/9)/(sigma D0)`
- `PASS: delta kappa_W^(21) - (eps/2) kappa1` / `PASS: delta kappa_W^(22) + eps kappa1`
- `PASS: delta gamma_W^(21) - (eps/2) gamma1` / `PASS: delta gamma_W^(22) + eps gamma1`
- `PASS: kappa signature: 22 = -20` / `PASS: gamma signature: 22 = -20`

**Output freshness:** confirmed re-generated post-fix. mtimes:
- sympy script 16:00:29 < sympy output 16:10:17.
- mathematica script 16:00:04 < mathematica output 16:13:29.
Both `.txt` outputs are newer than their scripts, and both transcripts already contain the new
Section 5 block and the corrected STAGE 170 banner, so the captures are post-fix. (Note: the
`redteam/MANIFEST.yaml` mtime fields for unit 170 are stale — they predate these runs — but the
prompt directs using the orchestrator-refreshed `.txt` files, which are the source of truth here.)

## Material-change assessment

`material_change`: true.

The edit adds new derived assertions (the Sec. 5 weak-axisymmetric signature and the κ1/γ1 closed
forms). It does not change any pre-existing derived result — Sections 1-4 outputs are byte-identical
to before, and the F2 rework reaches the same boxed targets by a different path. The new material is
the weak-axisymmetric (1,1/2,-1) deliverable that downstream units may reference. Per the directive's
user-resolution note, stages 171/173 verify the same lane *pattern* on other quantities
(obstructions / load) rather than re-using stage 170's κ_W/γ_W amplitudes, so no specific downstream
unit is expected to break. The orchestrator will mark units > 170 `upstream_stale: true`; I see no
narrow concern that warrants a targeted re-audit beyond that default.

## Side observations (non-blocking)

- SymPy Section 5 log labels (lines 44-49) print the literal token "lambda" for all three lanes, so
  the 20/21/22 amplitude-check lines are visually indistinguishable in the transcript even though
  the assertions differ. Purely cosmetic; the Mathematica side labels them distinctly. Not raised as
  a finding.
- `redteam/MANIFEST.yaml` unit-170 mtimes are stale relative to the refreshed outputs. Bookkeeping
  only; does not affect this verification.

## Verdict justification

All three findings are `resolved`: F1's missing Sec. 5 check is present in both engines, is
genuinely non-tautological (amplitude literals are independent of the maps, and the maps are
anchored to the Section-2 derivation, so a wrong κ1/γ1 would fail), and both logs show the 10 new
checks passing; F2's `.wl` now uses `D[...,eps]/.eps->0` with direct `Solve` and no `du2sym`/`dP0sym`
idiom, an independent path to identical boxed targets; F3's banner now reads STAGE 170 in both
transcripts. Exec logs are fresh, complete, and show no FAIL/error; both engines exit 0. The git
diff contains only the three named changes with no collateral edits. Verdict: verified.
