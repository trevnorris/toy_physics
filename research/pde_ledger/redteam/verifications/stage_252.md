---
unit_id: 252
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T20:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 252

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New independent-route `.wl` created at `mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl` (275 lines). Covers M1–M9 with 25 `expectZero`/`expectApprox` checks (28 PASS lines incl. the M1 triple), each aborting via `Exit[1]` on failure.

**Assessment:**
Genuinely independent, not a transliteration. Uses native primitives and a distinct decomposition: M1 derives fractions via `Cancel[Together[(Evac/Eexp) /. {I3->rV I2}]]` and compares to independently-written closed forms `fvacExpected`/`flatExpected`; M2 via `D[flat,rV]`+`Together`; M3 via `Limit` with `stripConditional` (project idiom applied); M4 directly `Integrate[D[Vin Exp[s t],{t,k}]^2,{t,0,T}]` with positivity assumptions; M6/M7 build `safeCombo` and apply the safe rule `safeCombo -> mueta(s0^2-sc^2)` rather than asserting a self-identity; M8 ties `Numerator[Cancel[Together[flat/.rV->sc^2 - 3/4]]]` to the surface (`threeOneFactor=1`, pinned in-script). Symbol assumptions match the directive (Γ's ≥0; I2,rV,s,T,Vin,sc,s0,mueta>0; `flat<1` not assumed — emerges). No partition/conservation identity asserted by construction: M1 partition sum is computed from independently-built fractions then reduced. M4/M6/M7 carry the integrals and safe rule through, not by construction.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py:141-148`. The old `assert simplify(gamma_safe_eq - mu_eta*(s0**2-sc**2)/sc)==0` (with `gamma_safe_eq` *defined* at L122 as exactly that RHS — a pure X−X) is replaced by: build `gamma_safe_eq_bridged = simplify(((G3*sc**3+G5*sc**5)/sc).subs(G3*sc**3+G5*sc**5, mu_eta*(s0**2-sc**2)))`, then `assert simplify(gamma_safe_eq_bridged - mu_eta*(s0**2-sc**2)/sc)==0`. L140 left unchanged.

**Assessment:**
Non-tautological and falsifiable. The new residual references the Γ coefficients AND mu_eta/s0/sc through the safe-equality substitution. The substitution provably fired: had `.subs` not matched, `gamma_safe_eq_bridged` would reduce to `G3*sc^2+G5*sc^4`, whose difference from `mu_eta*(s0^2-sc^2)/sc` (disjoint free symbols) is non-zero — the assert would have failed, but the script exits 0, so the safe equality is genuinely carried. Mathematica M7 corroborates independently: `safeRateRaw - (Gamma3tot sc^2+Gamma5tot sc^4)==0` AND `safeRateReduced - mueta(s0^2-sc^2)/sc==0` under `safeRule`.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py:171-178`. The old `assert simplify(split_surface) == G3l+G5l*sc**2-3*G3v-3*G5v*sc**2` (expand(X)==expand(X)) is replaced by: `surface = expand((G3l+G5l*sc**2)-3*(G3v+G5v*sc**2))`, `resid = together(fl_sc - Rational(3,4))`, `num = numer(cancel(resid))`, `assert simplify(num - surface)==0`. `split_surface` (L155) is retained only for the print, no longer load-bearing.

**Assessment:**
The iff is now genuinely tied. `fl_sc` is the actual partition fraction `f_lat(s_c)` (L124, `fl.subs(rV,sc**2)`, depending on G3l,G5l,G3v,G5v,sc); the numerator of `f_lat(s_c)-3/4` is asserted equal (factor 1) to the 3:1 surface. Removing the `-3*(...)` term would break it — not a self-expansion. Mathematica M8 reproduces the same route (`threeOneNumerator - 1*threeOneSurface == 0`).

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py:180-201`. Added the vacuum φ-family check `fv_phi = simplify(fv.subs({G3l:phi*G3T, G3v:(1-phi)*G3T, G5l:phi*G5T, G5v:(1-phi)*G5T}))` with `assert simplify(fv_phi-(1-phi))==0`, set `phi_val = Rational(3,4)`, and replaced `frac_v_num=0.25 / frac_l_num=0.75` with `frac_v_num=float(1-phi_val) / frac_l_num=float(phi_val)`.

**Assessment:**
The 0.25/0.75 weights are now derived from the proven φ=3/4 family rather than typed in. Benchmark inputs (`t_cross_num`, `s0_num`, `E_diss_num`) and all literal assert targets (L223-232) unchanged, so the calibration still closes on the same paper values. Mathematica M9 is even stronger: `fracVacNum = N[fvacPhi /. phi->phiVal]`, `fracLatNum = N[flatPhi /. phi->phiVal]` — sourced from the symbolic φ-family fractions themselves.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "All symbolic and numerical checks passed."; `f_vac + f_lat = 1`; `gamma_eq_safe = mu_eta*(s_0**2 - s_c**2)/s_c`; benchmark `E_vac=0.00258365`, `E_lat=0.0077509...` (unchanged from pre-fix targets).

**Mathematica:** exit=0. Notable lines: 28 PASS across M1–M9; `PASS: M7 safe rate after safe rule`; `PASS: M8 three-to-one numerator is exact positive factor times surface`; `fracVacNum = 0.25`, `fracLatNum = 0.75` derived from `phiVal = 3/4`; all M9 diffs ≤ tol (max residual 2.84e-14 vs 1e-9 tol on the prefactor).

**Output freshness:** confirmed. Both `.txt` outputs (mtime 2026-06-03 13:16:42) are newer than both scripts (`.py` 13:10:27, `.wl` 13:13:38). Outputs regenerated post-fix.

## Material-change assessment

`material_change`: false.

All numeric and symbolic results are unchanged. The edits replaced three weak/vacuous checks (F2 X−X, F3 self-expansion, F4 hardcoded weights) with genuine falsifiable ones and added an independent second engine. No derived constant consumed downstream changed; in particular `gamma_lat^eq(s_c)=Γ₃ˡs_c²+Γ₅ˡs_c⁴` (consumed by stage 253) is untouched. No downstream re-audit needed on the merits.

## Side observations (non-blocking)

- The stage's "cold survival" content is verified exclusively through equality theorems (safe-edge energy M6, safe-edge rate M7) and the benchmark; there is no strict-inequality survival-margin check (e.g. `s_c < s_0` or `f_lat > 0`) in either engine. The directive's M1–M9 manifest did not require one and the auditor raised no such finding, so this is not a gap to block on — flagging only because the instruction asked me to watch for a non-strict/sample-point inequality (none is present to be weakened).
- `scripts/...py:123` defines `fv_sc = simplify(fl.subs(rV, sc**2))` with a comment "intentionally compute lattice fraction here first" — `fv_sc` is the lattice fraction (mislabeled) and is unused; the real vacuum fraction is `fvac_sc` (L125). Harmless dead variable, pre-existing, not part of any finding.

## Verdict justification

All four findings are resolved. F2's X−X tautology is replaced by a substitution that genuinely carries the Stage-251 safe equality and is provably falsifiable (the substitution fired, else the assert would fail). F3's expand(X)==expand(X) is replaced by a check tying the real partition fraction `f_lat(s_c)−3/4`'s numerator to the 3:1 surface. F4's hardcoded 0.25/0.75 are derived from the proven φ=3/4 family. F1 adds a genuinely independent Mathematica route (native `Cancel`/`Together`/`D`/`Limit`/`Integrate`, safe-rule substitution, iff-via-numerator) with 28 PASS, no by-construction conservation and no weakened inequality. Both engines exit 0, outputs are fresh, and no downstream-consumed value changed (material_change false).
