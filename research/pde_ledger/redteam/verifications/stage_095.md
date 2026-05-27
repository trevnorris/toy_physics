---
unit_id: 095
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 095

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py:35-44` adds three closed-form residual asserts comparing `K0corr/K2corr/K4corr` to the notes Section 2 closed forms (`-M0**2*chi**2/G0`, `G2*M0**2*chi**2/G0**2`, `M0**2*chi**2*(G0*G4 - G2**2)/G0**3`), followed by three `PASS:` print lines. Pre-existing tautological `.has(chi**2)` asserts retained at lines 46-48 as cheap sanity flags (matches directive instruction).

**Assessment:**
Edit matches the directive's required-change block verbatim. The three new asserts pin SymPy's series-extracted coefficients to independent algebraic closed forms — non-tautological because a regression in `series(...).coeff(w, k)` would now fail. The SymPy exec log shows all three PASS lines printed. No collateral edits to logic; only the docstring carry-forward annotation (per F3 Applied block).

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl:33-44` inserts an independent Schur derivation: defines `Lq = (1/2) qSym^2 dQsym + (1/2) gSym^2 dGsym + chi m0 qSym gSym`, solves `D[Lq, gSym] == 0` for `gStar`, substitutes back, extracts `2*Coefficient[LqEff, qSym, 2]` as `dEffCoeff`, computes `corrDerived = dEffCoeff - dQsym`, and asserts via `expectZero` that `corrDerived - (-chi^2*m0^2/dGsym) == 0`. Block stands alone with symbolic `dQsym, dGsym` so does not interact with downstream series.

**Assessment:**
Matches the directive's required-change block. The derivation now uses Mathematica's own `Solve`/`Coefficient` machinery on the L action — a genuinely independent path to `D_eff` rather than a port of SymPy's choreography. Exec log line 6 shows `PASS: D_eff Schur derivation matches -chi^2 m0^2 / dG`. The existing series checks (M1–M3) and downstream computation remain intact. No collateral edits beyond the banner sweep.

### F3 — paper_misalignment

**Classification:** resolved

**What changed:**
Per the orchestrator's batch-IV1-paper-alignment Cluster A direction (a) resolution, a docstring annotation was added to `scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py:8-14` naming stage 094 as the upstream orthogonality anchor. No paper edit, no new script-side check.

**Assessment:**
Resolution path was user-selected via batch-IV1 resolution document; the script-side carry-forward annotation is in place and references the stage 094 angular integrals + Laplace eigenvalue check as the upstream proof. This is a documentation-only edit; no derived quantities change. No mathematical claim was added or removed; the annotation simply explains why the local 2-mode action already encodes orthogonality.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
```
K0corr = -M0**2*chi**2/G0
PASS: K0corr matches closed form
PASS: K2corr matches closed form
PASS: K4corr matches closed form
Stage 095 theorem verified: contamination begins at O(chi^2).
```
All three F1 PASS lines present. No failures or warnings.

**Mathematica:** exit=0. Notable lines:
```
D_eff Schur derivation matches -chi^2 m0^2 / dG = 0
PASS: D_eff Schur derivation matches -chi^2 m0^2 / dG
PASS: K0corr / chi^2 static factor
PASS: K2corr / chi^2 dynamic factor
PASS: K4corr / chi^2 dynamic factor
PASS: d c_pole / dchi |0
Stage 095 Mathematica audit passed.
```
F2 PASS line present; pre-existing closed-form residual PASS lines and `d c_pole / dchi |0` PASS retained. No failures.

**Output freshness:** confirmed. SymPy script mtime 1779902069, output mtime 1779913717 (output newer by ~11648s). Mathematica script mtime 1779902076, output mtime 1779913770. Both outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

No derived symbolic result was modified. F1 adds new asserts that pin existing extracted coefficients to closed forms already present in the notes — the coefficient values themselves are unchanged (the log still prints the same `K0corr/K2corr/K4corr` expressions). F2 inserts an independent derivation block that confirms (rather than redefines) `D_eff = D_q - chi^2 m0^2 / D_g`. F3 is documentation-only. Downstream units do not depend on any new quantity.

## Side observations (non-blocking)

- The SymPy script's lines 20-21 declare `G0` twice — first in the joint `sp.symbols('chi M0 G0 G2 G4 w', real=True)` and then again on line 21 with `nonzero=True, real=True`. The second declaration shadows the first; this is the script's pre-existing pattern (not introduced by this fix) and works correctly, but is somewhat fragile. Not blocking.
- The directive's F3 Applied block is annotated as "orchestrator-direct, post-user-resolution"; the Codex iteration only mechanically applied F1 and F2. This matches the verifier's expected scope (paper_misalignment findings are user-resolved).

## Verdict justification

All three findings are resolved. F1 introduces three non-tautological closed-form residual asserts in SymPy that mirror the Mathematica substantive checks, with all PASS lines confirmed in the exec log. F2 adds an independent Schur derivation in Mathematica using `Solve[D[Lq, gSym] == 0]` and asserts the result equals `-chi^2 m0^2 / dG`, with the corresponding PASS line confirmed. F3 was resolved via the orchestrator's batch-IV1 paper-alignment Cluster A direction (a) with a script-side docstring carry-forward annotation pointing to stage 094. Both exec logs exit 0, outputs are fresh, no regressions in diff, no material change to downstream-visible results.
