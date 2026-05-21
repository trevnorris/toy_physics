---
unit_id: 011
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 011

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
Codex created
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl`
(157 lines) and produced
`/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.txt`.
Note: the auditor and directive both name the path under `scripts/`, but the
orchestrator/engine convention places `.wl` files under `mathematica/` and
their transcripts under `mathematica/output/`; the orchestrator moved the
file from Codex's initial `scripts/` location into `mathematica/`. All other
stage NNN `.wl` audits in this repo live in `mathematica/`, so the final
location is consistent with the convention.

The new `.wl` defines bundle moments via
`momentD[n,e]`/`momentN[n,e]`, builds `u2Bundle/u4Bundle/p0Bundle/p2Bundle/p4Bundle`
parametrically in `e`, and extracts first-order shifts via
`firstShift[expr_, var_] := Coefficient[Normal[Series[expr - (expr /. var -> 0), {var, 0, 1}]], var, 1]`.
It enforces M1-M11 with `check[label, expr] := If[FullSimplify[expr] =!= 0, Print["FAIL..."]; Exit[1], ...]`.
Per-claim coverage:

- M1: line 59 (`deltaU2 - (D0 z2 - D2 z0)/D0^2`)
- M2: line 62 (`deltaP0 - (D0 n0 + N0 z0)/D0^2`)
- M3: lines 68-71 (one-pole linear coefficient against closed form)
- M4: lines 85-88 (Solve-based K elimination vs closed-form
  `(N0+eps n0)/Ptarget - 3(S+eps z2)^2/(T+eps z4)`)
- M5: lines 91-94 + line 95 (z0 independence)
- M6: lines 105-108 (transported-target K_norm equals `B0+Z0slot+eps z0+D0target`)
- M7: lines 112-115 + line 116 (z0 independence)
- M8: lines 129-133 (lane ratios via `ThreeJSymbol`; same-sign selection rule
  checked for m=1,2)
- M9: lines 141-142 (`xbar = x0`, `bx = 3 ax`)
- M10: line 146 (static Xi1 slope `na/N0 + za/D0`)
- M11: lines 150-153 (u2 lane slope `(D0 z2a - D2 za)/D0^2`)

**Assessment:**
- **Independence:** The directive explicitly demanded `Series`/`Coefficient`
  (not the SymPy `(expr_lin - expr_base)/eps` trick) and a direct
  `ThreeJSymbol` Gaunt construction (not `sympy.physics.wigner.gaunt`).
  Codex's `firstShift` uses Mathematica-native `Series` + `Normal` +
  `Coefficient[..., var, 1]` to pick off the eps^1 coefficient — a genuinely
  different engine path. (The `expr - (expr/.var->0)` subtraction inside is
  merely a baseline shift; the linear coefficient is still produced by
  `Coefficient`, not by `/eps + simplify` as in the SymPy script.) M8 uses
  `ThreeJSymbol` directly to build Gaunt and verifies the same-sign
  selection rule independently. M10/M11 use direct `D[..., ea]`
  differentiation rather than `sp.series().removeO()` then `sp.diff`. M4
  uses `Solve` and compares against the closed form `compatFixedClosed`
  written symbolically, not by mimicking SymPy intermediate variables.
  The Mathematica derivation is not a line-by-line port.
- **Assumptions:** `$Assumptions` declares the relevant symbols `Real` with
  `D0, N0, T, Ptarget, D0target, lam` all `!= 0`. This matches the SymPy
  `nonzero=True` declarations and does not over-assume positivity.
- **Non-tautology:** Each `check` compares an independently computed lhs
  (from `Series`/`Solve`/`D`/`ThreeJSymbol`) against a closed-form rhs
  written from the manifest, with `FullSimplify` collapsing the residual.
  No claim is checked by `lhs - lhs`. M5/M7 z0-independence checks are
  `D[compatFixedShift, z0]` and `D[compatTransportShift, z0]`, which would
  be nonzero if the K-elimination cancellation had failed.
- **Coverage:** All 11 manifest items M1-M11 are present. Plus two
  additional selection-rule checks for M8 (m=1, m=2 same-sign Gaunt = 0)
  and two `D[..., z0]` independence checks (M5, M7) that strengthen
  M4-M7. No regression checks were removed.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script was modified this iteration; the
existing SymPy audit was untouched and `redteam/exec_logs/stage_011_sympy.log`
is not present. This matches the directive (only a new `.wl` was required).

**Mathematica:** exit=0. The log
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_011_mathematica.log`
shows 18 `... residual = 0` lines followed by `STATUS: PASS` and
`# exit_code: 0`. Notable lines:

```
M1 delta u2 residual = 0
M4 fixed-target K-eliminated surface residual = 0
M8 same-sign selection m=1 residual = 0
M11 u2 projected-Maxwell lane slope residual = 0
STATUS: PASS
```

**Output freshness:** `.wl` mtime 2026-05-21 11:34:57; saved transcript
mtime 2026-05-21 11:51:41. Transcript is newer than the script, so the
output was regenerated post-creation.

The captured diff patch
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_011_diff.patch`
is empty (0 bytes). This is expected because the entire change is a new
untracked file pair — no in-place edit to any existing tracked file
appears in `git diff` for tracked content. Not a regression.

## Material-change assessment

`material_change`: false.

The new `.wl` only verifies (independently) the same 11 claims the SymPy
script already asserts. No derived quantity, definition, or downstream
result was altered. Downstream units should not need re-auditing on
account of this change.

## Side observations (non-blocking)

- File placement: the directive named the Mathematica path under
  `scripts/`, but the convention in this repo (all other stage NNN audits
  in `mathematica/`) places `.wl` files under `mathematica/`. The
  orchestrator's relocation aligns with the repo norm; the original
  directive path should probably be updated in the auditor template for
  future stages so Codex doesn't have to be moved.
- The `firstShift` helper subtracts the baseline `(expr /. var -> 0)`
  inside `Series` before pulling `Coefficient[..., var, 1]`. The
  baseline subtraction is redundant (since `Coefficient[..., var, 1]`
  already discards the var^0 term) but harmless. Not a blocker.
- The `Quiet[..., ClebschGordan::phy]` wrapper around `ThreeJSymbol`
  in `g[m1,m2,m3]` is appropriate — Mathematica's `ThreeJSymbol`
  occasionally emits a `phy` message on integer arguments; suppressing
  it does not change the numerical result.

## Verdict justification

Codex correctly created the missing Mathematica companion with 11
independent checks (plus four supporting consistency checks), derived the
linear-in-eps coefficients using Mathematica-native `Series`+`Coefficient`
rather than porting the SymPy `(expr_lin - expr_base)/eps` choreography,
and used `ThreeJSymbol` directly for the Y20 lane factors instead of
`sympy.physics.wigner.gaunt`. The exec log shows all 18 residuals equal to
zero, `STATUS: PASS`, and `exit_code: 0`. The transcript is fresher than
the script, the file lives where every other Mathematica audit in this
repo lives, and the diff is clean. Single finding resolved; verdict
`verified`.
