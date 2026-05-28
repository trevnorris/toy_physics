---
unit_id: 133
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 133

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The `.wl` file at
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage133_coupled_mouth_fixedpoint_mathematica_audit.wl:28-57`
no longer constructs `u` from hand-supplied `cCoeff`/`aCoeff`. Lines 35-41 now
hold

```
uSol = DSolveValue[
  {-uFun''[x] + kappa^2*uFun[x] == gSrc*sigma, uFun[0] == 0, uFun'[1] == 0},
  uFun[x],
  x
];
u = FullSimplify[uSol, Assumptions -> $Assumptions];
```

`grep -n "cCoeff\|aCoeff"` against the post-edit `.wl` returns no matches, and
`DSolveValue` is present at line 36 with both BCs supplied. The downstream
`residual` / `bc0` / `bc1` / `sKernel` / `sTarget` definitions and the four
`expectZero` calls (`ODE residual`, `u(0)`, `u'(1)`, `mouth derivative kernel`)
at lines 43-57 are unchanged from the directive's "After" block. The
static-shell limit block (lines 59-63), the general two-channel print
(lines 65-70), and `Exit[0]` (line 75) are untouched.

The orchestrator-direct cluster-A mass-fix (banner relabel `STAGE 116` →
`STAGE 133`) is also in place: `.wl` lines 26 and 65 read `STAGE 133`, `.py`
lines 31 and 64 read `STAGE 133`, notes H1 line 2 reads
`Stage 133: Full Coupled Mouth-Layer Fixed-Point Law`. No stray `STAGE 116`
or `Stage 235` references remain in any of the three files.

**Assessment:**
The edit matches the directive precisely. `DSolveValue` now derives `u(x)`
from the PDE and the D/N boundary data without any hand-imported
intermediates, so the Mathematica audit becomes a genuine independent
re-derivation that cross-checks the SymPy hand ansatz. The four `expectZero`
assertions are non-tautological because the source of `u` is now Mathematica's
own DSolve answer: if `DSolveValue` returns anything inconsistent with the
ODE, the BCs, or the paper's boxed `S(Pi, kappa)` kernel, those assertions
would fail. The exec log confirms all four still report `0` and `PASS`. No
collateral edits beyond the directive's "After" block and the cluster-A
banner relabel were introduced.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- `ODE residual = 0`
- `u(0) = 0`
- `u'(1) = 0`
- `mouth derivative kernel = 0`
- `S(Pi,0) = 1`

(SymPy script was untouched by this directive; output above is from the
post-fix rerun.)

**Mathematica:** exit=0. Notable lines:

- `u(x) = (gSrc*piM*(E^(piM + kappa*x)*kappa - ... ))/(...)`  (full DSolve
  closed form printed)
- `ODE residual = 0` then `PASS: ODE residual`
- `u(0) = 0` then `PASS: u(0)`
- `u'(1) = 0` then `PASS: u'(1)`
- `mouth derivative kernel = 0` then `PASS: mouth derivative kernel`
- `S(Pi,0) = 1`, `static-shell limit = 0`, `PASS: static-shell limit`
- `Stage 133 Mathematica audit passed.`

The printed `u(x)` is an exponential rearrangement of the hyperbolic SymPy
form, but the `mouth derivative kernel = 0` line confirms that
`FullSimplify[(D[u, x] /. x -> 0)/gSrc - sTarget] === 0`, so the two engines
agree symbolically on the kernel and the agreement now carries the
independent-derivation guarantee.

**Output freshness:** `.wl` mtime `May 27 17:22`, Mathematica `.txt` mtime
`May 27 17:47` — output is newer than script. `.py` mtime `May 27 17:20`,
SymPy `.txt` mtime `May 27 17:45` — also newer than script. Both outputs
were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The edit only swaps the construction path of `u(x)` inside the Mathematica
audit: from hand-typed ansatz to `DSolveValue` of the same ODE/BC system.
The closed-form result, the boxed kernel `S(Pi, kappa)`, the static-shell
limit, and every downstream identity printed by stage 133 are unchanged.
Banners now show `STAGE 133` instead of `STAGE 116`, which is a cosmetic
relabel with no impact on derived quantities. No downstream unit reads
stage 133's printed numerical values; downstream units consume the
algebraic kernel `S`, which is unaltered.

## Side observations (non-blocking)

- The DSolve-returned `u(x)` in the Mathematica output is written entirely
  in exponentials over the common base `E^((kappa + piM)*x)`, while the
  SymPy hand form uses `sinh`/`cosh`. They are algebraically identical, and
  the `mouth derivative kernel = 0` assertion is what binds them. Nothing
  to fix.
- The general two-channel fixed-point law remains "printed, not asserted,"
  as the original auditor noted (the assembly is trivial once `S` is
  verified per channel). The directive did not request a change here, and
  the auditor explicitly classified the partial coverage as acceptable.

## Verdict justification

The single finding F1 was addressed exactly as the directive specified:
`cCoeff` and `aCoeff` are gone, `DSolveValue` derives `u(x)` from the PDE
and both boundary conditions, and the four downstream `expectZero`
assertions still pass with the DSolve answer as their input — so the
Mathematica audit is now a genuine second-engine cross-check rather than a
transliteration. The cluster-A banner relabel landed in all three target
files. Both exec logs exit 0, outputs are fresh, and the derived kernel is
unchanged. Verdict: `verified`.
