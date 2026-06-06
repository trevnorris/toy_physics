---
unit_id: 118
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 118

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Exactly two lines, matching the captured `stage_118_diff.patch` and the working-tree
diff:

- `scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:82`
  - Before: `K_q = sp.simplify((Zq/mu0) * (sp.pi**2*c_s**2/(4*L_W**2)))`
  - After:  `K_q = sp.simplify((Zq/mu0) * c_s**2 * chi_grad)`
- `mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:90`
  - Before: `kQ = FullSimplify[(zQ/mu0)*(Pi^2*cSound^2/(4*lW^2)), Assumptions -> $Assumptions];`
  - After:  `kQ = FullSimplify[(zQ/mu0)*cSound^2*chiGrad, Assumptions -> $Assumptions];`

The closed-form assertions are UNCHANGED:
- SymPy L97: `expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))`
- Mathematica L105: `expectZero["K_q closed form", kQ - (zQ/mu0)*Pi^2*cSound^2/(4*lW^2)];`

**Assessment:**
The de-taut is correct and complete in both engines. `chi_grad` is the independently
computed gradient integral `∫(χ')²dz` evaluated at py L42 (`sp.integrate(sp.diff(chi,z)**2,
(z,0,L_W))`) and asserted equal to `π²/(4L_W²)` at the "D/N stiffness check" (py L50). Its
Mathematica counterpart `chiGrad` is computed at wl L45 (`Integrate[D[chi,zTube]^2, ...]`)
and asserted at wl L54. With `K_q` now built as `(Zq/mu0)·c_s²·chi_grad`, the unchanged
L97/L105 assertion `K_q − (Zq/mu0)π²c_s²/(4L_W²)` is no longer `X − X`: it now reduces to
`(Zq/mu0)c_s²·(chi_grad − π²/(4L_W²))`, i.e. it is load-bearing and would FAIL if the
gradient integral did not reduce to `π²/(4L_W²)`. Both `chi_grad` and `chiGrad` are computed
upstream of the V-section symbol scope (the `chi`/`zTube`/`lW` symbols persist), so the
reference resolves correctly — confirmed by the exec logs producing a finite, correct `K_q`
and a `0`/PASS residual.

The new assertion is non-tautological: it depends on a genuinely derived integral. No
collateral edits — the surrounding `g_q`, `J_s`, `g_s`, `I_q`, `λ` lines (L83-104 / L91-111)
are untouched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `∫ (chi')^2 dz = pi**2/(4*L_W**2)` (the load-bearing integral feeding `K_q`)
- `D/N stiffness check = 0` (asserts the integral reduced correctly)
- `K_q = pi**2*Zq*c_s**2/(4*L_W**2*mu0)` (printed value — UNCHANGED)
- `K_q closed form = 0`
- `All Stage 118 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `Integral (chi')^2 dz = Pi^2/(4*lW^2)`
- `PASS: D/N stiffness check`
- `K_q = (cSound^2*Pi^2*zQ)/(4*lW^2*mu0)` (printed value — UNCHANGED)
- `PASS: K_q closed form`
- `Stage 118 Mathematica audit passed.`

**Output freshness:** confirmed fresh. Both committed `.txt` outputs have mtime
2026-06-06 11:26:33, newer than both scripts (mtime 2026-06-06 11:08:00). The outputs were
re-generated post-fix.

## Material-change assessment

`material_change`: false.

This is a value-preserving de-taut. The printed `K_q` value is identical pre- and post-fix
(`pi**2*Zq*c_s**2/(4*L_W**2*mu0)` in SymPy, `(cSound^2*Pi^2*zQ)/(4*lW^2*mu0)` in Mathematica),
and no other emitted value changed. Only the internal construction of `K_q` (and hence the
provenance/load-bearing nature of the L97/L105 check) changed — `chi_grad` reduces to exactly
the former literal `π²/(4L_W²)`, so all downstream consumers see the same value. No downstream
unit depends on a moved result.

## Side observations (non-blocking)

- The working tree contains modifications to other files (stage117 .wl, stage116 outputs,
  paper/*.tex, MANIFEST.yaml) unrelated to this fix; these belong to other stages/orchestration
  and are outside the stage-118 scope. The stage-118 working-tree diff is exactly the two
  intended lines, byte-identical to the captured `stage_118_diff.patch`.

## Verdict justification

The sole low-severity finding (F1, the `K_q closed form` X−X tautology) is fully resolved in
both engines with exactly the directive's two-line change and no deviation. The closed-form
assertion is now load-bearing on the independently verified gradient integral `chi_grad`/
`chiGrad`, while the printed `K_q` value and every downstream value are unchanged. Both engines
exit 0 with all checks PASS, and the committed outputs are fresh. No regressions, no collateral
edits. Verdict: verified, material_change false.
