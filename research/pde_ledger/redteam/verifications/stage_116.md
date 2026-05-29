---
unit_id: 116
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T21:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 116

This unit took two Codex iterations: iteration 1 applied F1 (independent Mathematica
eigenvalue solve) and F2 (delete-and-print renormalization), but F1 regressed on re-run
(non-deterministic `Solve` over symbolic `lW` left `kWValue` unbound, exit 1). Iteration 2
applied delta-directive `stage_116_delta1.md` (`F1-fix`): solve for the product `u = kW*lW`
on `(0, Pi)` so no symbolic division is needed. I verify the FINAL post-iteration-2 state.

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`.wl` lines 34-49 (was a posited `Sin[k x]` transliteration). The final form:
- `gensol = DSolve[{q''[xv] + kSym^2 q[xv] == 0, q[0] == 0, q'[0] == 1}, q, xv]` (line 38) —
  the mode is now *derived* from the ODE+left-BC, not posited.
- `qGenExpr` is checked against the ODE (line 41) and `q(0)=0` (line 43).
- `charEq = D[qGenExpr, xv] /. xv -> lW` (line 44) forms the Neumann characteristic from the
  derived mode.
- `uRoot = u /. First[Solve[Cos[u] == 0 && 0 < u < Pi, u, Reals]]` (line 46) — eigenvalue
  obtained from a solver acting on `Cos[u]==0` over the bounded interval; `kWValue = uRoot/lW`
  (line 47). It is NOT typed as `Pi/(2 lW)`.
- Two checks: `Cos[kWValue*lW]==0` (line 48) and `kWValue - Pi/(2 lW)==0` (line 49).

**Assessment:**
Correct and matches both the main directive's F1 ("eigenvalue must be derived, not posited")
and the delta-directive's `F1-fix` ("solve on the product `u = kW*lW` over `(0, Pi)`; recover
`kWValue = u/lW`; deterministic"). The key load-bearing requirement is met: `kWValue` comes
from a `Solve` on the characteristic equation `Cos[u]==0`, not a hand-typed `Pi/(2 lW)` later
verified by substitution (which is what the SymPy script still does, and is allowed — only the
*second* engine must be made independent). The `Pi/(2 lW)` literal appears only inside the
verification target of line 49 and the comment, never as the definition of `kWValue`. The
regression from iteration 1 is fixed: the Mathematica exec log exits 0 and both eigenvalue
checks PASS (log lines 14-17). No collateral edits beyond the eigenvalue block; lines 51+
(`OmegaW`, `kappa0Derived`, tube-length) consume `kWValue` exactly as before. The SymPy script
was correctly left as posit-and-verify for the eigenvalue (per directive: "Do not change the
SymPy script for F1"), so the two engines now use structurally different eigenvalue derivations
— the point of the two-engine policy is restored.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy lines 64-86 (was 64-94): the round-trip assertion, the divide-by-own-scale
  `expect_zero("final kappa_c - 1/3"...)` / `("final gamma_c - 1/9"...)`, the `D_bare` ansatz,
  and both coefficient read-backs are DELETED. Replaced by labeled prints (lines 77-81) of
  `kappa0_bare` (now `kappa0_from_tube.subs(L_W, L_W_required)`), `gamma0_bare`, `kappa_c`,
  `gamma_c`. No `expect_zero` in this block.
- `.wl` lines 65-78 (was 62-81): mirrored — `kappa0BareGeom` round-trip assertion, the two
  `expectZero` divide-by-scale checks, the `dBare` ansatz and both `Coefficient` read-backs are
  DELETED. Replaced by labeled `Print` lines (74-78). No `expectZero` in this block.

**Assessment:**
This is exactly the agreed `f2_revised` resolution recorded in the directive front-matter
(orchestrator+Codex consult): the audit's originally-proposed `3*kappa0 == 9*gamma0` replacement
was itself tautological (both operands hardcoded literals), so the agreed fix is *deletion +
labeled prints*, because `gamma0` has no independent in-stage derivation and any
`kappa_c`/`gamma_c` assertion would be tautological. Per the verifier brief, the ABSENCE of
these assertions IS the fix and must not be flagged as `insufficient_verification`/`not_attempted`.

I independently confirm the directive's verification checklist:
1. No `kappa0_bare_geom`/`kappa0BareGeom` round-trip-from-`L_W_required` assertion remains
   (grep clean in both files).
2. No `D_bare`/`dBare` ansatz or `coeff`/`Coefficient` read-back remains (grep clean).
3. No `expect_zero`/`expectZero` on `kappa_c`/`gamma_c` and no `3*kappa0-9*gamma0` relation
   remains (grep clean — the only `expect_zero` mentions in the renorm region are inside the
   explanatory comments).
4. The genuine load-bearing checks survive and PASS in both logs: D/N eigenvalue, `kappa0`
   collapse (`kappa0_derived - 4 L_W^2/(pi^2 a^2)`), and the tube-length law
   (`L_W - pi a sqrt((1+r_c)/3)/2`). SymPy log lines 10-16 (all `= 0`), Mathematica log lines
   10-23 (all PASS).
5. Both scripts print the renormalization values and exit 0.

The retained `kappa0_bare = kappa0_from_tube.subs(L_W, L_W_required)` is now a *reported*
provenance-explicit value, not an assertion, so the earlier tautology concern no longer applies.
`gamma0_bare = (1+r_c)/9` is honestly labeled as an upstream-carried Stage-98 input, not an
in-stage derivation. The summary blocks (sympy 83-86, wl 80-84) still reference the surviving
symbols and remain intact. No regression.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `tube-length law: L_W = pi a sqrt((1+r_c)/3) / 2 = 0` (line 16) — genuine load-bearing check PASS.
- `kappa0 from D/N eigenvalue matches geometric expression = 0` (line 13) — PASS.
- `Renormalization (definitional consequence, not an independent check):` (line 17) followed by
  printed `kappa_c = 1/3`, `gamma_c = 1/9` (lines 20-21) — reported, not asserted, as designed.

**Mathematica:** exit=0. Notable lines:
- `PASS: D/N eigenvalue solves Cos[kW lW]==0` (line 15) and `PASS: D/N eigenvalue kW = Pi/(2 lW)`
  (line 17) — the F1 regression is fixed; eigenvalue now derived and both checks pass.
- `PASS: tube-length law` (line 23) and `PASS: kappa0 from D/N eigenvalue matches geometric
  expression` (line 19) — load-bearing checks intact.
- `Stage 116 Mathematica audit passed.` (line 35), exit 0.

**Output freshness:** The exec logs (`redteam/exec_logs/stage_116_sympy.log` @ 21:02:52,
`stage_116_mathematica.log` @ 20:59:30) are both NEWER than their scripts (`.py` @ 20:16:05,
`.wl` @ 20:58:12) and reflect the post-fix structure — these are the authoritative run records
this verifier is instructed to use. HOWEVER, the canonical saved outputs under
`scripts/output/...sympy_audit.txt` and `mathematica/output/...mathematica_audit.txt` are STALE
(both @ 2026-05-27, older than the scripts) and still contain the PRE-fix content
(`gamma0 extracted from D_bare`, `kappa0_bare extracted from D_bare`, `D/N trial satisfies...`).
The exec runs apparently wrote to the exec_logs but did not refresh the committed `.txt`
artifacts. This does not affect the verdict (the authoritative exec logs pass and match the
edited scripts) but the orchestrator should regenerate the saved `.txt` outputs before commit so
the repository artifact matches the verified scripts. Flagged below as a non-blocking side
observation.

## Material-change assessment

`material_change`: false.

F1 changes only *how* the Mathematica eigenvalue is obtained; the resulting value
`kW = Pi/(2 lW)` is unchanged and verified identical. F2 deletes tautological assertions and
converts them to prints; no derived value changes (`kappa0`, `gamma0`, `kappa_c=1/3`,
`gamma_c=1/9`, and the tube-length law are all identical to pre-fix). No constant, sign, or
target propagates to downstream units. The original auditor likewise found no stop-cold and no
downstream propagation.

## Side observations (non-blocking)

- Stale committed `.txt` outputs (see "Output freshness" above): the saved `scripts/output/` and
  `mathematica/output/` `.txt` files still show pre-fix content. The verified, fresh run is in
  `redteam/exec_logs/`. Recommend regenerating the saved outputs before the batch commit so the
  committed artifact matches the verified scripts. Non-blocking for verification because the
  authoritative exec logs pass.
- `.wl` line 44 computes `charEq` (the characteristic equation from the derived mode) but the
  variable is not subsequently consumed — the eigenvalue is obtained from the independent
  `Solve[Cos[u]==0...]` at line 46. This is harmless (it documents that the derived mode's
  Neumann residual is `~Cos[k lW]`) and the directive explicitly allowed Codex to adapt syntax;
  not a defect.
- The SymPy `z = sp.symbols("z", ...)` declaration (line 22) is now unused after removing the
  `D_bare` ansatz; the directive explicitly said this may remain ("do not refactor the symbol
  declarations"). Not a defect.

## Verdict justification

Both findings are `resolved`. F1: the Mathematica second engine now derives the D/N eigenvalue
from a `DSolve` general solution plus a deterministic `Solve` on the characteristic equation
`Cos[u]==0` over `(0, Pi)` (not a transliterated posited `Sin[k x]` and not a hand-typed
`Pi/(2 lW)`), the iteration-1 non-determinism regression is fixed, and both eigenvalue checks
PASS with exit 0. F2: the unfalsifiable round-trip / divide-by-own-scale / coefficient-read-back
assertions are deleted in BOTH scripts and replaced with honestly-labeled consequence prints —
the agreed `f2_revised` treatment — while the genuine load-bearing checks (eigenvalue, `kappa0`
collapse, tube-length law) remain and PASS. Both exec logs exit 0. No regressions in the diff; no
material change to any downstream-consumable value. The only blemish is stale committed `.txt`
outputs, which is a non-blocking artifact-refresh task for the orchestrator, not a verification
failure.
