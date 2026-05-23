---
unit_id: 054
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T18:05:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 054

## Per-finding outcomes

### F1 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:40` was rewritten from the hardcoded form `bExpr = FullSimplify[a Tan[k ell], Assumptions -> $Assumptions];` to

```
bExpr = FullSimplify[b /. First@Solve[(D[psi, s] /. s -> ell) == 0, b], Assumptions -> $Assumptions];
```

matching the directive's required change literally. (The original directive named line 34; after preserving the existing blank/`Clear`/`$Assumptions` lines the edit landed at line 40 of the current file — the substantive content is identical to what was specified.)

**Assessment:**
The Mathematica engine now genuinely derives `b` from the bottom Neumann condition `D[psi, s] /. s -> ell == 0` via `Solve`, mirroring the SymPy script's `sp.solve(sp.Eq(sp.diff(psi, s).subs(s, L), 0), B)[0]` at `scripts/.../stage054_..._sympy_audit.py:33`. The downstream assertion `expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]];` (now line 46) is no longer tautological — its passing depends on `Solve` returning the correct `b = a Tan[k ell]`. The exec output line 5 prints `B from Neumann bottom = ConditionalExpression[a*Tan[ell*k], ...]`, where the `ConditionalExpression` wrapper is the diagnostic fingerprint of an actual `Solve` invocation (Mathematica's `Solve` attaches such wrappers when auxiliary inequalities are involved). The `expectZero` helper's `ConditionalExpression[e_, _] :> e` strip (orchestrator's generic idiom fix) correctly collapses this to zero on the declared domain, so the assertion still passes. No collateral edits.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:84` was rewritten from `xFloor = FullSimplify[4 - 4/zetaReq, Assumptions -> zetaReq > 0];` to

```
xFloor = FullSimplify[x /. First@Solve[aKMax == zetaReq, x], Assumptions -> zetaReq > 0];
```

matching the directive's required change literally. (Directive cited line 78 of the original file; current file has it at line 84 after the F1 edit preserved earlier structure — substantive content identical.)

**Assessment:**
`xFloor` is now obtained by inverting `aKMax == zetaReq` for `x` via `Solve`, mirroring the SymPy `sp.solve(sp.Eq(AK_max, zeta_req), x)[0]` at `scripts/.../stage054_..._sympy_audit.py:83`. The assertion `expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)];` (now line 89) is no longer the trivial `(4 - 4/zetaReq) - (4 - 4/zetaReq) == 0`; it tests that `Solve` returns the expected root. Exec output line 29 prints `x floor at saturation = ConditionalExpression[4 - 4/zetaReq, ... zetaReq > 1]`, again carrying the `ConditionalExpression` fingerprint of a genuine `Solve` invocation (with the auxiliary `zetaReq > 1` constraint required for the inversion to lie in the valid domain). After the helper's wrapper strip, both `x floor = 4 - 4/zeta_req` (line 30) and `A_K,max(x_floor) - zeta_req` (line 32) print PASS. No collateral edits.

## Exec log assessment

**SymPy:** exit=n/a (no SymPy edit required by directive; output regenerated only to confirm freshness). Notable lines from the saved `.txt`:

```
B from Neumann bottom = A*tan(L*k)
Robin equation -> k tan(kL) - h = 0
x floor = 4 - 4/zeta_req = 0
```

SymPy output is unchanged in form and shows the same algebraic results as before.

**Mathematica:** exit=0. Notable lines:

```
B from Neumann bottom = ConditionalExpression[a*Tan[ell*k], (Element[C[1], Integers] && ...) || 2*ell*k < Pi]
Robin equation -> k tan(kL) - h = 0
PASS: Robin equation -> k tan(kL) - h
x floor at saturation = ConditionalExpression[4 - 4/zetaReq, ... zetaReq > 1]
PASS: x floor = 4 - 4/zeta_req
PASS: A_K,max(x_floor) - zeta_req
Stage 054 Mathematica audit passed.
```

All 9 Mathematica assertions PASS (Robin equation, dimensionless form, K_W identity, A_K x-form, DN limit, soft-mouth limit, x floor, A_K,max(x_floor)). The `ConditionalExpression` wrappers on the two `Solve`-derived expressions are the expected post-Solve form; they collapse to bare zero in the `expectZero` helper as designed.

**Output freshness:** confirmed. Script mtime `2026-05-22 17:38`; Mathematica output `.txt` mtime `2026-05-22 17:38` (same minute, post-edit); SymPy output `.txt` mtime `2026-05-22 17:38`. Both outputs were regenerated after Codex's edit (fix-batch log line 74 records `Stage 054: refresh mathematica -> ...` at 17:38:49).

## Material-change assessment

`material_change`: false.

Both edits are local algebraic re-derivations of values that previously were hardcoded to the same closed forms. The output values (`B from Neumann bottom = a*Tan[ell*k]`, `x floor at saturation = 4 - 4/zetaReq`) are unchanged after `ConditionalExpression` unwrapping. No expression that downstream units consume has changed numerically or symbolically — only the provenance inside this unit shifted from hand-stated to `Solve`-derived. No downstream re-audit needed on account of unit 054.

## Side observations (non-blocking)

- The orchestrator note records a preemptive batch-wide `ConditionalExpression[e_, _] :> e` strip patch added to the `expectZero` helper (visible at lines 20-30 of the `.wl`). For unit 054 specifically, this strip is exactly what the two new `Solve` calls require to produce a clean zero residual (Solve wraps both `bExpr` and `xFloor` in `ConditionalExpression` because of the auxiliary inequality conditions on `k, ell, zetaReq`). The strip is generic and does not weaken the assertions — it only removes a domain-conditional wrapper that the `$Assumptions` block already guarantees.
- Output line 17 `Limit::alimv: Warning: Assumptions that involve the limit variable are ignored.` is unchanged from pre-fix output and is unrelated to the two findings. Side observation only.
- The unit's banner text still says `STAGE 037 — ROBIN-COMPLIANCE SOFTENING` (mathematica line 32) / `STAGE 37 — ROBIN-COMPLIANCE SOFTENING` (sympy line 25). Cosmetic; out of scope for this verification.

## Verdict justification

Both findings F1 and F2 are fully resolved. The two hardcoded literals (`bExpr` for the Neumann coefficient and `xFloor` for the saturation root) have been replaced with genuine `Solve` invocations matching the SymPy counterpart, exactly as the directive specified. The Mathematica script exits 0, all 9 assertions PASS, the printed intermediate values are unchanged after `ConditionalExpression` unwrapping, and output mtimes confirm regeneration post-edit. No regressions, no collateral edits beyond the orchestrator's pre-applied generic helper patch. Verdict: verified.
