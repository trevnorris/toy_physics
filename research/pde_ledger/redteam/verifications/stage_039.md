---
unit_id: 039
batch: III.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T12:35:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 039

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:114-117` — replaced
  `expect_zero("z1/z0 - (kappa1/kappa0) R_U", sp.simplify(z1/z0 - (kappa1/kappa0)*R_U))` with
  `expect_zero("z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))", sp.simplify(z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))))`.
  `R_U` is no longer referenced in the assertion; the `print("R_U =", R_U)` informational line at 113 is retained.
- `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:95-98` — same replacement
  applied to the Mathematica mirror; the `Print["R_U = ", fmt[rU]]` informational line at 94 is retained.

**Assessment:**

The edits match the directive verbatim in both files. The new check references `z0` and `z1` against
their explicit `(1+rho0)` / `(1+rho0/(1+deltaU))` rho-structure rather than against the named symbol
`R_U`, satisfying the auditor's stated verification criterion (the auditor's own self-test notes at
report line 205 affirm this form is acceptable). The transcripts show
`z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU)) = 0` printed and PASS in both engines (sympy
output line 32; mma output lines 37-38). No collateral edits.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:129-139` — inserted flat-U
  baseline definitions `M_mix_flat`, `R_target_flat` immediately before `M_mix_split`,
  `R_target_split`; replaced `expect_zero("product law", product - 8*Lambda*(1-eps_W_split)/pi**2)`
  with two substitution checks
  `expect_zero("M_mix split is M_mix_flat under eps_W -> eps_W_split", M_mix_split - M_mix_flat.subs(eps_W, eps_W_split))`
  and the analogous `R_target` check. `print("product =", product)` at line 137 retained.
- `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:109-119` — inserted
  `mMixFlat`, `rTargetFlat` definitions immediately after the section banner; replaced
  `expectZero["product law", product - 8 lambda (1 - epsWSplit)/Pi^2]` with the two analogous
  substitution checks. `Print["product = ", fmt[product]]` at line 117 retained.

**Assessment:**

Edits match the directive verbatim in both files. The new checks are structural identities that
will fail under exponent perturbations in either `M_mix_split` or `R_target_split` (per the
auditor's stated verification criterion), so they are non-tautological in the sense the auditor
specified. The transcripts show both new checks pass with residual `0` in both engines (sympy
output lines 43-44; mma output lines 49-52). The old `product law` assertion is gone. No collateral
edits.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**

- `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:59-62` — `deltaSplit`
  is no longer a postulated closed form; instead the script introduces
  `a1Direct = FullSimplify[a1 /. cEtaU^2 -> epsEta kU kEtaEff]`, then
  `deltaSplitDerived = FullSimplify[a1Direct/a0Expected - 1]`,
  `deltaSplitPostulated = (delta0 + epsEta deltaU/(1+deltaU))/(1-epsEta)`, and finally
  `deltaSplit = deltaSplitDerived`. A new check at line 71
  `expectZero["delta_split derived matches postulated", deltaSplitDerived - deltaSplitPostulated]`
  asserts the Mathematica-derived form equals the SymPy postulate.
- `mathematica/...stage039...wl:77-79, 83` — `epsWSplit` similarly switched from postulated to
  derived: `epsWSplitDerived = FullSimplify[epsWDirect /. cUW^2 -> epsW kU kWEff/sigma]`,
  `epsWSplitPostulated = epsW(1 - (2/11) deltaU/(1+deltaU))`, `epsWSplit = epsWSplitDerived`. The
  old `expectZero["eps_W direct - split formula", ...]` at the prior line 77 was replaced by
  `expectZero["eps_W_split derived matches postulated", epsWSplitDerived - epsWSplitPostulated]`.
- `mathematica/...stage039...wl:100-105` — `dDir` left as `FullSimplify[kappa0 z1 - kappa1 z0]`;
  added `dDirDerived = dDir`, `dDirPostulated = FullSimplify[-kappa0 kappa1 gW rho0 deltaU/(1+deltaU)]`,
  `dDirExpected = dDirDerived`; replaced
  `expectZero["direction-splitting invariant", dDir - dDirExpected]` with
  `expectZero["direction-splitting invariant derived matches postulated", dDirDerived - dDirPostulated]`.

**Assessment:**

All three restructurings match the directive verbatim. The transcript confirms the three new
`derived matches postulated` assertions appear and print residual `0` (mma output lines 20-21,
28-29, 40-41), and the original `A0 direct - expected` / `A1 direct - expected` assertions remain
and still pass (mma output lines 16-19). The Mathematica engine now computes `deltaSplit`,
`epsWSplit`, `dDir` from the direct expressions and checks them against the SymPy-side postulates,
so a sign or coefficient typo in either side's postulate would surface as a non-zero residual
(`engine_disagreement` mode unlocked).

One mild code-quality note (non-blocking): at mma:60 the assignment of `deltaSplitDerived`
references `a0Expected` before `a0Expected` is bound at mma:63. Because Mathematica's `Set`
evaluates the RHS at definition time, this leaves the bare symbol `a0Expected` inside the stored
value of `deltaSplitDerived`; the symbol is resolved later when `deltaSplitDerived` is used in
the `Print` and the `expectZero` (by which time `a0Expected` is defined). The transcript
(`delta_split derived matches postulated = 0`, PASS) confirms the resolution happens correctly,
so the check is sound. The directive specified this exact ordering, so Codex followed it
verbatim — not a deviation. Listed as a side observation below for awareness but not blocking.

## Exec log assessment

**SymPy:** exit=0 (inferred from transcript; the script raises `AssertionError` on failure, the
transcript completes with the full theorem ledger printed). Notable lines:

- `A0 direct - expected = 0` and `A1 direct - expected = 0` (lines 16-17).
- `eps_W direct - split formula = 0` (line 24).
- `z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU)) = 0` (line 32) — new F1 label.
- `direction-splitting invariant = 0` (line 34).
- `M_mix split is M_mix_flat under eps_W -> eps_W_split = 0` and
  `R_target split is R_target_flat under eps_W -> eps_W_split = 0` (lines 43-44) — new F2 labels.
- No `Traceback` or `AssertionError` anywhere.

**Mathematica:** exit=0. Notable lines:

- `PASS: A0 direct - expected` and `PASS: A1 direct - expected` (output lines 17, 19).
- `PASS: delta_split derived matches postulated` (line 21) — new F3 label.
- `PASS: eps_W_split derived matches postulated` (line 29) — new F3 label.
- `PASS: z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))` (line 38) — new F1 label.
- `PASS: direction-splitting invariant derived matches postulated` (line 41) — new F3 label.
- `PASS: M_mix split is M_mix_flat under epsW -> epsWSplit` and
  `PASS: R_target split is R_target_flat under epsW -> epsWSplit` (lines 50, 52) — new F2 labels.
- Terminates with `Stage 039 Mathematica audit passed.`. No `FAIL`, `$Failed`, or stack traces.

**Output freshness:** confirmed. Both saved `.txt` outputs (mtime 12:26) are newer than the
corresponding script files (mtime 12:25); freshly regenerated after the codex edits.

## Material-change assessment

`material_change`: false.

The edits only strengthen the verification of Stage 039's existing claims — they replace
tautological assertions with structurally-informative ones (F1, F2) and convert postulated closed
forms in the Mathematica mirror into derived ones (F3). No derived constant, sign, or symbolic
form propagated downstream changed. The closed-form expressions for `delta_split`, `eps_W_split`,
`D_dir`, `M_mix_split`, `R_target_split`, and `product` printed in the transcripts are unchanged
from the pre-fix audit. Downstream units > 039 do not need a re-audit on substantive grounds.

## Side observations (non-blocking)

1. The Mathematica forward-reference ordering at mma:60 (`a0Expected` used before bound at mma:63)
   is fragile-looking but works correctly because of how `Set` and lazy symbol resolution interact;
   the directive specified exactly this ordering, so Codex was correct to apply it as-written. A
   future cleanup might reorder the `a0Expected` definition to precede the `deltaSplitDerived`
   line, but this is purely stylistic.
2. The F1 sympy/mma new assertion (`z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU)) = 0`) is
   technically still an algebraic identity following from the literal definitions of `z0`, `z1`,
   but it no longer references the named symbol `R_U`, and the auditor explicitly accepted this
   form in the report's self-test notes (line 205) as the intended non-tautological structure.
   Listed for completeness, not as a verification failure.

## Verdict justification

All three findings were applied verbatim per the directive, both engines exit 0, all new
assertions print residual `0` and PASS with the directive-specified labels, the old tautological
assertions are gone, and the saved outputs are newer than the scripts. The F3 restructuring
successfully unlocks engine-disagreement detection at the three key derived quantities
(`deltaSplit`, `epsWSplit`, `dDir`) by deriving them in Mathematica and asserting equality with
the SymPy-side postulates. No closed-form constants or signs changed, so downstream units are not
materially affected.
