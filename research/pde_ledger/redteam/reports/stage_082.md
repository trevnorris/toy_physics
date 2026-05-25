---
unit_id: 082
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 082 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.txt`

## What the script claims to verify

The docstring claims four results for the "STAGE 65 MASTER QUADRUPOLE RESIDUAL" (note: filename is stage082, banner is "STAGE 65"): (1) an exact inverse map between `zeta_req(Pi_tr)` and `Pi_tr = C_mix * Q(zeta)`; (2) the product thresholds `Pi_suff = C_mix Q(zeta_-)` and `Pi_fail = C_mix Q(zeta_+)` invert correctly; (3) the Family-1 strength identity `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w`; and (4) the master residual `R_quad = zeta_req - zeta_phys` vanishes at `(Pi_suff, zeta_-)` and `(Pi_fail, zeta_+)`. In practice all five assertions are different evaluations of the single identity `zeta_req(C_mix Q(z)) == z`, plus two arithmetic-only checks on the integer constants 37, 100, 1369, 136900.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 46-49 | `zeta_req(C_mix*Q(zeta)) - zeta == 0` | partial (one algebraic identity) |
| A2 | sympy | 59 | `zeta_req(Pi_suff) - zeta_- == 0` | no (A1 with `zeta -> zeta_-`) |
| A3 | sympy | 60 | `zeta_req(Pi_fail) - zeta_+ == 0` | no (A1 with `zeta -> zeta_+`) |
| A4 | sympy | 68-71 | `R_quad(Pi_suff, zeta_-) == 0` | no (definition `R_quad = zeta_req - zeta_phys` + A2) |
| A5 | sympy | 72-75 | `R_quad(Pi_fail, zeta_+) == 0` | no (definition + A3) |
| A6 | sympy | 87-90 | `100*Theta_w*37**2 - 136900*Theta_w == 0` | no (pure integer arithmetic, `100*1369 == 136900`) |
| A7 | sympy | 91-94 | `(100*Theta_w)*37**2 - 100*Theta_w*37**2 == 0` | no (substitution into a one-line algebraic definition; tautological) |
| A8 | mathematica | 39-42 | `zetaReq(cMix*qMap) - zeta == 0` | partial (mirror of A1) |
| A9 | mathematica | 50 | `zetaReq(piSuff) - zetaMinus == 0` | no (mirror of A2) |
| A10 | mathematica | 51 | `zetaReq(piFail) - zetaPlus == 0` | no (mirror of A3) |
| A11 | mathematica | 56-59 | `rQuad(piSuff, zetaMinus) == 0` | no (mirror of A4) |
| A12 | mathematica | 60-63 | `rQuad(piFail, zetaPlus) == 0` | no (mirror of A5) |
| A13 | mathematica | 71 | `100*thetaW*37^2 - 136900*thetaW == 0` | no (mirror of A6) |
| A14 | mathematica | 72-75 | `(100*thetaW)*37^2 - 100*thetaW*37^2 == 0` | no (mirror of A7) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:87-94`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:71-75`

**What's wrong:**
The two "Family-1 strength identity" assertions are pure arithmetic on hand-baked integer constants — no physics is exercised:

```
Lambda_ell = sp.Integer(37)
Xi_F1_from_Upsilon = sp.simplify(Upsilon_w * Lambda_ell**2)     # = 1369 * Upsilon_w
Xi_F1_from_Theta   = sp.simplify(100 * Theta_w * Lambda_ell**2) # = 136900 * Theta_w
expect_zero("Xi_F1(Theta_w) - 136900 Theta_w",
            Xi_F1_from_Theta - sp.Integer(136900) * Theta_w)     # 136900*Theta_w - 136900*Theta_w
expect_zero("Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w)",
            Xi_F1_from_Upsilon.subs(Upsilon_w, 100*Theta_w) - Xi_F1_from_Theta)
```

The first check is `100 * Theta_w * 37**2 - 136900 * Theta_w`, which is just SymPy confirming the integer arithmetic `100 * 1369 == 136900` — it would fail only if SymPy's integer multiplication were broken. The second check substitutes `Upsilon_w -> 100 Theta_w` into `Upsilon_w * 1369` and compares to `100 * Theta_w * 1369`, which is identically `0` by associativity. The Mathematica script at lines 65-75 contains the mirror of these checks with the same defect.

The docstring claim 3 ("Exact Family-1 strength identity `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w`") embeds three integer constants (37, 100, 1369) plus the ratio `Upsilon_w = 100 Theta_w`. None of these are derived or anchored to an upstream symbolic quantity in this script — they are written into the definitions then trivially re-confirmed.

**Why this matters:**
The audit reports "Xi_F1 identity verified" but the assertions cannot fail for any physics-related reason. If the underlying physical relationship that motivates `Lambda_ell = 37` or `Upsilon_w = 100 Theta_w` were wrong, these checks would still pass.

**Required change:**
Either (a) replace the two checks with a non-tautological derivation that builds `Xi_F1` from a physical input chain (e.g., an upstream construction that produces `Lambda_ell` and the ratio `Upsilon_w/Theta_w` independently, so the equality `Upsilon_w * Lambda_ell**2 == 100 * Theta_w * Lambda_ell**2` becomes a substantive consequence rather than a definitional restatement), or (b) demote the two assertions to plain `print` statements and document that this script does not verify the Family-1 strength constants — relegating that verification to whichever upstream stage produces 37 and the 100-ratio. Since this auditor was instructed not to expand scope or invent physical derivations, option (b) is the safe mechanical correction. Apply in both `.py` (lines 87-94) and `.wl` (lines 71-75).

**Verification:**
After the change, the SymPy script's output should no longer contain "Xi_F1(Theta_w) - 136900 Theta_w = 0" as an asserted check; either it is gone or it carries an inline disclaimer. The Mathematica script's `PASS:` lines for the same two checks should disappear or be replaced with `Print[...]` only.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:46-75`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:39-63`

**What's wrong:**
Five of the seven assertions in the SymPy script (A1-A5) test a single algebraic identity. The "inverse map" check `zeta_req(C_mix * Q(zeta)) - zeta == 0` (lines 46-49) is the parent identity. The "Pi_suff" check (line 59) is the same identity with `zeta -> zeta_-`. The "Pi_fail" check (line 60) is the same identity with `zeta -> zeta_+`. The two `R_quad` checks (lines 68-75) substitute `(Pi_tr, zeta_phys) = (Pi_suff, zeta_-)` and `(Pi_fail, zeta_+)` into `R_quad = zeta_req - zeta_phys`, which by `R_quad`'s own definition (line 64) immediately reduces to `zeta_req(Pi_suff) - zeta_-` and `zeta_req(Pi_fail) - zeta_+` — i.e., A2 and A3.

The docstring lists four distinct claims (inverse map, product thresholds, Family-1 strength, master residual). Claims 1, 2, and 4 collapse to a single algebraic identity once you trace the substitutions. The "product thresholds" claim has no separate content beyond renaming `zeta -> zeta_+/-`; the "master residual" claim has no separate content beyond defining `R_quad := zeta_req - zeta_phys` and evaluating it where `zeta_req = zeta_phys`.

The script never tests anything inequality-flavored about `Pi_suff` and `Pi_fail` (e.g., that one is strictly above and the other strictly below a threshold), nor does it verify that `R_quad` has the correct sign on either side of its root. The "Theorem ledger" printed at the end claims directional theorems ("guaranteed success", "guaranteed failure") that are not exercised by any assertion — those would require monotonicity or sign checks of `R_quad` in `zeta_phys` and `Pi_tr`, which are absent.

**Why this matters:**
The script publishes five PASS lines, suggesting five independent algebraic facts, but they are five evaluations of one identity. If the inverse-map identity were to fail, all five would fail together — there is no redundancy or coverage gain. The "guaranteed success / guaranteed failure" theorems in the ledger text are unverified by any assertion.

**Required change:**
Add one substantive assertion that exercises the directional content claimed in the theorem ledger. Concretely, in the SymPy script, after line 75 (and the corresponding location in the Mathematica script after line 63), add a sign check for `R_quad` on each side of its root, using the existing symbols and explicit assumptions that match the script's existing positivity declarations.

Specifically, add to the SymPy script after line 75:

```python
# Directional sign of R_quad away from its root, holding the other slot fixed.
# At Pi_tr = Pi_suff, R_quad is increasing in zeta_phys (since dR_quad/dzeta_phys = -1).
dR_dzeta_phys = sp.simplify(sp.diff(R_quad, zeta_phys))
expect_zero("dR_quad/dzeta_phys + 1", dR_dzeta_phys + 1)
# At zeta_phys = zeta_-, R_quad is the difference (zeta_req(Pi_tr) - zeta_-),
# so dR_quad/dPi_tr is the same as dzeta_req/dPi_tr; verify it equals the
# closed-form derivative of zeta_req w.r.t. Pi_tr.
dzeta_req_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
dR_dPi = sp.simplify(sp.diff(R_quad.subs(zeta_phys, zeta_minus), Pi_tr))
expect_zero("dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)",
            dR_dPi - dzeta_req_dPi)
```

(Mirror in `.wl` with `D[expr, var]` and `FullSimplify`.) This adds real, non-tautological content: it verifies the partial-derivative structure that underwrites the directional theorems printed in the ledger.

**Verification:**
The SymPy output should contain two new lines `dR_quad/dzeta_phys + 1 = 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-) = 0`. The Mathematica output should contain the corresponding new `PASS:` lines. Both scripts should still exit 0.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:26-81` (entire body)

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script, not an independent re-derivation. Corresponding sections:

- SymPy 38: `zeta_req = sp.simplify((Pi_tr - Cmix) / (Cmix - eps_blk * (2 * Cmix - Pi_tr)))`
- Mathematica 33: `zetaReq = FullSimplify[(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr)), Assumptions -> $Assumptions];`

- SymPy 39: `Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))`
- Mathematica 34: `qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];`

- SymPy 79-81: `Lambda_ell = sp.Integer(37); Xi_F1_from_Upsilon = sp.simplify(Upsilon_w * Lambda_ell**2); Xi_F1_from_Theta = sp.simplify(100 * Theta_w * Lambda_ell**2)`
- Mathematica 65-67: `lambdaEll = 37; xiF1FromUpsilon = FullSimplify[upsilonW*lambdaEll^2, ...]; xiF1FromTheta = FullSimplify[100*thetaW*lambdaEll^2, ...]`

- SymPy 97-99 (theorem ledger prints): identical English text appears verbatim at Mathematica 79-81.

Every assertion in `.wl` is the syntactic Mathematica rendering of the corresponding `.py` assertion at the same algebraic step. There is no alternate derivation path (e.g., the Mathematica script does not start from a different parametrization of `Q`, or build `zeta_req` by solving an equation in Mathematica's `Solve` and comparing). The intent of the second-engine policy — that an independent symbolic engine corroborate the result via a different algebraic route — is not met.

**Why this matters:**
A transliteration cannot catch implementation errors in the shared algebra; it can only catch engine-specific arithmetic bugs. If the SymPy author made an error in writing `(1 - eps_blk * zeta)` instead of, say, `(1 + eps_blk * zeta)`, the Mathematica script — having copied the same expression — would not reveal it.

**Required change:**
In the Mathematica script, replace the hand-supplied closed form for `zetaReq` (line 33) with a Solve-based derivation that re-discovers it from the Q-map. Specifically, replace lines 33-34 with:

```mathematica
qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
(* Independently solve PiTr = cMix*qMap for zeta, expressed as zetaReq(PiTr). *)
zetaReqDerived = zeta /. First[Solve[PiTr == cMix*(qMap /. zeta -> zetaSym), zetaSym]];
zetaReqDerived = FullSimplify[zetaReqDerived /. zetaSym -> zeta, Assumptions -> $Assumptions];
zetaReq = zetaReqDerived;
```

(Adjust the variable name `zetaSym` to whatever is free in `$Assumptions`.) This forces the Mathematica engine to derive `zetaReq` from `qMap` via `Solve`, rather than restating the SymPy expression in Mathematica syntax. After this change, the existing `expectZero` assertions will be testing whether SymPy's hand-supplied `zetaReq` agrees with Mathematica's `Solve`-derived inverse of `qMap` — that is a genuine cross-engine check.

**Verification:**
The Mathematica script's printed `zeta_req` line should still simplify to the same closed form, but the source of that expression in the script is now `Solve[...]` rather than a hand-coded fraction. The verifier inspects line 33's right-hand side and confirms it contains `Solve` (or `Reduce`) and not the literal `(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr))`.

### F4 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:65`

**What's wrong:**
The integer `Lambda_ell = 37` is dropped into the script as a literal, with no comment citing the upstream unit that derives it and no in-script construction. Likewise the factor `100` in `Xi_F1_from_Theta = sp.simplify(100 * Theta_w * Lambda_ell**2)` (line 81) is dropped in literal with no provenance. These numbers then feed the Family-1 assertions A6/A7, which become arithmetic identities on the chosen literals (see F1).

**Why this matters:**
A reader of the script alone cannot confirm whether `37` and `100` are the correct upstream constants for this unit; if they were copied wrong from an upstream stage, the assertions would still PASS because the same wrong numbers appear on both sides of each `expect_zero`. The captured output (`Xi_F1 = 1369*Upsilon_w`, `Xi_F1 = 136900*Theta_w`) is then propagated as "verified" without provenance.

**Required change:**
Add an inline comment immediately above line 79 of the SymPy script (and above line 65 of the Mathematica script) citing the upstream stage and/or paper section that establishes `Lambda_ell = 37` and the ratio `Upsilon_w = 100 * Theta_w`. Format:

```python
# Lambda_ell = 37 (carry-forward from stage NNN; see docstring of that script)
# Ratio Upsilon_w = 100 * Theta_w is the Family-1 weight convention from stage MMM.
Lambda_ell = sp.Integer(37)
```

If the upstream stage IDs are not known to Codex, leave a `TODO(provenance):` marker rather than guessing. Do not modify the numerical values.

**Verification:**
The SymPy script line immediately preceding `Lambda_ell = sp.Integer(37)` contains a `#`-prefixed comment naming an upstream stage (or `TODO(provenance):`). The Mathematica script line preceding `lambdaEll = 37` contains an analogous `(* ... *)` comment.

## Independent-derivation check (Mathematica)

See F3. The `.wl` is a syntactic transliteration of the `.py`. Both define `zetaReq` and `qMap` as identical closed-form fractions, both perform the same five substitutions in the same order, both hard-code `Lambda_ell = 37` and `100`, both print the same English "Theorem ledger" text at the end. No alternate derivation route (e.g., `Solve`, `Reduce`, series inversion, or formulation as a fixed-point equation) is used.

## Engine cross-check

Both engines agree on the printed forms:

- SymPy: `zeta_req(Pi_tr,C_mix,eps_blk) = (-C_mix + Pi_tr)/(C_mix - eps_blk*(2*C_mix - Pi_tr))`
- Mathematica: `zeta_req(Pi_tr,C_mix,eps_blk) = -((cMix - PiTr)/(cMix - 2*cMix*epsBlk + epsBlk*PiTr))`

These are algebraically identical (multiply numerator and denominator by `-1` and expand the `-eps_blk*(2*C_mix - Pi_tr)` term).

- SymPy: `Q = (zeta*(2*eps_blk - 1) - 1)/(eps_blk*zeta - 1)`
- Mathematica: `qMap = (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)`

Algebraically identical (multiply both num and denom by `-1`).

- Both produce `Xi_F1 = 1369*Upsilon_w` and `Xi_F1 = 136900*Theta_w`.

All five "physics" assertions and the two arithmetic assertions PASS in both engines. So `engines_agree = true`, but per F3 this agreement is weakened by the fact that the agreement is over identical inputs, not independent derivations.

## Verdict justification

Findings, not clean. The script does verify a single algebraic identity (the inverse-map relation between `zeta_req` and `Q`), and that identity does hold under attack — I traced the cancellation by hand and confirmed `zeta_req(C_mix*Q(z)) = z` reduces to `z*C_mix*(1-eps_blk)/(1-eps_blk*z) / [C_mix*(1-eps_blk)/(1-eps_blk*z)] = z`. The Mathematica engine concurs. What does not hold up: (a) the five "physics" PASSes are five instances of that one identity, so the verification coverage is overstated (F2); (b) the two "Family-1 strength" PASSes are pure integer arithmetic on hand-supplied literals 37, 100, 1369, 136900 (F1, F4); (c) the Mathematica script is a transliteration rather than an independent derivation (F3). None of these warrant `UNFIXABLE` or `CRITICAL_DOWNSTREAM` — the script's results are not wrong, they are under-verified. Codex can apply the F1, F2, F3, F4 corrections mechanically.

## Self-test notes

I walked through the proposed F2 derivative checks mentally: `R_quad = zeta_req(Pi_tr) - zeta_phys` depends on `Pi_tr, C_mix, eps_blk, zeta_phys`, so `sp.diff(R_quad, zeta_phys)` is well-defined and equals `-1` identically (variable independence: zeta_phys appears in R_quad). Substituting `zeta_phys -> zeta_minus` then differentiating w.r.t. `Pi_tr`: `zeta_minus` is independent of `Pi_tr`, so `dR/dPi_tr = dzeta_req/dPi_tr` (variable independence holds). The trivial-case check for the new assertions: at `eps_blk = 0`, `zeta_req = (Pi_tr - C_mix)/C_mix`, so `dzeta_req/dPi_tr = 1/C_mix`, which is nonzero — confirms the assert isn't trivially passing on zero. For F3 the `Solve`-based path: `PiTr == cMix*qMap` is linear in `zeta` after clearing the denominator, so `Solve` returns a single branch (no `ConditionalExpression` shenanigans expected, but the directive should mention `expectZero`-style stripping just in case). For F4 the change is comment-only and cannot break the assertions. Path check: SymPy script in `scripts/`, Mathematica script in `mathematica/`, both confirmed by the file listing.
