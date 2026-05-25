---
unit_id: 084
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: missing
  mathematica: present
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 084 red-team report

## Files reviewed

- sympy: (missing)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt`

Mathematica script mtime: 2026-05-11 11:56 (1778522213). Output mtime: 2026-05-11 13:00 (1778526051). Output is fresher than script; not stale.

Unit 084 has `is_status_only_candidate: true`, `is_checkpoint: false`. Per the manifest carve-out, missing SymPy is acceptable provided every result the script asserts is verified by an upstream unit.

Upstream verification check (positive):
- The five numeric `zeta` constants (`zeta_-^(chi) = 2.46622291347846`, `zeta_+^(chi) = 2.46752913273870`, `zeta_-^(J) = 2.44257571477179`, `zeta_+^(J) = 2.46752736855058`, `zeta_max^(F1) = 2.46752922945601`) all originate in `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:45-49` (and are reused in stages 086, 087, 090, 098).
- The inverse demand-map identity `zeta_req(C_mix*Q(zeta)) - zeta == 0` is independently verified in `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:47` and its Mathematica twin.
- The strength identity `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` is asserted in `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:80-89` (same construction).
- `kappa_F1 = 12321/5` and `eta_F1 = 37` are defined and used in stages 079, 080, 083, 089.

Conclusion: no `missing_sympy` finding for this unit.

## What the script claims to verify

The script consolidates the full reduced moving-throat PDE write-up. It restates: (1) the explicit physical-zeta formula `zeta_phys(Pe, eta; kappa)`, the required-zeta formula `zeta_req(Pi_tr, C_mix, eps_blk)`, and the demand-map `Q(zeta; eps_blk)`, and asserts the inverse round-trip `zeta_req(C_mix*Q(zeta)) - zeta == 0`; (2) the Family-1 strength identity `Xi_F1 = lambdaEll^2 * Upsilon_w = 1369 Upsilon_w` and `Xi_F1 = 100 lambdaEll^2 * Theta_w = 136900 Theta_w` at `lambdaEll = 37`; (3) the natural-window ordering gap `zeta_+^(chi) - zeta_-^(chi) ~ 1.306e-3` and the hard-ceiling gap `zeta_max^(F1) - zeta_+^(chi) ~ 9.67e-8`, given five upstream-imported numeric `zeta` values. The script docstring is the banner string, which mislabels the unit as "STAGE 067" (cosmetic / non-mathematical; not filed as a finding per the doc-alignment exclusion).

## Assertion inventory

| #  | Script        | Line | Form                                                                                             | Anchored to claim? |
|----|---------------|------|--------------------------------------------------------------------------------------------------|--------------------|
| A1 | mathematica   | 53   | `expectZero["inverse demand map", (zetaReq /. piTr -> cMix*qMap) - zeta]`                        | yes                |
| A2 | mathematica   | 65   | `expectZero["Xi_F1(Upsilon) - 1369 Upsilon_w", xiF1FromUpsilon - 1369*upsilonW]`                 | no (tautology)     |
| A3 | mathematica   | 66   | `expectZero["Xi_F1(Theta) - 136900 Theta_w", xiF1FromTheta - 136900*thetaW]`                     | no (tautology)     |
| A4 | mathematica   | 80   | `expectApprox["natural-window ordering gap", zetaPlusChi - zetaMinusChi, 0.00130621926024, 1e-12]` | no (hardcoded vs hardcoded) |
| A5 | mathematica   | 81   | `expectApprox["hard-ceiling gap", zetaMaxF1 - zetaPlusChi, 9.6717311e-8, 1e-12]`                 | no (hardcoded vs hardcoded) |

A1 is substantive: the symbolic round-trip really can fail if `zetaReq` and `qMap` are inconsistently defined; I worked the algebra by hand and confirmed the cancellation is non-trivial (numerator and denominator each pick up a `cMix*(1-eps)/(1-eps*zeta)` factor that cancels to leave `zeta`).

A2/A3 cannot fail: `xiF1FromUpsilon` is defined on line 58 as `lambdaEll^2*upsilonW` with `lambdaEll = 37` declared on line 55, so `xiF1FromUpsilon - 1369*upsilonW` is the closed-form `(37^2 - 1369)*upsilonW = 0` by integer arithmetic. Same for A3 with `100*lambdaEll^2*thetaW - 136900*thetaW`.

A4/A5 cannot fail either way it matters: both operands of the subtraction are arbitrary hardcoded high-precision floats (lines 68-72), and the expected difference target (lines 80-81) is also a hardcoded float that was clearly precomputed from the same numbers. The assertion checks IEEE arithmetic, not physics.

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:55-66`

**What's wrong:**
Lines 55-58 declare:
```
lambdaEll = 37;
kappaF1 = 12321/5;
etaF1 = 37;
xiF1FromUpsilon = lambdaEll^2*upsilonW;
xiF1FromTheta = 100*lambdaEll^2*thetaW;
```
Then lines 65-66 assert:
```
expectZero["Xi_F1(Upsilon) - 1369 Upsilon_w", xiF1FromUpsilon - 1369*upsilonW];
expectZero["Xi_F1(Theta) - 136900 Theta_w", xiF1FromTheta - 136900*thetaW];
```
Both residuals reduce to `(37^2 - 1369)*upsilonW` and `(100*37^2 - 136900)*thetaW` by direct substitution; they are 0 by construction regardless of any physics. The assertion cannot detect any incorrect upstream identity — it can only detect somebody changing the literal `37` on line 55 or the literal `1369`/`136900` on lines 65-66. (Note: the upstream verification of the substantive identity `Xi_F1 = Lambda_ell^2 * Upsilon_w` happens in `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:80-89` and is itself tautological in the same way; the meaningful physics is the *definition* `Lambda_ell = 37` from the underlying derivation, which is not exercised here.)

**Why this matters:**
A "PASS" line in the captured output (`PASS: Xi_F1(Upsilon) - 1369 Upsilon_w`) gives the false impression that the script verified a non-trivial identity. For a status-only consolidation that is supposed to demonstrate the full chain ties together, the only way to make this assertion non-trivial without re-deriving in this unit is to anchor `lambdaEll` against the upstream symbol it stands for (e.g., import the value from the upstream-defined ratio `100*Theta_w / Upsilon_w` or check `xiF1FromUpsilon - xiF1FromTheta.subs(thetaW -> upsilonW/100)`).

**Required change:**
Convert lines 65-66 from "subtract literal 1369/136900" into the substantive cross-check that the two `Xi_F1` formulas (the `Upsilon_w` route and the `Theta_w` route) agree under the upstream defining relation `Upsilon_w = 100*Theta_w`. Replace lines 65-66 with:
```
expectZero["Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)", (xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta];
```
This residual is non-trivially zero only when both routes share the same `lambdaEll`, and it can detect a sign flip or wrong power on either side.

**Verification:**
After the fix, line 65 of the script should show the new check name, and the captured output should print `Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta) = 0` followed by a `PASS:` line. The total number of `PASS:` lines in the output drops from 5 to 4 (line 66 is removed); the script must still exit 0.

### F2 — hardcoded_result

**Severity:** low
**Files:**
- `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:68-72,80-81`

**What's wrong:**
Lines 68-72 hardcode five `zeta` values as high-precision numeric literals:
```
zetaMinusChi = ToExpression["2.46622291347846`20"];
zetaPlusChi  = ToExpression["2.46752913273870`20"];
zetaMinusJ   = ToExpression["2.44257571477179`20"];
zetaPlusJ    = ToExpression["2.46752736855058`20"];
zetaMaxF1    = ToExpression["2.46752922945601`20"];
```
Lines 80-81 then assert:
```
expectApprox["natural-window ordering gap", zetaPlusChi - zetaMinusChi, ToExpression["0.00130621926024`20"], 10^-12];
expectApprox["hard-ceiling gap",            zetaMaxF1 - zetaPlusChi,    ToExpression["9.6717311`10*^-8"],    10^-12];
```
Both sides of each `expectApprox` are pure numeric literals; the "expected" target on the RHS is just the precomputed subtraction of the LHS literals. The assertion confirms that `2.46752913273870 - 2.46622291347846 == 0.00130621926024` and `2.46752922945601 - 2.46752913273870 == 9.6717311e-8`, both of which hold by basic IEEE arithmetic regardless of whether those `zeta` values mean anything physically. (The values themselves are legitimately carried forward from `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:45-49`, so the "hardcoded" finding is about the *form of the check*, not the provenance of the constants.)

**Why this matters:**
The `expectApprox` calls produce `PASS` lines that suggest the ordering and ceiling structure has been verified, but the only thing verified is float subtraction. The substantive claim — that `zeta_-^(chi) < zeta_+^(chi) < zeta_max^(F1)` and that the gaps are positive — can be expressed as direct inequalities, which would also catch a typo'd sign on any literal.

**Required change:**
Replace the two `expectApprox` calls on lines 80-81 with `expectZero` calls that exercise the ordering claim, not the precomputed difference. Insert before the closing `Print` block:
```
expectZero["chi-window ordering positive", If[TrueQ[zetaPlusChi > zetaMinusChi], 0, 1]];
expectZero["hard-ceiling gap positive",    If[TrueQ[zetaMaxF1   > zetaPlusChi  ], 0, 1]];
expectZero["J-window ordering positive",   If[TrueQ[zetaPlusJ   > zetaMinusJ   ], 0, 1]];
expectZero["chi-vs-J fail-side consistency", If[TrueQ[zetaPlusJ <= zetaPlusChi  ], 0, 1]];
```
Each `If` returns 0 when the inequality holds (so `expectZero` passes) and 1 when it fails (so `expectZero` fails with residual `1`). This exercises the ordering relations among the five carried-forward values rather than re-checking float subtraction. The unused `zetaMinusJ`/`zetaPlusJ` are also brought into the assertion set this way (currently they are printed but never asserted on, lines 74-78).

**Verification:**
After the fix, the captured output should show four new `PASS:` lines for the inequality checks and no more `natural-window ordering gap diff = ...` / `hard-ceiling gap diff = ...` lines. The script must still exit 0. If any of the upstream numeric values gets perturbed in a way that breaks an inequality, the corresponding assertion will report residual `1` and `Exit[1]`.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:39-46,55-57`

**What's wrong:**
Lines 39-46 define four symbolic expressions — `omegaPe`, `zetaPhys`, `zetaReq`, `qMap`, `rQuad` — that are then printed (lines 48-51) but never asserted against any structural property. Only `zetaReq`, `qMap`, and (implicitly) `zeta` are exercised by A1 on line 53. The freshly-defined `zetaPhys` (the physical-`zeta` formula in `pe`, `eta`, `kappa`) and the residual `rQuad = zetaReq - zetaPhys` are computed but the script never asserts anything about them.

Similarly, lines 55-57 define `kappaF1 = 12321/5` and `etaF1 = 37`, but neither value is used in any assertion — `zetaPhys` is built from generic symbols `pe`, `eta`, `kappa`, never evaluated at `(kappa -> kappaF1, eta -> etaF1)`. The output line 16 confirms `R_quad` is printed as a generic symbolic expression that the script does not pin to zero or to any structural form.

This is `insufficient_verification`: a status-only stage that defines five symbolic objects and two literal Family-1 parameters should at minimum substitute the Family-1 parameters into the physical formula and check that the resulting `zetaPhys` is consistent with the carried-forward numeric `zeta_max^(F1) = 2.46752922945601` at the upstream-fixed Pe (or that `rQuad` evaluates to a known closed form at `(kappa -> kappaF1, eta -> etaF1)`). As written, `kappaF1` and `etaF1` are dead code.

**Why this matters:**
The script presents itself as the "full reduced PDE write-up" consolidation, but the only substantive assertion (A1) is the algebraic round-trip of `zetaReq` against `Q(zeta)` — which is also verified upstream in stage 082. The Pe→zeta map (`zetaPhys`), which is the actual physical content tying the upstream numeric values to the demand-map framework, is never exercised. A reader of just the output would see five PASS lines and conclude the full chain is checked; in fact only one substantive algebraic identity is tested.

**Required change:**
Add a single substantive assertion that pins the carried-forward numeric `zetaMaxF1` to the physical formula `zetaPhys` evaluated at Family-1 parameters and the upstream-fixed Pe. Insert after line 51 (before the existing `expectZero["inverse demand map", ...]` on line 53):
```
peF1 = Pi/2;  (* upstream-fixed Pe at the Family-1 hard ceiling; mirror of scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py *)
zetaPhysF1Numeric = N[zetaPhys /. {pe -> peF1, kappa -> kappaF1, eta -> etaF1, y -> Sqrt[kappaF1]}, 20];
expectApprox["zetaPhys at (Pe=Pi/2, kappa_F1, eta_F1, y=sqrt(kappa_F1)) matches upstream zeta_max^(F1)",
  zetaPhysF1Numeric, zetaMaxF1, 10^-10];
```
BLOCKED IF: the actual upstream-fixed Pe for the Family-1 hard ceiling is not `Pi/2`, or the upstream definition of `y` at the maximizer is not `Sqrt[kappa_F1]`. In that case, append a `## Blocked: F3` block to the directive with a question rather than guessing — the correct binding must come from the same upstream source that produced `zetaMaxF1 = 2.46752922945601`.

**Verification:**
After the fix, the captured output should show a new line `zetaPhys at (Pe=Pi/2, kappa_F1, eta_F1, y=sqrt(kappa_F1)) matches upstream zeta_max^(F1) diff = <small number>` followed by `PASS:`. The script must still exit 0. If Codex cannot confirm the Pe / y values without reading docs, it should block this finding rather than apply a wrong binding.

## Independent-derivation check (Mathematica)

A SymPy counterpart for unit 084 does not exist (correctly, per the status-only carve-out); the upstream SymPy scripts (stages 081, 082, 086, 087, 098) provide the independent first-engine derivations. The Mathematica script's algebra is structurally similar to the SymPy bodies of stages 081-082 (same `zeta_req` formula, same `Q(zeta;eps_blk)` formula), but at this stage that is the *intended* behaviour — a status consolidation should mirror upstream definitions verbatim, not re-derive from scratch. Not flagged as `mathematica_transliteration`.

## Engine cross-check

n/a (no SymPy script). The carried-forward numeric values were independently produced by SymPy at stage 081 (`scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:45-49`) and the corresponding Mathematica twin at stage 081 — that is where engine agreement is established. Unit 084 just re-prints them; nothing for the two engines to disagree on here.

## Verdict justification

Three findings. None block: all are quality-of-verification issues rather than mathematical errors. F1 (tautological_check) catches two assertions (lines 65-66) that subtract the integer 1369 from `37^2`; F2 (hardcoded_result) catches two `expectApprox` calls (lines 80-81) where both operands and the target are hardcoded literals; F3 (insufficient_verification) notes that `zetaPhys`, `rQuad`, `kappaF1`, and `etaF1` are defined but never exercised in any assertion, and proposes a single substantive check that ties the upstream `zeta_max^(F1)` to the physical formula. Attacks that failed: I worked the inverse-demand-map algebra (A1) by hand and the cancellation is genuine, not symbolic accident; I checked that the carried-forward `zeta` constants are not invented locally but match upstream stage 081 to all 14 significant figures; I checked the integer arithmetic `37^2 = 1369` and `100*37^2 = 136900` to confirm F1's tautology claim (not the other way around). The script is honest about its status-only role, but its assertions could and should do more work.

## Self-test notes

- Variable independence trap: F3's proposed `zetaPhys /. {pe -> peF1, kappa -> kappaF1, eta -> etaF1, y -> Sqrt[kappaF1]}` substitutes into a `zetaPhys` that actually contains all four symbols (verified by reading line 43 of the existing script: `omegaPe^2*(kappa + Pi^2/4)/(kappa + y^2)` and `omegaPe` is defined on line 39-42 in terms of `pe`); the resulting numeric is well-defined.
- F2's `If[TrueQ[...], 0, 1]` idiom uses `TrueQ` so an undecidable comparison returns `False` (i.e., the assertion fails loudly) rather than a non-numeric residual — `expectZero` would see `1` and call `fail`.
- F3 is marked BLOCKED-IF on the upstream Pe / y binding because I cannot confirm `Pe = Pi/2` and `y = Sqrt[kappa_F1]` without reading the paper or notes/, which the audit rules forbid; if Codex can read the upstream sympy script (which is allowed under the directive) to confirm, it should apply; otherwise it should block.
