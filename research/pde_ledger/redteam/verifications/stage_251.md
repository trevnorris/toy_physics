---
unit_id: 251
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T13:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 251

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New independent Mathematica audit created at
`mathematica/moving_throat_pde_stage251_microscopic_damping_export_kernel_replacing_the_phenomenological_v_leg_envelope_law_mathematica_audit.wl`
(253 lines), covering claim manifest M1-M7 with repo-standard `expectZero`/`expectTrue`/`expectApprox`
helpers, `stripConditional`, and `Exit[1]` on failure.

**Assessment:**
Correct and genuinely independent route, not a transliteration of the `.py`:
- M1 (`.wl:83-90`): cubic coefficient via `Series[n0[omega],{omega,0,4}]` + `Coefficient[...,omega,2]` — a Series decomposition, whereas the SymPy script uses `Limit[N0/omega^2, omega->0]` (py:48). Different route, same target `eta0^2 omu0^4/delta0^2`.
- M5 (`.wl:162-186`): monotonicity AND positive-root uniqueness via `Resolve[ForAll[...]]` over the positive orthant, plus `F(0)<0`, `Limit[F,s->Infinity]==Infinity`, and a difference-quotient factorization identity — general positivity, not the SymPy "F' algebraic form + 2 sampled root counts." This is the Reduce/Resolve route the directive demanded.
- M6 (`.wl:190-201`): slowdown derived via `Series[fChar[s0+gBook*eps1] /. {kappaV->muEta s0^2, g3->gBook g3, g5->gBook g5}, {gBook,0,1}]`, solve O(g) coefficient, `expectZero` against the boxed shift.
- M7 (`.wl:211-248`): benchmark via exact rational inputs and `NSolve`; the F2 tautology is NOT re-introduced (M2 extracts s^5/s^3 coefficients from `kExpProjected` and checks `kProj^2` homogeneity, not a self-equality).
- Circularity warned by the prompt is ABSENT: neither script references the phenomenological `gamma_tot` law in any check (the only `gamma_safe` hits are docstring/banner prose at py:13, py:241 describing what the stage *replaces*).
Log shows 31 PASS, 0 FAIL, exit 0.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Old self-equality at py:81 (`Gamma5 - PiVm**2 * a**5 * beta0 * sminus/(27*cs**5*lamminus) == 0`,
identically zero by construction) is gone. Replaced (py:81-92) by:
structural coefficient extraction from `K_exp = Gamma3*s^3 + Gamma5*s^5` (`coeff_s5 - Gamma5 == 0`,
`coeff_s3 - Gamma3 == 0`); projection homogeneity `Gamma5.subs(PiVm,k*PiVm) - k**2*Gamma5 == 0`;
and three first-order microscopic-factor scalings of `P0_minus` in `beta0`, `sminus`, `1/lamminus`.

**Assessment:**
The X-built-from-then-asserted-against-the-same-factors tautology is removed. The new checks are
structural: the s^5/s^3 extraction verifies the quintic sits at s^5 with no even/cubic leakage (fails
if the kernel power structure is wrong); the `k**2` homogeneity fails if the `Pi_{V-}^2` projection
power is mistyped; the scaling checks fail if a microscopic-factor power slips. `Gamma5`'s symbolic
value (an imported placeholder per the card) is unchanged — confirmed identical to the original product
in the output. Matches the directive's specified acceptance criteria verbatim. The Mathematica mirror
(M2, `.wl:104-109`) carries the same structural checks. Note: `coeff_s5 - Gamma5` alone is weakly
discriminating since K_exp is built from Gamma5, but it is paired with the genuinely failing homogeneity
and scaling checks, exactly as the directive prescribed; this is acceptable, not a residual tautology.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`root_shift` (the boxed slowdown) is still hand-written and printed (py:137, py:153) as the directive
required. A genuine derivation was added (py:138-159): `F_weak = F.subs({kappaV: mu_eta*s0**2,
G3: g*G3, G5: g*G5, s: s0 + g*eps1})`, then `balance = series(F_weak,g,0,2).removeO().coeff(g,1)`,
`eps1_sol = solve(balance, eps1)[0]`, and the new assertion `assert sp.simplify(eps1_sol - root_shift) == 0`.

**Assessment:**
The shift is now derived from the characteristic polynomial F (via series in the bookkeeping kernel
strength g), not from the hand-written literal, and the two are asserted equal. The output confirms
"delta s (weak export)" and "delta s (derived)" print the identical expression
`s_0^2(-Gamma_3 - Gamma_5 s_0^2)/(2 mu_eta)`. The assertion would fail on a wrong coefficient (missing
factor of 2, wrong power of s0). The Mathematica mirror (M6, `.wl:190-201`) derives the same result by
an independent Series-in-gBook route. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- "delta s (weak export) = s_0**2*(-Gamma_3 - Gamma_5*s_0**2)/(2*mu_eta)" and "delta s (derived) = ..." identical (F3 derived==boxed).
- benchmark reproduces card: "safe half-plane rhs = 289.61004917557426", "G5hat_safe = 961.094295282802", "weight sc^2 = 0.3013336470698294".
- "All symbolic and numerical checks passed." / "# exit_code: 0".

**Mathematica:** exit=0. 31 PASS, 0 FAIL. Notable lines:
- "PASS: M1 omega^2 coefficient" (Series route), "PASS: M2 projection homogeneity", "PASS: M5 derivative positive for all s>0" / "PASS: M5 positive-root uniqueness quotient" (Resolve), "PASS: M6 root-shift coefficient", "PASS: M7 cubic/quintic positive root equals sc".
- "STAGE 251 MATHEMATICA AUDIT PASSED" / "# exit_code: 0".

**Output freshness:** confirmed. Output `.txt` mtimes (both 13:16:42) are newer than the SymPy script (13:13:32) and the new `.wl` (13:15:03). Outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. No symbolic form or numeric constant was changed — `Gamma5`, `root_shift`, and
all benchmark numbers are identical to the pre-fix script and the paper card. F2 and F3 only widened the
verification surface (replaced a tautology with structural checks; added a derivation+assertion of an
already-printed value), and F1 added a second engine. No derived result that downstream units depend on
changed.

## Side observations (non-blocking)

- The SymPy `N0_series` (py:47) is computed and printed but the actual coefficient used in the assertion comes from the independent `sp.limit` route (py:48); the series print is display-only. Not a defect.
- M5 in the `.wl` retains the algebraic `F'`-form check (`.wl:181`) in addition to the general Resolve positivity; harmless redundancy, strengthens rather than weakens.

## Verdict justification

All three findings are genuinely resolved. F2's identically-zero self-equality is replaced by structural
coefficient-extraction plus projection-homogeneity and microscopic-factor scaling checks that fail on a
wrong structure, with no value changed. F3's boxed slowdown is now derived from F via a series in kernel
strength and asserted equal to the printed expression, not merely restated. F1 adds a genuinely
independent Mathematica route (Series for the cubic, Resolve for monotonicity/uniqueness, Series-in-
strength for the slowdown), not a transliteration, and the prompt-warned circularity against the
phenomenological `gamma_tot` law is absent from every check in both engines. Both engines exit 0 with all
in-file assertions passing and outputs freshly regenerated. `material_change` is false.
