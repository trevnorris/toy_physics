---
unit_id: 065
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 065 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.txt`

## What the script claims to verify

The docstring lists six headline claims for the thin-wall confinement branch of the parent-wall potential V_conf(r;a)=V0 f((r-a)/ell): (1) g_phi = V0/ell; (2) the exact shell moment expansion I1 = 4*pi*ell*(a^2 J1 + 2 a ell J2 + ell^2 J3); (3) J2=0 for a centred symmetric wall layer; (4) the exact and thin-wall equilibrium gains and their O(ell) remainder; (5) the explicit V0_fail^2 and V0_suff^2 thresholds, with K_X cancelling once kappa = K_X L^2 / T_X is inserted; (6) the constant-compressibility reduction J1 = I_f / H_w. In practice the script *declares* (1), (2), (3), (6) as symbol definitions or substitutions and prints them, performs a manifestly-by-construction subtraction for the (4) remainder, and runs two genuinely non-trivial algebraic identities — the K_X-cancellation checks of (5).

## Assertion inventory

| #  | Script       | Line | Form                                                                                                          | Anchored to claim? |
|----|--------------|------|---------------------------------------------------------------------------------------------------------------|--------------------|
| A1 | sympy        | 81-84| `expect_zero("thin-wall remainder ...", G_eq_sym - G_eq_tw - 4*pi*V0**2*ell*J3/KX)`                            | no (tautological)  |
| A2 | sympy        | 104-107| `expect_zero("K_X cancellation in V0_fail^2", V0_fail_sq_k - TX*ell*Pe_req/(4*pi*a**2*L**2*J1*Deltainf))`     | partial            |
| A3 | sympy        | 108-111| `expect_zero("K_X cancellation in V0_suff^2", V0_suff_sq_k - TX*ell*Pe_req/(4*pi*a**2*L**2*J1*Delta0))`        | partial            |
| A4 | sympy        | 120-123| `expect_zero("constant-H fail threshold", V0_fail_const - Hw*TX*ell*Pe_req/(4*pi*a**2*L**2*If*Deltainf))`     | no (tautological)  |
| A5 | sympy        | 124-127| `expect_zero("constant-H suff threshold", V0_suff_const - Hw*TX*ell*Pe_req/(4*pi*a**2*L**2*If*Delta0))`        | no (tautological)  |
| A6 | mathematica  | 49-52 | `expectZero["thin-wall remainder ...", gEqSym - gEqTw - 4*Pi*v0^2*ell*j3/kx]`                                  | no (tautological)  |
| A7 | mathematica  | 73-76 | `expectZero["K_X cancellation in V0_fail^2", v0FailSqGeom - tx*ell*peReq/(4*Pi*a^2*len^2*j1*deltaInf)]`         | partial            |
| A8 | mathematica  | 77-80 | `expectZero["K_X cancellation in V0_suff^2", v0SuffSqGeom - tx*ell*peReq/(4*Pi*a^2*len^2*j1*delta0)]`            | partial            |
| A9 | mathematica  | 89-92 | `expectZero["constant-H fail threshold", v0FailConst - hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*deltaInf)]`        | no (tautological)  |
| A10| mathematica  | 93-96 | `expectZero["constant-H suff threshold", v0SuffConst - hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*delta0)]`           | no (tautological)  |

Anchoring summary:
- A2/A3 (and their Mathematica twins A7/A8) genuinely test the K_X exponent and the kappa substitution. A wrong sign in the exponent of K_X in V0_fail^2 would survive `solve` but fail the subtraction. These are real checks of claim (5).
- A1/A6 ("thin-wall remainder") are tautological. G_eq_sym is built from I1_sym = 4*pi*ell*(a^2 J1 + ell^2 J3), G_eq_tw is the J1 piece, and `4*pi*V0^2*ell*J3/KX` is the J3 piece. The residual is algebraically guaranteed zero by polynomial decomposition; it tests nothing about the underlying physics.
- A4/A5/A9/A10 ("constant-H ... threshold") are tautological substitutions. The "const" forms are produced by `.subs(J1, If/Hw)` on the geom expressions; the assertions then subtract the same algebraic re-arrangement. They confirm SymPy/Mathematica can do `1/J1 -> H_w/I_f`; they do not test the physical claim that h' is constant implies J1 = I_f/H_w.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:26-101`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:51-127`

**What's wrong:**
The `.wl` script is a line-for-line transliteration of the `.py` script. Three corresponding sections:

(a) Definition of g_phi and I1.
SymPy (l.61, l.66):
```
g_phi = sp.simplify(V0 / ell)
I1 = sp.simplify(4 * sp.pi * ell * (a**2 * J1 + 2 * a * ell * J2 + ell**2 * J3))
```
Mathematica (l.34-35):
```
gPhi = FullSimplify[v0/ell, ...];
i1 = FullSimplify[4*Pi*ell*(a^2*j1 + 2*a*ell*j2 + ell^2*j3), ...];
```
Identical algebraic form, identical variable choreography, only renamed.

(b) Thin-wall remainder residual.
SymPy (l.81-84):
```
expect_zero("thin-wall remainder ...", sp.simplify(G_eq_sym - G_eq_tw - 4 * sp.pi * V0**2 * ell * J3 / KX))
```
Mathematica (l.49-52):
```
expectZero["thin-wall remainder ...", gEqSym - gEqTw - 4*Pi*v0^2*ell*j3/kx];
```
Byte-for-byte the same residual expression.

(c) K_X-cancellation residual.
SymPy (l.104-107):
```
expect_zero("K_X cancellation in V0_fail^2", V0_fail_sq_k - TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Deltainf))
```
Mathematica (l.73-76):
```
expectZero["K_X cancellation in V0_fail^2", v0FailSqGeom - tx*ell*peReq/(4*Pi*a^2*len^2*j1*deltaInf)];
```
Same residual.

Both engines pull V0_fail^2 from the same `G_eq_tw == Pe_req/(kappa Delta_inf)` inversion and then subtract the same closed-form. Neither derives I1 from an actual ∫ (f')^2/h' (a+ell xi)^2 dxi over a chosen profile; neither derives J2=0 from a parity argument over a symmetric f; neither derives J1 = I_f/H_w from a constant-h' limit. So the engines are agreeing because they encode the same symbolic identity, not because they verify the physics from two independent routes.

**Why this matters:**
The second-engine policy requires Mathematica to provide an independent verification, not echo SymPy's algebra. A sign error or factor mistake propagated into the symbolic ansatz of either engine would survive both checks because the second engine literally copies the first engine's polynomial.

**Required change:**
Add at least one substantive independent derivation block to the `.wl` script. Concretely, derive I1 from the actual shell integral by choosing a concrete symmetric profile (Gaussian f' = (1/sqrt(2 Pi sigma^2)) Exp[-xi^2/(2 sigma^2)] with sigma = 1) and constant h' = h0, then:
  - Compute J1_num = Integrate[(f'[xi])^2/h0, {xi, -Infinity, Infinity}] and verify it matches I_f/H_w with I_f, H_w numerical.
  - Compute J2_num = Integrate[xi*(f'[xi])^2/h0, {xi, -Infinity, Infinity}] and assert J2_num == 0 (parity).
  - Compute the full shell weight (a + ell xi)^2 (f')^2 / h' integrated over xi (formally over (-Infinity, Infinity) which is the thin-wall limit of the layer), compare against 4*pi*ell*(a^2 J1_num + 2 a ell J2_num + ell^2 J3_num) with J3_num computed independently. This anchors claim (2) and claim (3) on Mathematica's side without echoing the SymPy polynomial.

Keep the existing K_X cancellation `expectZero`s in place — those remain valid as algebraic identity checks once the new independent shell integral is added.

**Verification:**
After patch, `mathematica/...stage065_...wl` contains a new block before line 54 that performs at least two `expectZero` calls anchored on concrete Gaussian f' integrals (J2 = 0 by parity, and the shell-integral expansion versus the polynomial form). The verifier confirms the script still exits 0 and the new `PASS:` lines appear in the output transcript.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:81-84`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:49-52`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:104-127` (constant-H block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:83-96`

**What's wrong:**
Two assertion families are tautological by construction.

(a) "thin-wall remainder after dropping O(ell/a) correction" — SymPy l.81-84 and Mathematica l.49-52.
`G_eq_sym` is defined as `g_phi^2 * I1_sym / KX = (V0/ell)^2 * 4*pi*ell*(a^2 J1 + ell^2 J3)/KX = 4*pi*V0^2*(a^2 J1)/(KX ell) + 4*pi*V0^2*ell*J3/KX`.
`G_eq_tw` is defined as `4*pi*a^2*V0^2*J1/(KX*ell)`, the first term.
The third term subtracted in the assertion, `4*pi*V0^2*ell*J3/KX`, is the second term.
So `G_eq_sym - G_eq_tw - 4*pi*V0^2*ell*J3/KX ≡ 0` as a polynomial identity in (J1, J3, V0, a, ell, KX). The assertion cannot fail under any physics; it only checks that SymPy can subtract a polynomial from itself.

(b) "constant-H fail/suff threshold" — SymPy l.120-127 and Mathematica l.89-96.
`V0_fail_const = V0_fail_sq_k.subs(J1, If/Hw) = TX*ell*Pe_req/(4*pi*a^2*L^2*(If/Hw)*Deltainf) = Hw*TX*ell*Pe_req/(4*pi*a^2*L^2*If*Deltainf)`.
The assertion subtracts `Hw*TX*ell*Pe_req/(4*pi*a^2*L^2*If*Deltainf)` — the exact result of the substitution — so the residual is identically zero by the substitution itself. The check tests that SymPy's `subs(J1, If/Hw)` is consistent with manual rewriting; it does not test the physical claim that constant h' implies J1 = I_f/H_w.

**Why this matters:**
Three of the five SymPy `expect_zero`s (and three of the five Mathematica `expectZero`s) are guaranteed by construction. They contribute zero attack surface — flipping the sign of J3 or of the constant-H reduction would not break them in the script-internal sense (it would just propagate consistently through both sides). The unit's substantive verification reduces to the two K_X-cancellation checks.

**Required change:**
Replace each tautological check with a check that actually exposes the underlying physics. Two concrete substitutions:

1. In both `.py` (l.81-84) and `.wl` (l.49-52): keep the residual print, but replace the tautological assertion with an independent test that the thin-wall expansion is the correct leading term. Use a concrete numeric profile (e.g. f'(xi) = Exp[-xi^2] with constant h' = 1) to compute the exact shell integral `∫ (a + ell xi)^2 (f'(xi))^2 dxi` for fixed numeric a, ell, and verify that
   exact_integral / (a^2 J1_num) → 1 as ell/a → 0 (e.g. take ell/a = 0.01 and 0.001, check the second is closer to 1). At minimum, replace the current `expect_zero` with a check that the ratio (G_eq_tw + 4*pi*V0^2*ell*J3/KX) / G_eq_sym evaluated under J2 -> 0 equals 1 *only at the moment level*, which is what we already have — so this finding's required change is: add the concrete-profile numeric check to make the assertion non-tautological.

   Minimal concrete patch: after the existing `expect_zero("thin-wall remainder ...")` (which can stay as an algebraic identity check), add a new `expect_zero` of the form: pick numeric J1_num, J3_num computed from a definite Gaussian profile, fix a/ell = 100 and a/ell = 1000, and assert that `(G_eq_sym - G_eq_tw)/G_eq_tw - (ell^2 J3_num)/(a^2 J1_num)` evaluates to zero (this tests the O(ell^2/a^2) scaling, which is non-trivial because the script's claim is that the dropped term scales like ell/a^2 relative to the leading 1/ell — i.e. the *ratio* of correction to leading term goes like (ell/a)^2 J3/J1).

2. In both `.py` (l.120-127) and `.wl` (l.89-96): the constant-h' reduction J1 = I_f/H_w is currently asserted by `.subs(J1, If/Hw)`. Replace this with an independent derivation: define I_f := ∫ (f'(xi))^2 dxi for a specific symmetric f (e.g. Gaussian), define J1_via_def := ∫ (f'(xi))^2 / H_w dxi (with constant H_w), and assert J1_via_def - I_f/H_w == 0 by direct integration in SymPy/Mathematica. This anchors the constant-h' claim on an actual integral identity.

**Verification:**
After patch, the .py and .wl each contain:
  - One additional `expect_zero`/`expectZero` that numerically/symbolically anchors the thin-wall ratio at two values of ell/a (or equivalently asserts the O((ell/a)^2) scaling using definite numeric J1, J3).
  - One additional `expect_zero`/`expectZero` deriving J1 = I_f/H_w from the integral definitions rather than substituting it.
Both scripts still exit 0; the output transcript shows the new PASS lines.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:51-84`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl:26-52`

**What's wrong:**
Of the six docstring claims, four (1, 2, 3, 6) are never tested by any assertion — they are introduced as definitions or substitutions and then printed.

  - Claim (1) "g_phi = V0/ell": l.61 sets `g_phi = sp.simplify(V0 / ell)` and only prints it. There is no derivation that the support loading amplitude of V0 f((r-a)/ell) under d/dr = (1/ell) d/dxi yields V0/ell, and no comparison against a competing form (e.g. V0 only, or V0/(2 ell)).
  - Claim (2) "I1 = 4*pi*ell*(a^2 J1 + 2 a ell J2 + ell^2 J3)": l.66 defines I1 as this polynomial directly. There is no integration ∫ chi_phi^2/h' (a+ell xi)^2 dxi performed against a concrete profile to confirm the (a^2, 2 a ell, ell^2) coefficients (factors of 2 are exactly the kind of thing this script claims to verify).
  - Claim (3) "J2 = 0 for symmetric layer": l.70 substitutes `J2 -> 0`. There is no parity argument — no check that ∫ xi (f'(xi))^2 dxi = 0 for an odd-symmetric f' (or even-symmetric f).
  - Claim (6) "J1 = I_f/H_w when h' constant": l.116 substitutes `J1 -> If/Hw`. (Covered by F2 too.)

Neither engine integrates a concrete f profile to exercise any of these.

**Why this matters:**
The script is presented as the verification of the thin-wall confinement branch. Its only non-tautological assertions test claim (5) (K_X cancellation in the inverted threshold expressions). Claims (1)-(3) and (6) — which form the algebraic chain leading into (5) — are taken on faith. A factor-of-2 error in claim (2)'s middle term, a sign error in claim (1)'s 1/ell, or a wrong parity argument in (3) would all survive both engines unchanged.

**Required change:**
Add concrete integral checks to both `.py` and `.wl`. Minimum acceptable set:

(a) g_phi from V_conf differentiation. Define V_conf(r) = V0 * f((r-a)/ell) with f a chosen symmetric symbolic profile (e.g. f(u) = Exp[-u^2]), compute (d V_conf / d r) explicitly with SymPy/Mathematica, and verify the amplitude of d V_conf/dr at the wall (r = a) equals V0/ell times f'(0). Since f'(0) is fixed by the profile, this anchors the 1/ell scaling.

(b) I1 polynomial coefficients. For a chosen f and constant h' = 1, compute
   I1_explicit = ∫_{-∞}^{∞} (f'(xi))^2 * (a + ell xi)^2 dxi   (after the 4*pi*ell shell factor)
and assert it expands to 4*pi*ell*(a^2 J1 + 2*a*ell*J2 + ell^2*J3) with J1, J2, J3 the corresponding moment integrals. The factor of 2 on the cross-term and the symmetric vanishing of J2 fall out of this single check.

(c) J2 = 0 by parity. With f symmetric (f(-xi) = f(xi)), assert ∫ xi (f'(xi))^2 dxi = 0 for the chosen profile — direct integration in SymPy/Mathematica.

These three additions transform claims (1), (2), (3) from declarations into checked derivations. Claim (6) is covered by F2 (b).

**Verification:**
After patch:
- `.py` and `.wl` each contain new blocks for (a)-(c) above with `expect_zero`/`expectZero` calls.
- Output transcripts show the new PASS lines.
- Both scripts still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` script does NOT derive the claim independently. Sections 26-52 mirror lines 51-84 of the `.py` exactly (same definitions of gPhi, i1, i1Sym, gEq, gEqSym, gEqTw, same residual expression for the "thin-wall remainder"). Sections 54-96 mirror the `.py`'s threshold and constant-H blocks identically. No integral over any concrete f profile appears. No parity argument for j2 appears. No alternative algebraic route (e.g. expand (a + ell xi)^2 first, then integrate term by term) is used. See F1 for quoted excerpts.

## Engine cross-check

The two engines produce identical results because they encode the same symbolic polynomials. Side by side (final V0_fail^2):
- SymPy: `Pe_req*T_X*ell/(4*pi*Delta_inf*J1*L**2*a**2)`
- Mathematica: `(ell*peReq*tx)/(4*a^2*deltaInf*j1*len^2*Pi)`
Same expression up to symbol renaming. All five SymPy `expect_zero` residuals print "0", and all five Mathematica `expectZero` residuals print "0". `engines_agree: true`, but the agreement is structural (both ran the same algebra) rather than corroborative.

## Verdict justification

Both scripts exist, run to exit 0, and agree on the printed forms — but the agreement is shallow because the .wl is a transliteration of the .py and four of the six docstring claims are encoded as definitions/substitutions rather than as checks. Two real algebraic identities are tested (K_X cancellation in V0_fail^2 and V0_suff^2); the other three "passes" per script are tautological. Verdict: `findings`, three findings, none stop-cold (the K_X cancellation results are correct and any downstream unit consuming V0_fail^2 = Pe_req*T_X*ell/(4*pi*Delta_inf*J1*L^2*a^2) will continue to see that expression — the fixes add verifications, they do not change the result).

## Self-test notes

- Variable independence: the proposed new `expect_zero`s involve integrals over `xi` of integrands that genuinely depend on `xi` (Gaussian f'(xi), and xi-weighted versions for J2); no degenerate `diff(EXPR, VAR)` with EXPR independent of VAR.
- Symmetry/parity: for the J2 = 0 check, integrand is xi * (f'(xi))^2 with f symmetric ⇒ f' antisymmetric ⇒ (f')^2 symmetric ⇒ xi * (f')^2 antisymmetric ⇒ integral over (-Infinity, Infinity) is zero. Verified by hand: parity is correct.
- Trivial-case pre-check: with f(u)=Exp[-u^2], f'(u)=-2 u Exp[-u^2], (f'(u))^2 = 4 u^2 Exp[-2 u^2]. Then J1 = ∫ 4 u^2 Exp[-2 u^2] du = 4 * sqrt(pi/2)/2/2 = sqrt(pi/2)/2 > 0 (nonzero, as required for J1 > 0 assumption in the script). J2 = ∫ u * 4 u^2 Exp[-2 u^2] du = 0 by parity (odd integrand). J3 = ∫ u^2 * 4 u^2 Exp[-2 u^2] du > 0 (even, positive integrand) — these all match the script's positivity assumptions on J1, J3.
- Path specifications: no `missing_verification_script` finding raised; both engines present.
