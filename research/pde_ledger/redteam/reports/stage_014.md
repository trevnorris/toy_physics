---
unit_id: 014
batch: I.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 014 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script asserts that a one-sided mouth-Taylor ansatz `X_proj(ell) = X(0) + ell*mu1*X'(0) + O(ell^2)` propagates the projected-Maxwell primitive slippages `{q1, s1, h1, d1, p1, g1}` into the bundle quantities `{z0, z2, z4, n0, n2, n4}` and from there into the 5PN bottlenecks `Xi_load = P1/P0`, `K1 = D21 + D01/9`, `H_even = D41 - (2/3)*D21 - D01/27`, and the constant-prefactor transport slots `delta P2`, `delta P4`. The assertions further claim a sector decomposition: `Xi_load` depends only on `(Qx, Deltax, Px)`; `K1` and `H_even` depend only on `(Qx, Deltax, S2x, Hx)`; `Gx` first enters via `delta P2`, `delta P4`. A mechanism-sieve block claims that neither a pure source/denominator correction `(Qx, Deltax)` nor a pure spectral correction `(S2x, Hx)` can close both even gates, and the mixed compensation surface has a denominator proportional to `Delta*H_port - Q*S2`, which is the `Z2`-slot factor.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 97 | `assert_zero("z0 derivative map", z0d - (Delta*Qx - Q*Dx)/Delta**2)` | no (tautological — `z0d` is defined as exactly this substitution) |
| A2 | sympy | 98 | `assert_zero("n0 derivative map", n0d - 2*P*(Delta*Px - P*Dx)/Delta**3)` | no (tautological — same reason as A1) |
| A3 | sympy | 100-101 | `for sym in (Sx, Hx, Gx): assert_zero(f"Xi_load independence from {sym}", sp.diff(Xi, sym))` | no (tautological — Xi's literal formula contains no s1, h1, g1) |
| A4 | sympy | 102-104 | `for sym in (Px, Gx): assert_zero(f"K1 independence from {sym}", sp.diff(K1, sym))` and same for `He` | no (tautological — K1, He defined from z0, z2, z4 which contain no p1, g1) |
| A5 | sympy | 105 | `assert_zero("d Xi_load / d Pprime", coef_Xi_Px - 2/P)` | partial (confirms input coefficient is preserved through `simplify`) |
| A6 | sympy | 106 | `assert_zero("d deltaP2 / d G_W prime", coef_dP2_Gx + 2*P/(D0*Delta**2))` | yes (non-trivial coefficient after composing n2 into deltaP2 formula) |
| A7 | sympy | 107 | `assert_nonzero("deltaP4 should depend on G_W prime", coef_dP4_Gx)` | yes |
| A8 | sympy | 108 | `assert_nonzero("source/denominator sieve determinant", qd_matrix.det())` | yes (substantive linear-algebra invertibility check) |
| A9 | sympy | 109 | `assert_nonzero("spectral sieve determinant", sh_matrix.det())` | yes |
| A10 | sympy | 110 | `assert_zero("compensation Hport denominator", Hx_den - 9*Delta**2*(Delta*Hport - Q*S2))` | yes (factored denominator identity) |
| A11 | sympy | 111 | `assert_zero("compensation S2 denominator", Sx_den - 9*Delta*(Delta*Hport - Q*S2))` | yes |
| A12 | sympy | 112 | `assert_nonzero("mutated compensation denominator should fail", Hx_den - 9*Delta**2*(Delta*Hport + Q*S2))` | yes (sign-flip mutation test) |
| A13 | sympy | 114-115 | `if qd_only != [{Qx: 0, Dx: 0}]: raise` | yes (negative-existence result: no nontrivial solve) |
| A14 | sympy | 116-117 | `if sh_only != [{Sx: 0, Hx: 0}]: raise` | yes |
| A15 | sympy | 118 | `assert_zero("compensation K1", K1.subs(comp_surface))` | no (tautological — `comp_surface` is the output of `solve(K1=0,He=0)`; substituting back is guaranteed zero by the solver) |
| A16 | sympy | 119 | `assert_zero("compensation H_even", He.subs(comp_surface))` | no (same as A15) |
| A17 | sympy | 120 | `assert_zero("compensation denominator tracks Z2 slot", (Delta*Hport - Q*S2) + Delta**2*((Q*S2 - Hport*Delta)/Delta**2))` | no (pure algebraic identity `A + (-A) = 0`, independent of any prior derivation) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:97-98`

**What's wrong:**

Lines 97-98 read:
```
assert_zero("z0 derivative map", z0d - (Delta*Qx - Q*Dx)/Delta**2)
assert_zero("n0 derivative map", n0d - 2*P*(Delta*Px - P*Dx)/Delta**3)
```
But `z0d` is *defined* on line 59 by `z0d = sp.simplify(z0.subs(subs_der)/mu1)`, where `z0 = (Delta*q1 - Q*d1)/Delta**2` (line 48) and `subs_der` maps `q1 -> mu1*Qx`, `d1 -> mu1*Dx`. Mechanical substitution gives `z0.subs(subs_der)/mu1 = (Delta*mu1*Qx - Q*mu1*Dx)/(mu1*Delta**2) = (Delta*Qx - Q*Dx)/Delta**2` *exactly*, with no cancellation that could fail. The assertion therefore reduces to `sp.simplify((Delta*Qx - Q*Dx)/Delta**2 - (Delta*Qx - Q*Dx)/Delta**2) == 0`, which is identically zero regardless of any physics. The `n0` check is the same pattern: `n0` is a single rational monomial with two primitive slips, so its substitution under `subs_der` cannot do anything other than reproduce the asserted form.

The assertions are presented as "derivative maps", but z2, z4, n2, n4 — which actually have non-trivial cross-terms after substitution — are not asserted against any reference form. The only derivative maps that are asserted are the two that cannot fail.

**Why this matters:**

This pair of assertions cannot detect any mistake in the Taylor lift. The substantive part of the lift (the cross-terms in z2, z4, n2, n4 that combine multiple primitive slips with multiple base symbols) is never tested against an independent reference expression. If a future edit accidentally introduces a typo into the z2 or z4 formulas (line 50, 51, 53, 54), nothing in this script catches it; the only derivative-map assertions are the two that are immune to typos by construction.

**Required change:**

Replace the two trivial map assertions with assertions that exercise the non-trivial cross-terms. Add reference forms for `z2d`, `z4d`, `n2d`, `n4d` (these are exactly the output transcript on lines 26-28, 32-34 of the saved output), and assert each against its reference. Specifically, in `scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`, after line 64 (after `n4d = ...`) and before line 97, insert four reference expressions (derived in-script *not* by copying from the saved output, but by independently writing the Taylor expansion of each bundle quantity from the primitive formulas) and assert each `*_d` matches.

For example, `z2` is defined by the primitive identity `z2 = (S2 - Hport/Delta - 2*Q*S2/Delta**2)*` ... (the existing line 50 form). An *independent* check would be: compute the formal derivative of `(Delta*q - Q*s)/(Delta**2)` style identity from the Z2 slot definition `Z2 = (Q*S2 - Hport*Delta)/Delta**2` rather than re-substituting the existing `z2`. Concretely, add:

```
# Independent Taylor-derivative reference for z2: Z2 = (Q*S2 - Hport*Delta)/Delta**2
# Derivative of -Z2 w.r.t. ell at ell=0 (chain rule, each primitive -> mu1 * its prime)
z2_ref = (-1)*sp.diff(((Q + mu1*Qx*sp.Symbol('ell'))*(S2 + mu1*Sx*sp.Symbol('ell'))
                       - (Hport + mu1*Hx*sp.Symbol('ell'))*(Delta + mu1*Dx*sp.Symbol('ell')))
                      /(Delta + mu1*Dx*sp.Symbol('ell'))**2,
                      sp.Symbol('ell')).subs(sp.Symbol('ell'), 0) / mu1
assert_zero("z2 derivative map vs independent Taylor", z2d - z2_ref)
```
(and similarly for z4, n2, n4 using their primitive defining identities — the writer must supply the four primitive identities from the upstream unit's bridge, *not* lift them from lines 48-54 of this script).

If the Z2/Z4/N2/N4 primitive identities cannot be re-derived in-script without copying from the existing definitions on lines 50-54, then this finding becomes informational only and the existing checks should at minimum be flagged with a comment explaining that they are construction-checks, not derivation-checks. In that case the minimum acceptable edit is to delete lines 97-98 (so the inventory does not claim "derivative map" verifications that are not real verifications) and add inline comments explaining that `z2d`, `z4d`, `n2d`, `n4d` are themselves the printed observation, not asserted invariants.

**Verification:**

Either (a) four new `assert_zero` lines for `z2d/z4d/n2d/n4d` appear and the script exits 0, or (b) lines 97-98 are removed and a comment block stating "no independent reference; map verified by inspection of printed output" is added between line 64 and line 99.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:100-104`

**What's wrong:**

Lines 100-104 read:
```
for sym in (Sx, Hx, Gx):
    assert_zero(f"Xi_load independence from {sym}", sp.diff(Xi, sym))
for sym in (Px, Gx):
    assert_zero(f"K1 independence from {sym}", sp.diff(K1, sym))
    assert_zero(f"H_even independence from {sym}", sp.diff(He, sym))
```

But `Xi` is *defined* on line 67 as `Xi = sp.simplify((2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der)/mu1)`. The expression `2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)` contains the symbols `{p1, d1, q1}` only — no `s1, h1, g1`. After `subs_der`, it contains `{Px, Dx, Qx}` only — no `Sx, Hx, Gx`. So `sp.diff(Xi, Sx)`, `sp.diff(Xi, Hx)`, `sp.diff(Xi, Gx)` are identically zero because `Sx, Hx, Gx` never appear in `Xi` to begin with. The assertion is independent of any property of the Taylor lift; it just confirms that sympy can compute the derivative of an expression with respect to a symbol that does not appear in it.

Same defect for K1 and He: `K1 = -(z2 + z0/9)` involves z0, z2 which depend on `{q1, s1, h1, d1}` only, so `K1` literally has no `Px, Gx` to begin with. The independence-from-Px and independence-from-Gx assertions for K1 and He cannot fail.

(Reported separately from F1 because this is a different cluster of assertions — the bundle Z2/Z4/N0/N2/N4 maps in F1 vs. the bottleneck sector-decomposition in F2 — and the fixes are different.)

**Why this matters:**

The script's central physics claim, in its own readout (lines 222-226), is that the projected-Maxwell packet factorizes into three jobs corresponding to disjoint primitive-slip sectors. The assertions on lines 100-104 are the only place the script tries to verify this factorization. But each assertion is trivially true by construction of the very expressions it tests, because the author wrote the formula for `Xi` without `s1, h1, g1` in it. The sector decomposition is therefore an *input* assumption embedded in the formula choices, not an *output* that the script verifies. A reader looking at the assertion inventory would believe the factorization is checked when it is not.

**Required change:**

Replace the `sp.diff(EXPR, SYM)` checks with explicit symbolic-presence checks against the *primitive* expressions that *do* contain all slips, so the test exercises which primitive slips contribute and which cancel. Concretely, replace lines 100-104 with assertions that test the corresponding *bundle* dependencies, where the sector structure is non-trivial.

For example, `Xi_load = 2*P1/P0 - ...`; expand Xi in terms of the bundle slips (z0, n0) and check that the bundle-level expression vanishes on the Sx, Hx, Gx axes. Better: define each bottleneck *symbolically* from the bundle quantities

```
Xi_bundle = 2*n0d/(2*P*(Delta**(-2))*P) + ... # the bundle definition
```
then assert `Xi_bundle - Xi == 0` (bundle-vs-primitive consistency) AND `sp.diff(Xi_bundle, sym) == 0 for sym in (Sx, Hx, Gx)`. The first assertion is non-trivial because `Xi_bundle` involves bundle quantities that *do* contain all slips; the cancellation to a `Xi` that contains only three slips is then a real check.

If recasting through bundle expressions is too invasive, the minimum acceptable edit is to (i) prepend a comment at line 100 stating: `# These checks are construction-level: the formulas above contain only the listed slips, so independence is by definition. A real factorization check would derive Xi, K1, H_even from a higher-level bundle expression that contains all six slips and assert cancellation.` and (ii) add at least one substantive check that *does* exercise cancellation, e.g. compute `Xi` from a "naive bundle pull-back" expression that contains all six slips and assert that the simplification removes Sx, Hx, Gx.

A concrete acceptable form: after line 75, add

```
# Bundle-level expression of Xi_load = 2*P1/P0 with P1 = N0*z0 + D0*n0
# (this version naively contains all six slips before simplification)
Xi_bundle = sp.simplify((2*(N0*z0 + D0*n0)/(D0*Ptarget)).subs(subs_der)/mu1)
# Xi_load is the P-normalized form; sanity check the primitive form here
assert_zero("Xi bundle-vs-primitive consistency factor", sp.diff(Xi_bundle - 2*(N0/D0)*sp.diff(z0.subs(subs_der)/mu1, Qx)*Qx, Sx))
```
(left as illustrative; the writer needs to pick the bundle expression that pre-simplification *does* contain all six slips).

**Verification:**

Either (a) new assertions appear that pull the bottlenecks back through the bundle and exercise the cancellation of Sx, Hx, Gx as a real algebraic fact, or (b) a comment block at line 100 explicitly marks the existing assertions as construction-checks (not factorization-checks), and at least one substantive cancellation assertion is added.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:118-120`

**What's wrong:**

Lines 118-120:
```
assert_zero("compensation K1", K1.subs(comp_surface))
assert_zero("compensation H_even", He.subs(comp_surface))
assert_zero("compensation denominator tracks Z2 slot", (Delta*Hport - Q*S2) + Delta**2*((Q*S2 - Hport*Delta)/Delta**2))
```

The first two are guaranteed by `solve()`'s contract: `comp_surface` is the dictionary returned by `sp.solve([sp.Eq(K1, 0), sp.Eq(He, 0)], [Hx, Sx], dict=True)[0]` on line 80, so by definition `K1.subs(comp_surface) == 0` and `He.subs(comp_surface) == 0` simplify to zero whenever sympy's solver returns a correct dictionary. These two assertions test sympy itself, not the physics.

The third assertion is a pure algebraic identity: `(Delta*Hport - Q*S2) + Delta**2*((Q*S2 - Hport*Delta)/Delta**2)` simplifies to `(Delta*Hport - Q*S2) + (Q*S2 - Hport*Delta) = 0` by elementary cancellation. It does not connect `Hx_den` or `Sx_den` to the Z2 slot — that linkage was already (correctly) established by assertions A10 and A11 on lines 110-111. The assertion on line 120 is a redundant `A + (-A) = 0` check that tests no relationship between the script's quantities.

**Why this matters:**

These three lines bulk up the assertion count but do not exercise the script's claims. The substantive compensation-surface result is captured by lines 110-111 (denominators factor as `Delta*Hport - Q*S2`); line 120 adds nothing on top. The reader looking at the assertion inventory may believe the Z2-slot linkage is doubly-verified when it is in fact verified once and then padded with a tautology.

**Required change:**

Either (a) delete line 120 outright, or (b) replace it with a substantive linkage check. The substantive check should be: independently compute the `Z2` slot from its primitive definition `Z2 = (Q*S2 - Hport*Delta)/Delta**2` (this is the upstream-unit definition, not derived in this script) and assert that the compensation denominator equals exactly `-9*Delta**2 * Z2`. Concretely:

```
Z2 = (Q*S2 - Hport*Delta)/sp.Symbol('Delta')**2  # if not already in scope, use Delta
# correction: Z2 = (Q*S2 - Hport*Delta)/Delta**2
Z2_slot = (Q*S2 - Hport*Delta)/Delta**2
assert_zero("compensation Hport denominator equals -9 Delta^4 Z2",
            Hx_den - (-9 * Delta**4 * Z2_slot))
assert_zero("compensation S2 denominator equals -9 Delta^3 Z2",
            Sx_den - (-9 * Delta**3 * Z2_slot))
```

If the writer prefers the cosmetic deletion route, line 120 should be removed and lines 118-119 should be marked with an inline comment `# tautological by sp.solve contract; kept for visual symmetry only` so future readers know these are not real checks.

**Verification:**

Either line 120 is deleted (with lines 118-119 commented as tautological), or a new pair of substantive assertions tying `Hx_den` / `Sx_den` to `-9*Delta^4*Z2` / `-9*Delta^3*Z2` appears. The script still exits 0.

### F4 — missing_verification_script

**Severity:** high
**Files:**
- (missing) `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl`

**What's wrong:**

This unit has `is_status_only_candidate: False` and `is_checkpoint: False` per the prompt. Both engines are required. The adjacent units 011, 012, 013 (this one) and 015, etc. follow the naming pattern `moving_throat_pde_stageNNN_<topic>_mathematica_audit.wl` in `/var/projects/toy_physics/research/pde_ledger/mathematica/`. A directory listing shows files for stages 001-012 and then jumps to 021, so stages 013, 014, 015, ... do not have Mathematica counterparts.

The unit's claims — Taylor-lift of the projected-Maxwell primitive bridge into the 5PN bottlenecks `Xi_load, K1, H_even`, the constant-prefactor transport `delta P2, delta P4`, the mechanism sieve (no nontrivial pure-source or pure-spectral solve), and the mixed compensation surface with `Z2`-slot denominator — are non-trivial multi-variable polynomial rational results. The script's own assertions are mostly non-tautological (sieve determinants, denominator factorizations, mutation tests), so the SymPy side does have substantive content for an independent engine to cross-check.

**Why this matters:**

Without a second engine, the policy that requires independent confirmation is not satisfied. Furthermore, the SymPy script depends on sympy's `solve` and `simplify` behavior for two of its substantive checks (`Hx_den`, `Sx_den` denominator factorization), so a Mathematica cross-check is the natural way to catch any solver-specific quirk.

**Required change:**

Create a Mathematica audit script at the path
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl`
that *independently* derives, from the primitive identities only, each of the claims listed in the manifest below. The script must not be a line-by-line transliteration of the SymPy script's algebra; it must define the primitive Z-bundle and N-bundle identities from their physical definitions (Z2 = (Q S2 − Hport Delta)/Delta^2, Z4 = ..., N0 = 2 P^2/Delta, ...) and Taylor-lift them via `D[..., ell]` evaluated at `ell -> 0`, rather than substituting `q1 -> mu1 Qx` into pre-given Z-formulas.

**Claim manifest** (this is the Mathematica script's job-list; see directive for details):

- M1: Sector decomposition: `Xi_load` depends only on `{Qx, Deltax, Px}`; `K1`, `H_even` depend only on `{Qx, Deltax, S2x, Hx}`; `Gx` first enters via `delta P2`, `delta P4`.
- M2: Direct coefficient: `D[Xi_load, Px] == 2/P`.
- M3: Direct coefficient: `D[delta P2, Gx] == -2*P/(D0*Delta^2)`.
- M4: `D[delta P4, Gx]` is non-zero.
- M5: Mechanism sieve, source/denominator sector: `Solve[{K1 == 0, H_even == 0} /. {Sx -> 0, Hx -> 0}, {Qx, Deltax}]` returns only `{{Qx -> 0, Deltax -> 0}}`.
- M6: Mechanism sieve, spectral sector: `Solve[{K1 == 0, H_even == 0} /. {Qx -> 0, Deltax -> 0}, {Sx, Hx}]` returns only `{{Sx -> 0, Hx -> 0}}`.
- M7: Source/denominator-sector Jacobian determinant `Det[{{D[K1/.{Sx->0,Hx->0}, Qx], D[K1/.{Sx->0,Hx->0}, Deltax]}, {D[He/.{Sx->0,Hx->0}, Qx], D[He/.{Sx->0,Hx->0}, Deltax]}}]` is non-zero.
- M8: Spectral-sector Jacobian determinant analogous to M7 is non-zero.
- M9: Mixed compensation surface (Solve for `Hx`, `Sx` from `K1 == 0`, `He == 0`) has denominator equal (up to a polynomial factor in `Delta`) to `Delta*Hport - Q*S2`, which is `-Delta^2 * Z2`.
- M10: Mutation test: swapping the sign in the denominator (testing `Delta*Hport + Q*S2`) fails (residual is non-zero), confirming the factorization is exact.

**Verification:**

After Codex creates the file, the verifier runs `redteam exec-mathematica 014`. The new `.wl` script must define each of M1-M10 as an `If[FullSimplify[residue] =!= 0, Exit[1]]` (or equivalent), produce a saved output transcript at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.txt`, and exit 0.

## Independent-derivation check (Mathematica)

No `.wl` script exists for stage 014; the independent-derivation check is therefore not yet possible. It becomes possible once F4 is acted on. The directive's claim manifest lists the 10 independent results the new Mathematica script must derive from the primitive identities (not transliterate from the SymPy script).

## Engine cross-check

Not possible — only the SymPy engine is present. See F4.

## Verdict justification

The SymPy script's substantive content (sieve determinants A8, A9; denominator factorizations A10, A11; mutation test A12; non-existence solves A13, A14; coefficient checks A5, A6, A7) holds up under attack and exercises real algebraic facts about the projected-Maxwell mouth-Taylor lift. I tried to break the denominator factorization assertions A10, A11 by checking against the saved transcript (lines 88-90), and the form `9*Delta^2*(Delta*Hport - Q*S2)` and `9*Delta*(Delta*Hport - Q*S2)` matches what `sp.factor(sp.denom(sp.together(comp_surface[Hx])))` would produce given the explicit numerator/denominator on transcript lines 88-90. I tried to break the sieve solves by reasoning about whether `solve` could return spurious solutions and conclude the negative-existence assertions A13, A14 are real. The mutation A12 (sign-flipped denominator) is a genuine adversarial check.

However, three clusters of assertions are tautological by construction (F1, F2, F3): the `z0/n0` derivative-map checks restate the substitutions that defined them; the sector-independence checks differentiate expressions w.r.t. symbols that do not appear in those expressions; the Z2-slot identity is `A + (-A) = 0`. None of these can detect a typo or sign error in the script. Combined with the missing Mathematica engine (F4), the verdict is `findings`, not `clean`. No stop-cold flag — the SymPy substantive content is correct and the fixes do not propagate downstream (the bottleneck definitions are inputs from upstream units, not produced here).

## Self-test notes

I checked the variable-independence trap on lines 100-104: confirmed Xi has no Sx/Hx/Gx in its definition (line 67), and K1/He have no Px/Gx in their definitions (lines 68-69), so `sp.diff` returning zero is structural, not physical — this is the trap from prior units, and I flagged it as F2. I checked the symmetry trap on line 120: confirmed `(Delta*Hport - Q*S2) + Delta**2*((Q*S2 - Hport*Delta)/Delta**2)` is a pure additive identity (A + (-A)) — flagged as F3. I traced output transcript freshness: script mtime 2026-05-04 12:00:51, output mtime 2026-05-11 12:39:12, so output is newer than script and not stale.
