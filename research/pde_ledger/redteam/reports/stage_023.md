---
unit_id: 023
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 023 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.txt`

## What the script claims to verify

The scripts verify a "full grouped bundle" bookkeeping layer for a real P2 system in the moving-throat PDE program. Concretely: (i) the weighted grouped metric `Ggrp = diag(1,2,2)` admits three Ggrp-orthogonal directions `ebar, ea, eb` with squared-norms 5, 20, 4, the corresponding projectors are idempotent and partition the identity, and any grouped vector decomposes into the three projected components; (ii) a "one-port" 2x2 BdG-like response function `(Q - H ω²)/(Δ - S ω² + ω⁴)` and the squared transfer `(P - g_W ω²)² / (Δ - S ω² + ω⁴)²` have closed-form Taylor coefficients Z_n, N_n in ω⁰, ω², ω⁴; (iii) the grouped decomposition is linear, so the three components of `D_{An} = K_A - B_{An} - Z_{An}` (etc.) read off from each lane's grouped pieces; (iv) on the isotropic branch the prefactor coefficients of `D₀·(N₀ + N₂ω² + N₄ω⁴)/Dcons²` collapse to `u₂ = -D₂/D₀`, `u₄ = (D₂² - D₀D₄)/D₀²`, `P₀ = N₀/D₀`, `P₂ = (D₀N₂ - 2D₂N₀)/D₀²`, `P₄ = (D₀²N₄ - 2D₀(D₂N₂ + D₄N₀) + 3D₂²N₀)/D₀³`; (v) the "constant-prefactor" branch conditions (`P₂ = P₄ = 0`) yield specific values for `N₂` and `N₄`; (vi) the Stage-4/5 outgoing radiation transfer coefficient is `Γ₅,port = a⁵/(27 c_s⁵)` and the required normalization `m̂² P₀ = 54 G c_s⁵ / (5 a⁵ c⁵)`; (vii) first-order anisotropy laws and the four monotonicity derivatives of `P₀(K,B₀,Z₀,N₀)`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 75-80 | norms/orthogonality of grouped basis | yes |
| A2 | sympy | 94-100 | projector idempotency, orthogonality, partition of unity | yes |
| A3 | sympy | 109-112 | exact decomposition `x = xbar ebar + ax ea + bx eb` | partial (Pbar*x by construction expands to xbar*ebar; still tests Matrix algebra) |
| A4 | sympy | 140-149 | Z₀, Z₂, Z₄ series coefficients of one-port denominator | yes |
| A5 | sympy | 159-175 | N₀, N₂, N₄ series coefficients of squared one-port | yes |
| A6 | sympy | 233-243 | grouped linearity of D_An decomposition | partial (follows from linearity of grouped_parts; still useful end-to-end) |
| A7 | sympy | 247-252 | `Nbar0 - (N020 + 2*N021 + 2*N022)/5 = 0`, `aN0 - (2*N020 - N021 - N022)/10 = 0`, `bN0 - (N021 - N022)/2 = 0` | **no — tautological (grouped_parts returns those exact formulas)** |
| A8 | sympy | 288-295 | exact `u₂, u₄, P₀, P₂, P₄` formulas | yes |
| A9 | sympy | 305-306 | `P2.subs(N2, N2_target) == 0`, `P4.subs(..., N4_target) == 0` | **no — N2_target/N4_target are solutions to P2=0/P4=0, so substitution is 0 by construction** |
| A10 | sympy | 311-314 | `P0 - N0/D0 == 0` (named "P0 normalization target") | **no — duplicate of A8's `P0 - N0/D0` at line 290** |
| A11 | sympy | 327 | `Gamma5_port - a^5/(27 c_s^5) == 0` | yes (substantive series claim about h₂ Hankel function) |
| A12 | sympy | 334-337 | `ratio_target.subs(mhat,1) - 54 G c_s^5 / (5 a^5 c^5) == 0` | yes |
| A13 | sympy | 370-371 | first-order anisotropy formulas `du2`, `dP0` | yes |
| A14 | sympy | 386-389 | grouped relabelings of A13 | partial (relabeling check) |
| A15 | sympy | 418-421 | monotonicity derivatives of `P0(K,B0,Z0,N0)` | yes |
| B1-B14 | mathematica | various | line-by-line mirror of A1-A15 with same algebraic recipe | yes (algebraically) but transliteration concern |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:247-252`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:130-132`

**What's wrong:**
The three "Nbar0 / aN0 / bN0 formula" assertions check `Nbar0 - (N020 + 2*N021 + 2*N022)/5 == 0`, `aN0 - (2*N020 - N021 - N022)/10 == 0`, `bN0 - (N021 - N022)/2 == 0`. But `(Nbar0, aN0, bN0) = grouped_parts(N020, N021, N022)`, and `grouped_parts` is literally defined (sympy lines 207-211, mathematica lines 31-35) to return exactly those three expressions. So the assertions reduce to `(N020 + 2*N021 + 2*N022)/5 - (N020 + 2*N021 + 2*N022)/5 == 0` and the two analogous identities. The check cannot fail for any input.

SymPy line 207-211:
```
def grouped_parts(x20, x21, x22):
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx
```
Then SymPy line 224 and 247-252:
```
Nbar0, aN0, bN0 = grouped_parts(N020, N021, N022)
...
expect_zero("Nbar0 formula", Nbar0 - (N020 + 2 * N021 + 2 * N022) / 5)
expect_zero("aN0 formula", aN0 - (2 * N020 - N021 - N022) / 10)
expect_zero("bN0 formula", bN0 - (N021 - N022) / 2)
```

**Why this matters:**
These assertions present as evidence that the outgoing-transfer N-bundle inherits the same grouped decomposition rules as the D-bundle (subbanner II.2). They don't — they only verify that the function returns what it returns. If the grouped formulas were wrong, both sides would shift together. The substantive content (D_An decomposition linearity) is already exercised at sympy lines 233-243 / mathematica lines 121-129, which use independently-defined `Dbar0 = grouped_parts(d020, d021, d022)[0]` versus `Kbar - Bbar0 - Zbar0` where each piece comes from a separate grouped_parts call — that has algebraic content. The N-side does no analogous independent comparison.

**Required change:**
Replace the three tautological assertions with non-tautological checks that exercise linearity of `grouped_parts` for the N-bundle in the same end-to-end form used for the D-bundle. The simplest fix: introduce N20 = K20-like microscopic split that's algebraically distinct, then compare `grouped_parts(<sum>)` to `<sum of grouped_parts>`. Concretely, since N_An is independent of any wall/BdG decomposition in this script, replace the three lines with checks that compare `grouped_parts(N020 + N220, N021 + N221, N022 + N222)` to `(Nbar0 + Nbar2, aN0 + aN2, bN0 + bN2)` componentwise — this verifies additivity rather than re-stating the formula. Alternative: simply delete the three tautological lines, since II.2's comment "Nothing to prove beyond linearity, but verify a representative identity" already concedes the section is informational; outright deleting honest-codes the situation. Either is acceptable.

**Verification:**
After the edit, either (a) the three "Nbar0 formula / aN0 formula / bN0 formula" lines no longer appear in the script and don't appear in the saved output, or (b) the new assertions compare two independent grouped_parts outputs and the saved output shows the new check names with value 0.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:308-314`

**What's wrong:**
The assertion at line 311-314 is named `"P0 normalization target"` but reads `expect_zero("P0 normalization target", P0 - N0/D0)`. This is byte-for-byte the same expression already verified at line 290 (`expect_zero("P0 - N0/D0", P0 - N0/D0)`). It is therefore both (a) a duplicate of a passing assertion and (b) misnamed: the surrounding subbanner is "III.3 — Universal normalization product", but this assertion does not exercise the normalization product `m̂² P₀ = 54 G c_s⁵ / (5 a⁵ c⁵)` at all. The actual normalization is tested at lines 334-337 via `ratio_target.subs(mhat,1) - 54 G c_s^5 / (5 a^5 c^5) == 0`. The intermediate `P0_target = sp.solve(sp.Eq(mhat**2 * P0, ...), N0)[0]` at line 309 is computed and then discarded; nothing tests it.

The Mathematica script does not contain an analogous misnamed assertion — `ratioTarget` is defined directly and tested once at line 178. SymPy has the duplicate but not Mathematica, which is a minor engine asymmetry; the Mathematica version is the cleaner of the two.

**Why this matters:**
A misnamed assertion creates the false impression that the universal normalization is tested in two places when in fact only line 334-337 tests it. A future reader reorganizing or pruning checks could remove line 334-337 thinking the line 311-314 "P0 normalization target" assertion still anchors the claim; the script would still pass but the actual normalization would no longer be verified.

**Required change:**
Either (a) delete lines 307-314 entirely (the `P0_target` computation that's unused and the duplicate assertion), or (b) replace the duplicate assertion with a substantive normalization check, e.g.
```
expect_zero("P0 normalization target", mhat**2 * (N0/D0) * Gamma5_port - gamma_GR).subs(N0, P0_target * D0)
```
or simpler:
```
expect_zero("P0_target satisfies mhat^2 P0 = 54 G c_s^5 / (5 a^5 c^5) at mhat=1",
            (mhat**2 * P0_target / D0).subs(mhat, 1) - 54*G*c_s**5/(5*a**5*c**5))
```
Option (a) is the safer mechanical fix; the rest of the section already covers the substantive content.

**Verification:**
After the edit, the saved output no longer contains a line `P0 normalization target = 0` that is byte-equal to the earlier `P0 - N0/D0 = 0` line, OR the new check exercises a non-trivial substitution involving Gamma5_port or gamma_GR and the print line shows a value reflecting that substitution (still 0 if correct).

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:305-306`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:157-158`

**What's wrong:**
`N2_target = sp.solve(sp.Eq(P2, 0), N2)[0]` and then `expect_zero("P2 under N2_target", P2.subs(N2, N2_target))` is checking that the solver returns a solution. By construction `P2.subs(N2, solve(P2==0, N2)[0])` simplifies to 0 — this verifies the SymPy solver's correctness, not the physics. The same is true for `P4 under N2_target,N4_target` at line 306 and for the Mathematica counterparts at lines 157-158.

The substantive content (i.e. what `N2_target` and `N4_target` actually equal in closed form) is already printed for human inspection at lines 300-303 / Mathematica lines 155-156, and the `pprint` output `N2 = 2*D2*N0/D0` and `N4 = N0*(2*D0*D4 + D2**2)/D0**2` are the actual reported values.

**Why this matters:**
Low severity because (a) the substitution does serve as a cheap sanity check on `solve`, and (b) the formulas are printed for visual inspection. But the "0 = 0" lines in the output transcript could be misread as independent verification that the closed-form solutions are correct, when they are merely solver consistency.

**Required change:**
Either (a) delete the two assertions (lines 305-306 in SymPy; 157-158 in Mathematica), keeping the `pprint`/Print of the closed-form solutions, or (b) replace with a non-tautological check: independently form the closed-form targets and verify they agree with the solver output, e.g.
```
N2_target_closed_form = 2 * D2 * N0 / D0
expect_zero("N2_target closed form", N2_target - N2_target_closed_form)
N4_target_closed_form = N0 * (2 * D0 * D4 + D2**2) / D0**2
expect_zero("N4_target closed form", N4_target - N4_target_closed_form)
```
This tests that the solver returns the documented closed form, not just that substituting its output back gives 0.

**Verification:**
After the edit, either (a) the two "P2 under N2_target" / "P4 under N2_target,N4_target" lines no longer appear in the saved output, or (b) new lines `N2_target closed form = 0` and `N4_target closed form = 0` appear, witnessing independent derivation of the closed forms.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl` (whole file)
- compared to `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`

**What's wrong:**
The Mathematica script is a near line-by-line port of the SymPy script: same variable choreography, same intermediate algebraic steps, same physical premises (including identical comments such as "Stage-003 carry-forward: B_{A0}, B_{A2}, B_{A4} are the stable-BdG Schur sums..."). Both engines:

1. Define the same `groupedParts`/`grouped_parts` function with the same three formulas `(x20+2x21+2x22)/5`, `(2x20-x21-x22)/10`, `(x21-x22)/2` and use it identically.
2. Define `Delta_expr = OmegaU²·OmegaW² - Rmix²`, `S_expr = OmegaU² + OmegaW²`, `Q_expr = gU²·OmegaW² + 2·gU·gW·Rmix + gW²·OmegaU²`, `H_expr = gU² + gW²`, `P_expr = OmegaU²·gW + Rmix·gU` — identical algebraic expressions.
3. Compute `Z_one_port = Series[(Q - H*omega²)/(Delta - S*omega² + omega⁴), {omega, 0, ...}]` and `N_one_port = Series[(P - gW*omega²)²/(Delta - S*omega² + omega⁴)², ...]` — identical series setup.
4. Compute `Lambda2 = omega·D[h2, z]/h2 /. z -> omega·a/cS` then series in omega to obtain `Gamma5_port = a⁵/(27 c_s⁵)` — identical Hankel-function reciprocal-DtN path.

Quoted correspondences:

SymPy lines 126-130:
```
Delta_expr = sp.simplify(OmegaU**2 * OmegaW**2 - Rmix**2)
S_expr = sp.simplify(OmegaU**2 + OmegaW**2)
Q_expr = sp.simplify(gU**2 * OmegaW**2 + 2 * gU * gW * Rmix + gW**2 * OmegaU**2)
H_expr = sp.simplify(gU**2 + gW**2)
P_expr = sp.simplify(OmegaU**2 * gW + Rmix * gU)
```
Mathematica lines 77-81:
```
deltaExpr = omegaU^2*omegaW^2 - rMix^2;
sExpr = omegaU^2 + omegaW^2;
qExpr = gU^2*omegaW^2 + 2*gU*gW*rMix + gW^2*omegaU^2;
hExpr = gU^2 + gW^2;
pExpr = omegaU^2*gW + rMix*gU;
```
SymPy lines 317-321:
```
j2 = (sp.Rational(3, 1) / z**3 - sp.Rational(1, 1) / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
y2 = -((sp.Rational(3, 1) / z**3 - sp.Rational(1, 1) / z) * sp.cos(z) + 3 * sp.sin(z) / z**2)
h2 = sp.simplify(j2 + I * y2)
Lambda2 = sp.simplify(omega * sp.diff(h2, z) / h2).subs(z, omega * a / c_s)
```
Mathematica lines 165-168:
```
j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
h2 = FullSimplify[j2 + I y2, Assumptions -> $Assumptions];
lambda2 = FullSimplify[(omega D[h2, z]/h2) /. z -> omega a/cS, Assumptions -> $Assumptions];
```

The second-engine policy requires the two engines to derive the result independently. Here both follow the same algebraic recipe and verify the same intermediates against the same closed forms. A genuine independent verification would, for example, (a) substitute concrete numerical values (e.g. `OmegaU=2, OmegaW=3, Rmix=1, gU=1, gW=2`) and verify the Z_n, N_n coefficients numerically; (b) verify `Gamma5_port = a^5/(27 c_s^5)` by computing the radiation reaction directly from `h2` via residue/contour or an alternate closed-form path; or (c) compute `P₀, P₂, P₄` by an alternate route (e.g. via Cauchy product of two series rather than `Series[D₀·(N₀+N₂ω²+N₄ω⁴)/Dcons², ...]`).

**Why this matters:**
If both engines share the same algebraic mistake (e.g., a sign convention, a wrong series order, an off-by-one in coefficient extraction), neither will catch it — both will report PASS. The whole point of running two engines is to catch errors via algebraically distinct paths.

**Required change:**
Restructure the Mathematica script to verify the Section II.0 (one-port Z_n, N_n) and Section III.3 (Gamma5_port and the normalization product) checks via an algebraically distinct route, while keeping the same final claims. Concrete suggested changes:

(i) For Section II.0 (mathematica lines 77-91): instead of `Coefficient[Series[(Q-Hω²)/(Δ-Sω²+ω⁴), {ω,0,4}], ω, n]` versus a closed-form, substitute a fixed numerical realization (e.g. `omegaU -> 2, omegaW -> 3, rMix -> 1, gU -> 1, gW -> 2`) into both the rational function `(Q-Hω²)/(Δ-Sω²+ω⁴)` and the closed-form Z_n, evaluate at, say, omega = 1/10 (small enough for the series to dominate) and verify they match to high precision. This breaks the structural correspondence with the SymPy script.

(ii) For Gamma5_port (mathematica lines 165-175): instead of computing `omega·D[h2,z]/h2`, `Series` in omega, then `Coefficient[..., omega, 5]/I`, compute the same coefficient via the small-z expansion of `j₂(z) + i·y₂(z)` directly using the known closed-form expansion of the spherical Bessel functions at small argument, and verify that the resulting 5th-order coefficient in `omega·a/c_s` equals `a⁵/(27 c_s⁵)`.

Either change suffices; both would be best.

**Verification:**
The Mathematica script's lines 77-91 and/or lines 165-175 use a distinct algebraic mechanic from the SymPy version. The saved output shows pass lines whose names indicate the alternate route (e.g. `Z0 numerical at omegaU=2,omegaW=3,...` or `Gamma5_port via Bessel small-z expansion`).

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. The two share: same variable layout, same `groupedParts` definition with identical formulas, same `Delta/S/Q/H/P` expressions, same `Series[..., {omega, 0, n}]` choreography, same `j2 + I·y2` Hankel-function construction, same `omega·D[h2,z]/h2` derivative-ratio path, same Cauchy-style normalization formulas. See F4 for quoted side-by-side excerpts.

## Engine cross-check

Both engines report PASS / EXIT_CODE 0. SymPy lists 47 numbered identity checks; Mathematica lists 44. The discrepancy is benign: SymPy emits four scalar/matrix expand-and-print lines that do not appear as named assertions in the Mathematica transcript (e.g. the four `du2, dP0, dP0/dN0, dP0/dB0` `pprint` entries), and SymPy includes the redundant "P0 normalization target" line discussed in F2 which Mathematica does not have. Every Mathematica assertion's named claim corresponds to a SymPy assertion on the same mathematical identity, and all simplify to 0 in both engines. Numeric/symbolic agreement at the named-identity level is complete.

## Verdict justification

Verdict is **findings**, not stop_cold. The mathematical content of the unit (projector calculus, one-port Z/N coefficients, isotropic prefactor formulas, monotonicity derivatives, Gamma5_port = a⁵/(27 c_s⁵), and the normalization product 54 G c_s⁵/(5 a⁵ c⁵)) is internally consistent and holds up under attack: the projectors are genuinely idempotent and Ggrp-orthogonal; the one-port series coefficients algebraically reduce to the asserted closed forms; the Hankel-function reciprocal-DtN expansion does give the stated 5th-order coefficient; and the normalization arithmetic (2/5 · 27 = 54/5) checks out. Attempted attacks: (a) tried to find a sign or factor-of-two error in the `54 G c_s⁵/(5 a⁵ c⁵)` constant by recomputing `2·27/5` — got 54/5; correct. (b) Tried to find a wrong denominator in the projectors (using 5/20/4 from norms) — they match `||·||_G²` correctly; correct. (c) Tried to break the constant-prefactor branch arithmetic by substituting `N2_target = 2 D2 N0/D0` into the full P4 formula and reducing — got `N₀·(2 D₀ D₄ + D₂²)/D₀²`, matching the printed solver output; correct. (d) Tried to find an order-of-series-too-low issue in Mathematica `Series[..., {omega, 0, 4}]` for the N_one_port (which involves `(... )²/(... )²`) — the order 4 is sufficient because we only extract coefficients up to ω⁴. Holds up. The four findings are real but bounded: three are tautological-check housekeeping (F1, F2, F3) and one is engine independence (F4). None invalidates the unit's claims; they only weaken the evidentiary weight of specific assertion lines.

## Self-test notes

(1) Variable independence: none of my required changes introduce new `sp.diff` or `D[...]` calls over variables not already in the script, so no risk of identically-zero derivative traps. (2) Symmetry/parity: no new integrals proposed; the existing series expansions are not over symmetric unbounded domains. (3) Trivial-case pre-check: for F1's suggested replacement (`grouped_parts(N020+N220, ...) - (grouped_parts(N020,...) + grouped_parts(N220,...))`), substituting concrete `N020=1, N021=2, N022=3, N220=4, N221=5, N222=6` gives both sides equal to `((1+4)+2(2+5)+2(3+6))/5 = (5+14+18)/5 = 37/5` and `(1+2·2+2·3)/5 + (4+2·5+2·6)/5 = 11/5 + 26/5 = 37/5` — non-trivially equal. For F2's suggested replacement `(mhat²·P0_target/D0).subs(mhat,1) - 54·G·c_s⁵/(5·a⁵·c⁵)`, P0_target = solve(mhat²·N0/D0 = 54·G·c_s⁵/(5·a⁵·c⁵), N0)[0] = 54·G·c_s⁵·D0/(5·a⁵·c⁵·mhat²); substituting back gives 0 — passes. For F3's suggested closed-form check `N2_target - 2·D2·N0/D0`, the solver returns exactly that expression (visible in the printed output), so substituting concrete `D2=1, N0=2, D0=3` gives `4/3 - 4/3 = 0` — passes. (4) Path specifications: all file paths in the directive are absolute and verified to exist (see `ls` at audit start).
