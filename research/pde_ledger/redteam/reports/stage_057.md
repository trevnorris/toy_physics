---
unit_id: 057
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md"]
  paper_appendix: present
---

# Audit unit 057 red-team report (v2 — paper-grounded)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_057.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 92 references stage 057)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.txt`

## What the paper claims

Stage card Output line: the **physical lowest-lane support ratio** is
`zeta_0(Pe, eta; kappa) = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` with `y(eta) tan y(eta) = eta` (eq:app-stage057-zeta), plus the zero-bias special case `zeta_0(0, eta; kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2)`. The card states "It is monotone in the constructive Peclet direction".

The notes elaborate six load-bearing deliverables: (1) `x = pi^2/(kappa + pi^2/4)`; (2) `A_K(eta; kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2)`; (3) the full `zeta_0^(Pe+R)` formula above; (4) **three monotonicity statements** — `partial_Pe zeta > 0`, `partial_eta zeta < 0`, `partial_kappa zeta < 0` on the constructive branch, with `0 < y < pi/2`; (5) closure ceiling `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa`; (6) reachability `kappa <= pi^4/[4(4 zeta_req - pi^2)]`. Notes section 7 adds threshold formulas `Omega_req^2`, `y_req^2`, `kappa_req`.

Part appendix row 92: "Physical support ratio `zeta_0(Pe, eta; kappa)`" — consistent with the stage card.

## What the script claims to verify

Both engines verify (1), (2), (3), the algebraic *form* of `partial_kappa zeta` and `partial_y zeta` (without sign), (5), (6), and the threshold formulas (7) — with significant tautological content in the threshold checks (see F2 and F3). Neither engine verifies (4a) `partial_Pe zeta > 0`; the Robin transcendental `y tan y = eta` is never introduced, and the upper bound `y < pi/2` is not declared.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `x = pi^2/(kappa+pi^2/4)` | sympy L45, mathematica L43 (printed only) | match |
| (2) `A_K = (kappa+pi^2/4)/(kappa+y^2)` | sympy L49-52 (A1), mathematica L50 (A9) | match |
| (3) `zeta_0 = Omega_Pe^2*A_K` | sympy L56-59 (A2), mathematica L51-54 (A10) | match |
| (4a) `partial_Pe zeta > 0` | (absent) | **missing** |
| (4b) `partial_eta zeta < 0` (via `dy/deta > 0`) | only `partial_y` algebraic form; no `y(eta)` link | partial |
| (4c) `partial_kappa zeta < 0` (needs `y < pi/2`) | sympy L66-69 / mathematica L60-63 — form only, no sign check, no `y < pi/2` | partial |
| (5) `zeta_max(kappa)` | sympy L78-81 / mathematica L74-77 | match |
| (6) `kappa_max` | sympy L84-89 via `solve` (real); mathematica L70+L78 (kappaMax literal, identity tautological) | partial |
| (7) threshold formulas | `kappa_req` checked via solve (sympy) + round-trip (mathematica); `y_req` only via structural-guarantee substitution | partial |
| `y tan y = eta` defining transcendental | (absent) | missing (but plausibly carry-forward from Stage 054) |

Paper alignment: **partial**.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 49-52 | `simplify(A_K_kappa - (kappa+pi^2/4)/(kappa+y^2)) == 0` | (2) | yes |
| A2 | sympy | 56-59 | `simplify(zeta_phys - Omega_Pe^2*(kappa+pi^2/4)/(kappa+y^2)) == 0` | (3) | yes (follows from A1) |
| A3 | sympy | 66-69 | `simplify(dkappa - Omega_Pe^2*(y^2-pi^2/4)/(kappa+y^2)^2) == 0` | (4c) form only | partial (no sign) |
| A4 | sympy | 70-73 | `simplify(dy + 2*Omega_Pe^2*y*(kappa+pi^2/4)/(kappa+y^2)^2) == 0` | (4b/c) form only | partial |
| A5 | sympy | 78-81 | `zeta_max - limit(Pe->oo, y->0+) zeta_phys == 0` | (5) | yes |
| A6 | sympy | 84-89 | `kappa_max - pi^4/(4(4*zeta_req-pi^2)) == 0` (kappa_max from `sp.solve`) | (6) | yes |
| A7 | sympy | 94-106 | `kappa_req - (Omega_Pe^2*pi^2/4 - zeta_req*y^2)/(zeta_req - Omega_Pe^2) == 0` (kappa_req from `sp.solve`) | (7) kappa_req | yes |
| A8 | sympy | 107-112 | `zeta_req - Omega_Pe^2*(kappa+pi^2/4)/(kappa + y_req_sq) == 0` | (7) y_req | **tautological** — y_req_sq defined to make this identity hold |
| A9 | mathematica | 50 | `aKKappa - (kappa+Pi^2/4)/(kappa+y^2) == 0` | (2) | yes |
| A10 | mathematica | 51-54 | `zetaPhys - omegaPe^2*(kappa+Pi^2/4)/(kappa+y^2) == 0` | (3) | yes |
| A11 | mathematica | 60-63 | `dKappa - omegaPe^2*(y^2-Pi^2/4)/(kappa+y^2)^2 == 0` | (4c) form | partial |
| A12 | mathematica | 64-67 | `dY + 2*omegaPe^2*y*(kappa+Pi^2/4)/(kappa+y^2)^2 == 0` | (4b/c) form | partial |
| A13 | mathematica | 74-77 | `zetaMax - Limit[Limit[zetaPhys, pe->oo], y->0+] == 0` | (5) | yes |
| A14 | mathematica | 70, 78 | `kappaMax - Pi^4/(4(4 zetaReq - Pi^2)) == 0` where kappaMax is the literal `Pi^4/(4(4 zetaReq - Pi^2))` | (6) | **tautological** |
| A15 | mathematica | 79 | `(zetaMax /. kappa -> kappaMax) - zetaReq == 0` | (6) | yes (real round-trip) |
| A16 | mathematica | 88-91 | `zetaReq - (omegaPe^2*(kappa+Pi^2/4)/(kappa+y^2)) /. kappa -> kappaReq == 0` (kappaReq is literal closed form) | (7) | **structurally guaranteed once kappaReq is literal** (see F2) |
| A17 | mathematica | 83, 92-95 | `kappaReq - (omegaPe^2*Pi^2/4 - zetaReq*y^2)/(zetaReq - omegaPe^2) == 0` (kappaReq = same literal) | (7) | **tautological** |
| A18 | mathematica | 82, 96-101 | `zetaReq - omegaPe^2*(kappa+Pi^2/4)/(kappa + yReqSq) == 0` (yReqSq = `(omegaPe^2/zetaReq)(kappa+Pi^2/4) - kappa`) | (7) y_req | **tautological** |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_057.tex:23`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:140-148,168-177,303`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:62-73`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:56-67`

**Subtype:** `script_missing_paper_claim`

**What's wrong:**

Paper card (stage_057.tex, line 23 immediately after the boxed `zeta_0` formula): "It is monotone in the constructive Peclet direction".

Notes section 4 (lines 140-148) enumerates this as deliverable (4a):
> Increasing in transport bias `Pe` — From Stage 39, `dOmega_Pe/dPe > 0` on the constructive branch `Pe >= 0`, so `partial_Pe zeta_0^(Pe+R) > 0.`

Notes section 8 (line 303) restates: "this map is monotone increasing in `Pe` and monotone decreasing in both `eta` and `kappa`."

Script side: neither script computes `sp.diff(zeta_phys, Pe)` / `D[zetaPhys, pe]`, and neither evaluates the sign. Only `partial_kappa` and `partial_y` are touched, and even those are checked only in algebraic form (see F3). The notes flag the Pe-monotonicity as carry-forward from Stage 056 ("Stage 39" in legacy numbering), but the scripts do not reference Stage 056 or include a confirming spot-check.

**Why this matters:**

The Pe-monotonicity is one of three monotonicity claims that together establish that the physical placement map is "completely ordered" (notes line 173-177). Two of three appear in the script (in form, not sign — see F3); the third does not appear at all. The threshold formulas in section 7 (Omega_req^2 threshold, etc.) depend on this ordering to make sense as a half-line in `Pe`.

This is a paper-vs-script asymmetry whose resolution direction is the user's call: add the check to the script, or annotate the paper/notes that Pe-monotonicity is established at Stage 056 and not re-verified here.

**Required change:**

`paper_misalignment` — Codex must NOT mechanically edit. See directive `## Resolve before fix_loop` block.

**Verification:**

After user resolution, either (a) a new `partial_Pe sign` assertion exists in both scripts (sweeping a sample of `Pe > 0` and confirming `D[zeta_phys, Pe] > 0`, or a symbolic refinement proof), or (b) a follow-up directive annotates the paper/notes to mark Pe-monotonicity as Stage 056 carry-forward.

---

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:70,78`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:83,92-95`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:82,96-101`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:93,107-112`

**What's wrong:**

Several "identity" / "defining equation" assertions are algebraic tautologies — the residual is zero by construction, independent of whether the closed-form expression actually solves the underlying equation.

(a) Mathematica `kappa_max identity` (L70, L78): `kappaMax = FullSimplify[Pi^4/(4 (4 zetaReq - Pi^2)), ...]` then asserts `kappaMax - Pi^4/(4 (4 zetaReq - Pi^2)) == 0`. This is `X - X == 0`. The SymPy counterpart (L84) derives `kappa_max` via `sp.solve(sp.Eq(zeta_req, zeta_max), kappa)`, which IS non-trivial; the Mathematica version skips the solve and hard-codes the literal.

(b) Mathematica `kappa_req identity` (L83, L92-95): `kappaReq = FullSimplify[(omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2), ...]` is the literal RHS; asserting `kappaReq - (omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2) == 0` is `X - X == 0`. (SymPy correctly derives kappa_req via `sp.solve`.)

(c) Mathematica/SymPy `y_req defining equation` (L96-101 / L107-112):
- `y_req_sq := (Omega_Pe^2 / zeta_req)(kappa + pi^2/4) - kappa`
- assertion: `zeta_req - Omega_Pe^2 (kappa + pi^2/4) / (kappa + y_req_sq) == 0`

Substituting: `kappa + y_req_sq = (Omega_Pe^2/zeta_req)(kappa + pi^2/4)`. Then `Omega_Pe^2(kappa+pi^2/4) / [(Omega_Pe^2/zeta_req)(kappa+pi^2/4)] = zeta_req`, so residual is `zeta_req - zeta_req = 0` purely by algebraic cancellation — independent of whether y_req² actually solves the defining equation. A wrong formula like `y_req_sq = (Omega_Pe^2/zeta_req)(kappa + pi^2/3) - kappa` would also pass (the wrong `(kappa + pi^2/3)` factor would simply cancel through the division). This is the same structural tautology v1 thought it had fixed.

(Mathematica `kappa_req defining equation` at L88-91 is the same: `zetaReq - omegaPe^2(kappa+Pi^2/4)/(kappa+y^2) /. kappa -> kappaReq`. Because `kappaReq = (omegaPe^2*Pi^2/4 - zetaReq*y^2)/(zetaReq - omegaPe^2)` is the unique linear-in-kappa solution, the substitution gives a residual that vanishes by linear algebra. Strictly this IS structurally guaranteed once kappaReq is the literal closed form, but it requires the user to verify the closed form by inspection — borderline. I do not flag it separately; the kappaReq closed form is genuinely derived by sympy via `sp.solve` at the same point, which provides the cross-check.)

**Why this matters:**

These checks fail to verify the physics they purport to verify. A sign error or factor-of-two slip in any of `kappa_max`, `kappa_req`, or `y_req^2` would pass undetected by both engines, because each engine's "identity" check is just `literal - literal == 0` or its substitution-cancellation cousin. The SymPy script's checks for `kappa_max` and `kappa_req` (via `sp.solve`) are the only real verifications of the closed forms.

**Required change:**

(a) Mathematica `kappa_max` derivation. At line 70 replace
```
kappaMax = FullSimplify[Pi^4/(4 (4 zetaReq - Pi^2)), Assumptions -> zetaReq > Pi^2/4];
```
with
```
kappaMaxSol = Solve[zetaReq == (Pi^2/4) (kappa + Pi^2/4)/kappa, kappa, Reals];
kappaMax = FullSimplify[kappa /. First[kappaMaxSol], Assumptions -> zetaReq > Pi^2/4];
```
Then the existing `kappa_max identity` check at L78 becomes a real solve-vs-closed-form comparison.

(b) Mathematica `kappa_req` derivation. At line 83 replace
```
kappaReq = FullSimplify[(omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2), Assumptions -> $Assumptions];
```
with
```
kappaReqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2), kappa, Reals];
kappaReq = FullSimplify[kappa /. First[kappaReqSol], Assumptions -> $Assumptions];
```
Then `kappa_req identity` at L92-95 becomes a real solve-vs-closed-form comparison.

(c) Replace `y_req defining equation` (both scripts) with a true solve-vs-closed-form check. The non-tautological assertion is: solve the defining equation for `y^2`, then compare to the closed form `y_req_sq`.

SymPy at L107-112 replace
```python
expect_zero(
    "y_req defining equation",
    zeta_req - sp.simplify(
        Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y_req_sq)
    ),
)
```
with
```python
y_req_sq_solved = sp.solve(
    sp.Eq(zeta_req, Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y**2)),
    y**2,
)[0]
expect_zero(
    "y_req identity",
    y_req_sq - y_req_sq_solved,
)
```

Mathematica at L96-101 replace
```
expectZero[
  "y_req defining equation",
  zetaReq - FullSimplify[
    omegaPe^2 (kappa + Pi^2/4)/(kappa + yReqSq),
    Assumptions -> $Assumptions
  ]
];
```
with
```
yReqSqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + ysq), ysq, Reals];
yReqSqSolved = FullSimplify[ysq /. First[yReqSqSol], Assumptions -> $Assumptions];
expectZero["y_req identity", yReqSq - yReqSqSolved];
```

**Verification:**

After the patch, the SymPy transcript shows `y_req identity = 0` (replacing `y_req defining equation`); the Mathematica transcript shows `PASS: y_req identity`, `PASS: kappa_max identity`, `PASS: kappa_req identity` where the kappaMax and kappaReq closed forms are now derived from `Solve`. The verifier confirms that all three `Solve` calls appear in the .wl source and the `sp.solve` for `y_req_sq_solved` appears in the .py source.

---

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:35,62-73`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:35-36,56-67`

**What's wrong:**

Notes deliverable (4c) says `partial_kappa zeta < 0` because `0 < y < pi/2` (lines 159-171). The scripts assert the algebraic form
`partial_kappa zeta - Omega_Pe^2 (y^2 - pi^2/4)/(kappa + y^2)^2 == 0`
but never assert the *sign*. With `sp.symbols("y", positive=True)` (SymPy L35) and `y > 0` only (Mathematica L36), the factor `y^2 - pi^2/4` has indeterminate sign; the bound `y < pi/2` (which follows from `y tan y = eta` with `eta < oo`) is not declared anywhere. A `simplify(sign(partial_kappa zeta))` would not reduce to a definite sign.

Without the `y < pi/2` constraint, the script does not rule out a hypothetical `y > pi/2` branch flipping the monotonicity. The form-only check is insufficient for the load-bearing physical statement.

Note: `partial_y zeta < 0` is implicit from the form `dy = -2 Omega_Pe^2 y (kappa+pi^2/4)/(kappa+y^2)^2` because the prefactor `2 Omega_Pe^2 y (kappa+pi^2/4)/(kappa+y^2)^2 > 0` follows from declared positivity, so deliverable (4b)'s `partial_eta zeta < 0` is verified *modulo* the unverified `dy/deta > 0` (claimed in notes line 152-154 but not exercised in the script). Deliverable (4c)'s sign is the harder gap.

**Why this matters:**

Deliverable (4c) is one of the headline physical results ("a softer baseline support branch (smaller `kappa`) always helps", notes line 171). A script that only verifies the algebraic shape of the derivative but not its sign leaves the physics unverified.

**Required change:**

(a) In SymPy at L35, replace
```python
Pe, kappa, y, zeta_req = sp.symbols("Pe kappa y zeta_req", positive=True, real=True)
```
with
```python
Pe, kappa, zeta_req = sp.symbols("Pe kappa zeta_req", positive=True, real=True)
y = sp.Symbol("y", positive=True, real=True)
```
(no functional change yet — y stays positive). After the existing `partial_y identity` block at L70-73, append:
```python
# Sign check on partial_kappa over 0 < y < pi/2 (from y tan y = eta, eta finite)
for y_val in (sp.pi / 8, sp.pi / 6, sp.pi / 4, sp.pi / 3, sp.Rational(7, 16) * sp.pi):
    val = sp.simplify(dkappa.subs({Pe: 1, kappa: 1, y: y_val}))
    if val >= 0:
        raise AssertionError(f"partial_kappa zeta sign failed at y={y_val}: {val}")
print("partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS")
```

(b) In Mathematica at L35-36, add `y < Pi/2` to `$Assumptions`:
```
$Assumptions =
  Element[{pe, kappa, y, zetaReq, x}, Reals] && pe > 0 && kappa > 0 && y > 0 && y < Pi/2 && zetaReq > 0 && x > 0;
```
After the existing `partial_y identity` block at L64-67, append:
```
Module[{yvals, kv, signOk},
  yvals = {Pi/8, Pi/6, Pi/4, Pi/3, 7 Pi/16};
  signOk = AllTrue[yvals, TrueQ[N[(D[zetaPhys, kappa] /. {pe -> 1, kappa -> 1, y -> #})] < 0] &];
  If[signOk, pass["partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)"],
    fail["partial_kappa zeta sign sweep"]]
];
```

**Verification:**

After the patch, the SymPy transcript prints `partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS`; the Mathematica transcript prints `PASS: partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)`. The Mathematica `$Assumptions` now contains `y < Pi/2`.

---

### F4 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:38-44`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:37-46`

**What's wrong:**

The Mathematica script ports the SymPy script line-for-line: same variable names (snake_case → camelCase: `Omega_Pe`→`omegaPe`, `A_K_x`→`aKX`, `x_sub`→`xSub`, `A_K_kappa`→`aKKappa`, `zeta_phys`→`zetaPhys`, `kappa_max`→`kappaMax`, `Omega_req_sq`→`omegaReqSq`, `y_req_sq`→`yReqSq`, `kappa_req`→`kappaReq`), same algorithmic order, same intermediate definitions in the same order.

SymPy L37-46:
```python
Omega_Pe = sp.simplify(sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)/((4 * Pe**2 + sp.pi**2)*(sp.exp(Pe) - 1)))
A_K_x = sp.simplify(1 / (1 - x / 4 + x * y**2 / sp.pi**2))
x_sub = sp.simplify(sp.pi**2 / (kappa + sp.pi**2 / 4))
A_K_kappa = sp.simplify(A_K_x.subs(x, x_sub))
```
Mathematica L38-44:
```
omegaPe = FullSimplify[Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1)), ...];
aKX = FullSimplify[1/(1 - x/4 + x y^2/Pi^2), ...];
xSub = FullSimplify[Pi^2/(kappa + Pi^2/4), ...];
aKKappa = FullSimplify[aKX /. x -> xSub, ...];
```

The Mathematica script never derives `A_K` from the physical operator quantities the notes describe (notes section 2: `K_W^(eff) = (T_X/L^2)(kappa + pi^2/4)`, `K_phi0^(eff) = (T_X/L^2)(kappa + y^2)`, `A_K = K_W^(eff)/K_phi0^(eff)`). It echoes SymPy's choreography.

**Why this matters:**

A transliteration cannot catch a SymPy bug — both engines walk the same algebra. The agreement between the engines is then partially vacuous; it confirms only that two CAS systems can both simplify the same substitution chain, not that an independent derivation path arrives at the same answer.

This is severity low because the algebra is short and visible; for a 4-step substitution chain, the second-engine policy is less load-bearing than in heavier units. I flag it but do not require a full restructure.

**Required change:**

In the Mathematica script at L38-44, replace the substitution chain with an independent derivation from the physical operator. Define `KW := KX + Pi^2 TX/(4 L^2)` and `Kphi0 := KX + TX y^2/L^2`, form `aKKappa := KW/Kphi0` symbolically, then introduce `kappa := KX L^2 / TX` and simplify to show `aKKappa == (kappa + Pi^2/4)/(kappa + y^2)`. Concretely insert before line 38:
```
Clear[KX, TX, LL];
KW = KX + Pi^2 TX/(4 LL^2);
Kphi0 = KX + TX y^2/LL^2;
aKPhys = FullSimplify[KW/Kphi0, Assumptions -> KX > 0 && TX > 0 && LL > 0 && y > 0];
aKKappaFromPhys = FullSimplify[aKPhys /. KX -> kappa TX/LL^2, Assumptions -> $Assumptions];
expectZero["A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)",
  aKKappaFromPhys - (kappa + Pi^2/4)/(kappa + y^2)];
```
Then the original Mathematica L42-44 (the `aKX`/`xSub`/`aKKappa` chain) may remain as a separate consistency check, but the load-bearing derivation is now from the physical operator.

**Verification:**

After the patch, the Mathematica transcript contains a new `PASS: A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)` line; the SymPy transcript is unchanged. The Mathematica source contains the symbols `KX`, `TX`, `LL` (or similar) that do not appear in the SymPy source.

## Independent-derivation check (Mathematica)

As detailed in F4: the Mathematica script mirrors SymPy's substitution path. The recommended remedy is to add a physical-operator derivation in the Mathematica side without removing the existing check.

## Engine cross-check

Both engines print agreeing closed forms (with surface-syntax differences that are algebraically identical):
- `x(kappa)`: sympy `4*pi**2/(4*kappa + pi**2)`, mathematica `(4*Pi^2)/(4*kappa + Pi^2)` — agree.
- `A_K`: sympy `(kappa + pi**2/4)/(kappa + y**2)`, mathematica `(4*kappa + Pi^2)/(4*(kappa + y^2))` — algebraically identical.
- `zeta_phys`: both show `pi^2 Pe^2 (kappa+pi^2/4)(2 Pe exp(Pe)+pi)^2 / [(exp(Pe)-1)^2 (4 Pe^2+pi^2)^2 (kappa+y^2)]` — agree.
- `partial_kappa zeta`, `partial_y zeta`, `zeta_max`, `kappa_max`: agree.
- All checked residuals reduce to 0 in both engines.

Engine agreement is genuine at the surface level, but qualified by F4 — both engines walk the same algorithm.

## Verdict justification

The scripts non-tautologically verify the algebraic content of paper deliverables (2), (3), (5), the SymPy `kappa_max` (via `sp.solve`), and the SymPy `kappa_req` (via `sp.solve`). Where the scripts fall short:

- F1 (paper_misalignment): the Pe-monotonicity (`partial_Pe zeta > 0`) listed in paper card and notes is not checked in either script and not visibly linked as Stage 056 carry-forward.
- F2 (tautological_check): the Mathematica `kappa_max`/`kappa_req` identity checks compare a literal to itself; the `y_req defining equation` checks in both engines are structurally guaranteed by the algebraic definition of `y_req_sq` (the v1 fix replaced one tautology with another). The kappa_max/kappa_req paths must derive the closed form via `Solve`; the y_req path must derive via `Solve[..., y^2]`.
- F3 (insufficient_verification): `partial_kappa zeta < 0` is asserted only in algebraic form, never in sign; `y < pi/2` is missing from the symbol declarations in both engines.
- F4 (mathematica_transliteration, low severity): the Mathematica script mirrors SymPy's substitution chain rather than deriving from physical operators.

Attacks tried that failed: the `Omega_Pe -> pi/2` limit as `Pe -> oo` is correct (the `2 Pe exp(Pe) + pi` numerator divided by `4 Pe^2 + pi^2 ~ (exp(Pe) - 1)` factor — leading-order `2 Pe exp(Pe)` in numerator and `4 Pe^2 (exp(Pe) - 1) ~ 4 Pe^2 exp(Pe)` in denominator gives `Pe * 2 / (4 Pe) = 1/2`, times `pi` gives `pi/2`). The kappa_max closed form is consistent with the algebra. The `kappa_req` closed form is correct (a linear-in-kappa equation has a unique solution). The `sp.solve` branch picks (`[0]` index) for kappa_max and kappa_req are safe because both equations are linear in kappa (single-rooted).

Verdict: **findings**. No `stop_cold`: F1 needs user direction but doesn't propagate downstream until a direction is chosen; F2/F3/F4 are local script-quality issues that Codex can remediate once F1 is resolved (or independently, since they don't touch the paper side).

## Self-test notes

(1) Variable independence: the proposed F2 SymPy patch uses `sp.solve(..., y**2)[0]` — `zeta_req` depends on `kappa`, `y`, `Pe`; the equation is linear in `y^2`, single-rooted, so `[0]` is safe. (2) Symmetry/parity: not applicable; no unbounded-domain integrals here. (3) Trivial-case pre-check: F3 sign-sweep at `y=pi/4, kappa=1, Pe=1` gives `(y^2 - pi^2/4) = -3 pi^2/16 < 0`; the numerator-squared prefactor is positive; residual is negative — PASS. (4) Path specs: `.py` in `scripts/`, `.wl` in `mathematica/` — correct. (5) Paper round-trip: F2/F3/F4 edits change scripts only; no constants change; no new paper_misalignment introduced. F1 is the only paper-touching finding and is explicitly user-routed.
