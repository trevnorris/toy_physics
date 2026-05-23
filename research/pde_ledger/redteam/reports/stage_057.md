---
unit_id: 057
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 057 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.txt`

## What the script claims to verify

The scripts build the physical lowest-lane closure family ζ_0^(Pe+R)(Pe,κ,y) = Ω_Pe² · A_K(η;κ), starting from a velocity factor Ω_Pe(Pe) = π Pe(2Pe e^Pe + π)/((4Pe²+π²)(e^Pe−1)) and a Robin softening A_K(x,y) = 1/(1 − x/4 + x y²/π²). They verify (i) the substitution x → π²/(κ+π²/4) collapses A_K to (κ+π²/4)/(κ+y²); (ii) closed-form monotonic derivatives ∂_κ ζ and ∂_y ζ; (iii) a closure ceiling ζ_max(κ) = (π²/4)(κ+π²/4)/κ matching lim_{Pe→∞} lim_{y→0+} ζ_phys; (iv) the inverted stiffness ceiling κ_max(ζ_req) = π⁴/(4(4ζ_req−π²)); and (v) closed-form parameter-threshold formulas for Ω_req², y_req², κ_req solving ζ_req = Ω_Pe²(κ+π²/4)/(κ+y²). The Mathematica script additionally checks that κ_req satisfies the defining equation (round-trip) and that ζ_max(κ_max) = ζ_req.

## Assertion inventory

| #  | Script | Line     | Form | Anchored to claim? |
|----|--------|----------|------|---------------------|
| A1 | sympy        | 49-52  | `A_K_kappa - (κ+π²/4)/(κ+y²) == 0` | yes |
| A2 | sympy        | 56-59  | `zeta_phys - Ω²·(κ+π²/4)/(κ+y²) == 0` | partial (follows from A1 by multiplication) |
| A3 | sympy        | 66-69  | `dκ ζ − Ω²(y²−π²/4)/(κ+y²)² == 0` | yes |
| A4 | sympy        | 70-73  | `dy ζ + 2Ω² y (κ+π²/4)/(κ+y²)² == 0` | yes |
| A5 | sympy        | 78-81  | `ζ_max − lim_{Pe→∞, y→0+} ζ_phys == 0` | yes |
| A6 | sympy        | 86-89  | `κ_max − π⁴/(4(4ζ_req−π²)) == 0` | yes |
| A7 | sympy        | 103-106 | `κ_req − closed-form == 0` | yes |
| A8 | sympy        | 107-110 | `y_req_sq − ((Ω²/ζ_req)(κ+π²/4) − κ) == 0` | **no — tautology** |
| B1 | mathematica  | 44      | `aKKappa − (κ+π²/4)/(κ+y²) == 0` | yes |
| B2 | mathematica  | 45-48   | `zetaPhys − Ω²·(κ+π²/4)/(κ+y²) == 0` | partial (follows from B1) |
| B3 | mathematica  | 54-57   | `dκ ζ identity == 0` | yes |
| B4 | mathematica  | 58-61   | `dy ζ identity == 0` | yes |
| B5 | mathematica  | 68-71   | `ζ_max − lim ζ_phys == 0` | yes |
| B6 | mathematica  | 72      | `κ_max identity == 0` | yes |
| B7 | mathematica  | 73      | `ζ_max(κ_max) − ζ_req == 0` (round-trip) | yes |
| B8 | mathematica  | 82-85   | `ζ_req − Ω²(κ+π²/4)/(κ+y²)|_{κ→κ_req} == 0` (round-trip) | yes |
| B9 | mathematica  | 86-89   | `κ_req identity == 0` | yes |
| B10| mathematica  | 90-93   | `y_req_sq − ((Ω²/ζ_req)(κ+π²/4) − κ) == 0` | **no — tautology** |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:92,107-110`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:76,90-93`

**What's wrong:**
The SymPy script defines

```python
y_req_sq = sp.simplify((Omega_Pe**2 / zeta_req) * (kappa + sp.pi**2 / 4) - kappa)
```

then asserts

```python
expect_zero(
    "y_req identity",
    y_req_sq - ((Omega_Pe**2 / zeta_req) * (kappa + sp.pi**2 / 4) - kappa),
)
```

The asserted residual is `y_req_sq − (the literal expression y_req_sq was defined as)`. This is `def x = E; assert x − E == 0` — algebraically guaranteed by construction, independent of any physics. Nothing about the threshold relation ζ_req = Ω_Pe²(κ+π²/4)/(κ+y²) is exercised.

The Mathematica script reproduces the same tautology verbatim (line 76 defines `yReqSq` as that expression; line 90-93 subtracts the same expression).

Contrast with the genuine check the Mathematica script does perform on κ_req (line 82-85, "kappa_req defining equation"), which substitutes the candidate back into the defining equation ζ_req = Ω²(κ+π²/4)/(κ+y²) and verifies the residual vanishes. No such round-trip exists for y_req² in either engine.

**Why this matters:**
y_req² is one of the three parameter-threshold formulas claimed by the unit. The current assertion adds no verification: a wrong derivation of y_req² (e.g. a sign error like `kappa - (Omega_Pe²/zeta_req)(kappa+π²/4)`, or a factor-of-two slip in the (κ+π²/4) term) would still pass the tautological check, because both sides of the comparison would carry the same mistake. The script's claim "exact parameter-threshold formulas" is therefore unverified for the y-threshold branch.

**Required change:**
Replace the tautological assertion with a round-trip check: substitute y² → y_req_sq into the defining equation Ω_Pe²(κ+π²/4)/(κ+y²) and verify the result equals ζ_req. Concretely, in SymPy at the line currently containing the `y_req identity` block (lines 107-110):

```python
expect_zero(
    "y_req defining equation",
    zeta_req - sp.simplify(
        (Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y_req_sq))
    ),
)
```

In Mathematica at lines 90-93, mirror with:

```mathematica
expectZero[
  "y_req defining equation",
  zetaReq - FullSimplify[
    (omegaPe^2 (kappa + Pi^2/4)/(kappa + yReqSq)),
    Assumptions -> $Assumptions
  ]
];
```

**Verification:**
After the patch, the SymPy output should print a new line `y_req defining equation = 0` (replacing `y_req identity = 0`), and Mathematica should print `PASS: y_req defining equation`. The check now genuinely depends on the formula y_req² = (Ω²/ζ_req)(κ+π²/4) − κ: any other expression for y_req² would not reduce the substituted equation to ζ_req.

## Independent-derivation check (Mathematica)

The Mathematica script is largely a structural mirror of the SymPy script — same variable order (Ω_Pe, A_K_x, x_sub, A_K_kappa, zeta_phys), same intermediate names (zetaMax, kappaMax, omegaReqSq, yReqSq, kappaReq), and the same assertion list in the same order. For instance:

SymPy lines 54-59:
```python
zeta_phys = sp.simplify(Omega_Pe**2 * A_K_kappa)
... expect_zero("zeta_phys - Omega_Pe^2*(kappa+pi^2/4)/(kappa+y^2)", zeta_phys - Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y**2))
```

Mathematica lines 39, 45-48:
```mathematica
zetaPhys = FullSimplify[omegaPe^2 aKKappa, ...];
... expectZero["zeta_phys - Omega_Pe^2*(kappa+Pi^2/4)/(kappa+y^2)", zetaPhys - omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2)];
```

However, this is not pure transliteration: Mathematica adds two independent checks absent in SymPy — `zeta_max(kappa_max) - zeta_req` (line 73) and `kappa_req defining equation` (line 82-85). Both are genuine round-trip verifications that constrain the closed forms beyond what SymPy's `sp.solve` autoguarantee provides. The shared tautology on y_req² is a co-bug, not a transliteration artifact (each engine independently defines y_req_sq as the algebraic expression and subtracts it from itself). I do not raise `mathematica_transliteration`; the engines verify the same identity backbone but Mathematica is strictly more aggressive on round-trips.

## Engine cross-check

The two engines print the same closed forms in different surface syntax:

- SymPy `x(kappa) = 4*pi**2/(4*kappa + pi**2)` vs Mathematica `x(kappa) = (4*Pi^2)/(4*kappa + Pi^2)` — agree.
- SymPy `A_K = (kappa + pi**2/4)/(kappa + y**2)` vs Mathematica `A_K = (4*kappa + Pi^2)/(4*(kappa + y^2))` — algebraically identical.
- SymPy `zeta_max(kappa) = pi**2/4 + pi**4/(16*kappa)` vs Mathematica `zeta_max(kappa) = (4*kappa*Pi^2 + Pi^4)/(16*kappa)` — algebraically identical (just unexpanded).
- SymPy `kappa_max = pi**4/(4*(4*zeta_req - pi**2))` vs Mathematica `kappa_max = -1/4*Pi^4/(Pi^2 - 4*zetaReq)` — agree (sign distributed differently).
- All assertions reduce to 0 in both engines.

The Limit::alimv warnings in the Mathematica output are nominal (Mathematica drops the assumption on the limit variable, but the limit is correct: Ω_Pe² → π²/4 as Pe→∞, then ζ_phys → (π²/4)(κ+π²/4)/(κ+y²); taking y→0+ gives ζ_max). Both engines arrive at residual 0. No engine_disagreement.

## Verdict justification

Eight (sympy) and ten (mathematica) assertions on a substantive set of physical-parameter identities reduce to zero in both engines, and they agree across engines. The unit holds up under attack except for one tautology on the y_req² threshold formula: in both scripts the assertion compares y_req_sq to the literal expression it was defined as, exercising no physics. I verified `sp.solve`-based assertions (κ_max, κ_req) by re-deriving the linear-in-κ solutions by hand and matching the printed closed forms. I checked the Pe→∞, y→0+ chain for ζ_max and confirmed Ω_Pe → π/2 in that limit. I checked the symbol assumptions (all positive reals — consistent with the physical setup) and the [0] branch picks (both solves are linear in κ, hence single-rooted — safe). The remaining check on y_req² needs to be a round-trip into the defining equation rather than a self-subtraction; that is the sole finding (F1). Verdict: findings.

## Self-test notes

Mental walk-through of the F1 patch on a trivial case: substitute y² → y_req_sq = (Ω²/ζ_req)(κ+π²/4) − κ into ζ_phys formula. Then κ+y² = κ + (Ω²/ζ_req)(κ+π²/4) − κ = (Ω²/ζ_req)(κ+π²/4). Substituting into Ω²(κ+π²/4)/(κ+y²) gives Ω²(κ+π²/4) / [(Ω²/ζ_req)(κ+π²/4)] = ζ_req — the residual `zeta_req − ζ_req` is identically zero. The check has no derivative variable so the independence-of-variables trap doesn't apply; no unbounded integrals so the parity trap doesn't apply; the trivial-case pre-check confirms residual is structurally zero (not just for specific numerics, but as an algebraic identity that depends on the correct form of y_req_sq). Path specs: SymPy edit lives in `scripts/`, Mathematica edit lives in `mathematica/` — both named explicitly with absolute paths in the directive.
