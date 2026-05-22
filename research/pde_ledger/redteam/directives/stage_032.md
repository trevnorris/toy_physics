---
unit_id: 032
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T17:24:21-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 032

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:169-179`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:134-140`

**Issue:**
In Stage 15.5, both scripts assert `Nprod - Nprod_target == 0` where `Nprod` is constructed by multiplying `s_minus/kappa0_sq` by `P0_from_stage13 = beta0*s_minus/lambda_minus`, and `Nprod_target` is written as `beta0*s_minus**2/(kappa0_sq*lambda_minus)`. These two expressions are *identical* by basic multiplicative associativity — the assertion cannot fail regardless of `s_minus`, `lambda_minus`, `beta0`, or `kappa0_sq` values. The stage's docstring claims to "eliminate the abstract source-map factor" but nothing physical is being checked.

**Required change:**

(a) **SymPy script.** Replace the contents of Stage 15.5 (lines 169-179) with a substantive evaluation. After the existing `banner("STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR")` (keep this banner), apply the natural-D/N substitution `subs_nat` from Stage 15.4 plus a `kappa0_sq -> kappa0**2` map, and verify two non-trivial properties of the resulting product:

Before applying, the relevant block looks like:
```python
banner("STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR")

beta0, kappa0_sq = sp.symbols("beta0 kappa0_sq", positive=True, real=True)
P0_minus = sp.symbols("P0_minus", real=True)

P0_from_stage13 = sp.simplify(beta0 * s_minus / lambda_minus)
Nprod = sp.simplify((s_minus / kappa0_sq) * P0_from_stage13)
Nprod_target = sp.simplify(beta0 * s_minus**2 / (kappa0_sq * lambda_minus))
expect_zero("mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)", Nprod - Nprod_target)
```

Replace with (keep the banner line and the `beta0` declaration; remove the unused `kappa0_sq`, `P0_minus` symbols):
```python
banner("STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR")

beta0 = sp.symbols("beta0", positive=True, real=True)

# Closed-form product under the natural D/N substitution.
P0_minus_nat = sp.simplify((beta0 * s_minus / lambda_minus).subs(subs_nat))
Nprod_nat = sp.simplify((s_minus / kappa0**2).subs(subs_nat) * P0_minus_nat)

# (i) At alpha0 = 0 the elimination must yield beta0 * kappa0^2 / A.
Nprod_alpha0 = sp.simplify(Nprod_nat.subs(alpha0, 0))
expect_zero(
    "Nprod(alpha=0) - beta0 * kappa0^2 / A",
    Nprod_alpha0 - beta0 * kappa0**2 / A,
)

# (ii) As alpha0 -> oo the product vanishes (lambda_minus diverges, s_minus is finite).
Nprod_inf = sp.limit(Nprod_nat, alpha0, sp.oo)
expect_zero("limit_{alpha->oo} Nprod_nat", Nprod_inf)

print("All Stage 15 checks passed.")
```

Notes for Codex:
- `subs_nat`, `s_minus`, `lambda_minus`, `kappa0`, `alpha0`, `A` are already defined earlier in the file (Stages 15.1, 15.4). Do not redeclare them.
- The fresh symbol `kappa0_sq` (current line 171) and the unused `P0_minus` symbol (current line 172) should be removed.
- Keep the final `print("All Stage 15 checks passed.")` line.

(b) **Mathematica script.** Replace the contents of Stage 15.5 (lines 134-138) with the analogous interior checks. The current block:
```mathematica
banner["STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR"];
p0Minus = FullSimplify[beta0*sMinus/lamMinus, Assumptions -> $Assumptions];
nProd = FullSimplify[(sMinus/kappa0^2)*p0Minus, Assumptions -> $Assumptions];
nProdTarget = FullSimplify[beta0*sMinus^2/(kappa0^2*lamMinus), Assumptions -> $Assumptions];
expectZero["mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)", nProd - nProdTarget];
```

Replace with:
```mathematica
banner["STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR"];
p0MinusNat = FullSimplify[(beta0*sMinus/lamMinus) /. subsNat, Assumptions -> $Assumptions];
nProdNat = FullSimplify[((sMinus/kappa0^2) /. subsNat)*p0MinusNat, Assumptions -> $Assumptions];

(* (i) alpha0 = 0 must give beta0 * kappa0^2 / a. *)
nProdAt0 = FullSimplify[nProdNat /. alpha0 -> 0, Assumptions -> $Assumptions];
expectZero["Nprod(alpha=0) - beta0 * kappa0^2 / a", nProdAt0 - beta0*kappa0^2/a];

(* (ii) alpha0 -> Infinity must give zero. *)
nProdInf = Limit[nProdNat, alpha0 -> Infinity, Assumptions -> $Assumptions];
expectZero["limit_{alpha->oo} Nprod_nat", nProdInf];
```

Notes for Codex:
- `subsNat`, `sMinus`, `lamMinus`, `kappa0`, `alpha0`, `a` are already defined earlier in the `.wl`. Do not redeclare.
- The `Limit::alimv` warning that currently appears for the alpha0→∞ limit is acceptable provided `expectZero` reports PASS.
- Keep the final `Print["All Stage 15 checks passed."]; Exit[0];` block.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 032` and `redteam exec-mathematica 032`. Both scripts must exit 0, and the captured outputs must contain the new assertion names `"Nprod(alpha=0) - beta0 * kappa0^2 / A"` (sympy) / `"Nprod(alpha=0) - beta0 * kappa0^2 / a"` (mathematica) and `"limit_{alpha->oo} Nprod_nat"`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Replaced the tautological Stage 15.5 product check with natural-branch alpha0 endpoint and infinite-limit checks.
- deviation: none

---

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:145-167`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:117-132`

**Issue:**
Stage 15.4 verifies the closed form `mhat_-^2 = s_minus/kappa0^2` only at the two endpoints `alpha0 = 0` (trivially `1`) and `alpha0 -> oo` (trivially `11/9`). The interior of the parameter space is not exercised, so a structural error in the cross term `(DK + alpha0*delta_kappa)*delta_kappa + 4*alpha0*Kprod` would not be caught. Additionally, the identity `delta_kappa^2 + 4*Kprod == sigma^2` (which holds for `sigma = kappa0^2+kappa1^2`, `delta_kappa = kappa0^2-kappa1^2`, `Kprod = kappa0^2 kappa1^2` and which is what makes Stage 15.4's limits clean) is asserted nowhere.

**Required change:**

(a) **SymPy script.** Add two new `expect_zero` calls between current lines 166 and 167 (i.e., after the `alpha=0` check and before the `alpha->oo` check). Specifically:

Insert immediately after current line 166 (the `mhat_-^2(alpha=0) - 1` assertion):
```python
# Non-trivial identity on the natural-D/N kappa products.
expect_zero(
    "delta_kappa^2 + 4*Kprod - sigma^2 (natural)",
    (delta_kappa**2 + 4 * Kprod - sigma_sym**2).subs(subs_nat),
)

# Interior consistency: s_minus_nat at alpha0 = 1, DK = 1 equals the
# closed form obtained using the natural-D/N identity above.
R_nat = sp.sqrt(DK**2 + 2 * alpha0 * DK * delta_kappa + alpha0**2 * sigma_sym**2)
s_minus_nat_simplified = sp.simplify(
    (sigma_sym + (DK * delta_kappa + alpha0 * sigma_sym**2) / R_nat) / 2
)
expect_zero(
    "s_minus_nat - s_minus_nat_simplified (interior identity)",
    sp.simplify((s_minus_nat - s_minus_nat_simplified.subs(subs_nat))),
)
expect_zero(
    "s_minus_nat at (alpha0=1, DK=1) interior point",
    sp.simplify(
        s_minus_nat.subs({alpha0: sp.Integer(1), DK: sp.Integer(1)})
        - s_minus_nat_simplified.subs({alpha0: sp.Integer(1), DK: sp.Integer(1)}).subs(subs_nat)
    ),
)
```

Notes for Codex:
- `delta_kappa`, `Kprod`, `sigma_sym`, `subs_nat`, `DK`, `alpha0`, `s_minus_nat` are already defined in Stage 15.4. Do not redeclare.
- The two `s_minus_nat_simplified` checks together exercise (i) the identity holds at the symbolic level under subs_nat, and (ii) it specifically holds at a concrete interior point. Both add interior coverage that the existing two-endpoint checks omit.

(b) **Mathematica script.** Insert the analogous block between current lines 131 and 132 (after the `mhat_-^2(alpha=0)` check and before the `Limit[mhatSq, alpha0 -> Infinity]` check):
```mathematica
(* Non-trivial identity on the natural-D/N kappa products. *)
expectZero[
  "delta_kappa^2 + 4*Kprod - sigma^2 (natural)",
  (deltaKappa^2 + 4*kappaProd - sigmaSym^2) /. subsNat
];

(* Interior consistency: derive s_minus_nat via the simplified R form. *)
rNat = Sqrt[dK^2 + 2*alpha0*dK*deltaKappa + alpha0^2*sigmaSym^2];
sMinusNatSimplified = FullSimplify[
  (sigmaSym + (dK*deltaKappa + alpha0*sigmaSym^2)/rNat)/2,
  Assumptions -> $Assumptions
];
expectZero[
  "s_minus_nat - s_minus_nat_simplified (interior identity)",
  FullSimplify[sMinusNat - (sMinusNatSimplified /. subsNat), Assumptions -> $Assumptions]
];
expectZero[
  "s_minus_nat at (alpha0=1, dK=1) interior point",
  FullSimplify[
    (sMinusNat /. {alpha0 -> 1, dK -> 1})
      - ((sMinusNatSimplified /. {alpha0 -> 1, dK -> 1}) /. subsNat),
    Assumptions -> $Assumptions
  ]
];
```

Notes for Codex:
- `deltaKappa`, `kappaProd`, `sigmaSym`, `subsNat`, `dK`, `alpha0`, `sMinusNat` are already defined. Do not redeclare.
- The `sMinusNat` symbol referenced is the existing assignment at line 127 (`sMinusNat = FullSimplify[sMinus /. subsNat, ...]`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 032` and `redteam exec-mathematica 032`. Both scripts must exit 0. The captured outputs must contain the new assertion names:
- `"delta_kappa^2 + 4*Kprod - sigma^2 (natural)"`
- `"s_minus_nat - s_minus_nat_simplified (interior identity)"`
- `"s_minus_nat at (alpha0=1, DK=1) interior point"` (sympy) / `"s_minus_nat at (alpha0=1, dK=1) interior point"` (mathematica)

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Added the natural kappa-product identity check and symbolic plus concrete interior checks for Stage 15.4.
- deviation: none

---

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:87-115` (Stage 15.3)

**Issue:**
The Mathematica script is, throughout, a line-by-line transliteration of the SymPy script (same variable choreography, same matrix layouts, same operand ordering in Stage 15.4 formulas). The most concentrated block of copying is Stage 15.3, where the Schur complement is computed via `Bmat.Inverse[kInt].Transpose[bMat]` with the same 4x4 matrix layout as SymPy's `Bmat * Kint.inv() * Bmat.T`. Replace the central Schur computation with an independent derivation via `LinearSolve`. The other stages should remain (rewriting them all is too large a refactor; this single block is the most damaging instance because of the matrix layout copying).

**Required change:**

After the current line 102 (the `bMat = {...}` definition), and before the current line 103 (`sigmaMat = FullSimplify[bMat.Inverse[kInt].Transpose[bMat], ...]`), add an independent-derivation block, then keep both computations and assert their equality. Specifically:

Replace the block at lines 91-115 with the following (preserving the banner and the existing variable definitions for `v`, `i2`, `kInt`, `bMat`):

```mathematica
banner["STAGE 15.3 — EXACT SCHUR-COMPLEMENT DECOMPOSITION"];
Clear[d0, d1, aPhi, aU, aW];
$Assumptions = Element[{d0, d1, aPhi, aU, aW, gU, gB, gW, gR}, Reals] && aU != 0 && aPhi != 0;

v = {{kappa0}, {kappa1}};
i2 = IdentityMatrix[2];
kInt = {
  {aU, 0, 0, -gR*kappa0},
  {0, aU, 0, -gR*kappa1},
  {0, 0, aPhi, 0},
  {-gR*kappa0, -gR*kappa1, 0, aW}
};
bMat = {{gU, 0, gB*kappa0, gW*kappa0}, {0, gU, gB*kappa1, gW*kappa1}};

(* Independent derivation: do not invert kInt directly.
   Solve kInt . y == Transpose[bMat] . z for y in terms of an arbitrary external z = {z1, z2}.
   Then Bmat . y is, by construction, Bmat . kInt^{-1} . Transpose[bMat] . z, i.e.
   the Schur image of z. We then read off sigmaMat from the linear coefficients of z1, z2. *)
zVec = {z1, z2};
rhs = Transpose[bMat] . zVec;
ySol = LinearSolve[kInt, rhs];
schurImage = FullSimplify[bMat . ySol, Assumptions -> $Assumptions];

sigmaMatViaSolve = FullSimplify[
  {{Coefficient[schurImage[[1]], z1], Coefficient[schurImage[[1]], z2]},
   {Coefficient[schurImage[[2]], z1], Coefficient[schurImage[[2]], z2]}},
  Assumptions -> $Assumptions
];

(* Original direct-inversion route, kept for cross-checking the two methods. *)
sigmaMat = FullSimplify[bMat.Inverse[kInt].Transpose[bMat], Assumptions -> $Assumptions];
expectZero["Schur via Inverse vs LinearSolve", sigmaMat - sigmaMatViaSolve];

sigma = FullSimplify[(Transpose[v].v)[[1, 1]], Assumptions -> $Assumptions];
xi = FullSimplify[gU^2/aU, Assumptions -> $Assumptions];
alphaCoeff = FullSimplify[
  gB^2/aPhi + (aU*gW + gR*gU)^2/(aU*(aU*aW - gR^2*sigma)),
  Assumptions -> $Assumptions
];
sigmaTarget = FullSimplify[xi*i2 + alphaCoeff*(v.Transpose[v]), Assumptions -> $Assumptions];

Print["Sigma = ", fmt[sigmaMat]];
Print["Xi = ", fmt[xi]];
Print["alpha = ", fmt[alphaCoeff]];
expectZero["Sigma - [Xi I + alpha vv^T]", sigmaMat - sigmaTarget];
expectZero["sigmaMatViaSolve - [Xi I + alpha vv^T]", sigmaMatViaSolve - sigmaTarget];
```

Notes for Codex:
- The `LinearSolve` route does not invert `kInt` as a 4x4 matrix; it solves for the four components of `y` given `z`. The intermediate variable choreography (introducing `zVec`, `rhs`, `ySol`, `schurImage`) is structurally different from the SymPy script's flow.
- The two routes (`sigmaMatViaSolve` vs `sigmaMat`) are then cross-checked, AND each is independently verified against the rank-1+identity target. This produces three substantive checks where there was one.
- Do NOT modify Stages 15.1, 15.2, 15.4 in this finding's edit — the transliteration is most damaging at 15.3; the other stages are addressed indirectly by F3's interior derivation in Stage 15.4 and F1's restructured Stage 15.5.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 032`. The script must exit 0, and the captured Mathematica output must contain three Stage-15.3 PASS lines for the assertion names `"Schur via Inverse vs LinearSolve"`, `"Sigma - [Xi I + alpha vv^T]"`, and `"sigmaMatViaSolve - [Xi I + alpha vv^T]"`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Added a LinearSolve-derived Schur image and cross-checked it against both the inverse route and the rank-one target.
- deviation: none
