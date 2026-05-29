---
unit_id: 032
batch: II.1
created_at: 2026-05-25T22:39:50-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-25T23:47:12-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 032

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification (add independent (v.e_-)^2 check)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:163-191` (insertion point: after current `mhat_sq` definition at line 163, before the `expect_zero("mhat_-^2(alpha=0) - 1", ...)` at line 166)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:148-177` (insertion point: after current `mhatSq` definition at line 148, before the `expectZero["mhat_-^2(alpha=0) - 1", ...]` at line 151)

**Issue:**
The paper's headline claim `mhat_- = (v . e_-)/kappa_0` (paper card lines 36-41) is verified only by the rename `mhat_sq := s_minus_nat / kappa_0^2` plus two endpoint checks. The selected eigenvector `e_-` is never constructed and `(v . e_-)^2` is never compared to the imported closed-form `s_-`. A bug in the upstream `s_-` formula would not be caught here.

**Required change:**

### SymPy (insert between lines 165 and 166)

Insert a block that:

1. Builds the loaded 2x2 wall operator in symbolic form using the same `(A, DK, alpha0)` symbols already in scope plus the `kappa0, kappa1` from Stage 15.1:
   ```python
   # Independent (v.e_-)^2 check: construct the loaded 2x2 wall operator
   # M(alpha0) = diag(A, A+DK) - alpha0 * v v^T and verify that
   # (v . e_-)^2 equals the closed-form s_minus_nat above.
   v2 = sp.Matrix([kappa0, kappa1])
   M_loaded = sp.diag(A, A + DK) - alpha0 * (v2 * v2.T)
   eigendata = M_loaded.eigenvects()
   # Sort by eigenvalue: pick the lower (lambda_-) branch.
   # eigenvects returns [(eigenvalue, multiplicity, [eigenvectors]), ...]
   # Identify the lower one by symbolic comparison at a probe point (alpha0=1, A=1, DK=1).
   probe = {alpha0: sp.Integer(1), A: sp.Integer(1), DK: sp.Integer(1)}
   ev_sorted = sorted(
       eigendata,
       key=lambda triple: float(sp.simplify(triple[0].subs(probe)))
   )
   lam_minus_sym, _, evecs_minus = ev_sorted[0]
   e_minus_raw = evecs_minus[0]
   # Normalize e_- to unit length.
   norm_sq = sp.simplify((e_minus_raw.T * e_minus_raw)[0])
   e_minus = sp.simplify(e_minus_raw / sp.sqrt(norm_sq))
   s_check = sp.simplify(((v2.T * e_minus)[0])**2)
   expect_zero(
       "s_check - s_minus_nat (independent (v.e_-)^2 construction)",
       sp.simplify(s_check - s_minus_nat),
   )
   expect_zero(
       "lam_minus_sym - lambda_minus (independent eigenvalue construction)",
       sp.simplify(lam_minus_sym - lambda_minus),
   )
   ```

   Notes for the implementer:
   - `A`, `DK`, `alpha0` are already declared at sympy line 147 with `positive=True, real=True`.
   - `kappa0`, `kappa1` are already in scope from Stage 15.1.
   - `lambda_minus` and `s_minus_nat` already exist (lines 154 and 162).
   - The `subs_nat` substitution is needed on `s_check`/`lam_minus_sym` only if they involve `sigma_sym`, `delta_kappa`, or `Kprod`; with this construction they do not (they are expressed directly in `kappa0`, `kappa1`, `A`, `DK`, `alpha0`). The comparison target `s_minus_nat` already has `subs_nat` baked in; `lambda_minus` does not, so wrap the second `expect_zero`'s difference with `.subs(subs_nat)` if SymPy struggles to verify equivalence symbolically.

   Adjusted second assertion if needed:
   ```python
   expect_zero(
       "lam_minus_sym - lambda_minus (independent eigenvalue construction)",
       sp.simplify(lam_minus_sym - lambda_minus.subs(subs_nat)),
   )
   ```

2. If SymPy's symbolic eigenvector simplification stalls (this is possible for the radical-laden 2x2), fall back to a numeric check at three probe points:
   ```python
   probe_points = [
       {alpha0: sp.Rational(1, 4), A: sp.Integer(1), DK: sp.Integer(1)},
       {alpha0: sp.Rational(1, 1), A: sp.Integer(1), DK: sp.Integer(2)},
       {alpha0: sp.Rational(4, 1), A: sp.Integer(2), DK: sp.Integer(1)},
   ]
   for i, pt in enumerate(probe_points):
       Mnum = M_loaded.subs(pt)
       eigs_num = Mnum.eigenvects()
       eigs_num_sorted = sorted(eigs_num, key=lambda t: float(sp.N(t[0])))
       _, _, evec_list = eigs_num_sorted[0]
       e_n = evec_list[0]
       e_n = e_n / sp.sqrt((e_n.T * e_n)[0])
       s_num = sp.simplify(((v2.T * e_n)[0])**2)
       s_target = sp.simplify(s_minus_nat.subs(pt))
       expect_zero(
           f"s_check_numeric_pt{i} - s_minus_nat_at_pt{i}",
           sp.simplify(s_num - s_target),
       )
   ```
   Prefer the symbolic version if it terminates; use the numeric version only as a backup if SymPy hangs or produces a non-zero residual that does not reduce.

### Mathematica (insert between lines 150 and 151)

Insert an analogous independent construction using `Eigensystem`, written in Mathematica idiom rather than mirroring the SymPy syntax:

```mathematica
(* Independent (v.e_-)^2 check: construct the loaded 2x2 wall operator
   M(alpha0) = DiagonalMatrix[{a, a+dK}] - alpha0 * v.Transpose[v]
   and verify that (v.e_-)^2 equals the closed-form sMinusNat. *)
v2 = {{kappa0}, {kappa1}};
mLoaded = DiagonalMatrix[{a, a + dK}] - alpha0 * (v2 . Transpose[v2]);
{eigVals, eigVecs} = Eigensystem[mLoaded];
(* Pick the lower-eigenvalue index by evaluating at a probe point. *)
probeRule = {alpha0 -> 1, a -> 1, dK -> 1};
eigValsProbe = eigVals /. probeRule;
lowerIdx = First[Ordering[N[eigValsProbe]]];
lamMinusIndep = FullSimplify[eigVals[[lowerIdx]], Assumptions -> $Assumptions];
eMinusRaw = eigVecs[[lowerIdx]];
eMinusNormSq = FullSimplify[eMinusRaw . eMinusRaw, Assumptions -> $Assumptions];
eMinus = FullSimplify[eMinusRaw/Sqrt[eMinusNormSq], Assumptions -> $Assumptions];
sCheck = FullSimplify[(Flatten[Transpose[v2]] . eMinus)^2, Assumptions -> $Assumptions];
expectZero[
  "s_check - s_minus_nat (independent (v.e_-)^2 construction)",
  FullSimplify[sCheck - sMinusNat, Assumptions -> $Assumptions]
];
expectZero[
  "lam_minus_indep - lamMinus (independent eigenvalue construction)",
  FullSimplify[lamMinusIndep - (lamMinus /. subsNat), Assumptions -> $Assumptions]
];
```

Notes for the implementer:
- `a`, `dK`, `alpha0`, `kappa0`, `kappa1`, `sMinusNat`, `lamMinus`, `subsNat` are already in scope (defined at .wl lines 139, 39-40, 144, 147, 146).
- If `FullSimplify` struggles, fall back to the numeric-probe strategy (mirror the SymPy fallback above using three probe points and verify residuals via `N[FullSimplify[...]] == 0` or a `Chop` tolerance).

**Claim manifest** (NEW checks the script must independently verify after this edit):

- M1: For the loaded 2x2 wall operator `M(alpha0) = diag(A, A+DK) - alpha0 * v v^T` with `v = (kappa_0, kappa_1)^T`, kappa_0 = 2 sqrt(2)/pi, kappa_1 = -4/(3 pi), the squared projection `(v . e_-)^2` of v onto the unit lower eigenvector `e_-` equals the closed-form expression `s_minus_nat` already in the script (with `sigma -> kappa_0^2 + kappa_1^2`, `delta_kappa -> kappa_0^2 - kappa_1^2`, `Kprod -> kappa_0^2 kappa_1^2`).
- M2: The lower eigenvalue from `Eigensystem`/`eigenvects` agrees with the closed-form `lambda_minus = (2A + DK - alpha_0 sigma - R)/2` under the same natural substitution.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 032` and `redteam exec-mathematica 032` and confirm:
- The new assertions `s_check - s_minus_nat (independent (v.e_-)^2 construction)` and `lam_minus_sym - lambda_minus (independent eigenvalue construction)` appear in both transcripts (or their numeric-probe variants).
- Both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Added independent lower-eigenvector projection and eigenvalue checks for the loaded 2x2 wall operator.
- deviation: none

## F2 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:143-176` (Stage 15.4 block)

**Issue:**
The Mathematica script's Stage 15.4 (and 15.5) is a transliteration of the SymPy script's Stage 15.4 (and 15.5): same `R` construction, same `lambda_minus`/`s_minus` formulas, same interior-identity scaffolding. This violates the second-engine policy for the headline-claim section of the unit.

**Required change:**

F1's prescribed `Eigensystem`-based block (inserted between .wl lines 150 and 151) is also the F2 fix for Stage 15.4 — it gives the .wl an independent algebraic path for the lower-eigenvalue branch that does not exist in the .py.

For Stage 15.5, also have the .wl compute `nProdNat` via a second route using the eigenvector-based `sCheck`:
```mathematica
nProdIndep = FullSimplify[
  (sCheck/kappa0^2) * (beta0 * sCheck / lamMinusIndep),
  Assumptions -> $Assumptions
];
expectZero[
  "nProdNat - nProdIndep (independent eigenvector path)",
  FullSimplify[nProdNat - nProdIndep, Assumptions -> $Assumptions]
];
```
Insert this in the Stage 15.5 block immediately after the existing `nProdNat = ...` assignment at .wl line 181. This gives the elimination check an independent derivation path in the .wl that the .py does not have.

If F1's symbolic `Eigensystem` route requires the numeric-probe fallback, the same numeric probes can be used here: compute `nProdNat` and `nProdIndep` at the three probe points and assert their difference is 0 numerically.

**Verification command:**
After Codex applies, the .wl transcript must contain at least one `Eigensystem` call (does not appear in the .py) AND must contain `nProdNat - nProdIndep` (or its numeric-probe variant) assertion passing. Both scripts exit 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Added the Mathematica Eigensystem path and an eigenvector-derived Stage 15.5 product cross-check.
- deviation: none

## F3 — insufficient_verification (drop or relabel algebraic-identity scaffolding)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:168-190` (Stage 15.4 "interior consistency" block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:153-176` (parallel block)

**Issue:**
The three assertions `delta_kappa^2 + 4 Kprod - sigma^2 = 0`, `s_minus_nat - s_minus_nat_simplified = 0`, and `s_minus_nat at (alpha0=1, DK=1) interior point = 0` are algebraic-rearrangement consistency checks; the comment labels them "interior consistency" but they test no interior physics. Once F1 lands, they are redundant.

**Required change:**

Apply F1 first. Then, if F1's check passes cleanly, delete the three assertions:

### SymPy
Delete lines 169-190 inclusive (the comment "Non-trivial identity on the natural-D/N kappa products." through and including the closing `)` of the `s_minus_nat at (alpha0=1, DK=1) interior point` block). Keep line 191 (`limit_{alpha->oo} mhat_-^2 - 11/9`).

Before deletion, lines 168-191 read:
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
expect_zero("limit_{alpha->oo} mhat_-^2 - 11/9", sp.simplify(sp.limit(mhat_sq, alpha0, sp.oo) - sp.Rational(11, 9)))
```

After deletion, this region should read:
```python
expect_zero("limit_{alpha->oo} mhat_-^2 - 11/9", sp.simplify(sp.limit(mhat_sq, alpha0, sp.oo) - sp.Rational(11, 9)))
```

### Mathematica
Delete .wl lines 153-176 inclusive (the `(* Non-trivial identity ... *)` comment block through and including the closing `];` of the `s_minus_nat at (alpha0=1, dK=1) interior point` `expectZero`). Keep line 177 (`expectZero["limit_{alpha->oo} mhat_-^2 - 11/9", ...]`).

If F1's symbolic route required the numeric-probe fallback (i.e., SymPy/Mathematica could not prove the symbolic equivalence cleanly), keep the algebraic-identity scaffolding (do not delete) and only update the comment to read `# Algebraic re-arrangement consistency: two forms of s_minus_nat agree because (a-b)^2 + 4ab = (a+b)^2.` — that documents what the block actually does without claiming it tests interior physics.

**Verification command:**
After Codex applies, re-run `redteam exec-sympy 032` and `redteam exec-mathematica 032`. The four removed assertions should no longer appear in the transcripts (if deleted), or the comment text should match (if relabeled). Both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Deleted the redundant natural-kappa identity and s_minus_nat algebraic rearrangement assertions from Stage 15.4.
- deviation: none
