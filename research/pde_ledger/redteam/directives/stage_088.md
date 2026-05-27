---
unit_id: 088
batch: III.5
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
verification_status: scripts_pass_pending_verifier
needs_user_resolution: false
---

# Codex directive — unit 088

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT run python or mathematica. Only edit files. Do NOT touch paper.tex, notes/, or any prose documents.

The numeric paper values that must remain untouched in the asserted outputs: `c0 = 3/4`, `c1 = 1/4`, `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pi_tr = (4/3) C_mix`. Any new derivation path must arrive at these same numbers.

## F1 — tautological_check (Pi_tr identity must invoke rho_alpha, not literal)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:81-87`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:58-67`

**Issue:**
Both scripts hardcode `Pi_tr = (4/3) C_mix` and then "verify" `Pi_tr/C_mix - 4/3 == 0`, which is identically zero with no physics in play. The paper / notes specify that this identity comes from Stage 68's product identity `Pi_tr / C_mix = alpha_req / alpha_mix = rho_alpha`. The check must exercise that chain (`rho_min` derived from `c0=3/4`, then `Pi_tr = rho_min * C_mix`, then compare to `(4/3) C_mix`) so that a wrong upstream `rho_min` propagates to a nonzero residual.

**Required change (sympy, replace lines 81-87):**

Before:
```python
# Product-language regime
Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)
expect_zero("Pi_tr/C_mix - 4/3", sp.simplify((sp.Rational(4,3)*C_mix)/C_mix - sp.Rational(4,3)))

print("\nRegime classification:")
print("  Pi_tr = (4/3) C_mix")
print("  therefore C_mix < Pi_tr < 2 C_mix")
print("  and zeta_req = 1/3 < 1, so the symmetric lowest twin suffices.")
```

After:
```python
# Product-language regime: Pi_tr = rho_alpha * C_mix is the Stage 68 identity.
# Substitute rho_min (derived above from c0 = 3/4) and confirm Pi_tr = (4/3) C_mix.
C_mix = sp.symbols("C_mix", positive=True, real=True)
Pi_tr_from_rho = rho_min * C_mix
expect_zero("Pi_tr_from_rho - (4/3) C_mix", Pi_tr_from_rho - sp.Rational(4,3)*C_mix)

# Regime gate: with rho_min = 4/3, both bounds are exercised as Rational inequalities.
assert rho_min > 1, "rho_min must exceed 1 for the mixed-baseline-not-enough regime"
assert rho_min < 2, "rho_min must lie below 2 to stay in the symmetric-lowest-twin regime"

print("\nRegime classification:")
print("  Pi_tr = rho_min * C_mix = (4/3) C_mix")
print("  therefore C_mix < Pi_tr < 2 C_mix")
print("  and zeta_req = 1/3 < 1, so the symmetric lowest twin suffices.")
```

**Required change (mathematica, replace lines 58-67):**

Before:
```
rhoMin = FullSimplify[rhoFromC0 /. c0 -> 3/4];
zetaMin = FullSimplify[zetaFromC /. {c0 -> 3/4, c1 -> 1/4}];
piMin = FullSimplify[(4/3)*cMix];

Print["rho_alpha(minimal isotropic module) = ", fmt[rhoMin]];
Print["zeta_req(minimal isotropic module) = ", fmt[zetaMin]];
expectZero["rho_min - 4/3", rhoMin - 4/3];
expectZero["zeta_min - 1/3", zetaMin - 1/3];
expectZero["Pi_tr/C_mix - 4/3", piMin/cMix - 4/3];
expectTrue["C_mix < Pi_tr < 2 C_mix", cMix < piMin < 2*cMix];
```

After:
```
rhoMin = FullSimplify[rhoFromC0 /. c0 -> 3/4];
zetaMin = FullSimplify[zetaFromC /. {c0 -> 3/4, c1 -> 1/4}];
(* Stage 68 product identity: Pi_tr = rho_alpha * C_mix; substitute rhoMin. *)
piFromRho = FullSimplify[rhoMin*cMix];

Print["rho_alpha(minimal isotropic module) = ", fmt[rhoMin]];
Print["zeta_req(minimal isotropic module) = ", fmt[zetaMin]];
expectZero["rho_min - 4/3", rhoMin - 4/3];
expectZero["zeta_min - 1/3", zetaMin - 1/3];
expectZero["Pi_tr_from_rho - (4/3) C_mix", piFromRho - (4/3)*cMix];
expectTrue["1 < rho_min < 2 (symmetric-lowest-twin regime)", 1 < rhoMin < 2];
```

**Verification:**
After edits, residual for `Pi_tr_from_rho - (4/3) C_mix` is `(4/3)*C_mix - (4/3)*C_mix = 0` ✓; substituting `c0 = 1/2` upstream (mutation check) would yield `rho_min = 2`, residual `2 C_mix - (4/3) C_mix = (2/3) C_mix ≠ 0` — so the check is non-tautological.

## F2 — tautological_check (contact-plus-pole identity derived rather than restated)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:46-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:46-56`

**Issue:**
`c0_from_rho := 1/rho_alpha` and `c1_from_rho := (rho_alpha-1)/rho_alpha` are introduced as definitions, then asserted equal to those very pieces of `Y_rho` — tautological. To make the check exercise something, the constants must be recovered from `Y_rho` (or the paper's `Y_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)`) by an independent procedure (partial fractions / static-limit / pole residue), not by definition.

**Required change (sympy):**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`, replace the block at lines 46-69 with the following. The block extracts c0 and c1 from `Y_rho` by static-limit (omega -> 0) and pole-residue identification at `omega^2 = Omega_Q^2`, then runs the round-trip checks.

```python
# Recover (c0, c1) from Y_rho by independent paths, rather than defining them.
# c0 is the static contact (omega -> 0); c1 is the residue at the pole 1/(1 - u),
# where u = omega^2/Omega_Q^2.  We extract via Apart in the variable u to avoid
# the tautological "define-then-assert" pattern.
u = sp.symbols("u", real=True)
Y_rho_u = Y_rho.subs(omega**2/Omega_Q**2, u)
Y_rho_u = sp.together(Y_rho_u)
# Apart in u splits into constant part + pole part 1/(1-u).
Y_rho_apart = sp.apart(Y_rho_u, u)
c0_extracted = sp.simplify(Y_rho_apart.subs(u, 0))            # static limit gives c0
# pole residue at u = 1: coefficient of 1/(1-u) term
c1_extracted = sp.simplify(sp.limit((1 - u) * Y_rho_u, u, 1))

print("c0 extracted from Y_rho (static limit) =", c0_extracted)
print("c1 extracted from Y_rho (pole residue) =", c1_extracted)

expect_zero("c0_extracted - 1/rho_alpha", c0_extracted - 1/rho_alpha)
expect_zero("c1_extracted - (rho_alpha-1)/rho_alpha",
            c1_extracted - (rho_alpha - 1)/rho_alpha)

# With c0, c1 now derived (not defined), the contact-plus-pole reconstruction
# is a real check.
expect_zero(
    "contact-plus-pole reconstruction",
    Y_rho - (c0_extracted + c1_extracted/(1 - omega**2/Omega_Q**2))
)
expect_zero("c0 + c1 - 1", c0_extracted + c1_extracted - 1)

rho_from_c0 = sp.simplify(1/c0)
rho_from_c1 = sp.simplify(1/(1-c1))
zeta_from_c = sp.simplify(c1/c0)

print("rho_alpha from c0 =", rho_from_c0)
print("rho_alpha from c1 =", rho_from_c1)
print("zeta_req from (c0,c1) =", zeta_from_c)

expect_zero("rho(c0(rho)) - rho", rho_from_c0.subs(c0, c0_extracted) - rho_alpha)
expect_zero("rho(c1(rho)) - rho", rho_from_c1.subs(c1, c1_extracted) - rho_alpha)
expect_zero("zeta(c(rho)) - (rho-1)",
            zeta_from_c.subs({c0: c0_extracted, c1: c1_extracted}) - (rho_alpha - 1))
```

Then in the "Minimal isotropic quadrupole module" block (lines 71-79), replace the bare substitution-based `rho_min`/`zeta_min` extraction with one that goes through the paper-form precursor:

```python
# Minimal isotropic quadrupole module: paper Y_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)
Y_Q_paper = sp.Rational(3, 4) + sp.Rational(1, 4)/(1 - omega**2/Omega_Q**2)
Y_Q_paper_u = Y_Q_paper.subs(omega**2/Omega_Q**2, u)
c0_paper = sp.simplify(sp.apart(Y_Q_paper_u, u).subs(u, 0))
c1_paper = sp.simplify(sp.limit((1 - u) * Y_Q_paper_u, u, 1))
expect_zero("c0_paper - 3/4", c0_paper - sp.Rational(3, 4))
expect_zero("c1_paper - 1/4", c1_paper - sp.Rational(1, 4))

rho_min = sp.simplify(rho_from_c0.subs(c0, c0_paper))
zeta_min = sp.simplify(zeta_from_c.subs({c0: c0_paper, c1: c1_paper}))

print("rho_alpha(minimal isotropic module) =", rho_min)
print("zeta_req(minimal isotropic module) =", zeta_min)

expect_zero("rho_min - 4/3", rho_min - sp.Rational(4, 3))
expect_zero("zeta_min - 1/3", zeta_min - sp.Rational(1, 3))
```

**Required change (mathematica):** see F3 below — F3 already prescribes an Apart/Residue-based independent path for the Mathematica file that subsumes F2's mathematica fix.

**Verification:**
- `sp.apart((1 + (rho_alpha - 1)/(1-u))/rho_alpha, u)` should yield `1/rho_alpha + (rho_alpha-1)/rho_alpha/(1-u)` (or equivalent), and the static-limit substitution `u -> 0` returns `1/rho_alpha`. The pole-residue limit returns `(rho_alpha-1)/rho_alpha`. So `c0_extracted - 1/rho_alpha = 0` and `c1_extracted - (rho_alpha-1)/rho_alpha = 0` both reduce to literal zero ✓.
- Substituting `c0_paper = 3/4` into `1/c0` gives `4/3` ✓; substituting `(c0, c1) = (3/4, 1/4)` into `c1/c0` gives `1/3` ✓.

## F3 — mathematica_transliteration (independent derivation path)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage088_loading_ratio_from_minimal_module_mathematica_audit.wl:32-72`

**Issue:**
The Mathematica file mirrors the SymPy file step-for-step (same variable construction, same assertion order, same labels). The second engine must derive the claim by an independent algebraic route.

**Required change:**
Replace the Mathematica audit body (lines 32-70, keeping the existing header `banner`/`expectZero`/`expectTrue`/`pass`/`fail`/`fmt` helpers at lines 1-30 unchanged) with the following Apart/Limit/Residue-based derivation. Note: keep the banner title `"STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE"` (fix the stale 071 label while we're here).

```
banner["STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE"];

Clear[rhoAlpha, omega, omegaQ, u, cMix];
$Assumptions =
  Element[{rhoAlpha, omega, omegaQ, u, cMix}, Reals] &&
  rhoAlpha > 0 && omegaQ > 0 && cMix > 0 && -1 < u < 1;

(* Paper Input: Y_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2). *)
yQpaper = 3/4 + (1/4)/(1 - omega^2/omegaQ^2);
yQpaperU = yQpaper /. omega^2/omegaQ^2 -> u;

(* Independent extraction of c0 and c1 from Y_Q^cons via Apart in u. *)
yQpartial = Apart[yQpaperU, u];
c0Paper = FullSimplify[Limit[yQpaperU, u -> 0]];
c1Paper = FullSimplify[Limit[(1 - u)*yQpaperU, u -> 1]];

Print["Y_Q^cons (paper)   = ", fmt[yQpaper]];
Print["Y_Q^cons partial   = ", fmt[yQpartial]];
Print["c0 (static limit)  = ", fmt[c0Paper]];
Print["c1 (pole residue)  = ", fmt[c1Paper]];

expectZero["c0_paper - 3/4", c0Paper - 3/4];
expectZero["c1_paper - 1/4", c1Paper - 1/4];
expectZero["c0_paper + c1_paper - 1", c0Paper + c1Paper - 1];

(* Loading-ratio extraction from coefficients: rho_alpha = 1/c0, zeta = c1/c0. *)
rhoMin = FullSimplify[1/c0Paper];
zetaMin = FullSimplify[c1Paper/c0Paper];

Print["rho_alpha (= 1/c0) = ", fmt[rhoMin]];
Print["zeta_req (= c1/c0) = ", fmt[zetaMin]];

expectZero["rho_min - 4/3", rhoMin - 4/3];
expectZero["zeta_min - 1/3", zetaMin - 1/3];

(* Reconstruct the contact-plus-pole form from extracted (c0, c1) and confirm
   it matches the paper precursor.  This is the actual coefficient-matching
   claim of the stage. *)
yRhoFromCoeffs = c0Paper + c1Paper/(1 - omega^2/omegaQ^2);
expectZero["paper form - reconstruction from extracted (c0, c1)",
           yQpaper - yRhoFromCoeffs];

(* Also confirm the general rho-parameterized form rebuilds yQpaper after the
   substitution rhoAlpha -> rhoMin.  This exercises the rho parameterization. *)
yRhoParam = 1/rhoAlpha + ((rhoAlpha - 1)/rhoAlpha)/(1 - omega^2/omegaQ^2);
expectZero["rho-parameterized form (rhoAlpha -> rhoMin) - paper form",
           (yRhoParam /. rhoAlpha -> rhoMin) - yQpaper];

(* Stage 68 product identity: Pi_tr = rho_alpha * C_mix. *)
piFromRho = FullSimplify[rhoMin*cMix];
Print["Pi_tr (= rho_min * C_mix) = ", fmt[piFromRho]];
expectZero["Pi_tr_from_rho - (4/3) C_mix", piFromRho - (4/3)*cMix];
expectTrue["1 < rho_min < 2 (symmetric-lowest-twin regime)", 1 < rhoMin < 2];

Print[""];
Print["Stage 088 Mathematica audit passed."];

Exit[0];
```

**Verification:**
- `Limit[yQpaperU, u -> 0]` for `yQpaperU = 3/4 + (1/4)/(1 - u)` returns `3/4 + 1/4 = 1`? — NO. Self-test caught this. The static-contact reading of "c0" is *not* `Limit[yQpaperU, u -> 0]` of the full expression (that gives Y_Q^cons evaluated at omega=0 which is 1). It is the constant part of the partial-fraction expansion in `1/(1-u)`. Use `Apart` to split, then extract the part without `1/(1-u)`.

  Corrected extraction:
  ```
  yQpartial = Apart[yQpaperU, u];          (* should yield 3/4 + (1/4)/(1 - u) *)
  c0Paper = FullSimplify[yQpartial /. u -> Infinity];   (* limit as u -> Infinity gives 3/4 *)
  ```
  Or, more robust:
  ```
  c0Paper = FullSimplify[Limit[yQpaperU - (1/4)/(1 - u) /. {}, u -> 0]];
  ```
  Cleanest: use Coefficient-style extraction via series at u = 1.
  ```
  c0Paper = FullSimplify[yQpaperU /. u -> 0] - FullSimplify[(1/4)/(1 - 0)];  (* NO — still circular *)
  ```

  **Use this clean form instead, which is robust and non-circular:**
  ```
  (* c0 = constant part of Apart[yQpaperU, u]; c1 = pole residue at u = 1. *)
  yQpartial = Apart[yQpaperU, u];
  (* The pole part is the only term containing (1 - u) in the denominator.
     Strip it by setting (1/(1-u)) -> 0 in yQpartial to read off c0. *)
  c0Paper = FullSimplify[yQpartial /. (a_./(1 - u)) -> 0];
  c1Paper = FullSimplify[Residue[yQpaperU, {u, 1}] * (-1)];  (* Residue at u=1 of f(u)/(1-u)-style: 1/(1-u) has residue -1 *)
  ```

  Sympy-side mirror: keep the `apart` + static-limit + pole-residue approach but with the correction noted.

  Implementation note for Codex: the cleanest cross-engine pattern is `series in 1/(1-u)`:
  - c0 = `lim_{u->-inf} yQpaperU` (since `1/(1-u) -> 0` as `u -> -inf` if domain real-valued; but Mathematica may not evaluate cleanly).
  - **Safer**: use `Coefficient` on `(1-u)*yQpaperU` expanded around u=1: `Series[(1 - u)*yQpaperU, {u, 1, 0}]` first term is c1; then `c0 = yQpaperU - c1/(1-u)` evaluated at any safe u (say u=0): `c0 = (yQpaperU /. u -> 0) - c1/(1-0)` = `1 - c1 = 1 - 1/4 = 3/4`. This gives c0 by subtraction-after-residue-extraction, no Apart-pattern-matching needed.

  **Codex, apply this corrected extraction in both engines (sympy and mathematica). Verbatim form for sympy:**

  ```python
  # Extract c1 first via pole residue at u = 1, then c0 by subtraction at u = 0.
  c1_extracted = sp.simplify(sp.limit((1 - u) * Y_rho_u, u, 1))
  c0_extracted = sp.simplify(Y_rho_u.subs(u, 0) - c1_extracted / (1 - 0))
  ```

  Mathematica:
  ```
  c1Paper = FullSimplify[Limit[(1 - u)*yQpaperU, u -> 1]];
  c0Paper = FullSimplify[(yQpaperU /. u -> 0) - c1Paper/(1 - 0)];
  ```

  Trivial-case verification:
  - `Y_rho_u = (1/rho_alpha) + ((rho_alpha-1)/rho_alpha)/(1-u)`.
  - `Limit[(1-u)*Y_rho_u, u->1] = (rho_alpha-1)/rho_alpha` ✓ (c1).
  - `Y_rho_u.subs(u,0) = 1/rho_alpha + (rho_alpha-1)/rho_alpha = 1`; subtract `c1/(1-0) = (rho_alpha-1)/rho_alpha` gives `1 - (rho_alpha-1)/rho_alpha = 1/rho_alpha` ✓ (c0).
  - For paper input `Y_Q^cons = 3/4 + (1/4)/(1-u)`: `c1 = 1/4`, `c0 = (3/4 + 1/4) - 1/4 = 3/4` ✓.

  **Codex must use these corrected extraction lines in both engine fixes for F2 (sympy) and F3 (mathematica). The other Apart-based fragments above were exploratory; ignore them. The two-line extract-c1-then-c0-by-subtraction pattern is the one to apply.**

**Verification (post-fix):**
After Codex applies F1+F2+F3, the verifier will run `redteam exec-sympy 088` and `redteam exec-mathematica 088` and confirm:
1. Both exit 0.
2. Final printed values: `rho_alpha = 4/3`, `zeta_req = 1/3`, `Pi_tr_from_rho = (4/3) C_mix`.
3. The Mathematica script's intermediate symbol sequence (Apart/Limit/Residue/Series) no longer mirrors the SymPy intermediate sequence line-by-line.
4. Mutating the paper Input (e.g. `Y_Q^cons -> 1/2 + (1/2)/(1-omega^2/Omega_Q^2)`) propagates to `c0=1/2, rho_min=2, Pi_tr=2*C_mix` — exercising the chain rather than restating constants.
