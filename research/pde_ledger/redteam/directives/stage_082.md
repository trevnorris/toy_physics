---
unit_id: 082
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-25T00:23:28-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 082

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:87-94`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:71-75`

**Issue:**
The two `expect_zero` checks at lines 87-94 of the SymPy script (and lines 71-75 of the Mathematica mirror) are arithmetic identities on hand-baked integer constants — `100 * 37**2 == 136900` and the substitution `Upsilon_w -> 100*Theta_w` into `Upsilon_w * 1369`. They can fail only if the engine's integer multiplication is broken, not if any physics is wrong. Replace the `expect_zero` (assertion) form with plain `print` (display-only) form, and add a comment stating that this script does not verify the Family-1 strength constants.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`, replace lines 87-94 (the two `expect_zero(...)` blocks for Xi_F1) with:

```python
# Note: the two equalities below are arithmetic on the hand-supplied integers
# 37, 100, 1369, 136900. They are not independent verifications of the Family-1
# strength identity — the upstream stage that derives Lambda_ell = 37 and the
# convention Upsilon_w = 100 * Theta_w is responsible for those facts. Here we
# only display the resulting Xi_F1 forms for downstream readers.
print(f"Xi_F1(Theta_w) - 136900 Theta_w = {sp.simplify(Xi_F1_from_Theta - sp.Integer(136900) * Theta_w)}  (display only)")
print(f"Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w) = {sp.simplify(Xi_F1_from_Upsilon.subs(Upsilon_w, 100 * Theta_w) - Xi_F1_from_Theta)}  (display only)")
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, replace lines 71-75 (the two `expectZero[...]` blocks for Xi_F1) with:

```mathematica
(* Note: the two equalities below are arithmetic on the hand-supplied integers
   37, 100, 1369, 136900. They are not independent verifications of the
   Family-1 strength identity — the upstream stage that derives lambdaEll = 37
   and the convention upsilonW = 100 * thetaW is responsible for those facts.
   Here we only display the resulting Xi_F1 forms for downstream readers. *)
Print["Xi_F1(Theta_w) - 136900 Theta_w = ", fmt[FullSimplify[xiF1FromTheta - 136900*thetaW]], "  (display only)"];
Print["Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w) = ",
      fmt[FullSimplify[(xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta, Assumptions -> thetaW > 0]],
      "  (display only)"];
```

Do not change lines 78-94 of the SymPy script that *define* `Lambda_ell`, `Xi_F1_from_Upsilon`, `Xi_F1_from_Theta`. Only the two `expect_zero` calls become `print`. Same for Mathematica lines 65-70.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 082` and `redteam exec-mathematica 082` and confirm:
- the output no longer contains an assertion line `Xi_F1(Theta_w) - 136900 Theta_w = 0` followed by `PASS:` (Mathematica) or a successful `expect_zero` print (SymPy);
- the output instead contains `(display only)` suffix on those two lines;
- both scripts still exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- summary: Replaced the Family-1 strength assertions with display-only prints and comments explaining the constants are not independently verified here.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:75` (insert after this line)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:63` (insert after this line)

**Issue:**
Five of the seven assertions test the single algebraic identity `zeta_req(C_mix * Q(z)) == z` at different symbolic z. The directional content advertised in the "Theorem ledger" (`Pi_tr <= ...` → success, `Pi_tr >= ...` → failure) is not exercised. Add two new derivative-based assertions that exercise the partial-derivative structure underwriting those directional theorems.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`, insert a new block immediately after line 75 (after the `expect_zero("R_quad(Pi_fail, ...)", ...)` call's closing parenthesis) and before line 77's blank line / Family-1 section:

```python
# Directional content of R_quad: verify the partial derivatives that
# underwrite the "guaranteed success / guaranteed failure" theorems.
dR_dzeta_phys = sp.simplify(sp.diff(R_quad, zeta_phys))
expect_zero("dR_quad/dzeta_phys + 1", dR_dzeta_phys + 1)

dzeta_req_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
dR_dPi_at_zeta_minus = sp.simplify(sp.diff(R_quad.subs(zeta_phys, zeta_minus), Pi_tr))
expect_zero(
    "dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)",
    dR_dPi_at_zeta_minus - dzeta_req_dPi,
)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, insert immediately after line 63 (after the closing `];` of the `expectZero["R_quad(Pi_fail, ...)", ...]`) and before line 65's `lambdaEll = 37;`:

```mathematica
(* Directional content of R_quad: verify the partial derivatives that
   underwrite the "guaranteed success / guaranteed failure" theorems. *)
dRDzetaPhys = FullSimplify[D[rQuad, zetaPhys], Assumptions -> $Assumptions];
expectZero["dR_quad/dzeta_phys + 1", dRDzetaPhys + 1];

dZetaReqDPi = FullSimplify[D[zetaReq, PiTr], Assumptions -> $Assumptions];
dRDPiAtZetaMinus = FullSimplify[D[rQuad /. zetaPhys -> zetaMinus, PiTr], Assumptions -> $Assumptions];
expectZero[
  "dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)",
  FullSimplify[dRDPiAtZetaMinus - dZetaReqDPi, Assumptions -> $Assumptions]
];
```

Pre-check (already done in self-test): `zeta_phys` appears in `R_quad` (line 64 of SymPy defines `R_quad = zeta_req - zeta_phys`), so `diff(R_quad, zeta_phys) = -1` non-trivially. `Pi_tr` appears in `zeta_req` (line 38), and after `zeta_phys -> zeta_minus` the `R_quad` still depends on `Pi_tr`, so `diff(R_quad.subs(zeta_phys, zeta_minus), Pi_tr) = diff(zeta_req, Pi_tr)`, which is nonzero (equals `(C_mix*(1-eps_blk))/(C_mix - eps_blk*(2*C_mix - Pi_tr))**2` after simplification). Variable-independence trap avoided.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 082` and `redteam exec-mathematica 082` and confirm:
- SymPy output contains the new lines `dR_quad/dzeta_phys + 1 = 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-) = 0`;
- Mathematica output contains the corresponding `PASS:` lines;
- both scripts still exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- summary: Added derivative-based checks for the quadrupole residual's zeta_phys and Pi_tr directional structure.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:33` (replace this line and add supporting setup)

**Issue:**
The `.wl` script's `zetaReq` (line 33) is a syntactic transliteration of the SymPy script's `zeta_req` (line 38). The "inverse map" assertion at line 39 then verifies that the same closed-form expression (rendered in two engine syntaxes) is the inverse of `qMap` — agreement is guaranteed because both sides are copies. Replace the hand-supplied `zetaReq` with one derived from `qMap` via Mathematica's `Solve`, so the engine independently re-derives the inverse.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, replace line 33 (currently `zetaReq = FullSimplify[(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr)), Assumptions -> $Assumptions];`) with:

```mathematica
(* Independent re-derivation: solve PiTr == cMix * qMap(zetaSym) for zetaSym,
   then rename to zeta. This forces Mathematica to find the inverse of qMap
   rather than restating the SymPy-side closed form. *)
qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
zetaSym = Unique["zetaSym"];
zetaReqSolList = Solve[PiTr == cMix*((1 + (1 - 2*epsBlk)*zetaSym)/(1 - epsBlk*zetaSym)), zetaSym];
zetaReq = FullSimplify[(zetaSym /. First[zetaReqSolList]) /. ConditionalExpression[x_, _] :> x,
                       Assumptions -> $Assumptions];
```

Then DELETE line 34 (the now-duplicate `qMap = ...` definition), since the new block above already defines `qMap` before using it. (If Codex finds this rearrangement awkward, alternative: keep line 34 in its current location, and remove the duplicate `qMap = ...` inside the new block above — the key requirement is that `zetaReq` be set via `Solve`, not via the literal fraction.)

The `ConditionalExpression` strip is a defensive pattern from project memory (`feedback_mathematica_script_idioms`) — `Solve` may wrap branches in `ConditionalExpression`, which `FullSimplify` would otherwise carry forward and break subsequent `=== 0` checks.

Pre-check (self-test): `PiTr == cMix * qMap(zetaSym)` is linear in `zetaSym` once you clear the denominator `(1 - epsBlk*zetaSym)`, so `Solve` returns exactly one branch. Substituting back, the derived `zetaReq` must algebraically equal the closed form `(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr))`. All downstream `expectZero` assertions (lines 39-51 and the F2-added derivative checks) continue to hold because they only reference `zetaReq` and `qMap`, not the literal expression.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 082` and confirm:
- the `.wl` file's source no longer contains the literal expression `(PiTr - cMix)/(cMix - epsBlk*(2*cMix - PiTr))` outside a comment;
- the `.wl` file contains a `Solve[PiTr == cMix*...` call setting `zetaReq` (or a temporary that flows into `zetaReq`);
- the script still exits 0 and the printed `zeta_req` line matches the previous output up to algebraic equivalence.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- summary: Re-derived zetaReq in Mathematica via Solve from the qMap relation instead of assigning the closed-form inverse directly.
- deviation: none

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:65`

**Issue:**
The integer `Lambda_ell = 37` and the ratio factor `100` (between `Upsilon_w` and `Theta_w`) are dropped into the script as literals with no provenance comment. Add a comment naming the upstream stage that establishes these constants. If the upstream stage is unknown to Codex, use a `TODO(provenance):` marker — do not invent a stage number.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`, immediately before line 79 (`Lambda_ell = sp.Integer(37)`), insert:

```python
# TODO(provenance): Lambda_ell = 37 and the convention Upsilon_w = 100 * Theta_w
# are carry-forward constants. Cite the upstream stage's script (likely an
# earlier moving_throat_pde_stage*_sympy_audit.py) that derives them. Until
# then, this stage treats them as inputs and only displays their consequences.
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, immediately before line 65 (`lambdaEll = 37;`), insert:

```mathematica
(* TODO(provenance): lambdaEll = 37 and the convention upsilonW = 100 * thetaW
   are carry-forward constants. Cite the upstream stage's script (likely an
   earlier moving_throat_pde_stage*_mathematica_audit.wl) that derives them.
   Until then, this stage treats them as inputs and only displays their
   consequences. *)
```

Do not change any numeric values. Comment-only edit.

**Verification command:**
After Codex applies, the verifier will inspect the two files (no execution needed) and confirm:
- the SymPy script contains a `# TODO(provenance):` comment immediately preceding `Lambda_ell = sp.Integer(37)`;
- the Mathematica script contains a `(* TODO(provenance): ... *)` comment immediately preceding `lambdaEll = 37;`;
- both scripts still parse and exit 0 when re-run.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- summary: Added provenance TODO comments immediately before the Family-1 carry-forward constants in both scripts.
- deviation: none
