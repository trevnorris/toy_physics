---
unit_id: 189
batch: V.3
created_at: 2026-06-01T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-01T11:48:58-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
iteration: 2
---

# Codex directive — unit 189

Apply F1 and F2 (script-side, non-paper_misalignment) in order. After applying each, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F3 is a `paper_misalignment` (stale stage-number label) — now RESOLVED via the settled canonical-stage-number convention (see `## RESOLVED — F3` below). Apply the authorized banner relabel as part of this directive. Do NOT edit paper.tex or notes/.

If a required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the script (`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:78-89`

**Issue:** The Section II selected-branch identity check is a self-inverted round-trip. `Rtarget_formula` (line 83) is defined as exactly `Lambda0*(1-epseta)/T2_direct`, so `T2_selected = Lambda0*(1-epseta)/Rtarget_formula` (line 84) is forced to equal `T2_direct`, and the check at line 89 (`T2_selected - T2_direct`) is identically zero regardless of any input. The check verifies nothing about the physics of `R_target`, and the paper's load-bearing constant `Lambda0 = 27 pi^2 G c_s^5 / (20 a^5 c^5)` (notes line 144) is never instantiated — `Lambda0` stays a free symbol.

**Required change:**
Make the selected-branch identity a genuine cancellation by treating `R_target` as an INDEPENDENT carried symbol (not a back-inversion of `T2_direct`) and checking the boxed product relation `R_target * T^2 = Lambda_0 (1 - epseta)` directly, AND instantiate the `Lambda_0` value.

Replace the block at lines 81-89:

Before (lines 81-89):
```python
Rtarget = sp.symbols("R_target", positive=True, real=True)
T2_direct = sp.simplify(Zw * (1 + rho) ** 2 / (OmW2 * (1 - epsW) ** 2))
Rtarget_formula = sp.simplify(Lambda0 * OmW2 * (1 - epseta) * (1 - epsW) ** 2 / (Zw * (1 + rho) ** 2))
T2_selected = sp.simplify(Lambda0 * (1 - epseta) / Rtarget_formula)
print("T_A^2 (direct continuum form) =")
sp.pprint(T2_direct)
print("R_target (selected-branch form) =")
sp.pprint(Rtarget_formula)
expect_zero("Lambda_0 (1-epseta) / R_target - T_A^2", T2_selected - T2_direct)
```

After:
```python
# G, c, a, cs are needed for the explicit Lambda_0 value; declare here so the
# value check is available in this section (they are also declared in VI; if a
# NameError results from ordering, hoist these symbol declarations above Section II).
G_II, c_II, a_II, cs_II = sp.symbols("G c a c_s", positive=True, real=True)
Lambda0_val = sp.Rational(27, 20) * sp.pi**2 * G_II * cs_II**5 / (a_II**5 * c_II**5)

T2_direct = sp.simplify(Zw * (1 + rho) ** 2 / (OmW2 * (1 - epsW) ** 2))
# R_target is an INDEPENDENT continuum-kernel quantity; its one-port realization
# (carried from the continuum-kernel stages, notes 242) is the selected-branch
# demand. Verify the boxed selected-branch identity R_target * T^2 = Lambda_0 (1-epseta)
# as a genuine relation, not by inverting T2_direct.
Rtarget_oneport = sp.simplify(Lambda0 * (1 - epseta) / T2_direct)  # the demanded R_target
print("T_A^2 (direct continuum form) =")
sp.pprint(T2_direct)
print("R_target (one-port selected-branch demand) =")
sp.pprint(Rtarget_oneport)
# Genuine product identity (boxed, notes 139): this is a real cancellation, and it
# is the relation the stage claims, with Lambda_0 carried as the front factor.
expect_zero("R_target * T_A^2 - Lambda_0 (1-epseta)", Rtarget_oneport * T2_direct - Lambda0 * (1 - epseta))
# Confirm the load-bearing front-factor value stated in the notes (line 144).
expect_zero("Lambda_0 - 27 pi^2 G c_s^5 / (20 a^5 c^5)", Lambda0_val - sp.Rational(27, 20) * sp.pi**2 * G_II * cs_II**5 / (a_II**5 * c_II**5))
```

NOTE on the `Lambda_0` value check: as written above the `Lambda_0 - ...` check is itself trivially zero (value minus its own definition). That is intentional only as a placeholder — the substantive content is that `Lambda0` (the free symbol used in the identity) is documented as equal to `Lambda0_val`. If you can wire `Lambda0` to `Lambda0_val` upstream (substitute `Lambda0 -> Lambda0_val` in the product identity and confirm it still cancels), do that instead and drop the trivial value line. Prefer:
```python
expect_zero("R_target * T_A^2 - Lambda_0 (1-epseta) [with Lambda_0 value]",
            (Rtarget_oneport * T2_direct - Lambda0 * (1 - epseta)).subs(Lambda0, Lambda0_val))
```
The key requirement: the surviving Section II check must NOT be `T2_selected - T2_direct` where `T2_selected` is `Lambda0*(1-epseta)/(Lambda0*(1-epseta)/T2_direct)`. It must be the product form `Rtarget * T2 - Lambda0*(1-epseta)`.

**Self-test (already done by auditor):** With `Rtarget_oneport = Lambda0*(1-epseta)/T2_direct`, the product `Rtarget_oneport*T2_direct - Lambda0*(1-epseta) = Lambda0*(1-epseta) - Lambda0*(1-epseta) = 0` — still zero, BUT now it is a cancellation in the product (an additive identity), not a `x/x` round-trip; it would catch a wrong factor in either `R_target` or `T^2`. This is the faithful representation of the boxed claim `R_target T^2 = Lambda_0(1-epseta)`.

**Verification command:** The verifier runs `redteam exec-sympy 189` and confirms (a) the line `R_target * T_A^2 - Lambda_0 (1-epseta) = 0` (or the `[with Lambda_0 value]` variant) appears, (b) the old `Lambda_0 (1-epseta) / R_target - T_A^2` line is gone, and (c) the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
- summary: Replaced the selected-branch quotient round-trip with the product identity and checked it after substituting the explicit Lambda_0 value.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:40-72`

**Issue:** The Section I defect-notation checks (lines 70-72) and the compatibility-row check (line 63) re-assert their own definitions and cannot fail. `Xi := dlnT2 := dn - Bstar*dr` (lines 45, 68), so `Xi - (dn - Bstar*dr) == 0` identically (line 70). `Rcal := dlnRtarget` (lines 47, 69), so `Rcal - dlnRtarget == 0` identically (line 71). `dlnRtarget := dlnOneMinus - dlnT2` (line 47), so `(Rcal+Xi)-dlnOneMinus == 0` (line 72) and `trf[2]+trf[0]-trf[1] == 0` (line 63) are both forced. The notes present the compatibility relation `δln R_target + δln T^2 - δln(1-epseta) = 0` (notes 217-224) as a DERIVED rank-2 consequence of the three transfer identities, not as a definition.

**Required change:**
Introduce the Stage-239 (notes "Stage 239") observable relations as INDEPENDENT defect symbols and derive the transfer identities from them, so each check is a real cancellation. The matrix check at line 62 (A1) is already non-tautological — keep it.

Add, after the existing matrix/compatibility prints (i.e., replacing lines 66-72), independent defect symbols and substitution-based checks:

Before (lines 66-72):
```python
# Match to defect notation from Stage 188.
Theta = dr
Xi = dlnT2
Rcal = dlnRtarget
expect_zero("Xi - (dln N_* - B_* dln R_tr)", Xi - (dn - Bstar * dr))
expect_zero("R_1 - dln R_target", Rcal - dlnRtarget)
expect_zero("(R_1 + Xi_1) - dln(1-epseta)", (Rcal + Xi) - dlnOneMinus)
```

After:
```python
# Match to defect notation from the source branch-observable stage (notes "Stage 239").
# Introduce the upstream defect scalars as INDEPENDENT symbols and impose the
# Stage-239 observable definitions (notes lines 120-122):
#   dln R_tr   = Theta1
#   dln N_*    = Xi1 + B_* Theta1     (since N_* := T^2 R_tr^{B_*})
#   dln epseta = Sigma_eta
# Then the transfer identities dln T^2 = Xi1, dln R_target = R1, dln(1-epseta)=R1+Xi1
# must FOLLOW (not be assumed).
Theta1, Xi1, Sigma_eta = sp.symbols("Theta_1 Xi_1 Sigma_eta", real=True)
obs_sub = {dr: Theta1, dn: Xi1 + Bstar * Theta1, de: Sigma_eta}
# dln T^2 computed from the observable packet must equal the independent Xi1:
expect_zero("dln T^2 - Xi_1", dlnT2.subs(obs_sub) - Xi1)
# Define R1 from the selected-branch identity dln R_target (notes 172):
R1 = dlnRtarget.subs(obs_sub)
# Complementary dressing relation dln(1-epseta) = R1 + Xi1 must follow:
expect_zero("dln(1-epseta) - (R_1 + Xi_1)", dlnOneMinus.subs(obs_sub) - (R1 + dlnT2.subs(obs_sub)))
# Rank-2 compatibility as a DERIVED relation among the three identities
# (computed independently from R1, Xi1, and the dressing relation), not from
# the definition dlnRtarget := dlnOneMinus - dlnT2:
expect_zero(
    "compatibility: dln R_target + dln T^2 - dln(1-epseta)",
    R1 + dlnT2.subs(obs_sub) - dlnOneMinus.subs(obs_sub),
)
```

Also REMOVE the now-redundant tautological compatibility-row check at line 63 (`expect_zero("selected-branch compatibility row", trf[2] + trf[0] - trf[1])`), OR keep it but understand it is subsumed; preferred: remove line 63 since the new derived-compatibility check above replaces it. Keep the matrix check at line 62 (it is genuine).

**Self-test (already done by auditor):**
- `dlnT2.subs(obs_sub) - Xi1 = (dn - Bstar*dr).subs{dr:Theta1, dn:Xi1+Bstar*Theta1} - Xi1 = (Xi1 + Bstar*Theta1 - Bstar*Theta1) - Xi1 = 0` — genuine cancellation that depends on the observable substitution; it would FAIL if `N_*`'s `B_*` weight were wrong. Non-vacuous. Good.
- `dlnOneMinus.subs(obs_sub) - (R1 + dlnT2.subs(obs_sub))`: `dlnOneMinus = -epsetas/(1-epsetas)*Sigma_eta` after sub; `R1 = (dlnOneMinus - dlnT2).subs(obs_sub) = -epsetas/(1-epsetas)*Sigma_eta - Xi1`; `R1 + dlnT2.subs = (-epsetas/(1-epsetas)*Sigma_eta - Xi1) + Xi1 = -epsetas/(1-epsetas)*Sigma_eta = dlnOneMinus.subs`. Residual 0 — genuine. Good.
- The compatibility check: `R1 + dlnT2.subs - dlnOneMinus.subs = (dlnOneMinus.subs - dlnT2.subs) + dlnT2.subs - dlnOneMinus.subs = 0` — this remains an additive cancellation through the substituted (physical) quantities rather than through the bare definition; it is the derived rank-2 condition. Acceptable. (Note: because `R1` is itself defined via `dlnRtarget = dlnOneMinus - dlnT2`, this last check is still partly definitional. If you want it fully independent, define `R1` directly as `dlnOneMinus.subs(obs_sub) - Xi1` from the dressing relation `dln(1-epseta)=R1+Xi1` and then check `R1 == dlnRtarget.subs(obs_sub)` as the genuine consistency. Prefer that stronger form if it runs cleanly.)

**Verification command:** The verifier runs `redteam exec-sympy 189` and confirms the new lines `dln T^2 - Xi_1 = 0`, `dln(1-epseta) - (R_1 + Xi_1) = 0`, and `compatibility: ... = 0` appear, the old lines 70-72 and the old line-63 compatibility-row check are gone, and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
- summary: Removed the tautological compatibility-row check and derived the transfer identities from independent observable defect substitutions.
- deviation: none

## F3 — paper_misalignment (subtype: notes_contradicts_script)

DO NOT APPLY. Hold for user resolution.

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_189.tex:1` quote: "Stage 189: Transfer-shape / outgoing-prefactor compiler" (this is Stage 189; filename is `stage189`).
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage189_transfer_shape_prefactor_compiler.md:96` quote: "So Stage 240 is the exact compiler ..." (notes body calls it Stage 240, source = Stage 239).

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:35` quote: `banner("STAGE 172 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:182` quote: `banner("STAGE 172 LEDGER")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py:38,66,183` quote: "Stage 188 observable packet" / "from Stage 188"

## Resolve before fix_loop

Four different stage numbers (172 / 188 / 189 / 240) appear across the script banner, script comments, paper card, and notes for this one stage. The restored card and filename say 189; the notes body says 240 with upstream source 239; the script banner/ledger say 172 with upstream "Stage 188". The verified MATH is correct and aligned regardless — this is purely a label/traceability inconsistency. Which display number should the script banner/ledger/comments use, and which upstream-source number is correct (188 or 239)?

Possible directions (the user picks one):
- (a) Use restored numbering → change script banner (line 35), ledger (line 182) to "STAGE 189", and the upstream "Stage 188" references (lines 38, 66, 183) to match the restored upstream stage number for the branch-observable packet. No assertion changes.
- (b) Use original numbering → change banner/ledger to "STAGE 240" and upstream to "Stage 239" to match the notes body. No assertion changes.
- (c) The "172"/"188" labels point to a genuinely different upstream source than 239 → flag for deeper review of which observable packet is actually being compiled before any relabel.

F1 and F2 (above) are independent of F3 and may be applied now.

## RESOLVED — F3 (settled canonical-stage-number convention)

Resolution = direction (a). Per the project's settled convention (stale script "STAGE NN" labels → the canonical internal stage number; Codex-CONCUR, established in the remediation batches), the canonical number for this unit is **189** (paper card + filename). Authorized script-side edit — apply it:
- Relabel the two banner string literals only: line 35 `banner("STAGE 172 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER")` → `banner("STAGE 189 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER")`; line 182 `banner("STAGE 172 LEDGER")` → `banner("STAGE 189 LEDGER")`. String-only; no math change.
- Leave the "Stage 188" upstream references (lines 38, 66, 183) AS-IS — 188 is the canonical branch-observable source stage, already consistent with restored numbering.
- The notes-prose legacy numbers (240 / source 239) are OUT OF red-team scope (notes are paper files) → to be logged to PAPER_CLEANUP_TRACKER, not edited here.

`needs_user_resolution` is cleared.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
- summary: Relabeled only the two authorized script banner string literals from STAGE 172 to STAGE 189.
- deviation: none

## F4 — missing_mathematica

**Issue:** Stage 189 is dual-engine-capable (exact rational/series identities and a prefactor compiler) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the CORRECTED SymPy script (after F1/F2 above) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_189.tex` and the notes file are the source of truth. Mirror the CORRECTED checks (the product-form selected-branch identity `R_target*T^2 - Lambda_0(1-epseta)` and the observable-symbol-derived transfer identities), NOT the broken tautological originals.
- Use Mathematica-NATIVE primitives (`Series`+`Coefficient`, `D[...]/.eps->0`, `Solve`/`Reduce` for the compiler inversion) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (`expectZero`, `$Assumptions` positivity, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route. RUN it (`timeout 600 math -script <path>`) and iterate to exit 0; a timeout (124) is a failure — reformulate, don't raise the cap.

**Verification command:** the verifier runs `redteam exec-mathematica 189`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route whose conclusions agree with the SymPy engine.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.wl`
- summary: Added an independent Mathematica audit for the corrected transfer-packet, selected-branch, prefactor, outgoing, normalization, and constant-prefactor checks.
- deviation: none

## Iteration 2 (delta) — F1 re-fix + F4 mirror (consult 019e843e, CONCUR)

The iter-1 F1 fix did NOT de-tautologize. `Rtarget_oneport = Lambda0*(1-epseta)/T2_direct` (line ~105) is a back-definition, so the line-~114 check `Rtarget_oneport*T2_direct - Lambda0*(1-epseta)` is identically 0 for ANY `T2_direct`/`Lambda0` (the `.subs(Lambda0, Lambda0_val)` cancels too). The selected-branch identity `R_target*T^2 = Lambda_0(1-epseta)` is DEFINITIONAL (rank-2 compatibility, notes 215-226) — it defines R_target, so it cannot serve as a verification.

**Re-fix — SymPy (Section II, ~lines 96-114):**
1. DEMOTE the selected-branch identity to a printed DEFINITION of R_target (a comment + `print` showing `R_target := Lambda_0(1-epseta)/T^2` and `Lambda_0`'s value). REMOVE the tautological `expect_zero("R_target * T^2 - Lambda_0(1-epseta)", ...)` entirely.
2. ADD the genuine direct-slope bridge (notes §2.1-2.3, eq line 291 `δln T_A^2 = ε λ_A Ξ_1`): perturb the CONCRETE coherent continuum transfer shape `T2_coh = Zw(1+chi0)^2/(OmW2(1-eps)^2)` (notes §2.2 line 280; variables already in Section II) to first order in the perturbation parameter, compute the first-order log-slope `Series[log(T2_pert/T2_0)]` (or the `D[log(...)]|_{0}` equivalent), and assert it equals `eps*lambda_A*Xi1_closed` — where `Xi1_closed` is the INDEPENDENTLY-defined first grouped defect scalar, copied from the upstream defect law, NOT defined from this derivative. Follow the stage-181 pattern: `scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:76` (eps split + Xi1 definition) then `:110/:114` (the Xi1-derived-from-T^2 derivative check). This check FAILS if the continuum closed form or Xi1's definition is wrong.

**Re-fix — Mathematica `.wl` (F4 mirror):** apply the SAME correction. NO `Rtarget = Lambda0(1-epseta)/T2; ... Rtarget*T2 - ...` assertion (same tautology). Use a printed R_target DEFINITION + the non-tautological first-order slope check (perturb `T2_coh`; `Series[Log[...]]` or `D[Log[...],eps]/.eps->0`; compare to `eps*lambda_A*Xi1` with Xi1 supplied independently).

Run BOTH engines to exit 0 and iterate. The verifier will confirm: (a) the Section II selected-branch product `expect_zero` is GONE in both engines, (b) the new slope-bridge check against an INDEPENDENTLY-defined Xi1 is present in both engines and passes, (c) exit 0 both engines.

## Applied: Iteration 2

- files_changed:
  - `scripts/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage189_transfer_shape_prefactor_compiler_mathematica_audit.wl`
- summary: Demoted the selected-branch product relation to a printed definition and added the non-tautological coherent transfer-shape slope bridge in both engines.
- deviation: none
