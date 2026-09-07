# Independent review — S11c-b carrier ablation-harness knife-list (orchestrator-written directive)

You are one of two independent review legs. Working dir: `/var/projects/toy_physics`. You have repo read access;
you may write and run small check scripts under `/tmp`. ⛔ Do not modify the working tree. ⛔ You are not handed
the other leg's output or any prior review of this artifact — assess independently.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_carrier_ablation_harness_directive.md`

This is an **orchestrator-owned knife-list**: it tells a blind builder (`gpt-6-astra`) exactly which construction
sites to ablate, with which FORM perturbation, in two committed ablation harnesses (one per engine) that certify
the **pressure-slot carrier** of the S11c-b slab operator is genuinely *computed* (not hand-typed / frozen /
tautological). It is a **re-author**: a prior draft was found to target the wrong object, so two knives were
dropped and the face-law channels rewritten. Your job is to judge whether **this** list targets the **right
object**, at the **right sites**, with **genuine FORM** perturbations and **correct impact cones**, and is
**implementable** as written.

## The engines the directive governs (verify every site claim against these)
- SymPy: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
- WL:    `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- Reference for the carrier's atoms/reduction (⛔ the harness must NOT import it, but you may read it to check the
  claim): `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py`
  (`pressure_coefficients`, `pressure_slots`).

## What to check — derive/verify independently, cite line numbers for every claim
Report a finding only if it catches a way the **harness would certify the wrong thing** (or fail to run). For each
item, say whether the directive is correct, and if not, the exact site + the fix.

1. **Is the target object right?** The directive claims the carrier `∂(emitted slab rows)/∂(native pressure
   atoms)|_{P=0}` is **face-law-sourced** — pressure enters the rows ONLY through the permeable-face laws (Λ_A
   response → θ row; traction/virtual-work → U/e_W rows) — and that the **§3a energy basis** and the
   **virtual-constraint fold** are **pressure-independent** (so the two dropped knives, old K1 energy-spurion /
   old K2 constraint-fold, would have moved objects NOT in the carrier). **Verify this independently** against the
   engines (e.g. does the energy density / constraint source carry any pressure atom? do the face laws?). If the
   drop is wrong — if energy or constraint *does* carry a pressure atom that reaches the emitted rows — that is a
   top-severity finding.

2. **Site accuracy.** For each knife (K_A, K_T, K_W), do the cited functions/lines in BOTH engines carry what the
   directive says? Specifically: K_A — the SymPy closure response-flux fold (`closure_residuals` /
   `closure_residual_sum` / the `mass_balance_source - closure_residual_sum` fold) and the WL `flux =
   lambdaAResponse affinity + lambdaVResponse normalVelocity`; K_T — the SymPy `face_generalized_force_rows` +
   its addition into `U_BODY_BALANCE`/`E_W_BALANCE`, and the WL `tractionPressure = pressureField[sign] + …` →
   `virtualWork` → `faceGeneralizedRows`; K_W — the SymPy `delta_p_*` registration (claimed theater) vs the
   consumed row expressions, and the WL `pressureField[±1] = pressureUpper/pressureLower`.

3. **Genuine FORM, and disjointness.** Is each FORM a real **structural** change (leaves the variable-coefficient
   family), not a coefficient rescale or a no-op? Check the disjointness claim: K_A and K_T patch **disjoint**
   downstream sites and **neither** patches the shared `affinity` (WL `:1079`) / shared source — so K_A moves only
   the θ-carrier and K_T only the U/e_W-carrier. If patching K_A's site also moves the U/e_W carrier (or vice
   versa), the channels are not disjoint and the cone is wrong. Check the K_W claim that the value-slot ↔ ∂_w-slot
   swap is **dimensionally invalid** (`delta_p` vs `d_w_delta_p`) and that only a same-dimension **face** swap is
   admissible.

4. **Impact cones.** Are the must-MOVE / must-NOT-move sets correct for each knife? Is there a **coverage hole** —
   a pressure channel that reaches an emitted row but has **no** knife? Is there a knife whose must-MOVE object is
   actually pressure-independent (a wrong-object knife, the prior draft's defect)?

5. **Observation + entrypoint implementability.** Is the pinned observation correct and runnable: SymPy — is
   `import` truly side-effect-free (4-case emit under `if __name__=="__main__"`) so `build_operator("LAB_HELD",
   "RHO4_CONSTANT","EULERIAN")` yields the single-case operator, and are the atoms `delta_p_{plus,minus}` /
   `d_w_delta_p_{plus,minus}` the right ones (match N6)? WL — does `S11CB_PRIMARIES_ONLY` + restricting
   `branches`/`densities` to the single case actually yield `slabOperatorPayload` for that case without the
   ≥64 GB `COUPLING_KERNEL`/tower path, and do the emitted rows carry `pressureUpper/pressureLower` so `∂/∂`
   them has bite? Is the `→0` (P=0) reduction the right one? Does the directive correctly forbid reimplementing
   N6's row re-derivation while allowing the trivial generic differentiation?

6. **Leaks + discipline.** Does the directive leak any expected **value / count / ratio** or a **pass-condition**
   (anything a builder could iterate toward)? Are the three PRINT-not-assert clauses intact? Is exactly ONE
   function + ONE FORM pinned per knife (no "or"/"e.g." builder discretion)? Is the self-ablation spec (dead-path
   + carrier-extractor mutations → must-MOVE diff identically zero; live rescale → nonzero) sound and sufficient
   to catch a self-reporting harness?

## Output
A per-item verdict with cited lines, then a final line: **SOUND** (ready to build) or **NOT-SOUND** (with the
exact site + fix for every blocking issue). Prose assertions about the physics are discounted unless grounded in
a cited engine line or a small check script whose path + literal stdout you report.
