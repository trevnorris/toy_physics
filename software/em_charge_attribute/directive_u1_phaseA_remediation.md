# Directive — U1 Phase-A remediation round (post four-leg verification)

**Status:** the Phase-A REBUILD's physics verdict — **`U1_BASE_OK`** (normalizable translational zero mode; 4 power-law channels ν=2, 3 exponential gaps, 1 algebraic gap) — was **VERIFIED GENUINE** by an independent adversarial recomputation (all 8 channels re-derived from the whitelisted action and matched to ≤1e-8; verdict unflippable by planted tokens; flips only with physics changes, with computed leaf names; d=5/d=6 generality probes clean). Arbiter reproduced deterministically; runner honesty confirmed. **Do NOT change the physics approach.** This round fixes the remaining acceptance defects found by the four-leg review (each confirmed by ≥2 independent reviewers). Governing contracts: `directive_u1_body_dynamics.md` (v5) + `directive_u1_body_dynamics_rebuild.md`. Effort xhigh. Files: `u1_body_dynamics_{sympy.py,dual.wl,compare.py,inputs.yaml,fixtures.yaml}`, `run_u1_body_dynamics.sh`.

## Fix list

- **R1 (PHYSICS-RELEVANT, priority).** The channel table is hand-typed; the "assembled action" is decorative — deleting `quantum_pressure`/`Pu_coupling`/`brane_shear_gradient` from inputs changes NOTHING (demonstrated). Make the linearized channel operator **DERIVED from the action data** (consume the action-term expressions; a deleted term must change the derived operator and fail the source-completeness gate). Include the omitted cross-couplings — `lambdaPu` (P–u; couples the two massless channels), density–phase via the drain, EOS/well curvature in `K_pp` — and either solve the coupled radial system or **PROVE by computed indicial/degree analysis** (executed in SymPy + WL, not prose) that they do not shift the leading exponents/gaps. If the coupled analysis changes any channel class, report the new verdict honestly.
- **R2.** Co-moving reductions must be an explicit computed change of variables on symbolic `f(x−Vt)` fields (chain rule executed; residual an output) — not the `(map_sign+1)·V∇` toggle.
- **R3.** RECONSTRUCTION must independently substitute the rigid embedding into the assembled action and compare to the claimed `L_eff` (no `copy.copy` of the same term list).
- **R4.** ANCESTRY/NATIVE_PADDING must actually remove the term and RE-DERIVE (classification/monomial change computed); the provenance graph must be **derived by dependency traversal of the actual expressions** (the comparator's free-symbol extraction is the model), with the forbidden-node check running on the derived graph.
- **R5.** E1–E4 endpoint trace systems are identity solves (only E5 genuine) — make them real solves of distinct declared BC/constraint systems, or downgrade the report's claim honestly.
- **R6 (harness hole, demonstrated).** The engines' "mutation did not reach its own assert" sentinel emits the same `ASSERT_FAIL:<tooth>:` prefix as a genuine assert, so the runner's awk acceptance cannot tell them apart — a warehoused solver that ignores an attack passes the TAIL_ODE ablation gate. Give the noop sentinel a distinct tag (`ASSERT_FAIL:MUTATION_NOOP:`) and make the runner treat it as gate FAILURE. Fix in sympy, .wl, and comparator.
- **R7.** The .wl G1 gate has a duplicated condition making its `NEGATIVE_COMPUTED` branch dead — fix.

## Acceptance
Full runner exits 0 (or lands an honest changed verdict); all ablation teeth per-tooth verified under the R6-fixed acceptance; report + results.yaml regenerated (no verdict literals anywhere; summary parsed from artifacts). HALT after Phase A; 10-line summary. Working notes in `_scratch/`.

## After the build
Focused re-verification: fresh adversarial spot-check on R1 (the coupled/derived channel operator) + R2–R7 fixes; arbiter re-run. Then Phase A is closed and the U1 stack goes to commit.
