# Review (round 2) — S11c-c2 N6 ablation-harness knife-list directive (re-authored by a code engine)

You are one of two independent review legs (a fresh Claude agent + Grok) for a **physics-bearing build directive**
`directives/S11c_c2_N6_ablation_harness_directive.md` that specifies four committed ablation harnesses (one per N6
engine). It was just re-authored by a code engine to fix a prior round of review findings, and it now contains
concrete code PATCHES (Python fragments and a WL edit) the builder will apply verbatim. Working dir
`/var/projects/toy_physics`. Repo read access; scratch under `/tmp` only; ⛔ do NOT modify the working tree.
**Ground every finding in a cited engine line (read the file) or a small check script whose absolute path + literal
stdout you report — a prose claim is discounted.** This is review-until-clear: catch anything that would make a
harness certify the WRONG object, fail to bite, self-report, leak an expected value, or CRASH — before the build.

## What you are handed
- The directive: `directives/S11c_c2_N6_ablation_harness_directive.md`.
- The four engines + siblings it wraps/patches:
  - `scripts/S11c_c2_N6_diagnostic_sympy.py`, `scripts/S11c_c2_N6_reconcile_sympy.py`,
    `scripts/S11c_c2_N6_covariance_sympy.py`, `mathematica/S11c_c2_N6_mathematica_audit.wl`
  - siblings: `scripts/S11c_a_interface_geometry_sympy_audit.py`, `scripts/S11c_b_brane_operator_sympy_audit.py`
- Clearance records (the certified per-engine results being made durable):
  `_measurements/S11c_c2_N6_{build_clearance,reconcile_build_clearance,covariance_build_clearance,wl_rebuild_leg1_review}.md`.
- Template: `scripts/S11c_b_carrier_ablation_harness_sympy.py`, `directives/S11c_b_carrier_ablation_harness_directive.md`.
- ⛔ You are NOT handed an expected-cone/answer key. Confirm each knife hits the right object by READING the engine.

## Check, per knife (H1 WL: K_carrier, K_junk, K_split_route; H2 covariance: K_junk, K_circular, K_rank; H3 diagnostic: K_EW_rowdrop, K_slotdrop; H4 reconcile: K_normal, K_source_route, K_operand_swap)
1. **Site exists + unique + is the object's construction.** Does the cited line/scope hold the exact old fragment
   the directive quotes, occurring exactly once in that scope? (Line numbers are anchors; the fragment + scope are
   authoritative — verify the fragment.)
2. **The concrete PATCH is correct code.** For every Python/WL replacement fragment the directive gives: do the
   symbols/APIs it uses actually exist (e.g. `a.grad_W`, `a.dot`, `n.replace`, `a.build_material_face_source`,
   `a._FACE_CACHE`, `globals()['pressure_coefficients']`, `paths[atom][1]`, `total_derivative(..., background_depth=)`)?
   Would the patched engine PARSE and RUN, or does it crash / raise (e.g. the `K_EW_rowdrop` `rows.pop('E_W')` vs
   `build_increment`'s fixed output-domain completion; the `K_rank` coverage check `if uncovered: raise`)?
3. **FORM vs COEFFICIENT is correctly classified.** Is each FORM knife a genuine STRUCTURAL change (leaves the
   family), and each ×2 / sign-flip companion correctly a COEFFICIENT? (`K_EW_rowdrop` must be a real row-family
   removal, not a scalar; `K_operand_swap` must be the `es−es` collapse, not the ms↔es ×−1 swap; `K_split_route`
   must actually break the affine-split identity.)
4. **Right certified object / correct one-sided cone.** Trace the data flow: does each FORM knife actually reach the
   engine's load-bearing certified object (`R_N6`, `SPLIT_CHECK`, `R_cov`, the carrier/source bridges) and leave its
   named DEAD objects genuinely independent? Is any DEAD object actually dependent on the knife (a wrong DEAD claim)?
   The prior round caught: K_carrier's DEAD wrongly including `R_COV_INCREMENT`/`ACTUAL_CONTROL_PARAMETERS`; a knife
   that is 0 at the pinned case (the dropped advection knife) — check the pinned case `LAB_HELD × RHOBR_CONSTANT`
   actually exercises each knife.
5. **Reconcile K_normal is carrier-only.** Confirm the scoped `build_material_face_source` wrapper (restored in
   `finally`) corrupts only the carrier path and does NOT move the source bridge via `face_velocity_raw`
   (S11c-a:890-891) — i.e. `build_material_velocity` (rec:82-92) has already supplied `m_v` before the wrapper is
   installed. Confirm it does not globally patch S11c-a, does not touch the Eulerian normal (a:850) or
   `material_inverse_transpose` (a:694).
6. **PRINT-not-PASS + no leak.** No PASS/FAIL/verdict/expected-value/"must vanish"/"non-trivial" anywhere in the
   builder-facing directive; only site + FORM + "print these named objects". Builder correctly bounded
   (build→run→report→stop, no other AI, no commit).
7. **Certified-object completeness.** Does each harness certify its engine's real load-bearing claim (WL blind
   reproduction of R_N6/SPLIT_CHECK/R_cov under its own primes; reconcile two-channel separation; covariance
   non-circularity + rank-2 Φ incl. FROZEN_PHI/PHI_DOMAIN_CENSUS; diagnostic R_N6 + the TILT/N4 control triples)?

## Physics filter
Report a finding only if it would make a harness certify the wrong thing, fail to bite, self-report, leak, or crash.
⛔ Not style; ⛔ not "wrong on a different engine".

## Output
Per-item findings with cited engine lines (+ any check-script path + literal stdout), then a final line:
**SOUND** (build-clear as written) or **NOT-SOUND** (exact section + knife + the fix for every blocking issue).
