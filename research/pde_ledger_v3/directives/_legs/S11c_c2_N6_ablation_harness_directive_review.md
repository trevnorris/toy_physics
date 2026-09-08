# Decision review — S11c-c2 N6 ablation-harness knife-list directive (orchestrator-written)

You are one of two independent decision legs (Codex + Grok) for an **orchestrator-written, physics-bearing build
directive**: `directives/S11c_c2_N6_ablation_harness_directive.md`. It specifies four committed ablation harnesses
(one per N6 instrument) that must make each engine's load-bearing checks BITE under a structural (FORM) corruption —
the durable per-engine cert. Working dir `/var/projects/toy_physics`. Repo read access; scratch under `/tmp` only;
⛔ do NOT modify the working tree. **A prose claim is worth nothing — ground every finding in a cited engine line
(read the file) or a small check script whose absolute path + literal stdout you report.**

This is a **physics-bearing directive → review-until-clear** (a wrong-object or wrong-site knife-list is a physics
defect, ⛔ not a one-pass decision list). Your job: catch anything that would make a harness certify the WRONG object,
fail to bite, or self-report — BEFORE the build. Executable harness-behavior tests belong to the later build review;
here, verify the DESIGN against the real engine code.

## What you are handed
- The directive: `directives/S11c_c2_N6_ablation_harness_directive.md`.
- The four engines the harnesses wrap + certify:
  - `scripts/S11c_c2_N6_diagnostic_sympy.py` (918 ln) — clearance `_measurements/S11c_c2_N6_build_clearance.md`
  - `scripts/S11c_c2_N6_reconcile_sympy.py` (300 ln) — `_measurements/S11c_c2_N6_reconcile_{build_clearance,adjudication}.md`
  - `scripts/S11c_c2_N6_covariance_sympy.py` (277 ln) — `_measurements/S11c_c2_N6_covariance_build_clearance.md`
  - `mathematica/S11c_c2_N6_mathematica_audit.wl` (1029 ln) — `_measurements/S11c_c2_N6_wl_{build_review,rebuild_leg1_review}.md`
    (⚠ the `wl_build_review` is PRE-REBUILD/stale on line numbers; the `wl_rebuild_leg1_review` matches the current tree)
  - siblings the SymPy engines import: `scripts/S11c_a_interface_geometry_sympy_audit.py`, `scripts/S11c_b_brane_operator_sympy_audit.py`
- The template being imitated: `scripts/S11c_b_carrier_ablation_harness_sympy.py`, `mathematica/S11c_b_carrier_ablation_harness.wl`,
  directive `directives/S11c_b_carrier_ablation_harness_directive.md`.
- ⛔ You are NOT handed an expected-cone/adjudication key — you must confirm each knife hits the right object by
  READING the engine, not by being told the answer.

## Check, per knife in each of the 4 harness sections (cite lines)
1. **Site exists + is the right construction site.** Does the cited line in the named engine (or sibling) hold the
   construction the directive claims? Is it the object's CONSTRUCTION, not a downstream read? Confirm the exact-string
   patch would be unique (occurs once) at that site.
2. **FORM vs COEFFICIENT.** Is each knife labeled FORM actually a structural change that leaves the family (not a
   scalar rescale)? Is each COEFFICIENT/RESCALE companion correctly labeled? A FORM knife mislabeled (or a rescale
   passed off as FORM) is a finding — only a FORM change tests physics.
3. **Right certified object / correct cone.** Does the knife actually reach the engine's LOAD-BEARING certified
   object (`R_N6`, `SPLIT_CHECK`, `R_cov`, the carrier/source bridges, the guards) — the S11c-b failure mode was a
   knife on a P-INDEPENDENT object that could not move the carrier. Verify by reading the data flow: does K_ewsign /
   K_slotdrop actually reach `R_N6`? does K_rank / K_circular actually reach `R_cov`? does K_normal reach the carrier
   bridge and NOT the source bridge?
4. **The reconcile K_normal cross-file hazard.** Confirm the directive sites K_normal at `normal_exact` in
   `S11c_a_interface_geometry_sympy_audit.py` and ⛔ explicitly AVOIDS `material_inverse_transpose` (S11c-a:694,
   confirmed false-negative). Confirm the subprocess-isolation + patched-sibling-import mechanism is sound.
5. **Self-tests genuine.** For each knife: is a NO-OP (identity patch ⇒ diff identically 0) specified? a DEAD /
   one-sided separation? a RESCALE contrast? Are the DEAD channel-separation claims (which objects a knife must
   leave byte-identical) CORRECT per the engine's data flow — or does the directive assert a separation the code
   does not have?
6. **PRINT-not-PASS.** Does the directive forbid any PASS/FAIL/verdict/expected-value payload and require
   baseline/corrupted/diff triples? Is the builder correctly bounded (build→run→report→stop, no self-review, no
   calling other AI)?
7. **Leakage.** Does the directive leak an expected residual VALUE or an acceptance criterion a builder could iterate
   toward? (Knife knob input values like `materialNormalKnife 0→1` are inputs, not leaks; an expected residual
   tally or "must vanish" would be.)
8. **Certified-object completeness.** Does each harness certify its engine's actual load-bearing claim (the WL
   engine's blind reproduction of `R_N6`/`SPLIT_CHECK`/`R_cov` under DISJOINT primes; the reconcile's two-channel
   separation; the covariance's non-circularity + rank-2 Φ; the diagnostic's `R_N6` + control-independence)? Flag any
   engine whose real load-bearing control the directive omits or mis-targets.

## Physics filter
Report a finding only if it would make a harness certify the wrong thing, fail to bite, self-report, or leak. ⛔ Do
not report style, and ⛔ do not report "it would be wrong on a different engine".

## Output
Per-item findings with cited engine lines (and any check-script path + literal stdout), then a final line:
**SOUND** (the knife-list is correct and buildable as written) or **NOT-SOUND** (exact section + knife + the fix for
every blocking issue).
