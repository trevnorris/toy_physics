# Independent review — S11c-a cross-engine comparator BUILD DIRECTIVE rev 3 (decision list, pre-build)

## Your task
Review an orchestrator-written build directive BEFORE a builder is launched. It specifies a measurement
instrument (a cross-engine comparator) diffing two independent CAS engines across 39 tag families. A defect
here becomes a defect in the instrument that judges both engines — physics-bearing. Find what is wrong,
ambiguous, or missing; a review that finds nothing is weak evidence. Point at the exact artifact/line; cite
payload evidence; do not offer prose-only fixes.

This is **rev 3**. rev-1 (~14 defects) and rev-2 (~8) both failed independent-leg review; rev-3 folds those
verified corrections plus a two-advisor consult. The two known-hard failure modes to re-check hardest:
(i) case-keying that joins ~0 in practice, and (ii) a perturbation-current fold that manufactures FALSE
cross-engine agreement. Do NOT assume rev-3 fixed them — verify against the real payloads. Do NOT re-report a
rev-2 defect that rev-3 genuinely fixed.

## Read these (sources of truth FIRST; form your own view before judging)
- Directive: `research/pde_ledger_v3/directives/S11c_a_comparator_build_directive.md`
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — §1a, §1b (4D current), §3a, §3c, §5–§8.
- Committed transcripts (the real payloads): PY
  `research/pde_ledger_v3/scripts/out/S11c_a_interface_geometry_sympy_audit.out` (47 tags); WL
  `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (40 tags, multi-line).
- Reference implementations the directive re-expresses: `~/.s11_build/S11c_a_reconcile_fixed.py`,
  `S11c_a_cov_all.py`, `S11c_a_run_projection.py` (projection Dwin/integral bridge), `S11c_a_scratch_loader.py`
  (`load_py`/`load_wl`/`split_top`/`arrow`/`py_cases`/`canon_key`).
- Frozen contract for reuse: `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`.

## The questions that matter most
1. **CORRECTION 1 — axis-typed keying.** Will the specified keying actually JOIN all 39 families against the
   real keys? Check the hard cases against the payloads: VIRTUAL_WORK (PY positional `(BRANCH,dof,vdof,DENSITY)`
   vs WL `DOF_*`/`VIRTUAL_DOF_*`; order/role must be preserved, NOT a frozenset), CONTROL_FORM (PY carries a
   face `±1` AND a direction `±1`; must not conflate), FACE_SHIFT (PY 16 aggregate cases vs WL 80 `FIELD_*`),
   HOMOGENEITY (object-name rename). Is the PY positional→typed decode specified well enough to build, or will
   the builder guess? Does the accounting `{join, py_only, wl_only, duplicate_key, parse_failed}` make a
   0-join or silent-collapse impossible to hide?
2. **The perturbation-current fold (the #3 trap).** Is the hardened treatment sufficient? It must be a literal
   AST-head/symbol-base rename preserving every arg/arity, explicitly NOT the FIELD `AppliedUndef→bare Symbol`
   pattern; it emits (A) a raw residual and (B) a separately-labeled held-context diagnostic that must NOT be
   in the closed map and whose vanishing must NOT be counted as agreement. Verify against payloads: PY normal
   current is `Function('delta_j_bulk_4')(w)` on projection but bare `Symbol('delta_j_bulk_4')` on FACE_SHIFT;
   WL is `currentWPerturbation[x1,x2,x3,{w,t}]`; in-plane PY `delta_j_bulk_{1,2,3}` bare vs WL
   `currentXPerturbation{1,2,3}` fields. Can this fold, as specified, (a) never strip an arg to force
   agreement, and (b) never manufacture a false disagreement? Is the (B) "held-context" construction sound, or
   does it smuggle the inert-spectator assumption the spec has not authorized?
3. **CORRECTION 2 — integral canonicalization.** Is "retain bounds + binder identity; bound-equality is part
   of the compared object; pull out only integration-variable-independent factors; combine only over identical
   limits" correct and sufficient to (a) cancel the genuine in-plane-divergence linearity identity and (b) NOT
   silently pull a w-dependent factor out? Is the window `Dwin`/`O_window` bridge (fold 10) correctly scoped
   to the window AST so it cannot absorb a current/density factor?
4. **Per-family extractor leaves.** Cross the stated leaf paths against the transcripts: simple shape-derivs
   (PY `VALUE` ε¹ / WL `SHAPE_DERIVATIVE.EXPRESSION`); triples take PY `VALUE[2]`; KINEMATIC
   `OPERAND_A/B_SHAPE_DERIVATIVE`; VIRTUAL_CONSTRAINT `NORMALIZED_VIRTUAL_MASS_VARIATION`; VIRTUAL_WORK
   `UPPER/LOWER`; projection PY `VALUE`/WL `SHAPE_DERIVATIVE.EXPRESSION`|`EXPRESSION`; origins partition-sums
   (evolution named-map vs WL cases; projection DYNAMIC+STATIC two sums); FACE_SHIFT flatten
   `VALUE[0]→PRESSURE, [1][0:4]→VELOCITY, [2]→DENSITY, [3][0:4]→CURRENT`. Any wrong leaf?
5. **Closure & smuggles.** Are MUMAP (→ per-branch registry) and the CONORMAL Taylor phantom correctly
   excluded? Is the map genuinely closed — could a builder "complete" it (add a current arg-strip, a mu_theta
   collapse, an object coercion, a waveOrder omission)?
6. **Rule-5 leak.** Any zero/nonzero target or expected value for any MEASURED family? Flag it.
7. **Script discipline & reuse.** Prints A,B,A−B before any guard; asserts nothing; exit 0 on disagreement;
   emission not value-conditional; synthetic-fixture self-tests only; reuses S11b parse/residual utilities but
   NOT its `compare_records`/`render_*`/`main` run path (which would restore classify-first +
   `FINAL_OPERATIONAL_STATUS`).
8. **Anything missing** that would breed a build defect.

## Output
A ranked list of concrete defects (most severe first), each with the directive section/line, what is wrong,
and the spec/payload evidence. Rank a fold that could manufacture false cross-engine agreement, or a keying
that joins ~0, above coverage/cosmetic issues. If you judge an axis sound, say so and state exactly what you
checked (with the payload token).
