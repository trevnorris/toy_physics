# Plan review — S11c-a comparator SHARED-SCHEMA RE-EMIT roadmap

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_comparator_reemit_plan.md`

This is the orchestrator's roadmap for building the S11c-a cross-engine comparator by **re-emitting both
engines to a shared §7 emit schema** (instead of reconciling incompatible schemas in the comparator). Two
prior directive legs rejected the comparator-side approach as unsafe; this plan is the chosen replacement.
Your job: is the plan **accurate, sound, and complete enough to execute after a context compaction** — or are
there gaps/risks that would waste a cycle if we discover them mid-execution?

## What you are handed
- The plan (above). Its measurements twin: `_measurements/S11c_a_comparator_reemit_plan.md`.
- The two real engines + tag streams:
  - SymPy engine `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` (`9b6438fa`); its
    tag stream regenerate with `python3 …` (⚠ slow; rewrites `S11c_a_exports.py` in Dummy-index counters —
    `git checkout` after) or inspect the emit path.
  - Wolfram engine `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`
    (`ddecdbc2`); committed tag stream `mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`.
- The closed spec `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`) — §7 is the tag/
  emit grammar the plan proposes to amend.
- The frozen comparator precedent `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`.

## Check, concretely
1. **Divergence map accuracy.** Spot-check the plan's "divergence map" against the real `.out`/source: the
   case-key encodings (PY positional tuple w/ integer faces vs WL pipe-string; the extra axes DIRECTION/FIELD/
   ORIGIN/VIRTUAL_DOF; FACE_MAP `PLUS`/`MINUS`), the field envelope (PY flat `VALUE` vs WL nested
   `SHAPE_DERIVATIVE`/`EXACT_SOURCE`), and the CAS reps (2-arg window both sides; `Inactive[Integrate]` vs
   `Integral`; nested derivatives; integer Association keys). Is anything wrong or missing?
2. **Is the shared §7 schema actually achievable by an EMIT-LAYER change?** The plan claims both engines can
   be patched at the serialization layer only (physics untouched), verified by byte-comparing payload VALUES
   modulo relabeling. Test that claim against reality: (a) can WL flatten its `SHAPE_DERIVATIVE`/`EXACT_SOURCE`
   nesting so the computed result sits at `EXPRESSION` at case level WITHOUT recomputing physics? (b) does the
   PY engine already HOLD every field the shared schema needs (e.g. is there a stage the shared `EXPRESSION`
   maps to), or would matching WL force PY to expose an intermediate it does not currently compute? (c) are
   there objects where the two engines' case-DECOMPOSITION (which axes) genuinely differs, so no relabeling
   makes them match (a real finding, not serialization)? Name any object where the emit-only claim breaks.
3. **Is the mechanical physics-preservation gate sufficient?** "Extract each object's payload VALUE
   before/after and byte-compare modulo relabeling." Is that a well-defined, runnable check given the nesting/
   key changes — or does the relabeling itself make "the value" ambiguous to extract? Suggest the concrete
   form of the gate if the plan's is underspecified.
3b. **Face/direction convention.** The plan says pin ONE face convention. Given PY uses integer `1/-1` for
   face AND `1/2/3` for direction, and WL uses `FACE_PLUS/FACE_MINUS` and `DIRECTION_n` and `PLUS/MINUS`
   (FACE_MAP), is the proposed single convention actually consistent across every object, or will some object
   collide?
4. **Workflow soundness + ordering.** Is the sequence right (schema enumeration → §7 amendment (2 legs) →
   per-engine emit patch (2 legs each) → re-run → light re-review (2 legs each) → trivial comparator (2 legs)
   → freeze+run)? Is reopening the closed spec `2926c71c` for §7 the correct move, or is there a lighter way?
   Does "light re-review" adequately cover the risk that an emit patch silently drops/reorders a payload?
5. **Completeness / gaps.** Anything the plan omits that will surface mid-execution: local (`_LOCAL_`) tags,
   the `DIMENSIONS` whole-run object, the `HOMOGENEITY` controls, dimension vectors, the export/MANIFEST
   coupling, or a shape neither leg nor the plan mentions.
6. **Rule-2 / rule-5.** Does the plan (or its twin) assert any artifact fact without a command behind it? Does
   it leak an expected cross-engine RESULT (which objects will agree)? It must not.

## Output
A short verdict: is the plan sound and safe to execute post-compact as-is, or the exact gaps/edits to fix
first. Quote plan line numbers and cite `.out`/source evidence for any correction. Focus on point 2 — whether
the emit-only re-emit is actually achievable — since that is the load-bearing assumption of the whole plan.
