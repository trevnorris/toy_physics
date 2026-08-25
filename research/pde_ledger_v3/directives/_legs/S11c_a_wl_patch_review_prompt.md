# Independent physics review — S11c-a blind WL engine, post-patch

## Artifact
`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (working tree, uncommitted
patch). A blind Wolfram engine that re-derives S11c-a interface shape-derivatives from the spec and emits
`WL_S11CA_*` tags. It was just patched for SCHEMA + CONTROL-COVERAGE only.

## What you are handed
- The artifact above.
- The CLOSED spec `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — your source of truth.
- NOTHING ELSE. ⛔ Do not seek out, read, or reconstruct any sibling engine (SymPy) or any adjudication/verdicts
  document. This engine is BLIND and so is your review: derive every object yourself from the spec. If you find
  another engine's output, ⛔ do not open it.

## What to check (verify each against the spec — ⛔ I am giving you NO expected values; a residual is a
## measurement, never a target)
The patch made three changes plus a preservation claim. Verify each is physically correct per the spec:
1. **Projection density representative (spec §1b lines 74-83, §2b 228-230, §3a 312-318, §4 397-400/424-429).**
   The five T-f projection objects (`PROJECTION_SHAPE_DERIV`, `_STATIC_OPERAND`, `_DYNAMIC_OPERAND`,
   `_RESIDUAL`, `_TERM_ORIGINS`) now carry a density-representative axis `{RHO4_CONSTANT, RHOBR_CONSTANT}`.
   The projection integrates the 4D current `j = ρ_4D v_bulk`; `ρ_4D` is representative-dependent. ⭐⭐ CRUX:
   confirm each representative's payload is **REACHED BY COMPUTATION from that representative's own density
   source, not a hand-copied duplicate**. Test it: perturb ONE representative's density source only and confirm
   only that representative's emitted projection payload moves. A copied second representative would be
   byte-identical by construction, or both would move under a single perturbation. Report the literal before/
   after for the perturbed and the untouched representative.
2. **Kinematic/flux carry NO representative axis (spec §3b 349-353).** `KINEMATIC_BALANCE` and `RELATIVE_FLUX`
   are built through `J_s ≡ ρ_m(v_bulk − v_face)·n̂`. Decide from the spec whether these objects depend on the
   density REPRESENTATIVE (ρ_4D/ρ_br) at all. Confirm the physics of dropping the axis: is `ρ_m` a bound
   constant distinct from the representatives? Would emitting a representative axis here be a spurious label?
3. **Control coverage (spec §5b 487-492 "every T-object", §5c 494-499 "each S11c-a object").** `CONTROL_FORM_*`
   and `UNIFORM_LIMIT_*` now cover `EVOLUTION_TERM_ORIGINS` and the four projection operands. Confirm these
   added controls are GENUINE: the form ablation independently recomputes an ablated operand from a source-form
   change (⛔ not a coefficient rescale, ⛔ not a reuse of the already-emitted payload), and its residual can
   actually move. Ablate one and show it moving.
4. **Preserved physics unchanged.** The §5a rep-invariance / independence material-coordinate routes, T-h
   sourced mass balance, T-c′ two-operand kinematic, T-0 grading, and T-f window-inside-integral must be
   physically unchanged by the patch. Independently re-derive at least the ones you can and confirm.
5. **§5a REP_INVARIANCE material route — verify it genuinely computes and the control CAN FAIL.** For every
   compared T-object — especially `RELATIVE_FLUX`, `TRACTION`, `CLOSURE_SHAPE_DERIV` — the
   `REP_INVARIANCE_MATERIAL_OPERAND` must be a computed material-coordinate (route-2) expression, ⛔ NOT an
   inert/unevaluated symbol. Grep the WHOLE emitted stream for any inert `materialShape[` token — there must be
   NONE. Then confirm the control can fail: a one-sided corruption of the Eulerian route ALONE must move
   `REP_INVARIANCE_RESIDUAL`, and a one-sided corruption of the material route ALONE must move it (§5a lines
   466-474 give the corruptions). ⛔ Do not state what the residual equals; show it moves and PRINT it. Also
   note whether the `RELATIVE_FLUX` representative axis inside these §5a control tags is a harmless redundant
   label or a spurious test, and say which.

## Required method — SCRIPT review, independent derivation MANDATORY
- ⛔ **Write your OWN derivation script BEFORE opening the artifact**, and save BOTH the script and its literal
  stdout to named absolute paths. Without these, your derivation claims are discarded (a prose re-derivation is
  worth nothing).
- ⛔⛔ **A FORM ABLATION IS MANDATORY, not optional — it is the only control that has ever caught the worst
  defect.** Change the STRUCTURE of a load-bearing object (flip a sign AND an off-diagonal, collapse two
  independent symbols into one, break a normal/Jacobian/traction form) on a /tmp COPY, re-run, and report the
  LITERAL tag diff. A coefficient rescale tests arithmetic; only a FORM change tests physics.
- For every claim, ask **WHICH LINE COMPUTED THIS** and give the line number, or report it as uncomputed.
- Report any `assert`/hard-stop that PRECEDES the value it guards (it hides a form defect behind PASS-or-crash).
- Probe for: a value verified using the predicate that produced it; a conclusion emitted as an unconditional
  literal; a payload present under two representative keys but identical because it was COPIED (change 1 crux).

## Physics filter
Report a finding only if it catches a way the physics could be wrong (or a control that cannot fail, or a
representative axis that is fake/copied). ⛔ Do not report "the script would be wrong on a different input".

## Ablation sandbox + Mathematica operational constraints (BOTH legs receive these identically)
- ⛔ Copy the artifact to `/tmp` and ablate the COPY. ⛔ Never modify the working tree.
- ⛔ Wrap EVERY kernel run in `timeout 1800`. A 1800s hit is a FAILED ablation — report it and move on.
  ⛔ NEVER raise the timeout. ⛔ Never run more than one kernel at a time (the licence has TWO seats).
- ⚠ The full run is ~15-20 min and peaks ~9-10 GB RSS. Watch memory; a run over ~15 GB or a leftover kernel
  after your job ends is itself a finding — kill it by PID and note it.
- ⭐ Save every ablation script AND its literal stdout to named absolute paths, and report those paths.
- Use a distinct `/tmp` subdirectory for your scripts so the two legs never share files.
