# Independent physics review — S11c-a SymPy engine, post-patch

## Artifact
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` (working tree, uncommitted patch).
A SymPy engine that derives S11c-a interface shape-derivatives from the spec and emits `PY_S11CA_*` tags,
and feeds `scripts/S11c_a_exports.py`. It was just patched for CASE-STRUCTURE (three objects).

## What you are handed
- The artifact above.
- The CLOSED spec `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — your source of truth.
- NOTHING ELSE. ⛔ Do not read, search for, or reconstruct the sibling Wolfram engine or any adjudication/
  verdicts document. Derive every object yourself from the spec; if you find another engine's output, do not
  open it. (The two engines are compared later by a separate comparator — your independence protects that.)

## What to check (verify each against the spec — ⛔ NO expected values given; a residual is a measurement,
## never a target; ⛔ do not compare to another engine)
1. **Virtual-work full physical×virtual grid (spec §4 T-d, lines 417-419: "Which virtual-displacement
   pairings occur is part of the computation").** `virtual_work_cases` (line 854) now takes a PHYSICAL dof
   and a VIRTUAL dof and keys `(branch, physical_dof, virtual_dof, representative)`. Verify the decoupling is
   genuine: the physical DOF must drive the shape-derivative and traction (face geometry, normal, measure,
   pressure), and the virtual DOF must drive the virtual displacement δ_v contracted against that traction.
   ⭐ CRUX: are the OFF-diagonal (physical≠virtual) cases REACHED BY COMPUTATION, not the diagonal copied?
   Test it — e.g. perturb the physical source alone vs the virtual source alone and confirm each moves the
   expected factor; or confirm an off-diagonal case is built from the physical DOF's traction and the virtual
   DOF's δ_v. ⛔ Do not assert which pairings "should" survive — compute and report what the engine emits.
   Also confirm the §5 controls that consume virtual work (form ablation, uniform limit, rep-invariance/
   independence) carry the same full grid, and the uniform-limit INDEPENDENT reference (built inline in
   `uniform_reference_geometry`) is a genuine independent construction that matches the engine's virtual-work
   at the (η,σ_W)→0 limit (residual → 0), not a copy of the engine's own object.
2. **BACKGROUND_STATE boundary loads (spec §2d, lines 246-272: `𝔅⁰ ≡ {…, boundary loads}`,
   `𝒮_hold⁰ ≡ {f_hold⁰, t_hold,s⁰}`).** The state now appends `support_body`(=f_hold_0),
   `support_face_plus`(=t_hold_plus_0), `support_face_minus`(=t_hold_minus_0) after the existing θ⁰/V/J/A=0
   zeros. Confirm: the zeros are NOT duplicated; the loads are SUPPLIED premise symbols (not computed); and the
   BACKGROUND_STATE dimension tuple (task_homogeneity) matches the new field count with physically correct
   dimensions (body-force density vs face traction). Check dimensional consistency with units restored.
3. **BACKGROUND_DENSITY_MAP branch axis dropped (spec §2b, lines 228-230: "Emit the two computed maps" = per
   representative, on pre-anchoring y).** `build_background_density_raw` (line 897) now keys `(representative,)`
   only. Confirm the map `Σ_E⁰ = ρ_4D·W_bg` is genuinely branch-independent (does not depend on the anchoring
   branch), so dropping the branch axis loses no information; and the uniform reference matches.
4. **Preserved physics unchanged.** The shape-derivative machinery (`build_face_source`, `shape()`,
   `build_material_face_source`, face/evolution/projection/closure builders), the §5a routes, T-h, T-c′, T-0,
   T-f must be behaviourally unchanged. Independently re-derive what you can and confirm.

## Required method — SCRIPT review, independent derivation MANDATORY
- ⛔ Write your OWN derivation script BEFORE opening the artifact; save the script AND its literal stdout to
  named absolute paths. Prose-only derivation claims are discarded.
- ⛔⛔ A FORM ABLATION IS MANDATORY: on a /tmp COPY, change the STRUCTURE of a load-bearing object (flip a sign
  and an off-diagonal in a normal/traction/measure), re-run, report the LITERAL `PY_S11CA_*` diff. A
  coefficient rescale tests arithmetic; only a FORM change tests physics.
- For every claim ask WHICH LINE COMPUTED THIS and give the line number, or report it as uncomputed.
- Report any `assert`/hard-stop that PRECEDES the value it guards. Probe for: a value verified with the
  predicate that produced it; a conclusion emitted as an unconditional literal; a case present under two keys
  but identical because it was COPIED (the virtual-work off-diagonal crux).
- Run the engine with `python3 research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
  (~3 min). ⚠ Running it regenerates `scripts/S11c_a_exports.py` — that is a DELIBERATE side effect of this
  patch (VIRTUAL_WORK + BACKGROUND_DENSITY_MAP are export objects); ablate a /tmp COPY so you never touch the
  working tree.

## Physics filter
Report a finding only if it catches a way the physics could be wrong (or a control that cannot fail, or a
copied/fake case). ⛔ Do not report "the script would be wrong on a different input".

## Sandbox
⛔ Copy the artifact to /tmp and ablate the COPY. ⛔ Never modify the working tree. Save every ablation script
and its literal stdout to named absolute paths and report them. Use a distinct /tmp subdirectory.
