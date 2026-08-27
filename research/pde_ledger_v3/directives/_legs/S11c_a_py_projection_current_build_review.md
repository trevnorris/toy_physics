# Independent physics review — S11c-a PY projection-current fix (SCRIPT)

## Artifact
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` (working tree, modified) — the
functions `projection_terms` (~lines 1114–1171) and `uniform_projection_reference` (~lines 1660–1682), which
emit `PY_S11CA_PROJECTION_{SHAPE_DERIV,STATIC_OPERAND,DYNAMIC_OPERAND,RESIDUAL,TERM_ORIGINS}`.

## What to check (the physics claim)
A fix was applied so the dynamic-window projection carries the perturbation current's **normal variation
`∂_w δj_w`**, which was previously frozen (the normal term used a `w`-constant symbol `j_bulk[3]`). Verify two
things, by computation:
1. **The normal variation genuinely enters** the projected conservation law now (it is no longer identically
   independent of the current's `w`-profile).
2. **It is not double-counted.** The projection is formed by integrating `∂_tρ + ∇₄·δj` against `Ω` with
   integration by parts in `w`; the normal divergence `∂_w δj_w` must contribute **exactly once** (through the
   single post-IBP window-normal term `-∫ δj_w·∂_wΩ`), NOT also through a separate `+∫ Ω·∂_w δj_w` channel.
   Adding both would count the same contribution twice and silently corrupt the conservation projection.

## What you are handed
- The artifact above (read the raw file; the fix is in the working tree, uncommitted).
- The spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — §1b (the current law, conservation,
  IBP-in-`w`), §3c (background current zero; traced-field scope), §4 T-f, §5c (the uniform-limit reference),
  §6 (script discipline), §7 (tag grammar).
- The current declarations (~lines 176–181) and the trace construction (~line 651).

## Required method — DERIVE INDEPENDENTLY, then ABLATE (this is a script)
Write your own derivation **before** trusting the artifact, and save the script AND its literal stdout to
named absolute paths; a prose-only derivation is discarded.
- From §1b, derive what the projected normal-current contribution should be for a **generic `w`-dependent**
  perturbation current `δj_w(w)`, with the supplied dynamic window, IBP in `w`. State the boundary property
  you use.
- **FORM ABLATION (mandatory).** Import `projection_terms` (and `uniform_projection_reference`) and exercise
  them; do not just read. Probe the current's `w`-dependence:
  - Replace the projection's perturbation current by a **`w`-dependent probe** whose normal profile you
    control, recompute the projection operand, and confirm the normal-current contribution **tracks the
    probe's normal variation** (a genuine `∂_w δj_w` dependence). Then collapse the probe to a `w`-constant
    and confirm the contribution reduces to the old frozen form. Report the literal diff both ways.
  - **Double-count test:** count how many distinct terms in the assembled operand carry `∂_w δj_w` (or,
    equivalently, the window-normal derivative against the current). Confirm it is exactly one. A genuinely
    correct fix has the normal divergence entering once; a double-count shows up as two window-normal-current
    terms or as a residual that is twice the expected size. Show the operand's term structure, not a verdict.
- Confirm the in-plane divergence and time terms are unchanged, the five `PROJECTION_*` tags are all still
  emitted, and each still has its eight `(branch, dof, representative)` cases (no face axis introduced).
- Confirm `uniform_projection_reference` received the same correction and remains an *independently* built
  reference (not a copy of `projection_terms`).

## Physics filter
Report a finding only if it catches a way the physics could be wrong (normal variation absent, double-counted,
mis-placed, a wrong boundary assumption, a broken tag/case set, or the uniform reference left frozen). Do not
report "would be wrong on a different input."

## Ablation sandbox / ops
- Copy the artifact to `/tmp` and ablate the **copy**; never modify the working tree.
- The full engine audit takes ~5 minutes; do **not** run it repeatedly. Exercise `projection_terms` /
  `uniform_projection_reference` in **isolation** (import and call them with small inputs). Wrap any run that
  could hang in `timeout 600`; a 600 s hit is a failed ablation — report and move on.
- ⛔ Running the SymPy engine's full audit rewrites the committed `S11c_a_exports.py`; if you must run it, do
  it on a `/tmp` copy, never in the working tree.
- Save every ablation script and its literal stdout to named absolute paths and report those paths.
