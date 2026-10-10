# Survey: light's polarization, what experiments show, and where the toy model stands

## Why
The user, a programmer and not a physicist, wants to understand polarization well enough to reason about it in their
toy model. In the model, light is a transverse shear wave of a brane: an ordered, finite-thickness slab, centred at
`w = 0` in a four-dimensional superfluid bulk. The user asked: what kinds of polarization exist, what experiments show
about them, and what each fact means for the model.

## What to produce
Write one markdown file, `/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md`, in four parts.

1. **What polarization is.** Plain language first, for a reader who knows programming but little physics. Then the
   precise version. Cover:
   - linear, circular and elliptical polarization, and unpolarized and partially polarized light;
   - how polarization is described (Jones vectors, Stokes parameters, the Poincaré sphere);
   - how many independent polarization states light has, and why;
   - spin angular momentum, helicity and handedness;
   - orbital angular momentum;
   - how a single photon behaves under polarization measurement;
   - polarization entanglement.

   Say what is a classical-wave fact and what is a quantum fact.
2. **What experiments show.** A table of the experimental record. For each entry give:
   - what was measured;
   - who, when, and to what precision or bound;
   - a primary citation with a DOI or URL.

   Include at least:
   - the number of polarization states (transversality), and bounds on a longitudinal mode or a photon mass;
   - the angular momentum carried by circularly polarized light;
   - orbital angular momentum beams;
   - single-photon polarization statistics;
   - Bell tests that use polarization;
   - whether the speed of light in vacuum depends on polarization or handedness (astrophysical and cosmological
     birefringence bounds, and any reported signals);
   - polarization in strong magnetic or electric fields (Faraday rotation, vacuum birefringence searches and
     observations);
   - any dependence of gravitational effects on polarization or handedness (deflection, delay, spin Hall effects);
   - polarization produced by scattering and reflection.

   Add anything else you judge load-bearing. Mark any reported signal that is not established as such, and say what
   the evidence for it is.
3. **What the repository says.** Inventory where this repository discusses light's polarization, the number of
   transverse modes, the displacement directions available to the brane (including the `w` direction), handedness,
   helicity and birefringence. Quote each source with `file:line`, and give its status as the source states it:
   derived, postulated, conditional, ruled out, stale, or other. The current authority is the v3 ledger in
   `research/pde_ledger_v3/`: its step plan `V3_STEP_PLAN.md`, its step records in `steps/`, `DEFECT_REGISTER.md` and
   `SUBSTRATE_REQUIREMENTS.md`. Older material (`research/pde_ledger_v2/`, `docs/`, `research/*/paper/`) is a lead,
   not an authority. `docs/model_map.md` is known to be stale.
4. **Classification.** One row per experimental fact from part 2. Classify it against the model as:
   - reproduced, naming the step and the source that shows it;
   - required but with no mechanism in the model yet;
   - in apparent conflict, quoting both sides;
   - not addressed by the model;
   - a place where the model could make a testable statement.

   Classify only. Do not propose new mechanisms, and do not resolve conflicts. Where the model's status is unclear,
   say so and give the source.

## Method
- Every claim about a repository file carries the command you ran (`grep`, `sed -n`, or similar) and its literal
  output. Put these in an appendix of the same file.
- Every external fact carries a citation. Prefer primary papers. Say which citations you opened and which you only
  found referenced elsewhere.
- No CAS runs are needed.
- Do not read `~/.claude/projects/`. Do not modify, create or delete any file other than the output file. Do not
  commit, push or touch git state.

## Length
Part 1: under 1,500 words. Parts 2–4: as long as the material needs, in tables where possible.

## Bounds
Write the file and exit. Your final message is a short summary: the file's path, the number of experimental entries,
and the classification counts.
