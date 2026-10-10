# Repair 3: the polarization survey (you are a new author)

## The document
`/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md` is a survey of light's polarization for the
user, a programmer and not a physicist. In the user's toy model, light is a transverse shear wave of a brane: an
ordered, finite-thickness slab centred at `w = 0` in a four-dimensional superfluid bulk. The survey explains
polarization (Part 1), lists what experiments show (Part 2), inventories what this repository says (Part 3), and
classifies each experimental fact against the model (Part 4). Another author wrote it from the task in
`survey_prompt_r0.md` (same directory). Read that task for the method; it still applies.

Three review rounds have checked it. Every external value and repository quotation now checks out. Five items
remain. Fix exactly these, and change nothing else. Do not add experiments or sources beyond the ones named below.

## The classes in force
These replace the task's class list. Apply them in this order; a row takes the first class whose definition it
meets, except as stated in rule 4. State each row's conditions in its cell.
1. **Reproduced.** A v3 step computes the measured fact from inputs that do not themselves state it.
2. **In apparent conflict.** A model source states something opposed to the measured fact, or the ledger itself
   records a departure from it. An open condition, such as an unresolved coupling, makes the conflict conditional;
   state the condition. Against a reported signal that is not established, the conflict is conditional on the signal
   too.
3. **Testable.** A model source supplies or computes the object the measurement constrains (for example, the spectrum
   whose splitting a birefringence bound limits), but no step has compared the two.
4. **Required, no mechanism yet.** The fact concerns light within the model's current scope (propagation in the brane,
   and interaction with the brane's own structures and defects), and no model source supplies the object it
   constrains. **A fact that follows only from a supplied input that states it belongs here, and this takes
   precedence over rule 3**; name that input.
5. **Not addressed.** The fact depends on physics the model has not defined, such as the response of ordinary matter,
   material interfaces, or emission processes, and no model source bears on it.

## Must be true after the repair
1. **The class rule.** The survey's class-rule paragraph states these definitions as written, including rule 4's
   precedence over rule 3. Quote rule 5 as defined ("no model source bears on it"). Update the counts and the summary.
2. **E12.** E12 (anisotropic CMB polarization rotation) is classed by the same reasoning as E10, E11 and E13–E15,
   which rest on the same model sources and need the same CMB source and statistics. State its conditions.
3. **S11 on parity.** Part 3 and E10–E15 state what
   `research/pde_ledger_v3/steps/S11_stray_longitudinal.md:52–64` finds at `D = 3`, as opposed to `D = 2` and `D = 4`,
   for its class of terms. Locate the open parity/chirality condition outside that class, citing the sources the
   survey already gives. Do not resolve it.
4. **E01 and E22.** Their stated reasons follow the definitions above, including rule 4's precedence, and do not
   contradict `research/pde_ledger_v3/steps/S10_two_transverse_photons.md:73–76`.
5. **Additions, and only these:**
   - Part 2 and Part 4: the cosmological count of relativistic species (`N_eff`) from the Planck 2018 cosmological
     parameters, with its assumptions, classed under the same rule as E32–E34. In Part 2, mention the nonzero
     fast-radio-burst photon-mass value in the Particle Data Group's photon listing, marked as not established, as the
     listing marks it.
   - Part 3: `research/pde_ledger_v3/steps/S11bB_interface_assembly.md:74–76` and
     `research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:425` (entry R-S8-06), each with its status as the source states it.

## Method and bounds
- Every repository claim carries its command and literal output in the appendix. Every external fact carries a
  citation, and you say which sources you opened.
- Do not read `~/.claude/projects/`. Do not read the review reports or dispositions in the survey's directory. Modify
  only `POLARIZATION_SURVEY.md`. Do not commit, push or touch git state.
- Write the file and exit. Your final message: the file's path, the number of experimental entries, the
  classification counts, and one line per item above saying what you changed.
