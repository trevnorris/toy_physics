# Repair 1: the polarization survey

Two independent reviews of your survey, `/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md`,
found the problems below. Repair the file in place so that each statement under "Must be true" holds. Keep your
original task (`survey_prompt_r0.md` in the same directory) and its method. Change nothing the findings do not
touch.

## Must be true after the repair

1. **Classes.** The task's class definitions apply as written. Remove the rule that an apparent conflict requires an
   experimentally comparable model prediction.
   - "Reproduced" is used only where a step computes the fact from inputs that do not themselves state it. Where a
     fact follows from a supplied input that states it, the row names that input and is not "reproduced".
   - Apply one rule to every row, and state each row's conditions in its cell.
   - E01: the count of two follows from `D = 3`, which goes in as an input (`SUBSTRATE_REQUIREMENTS.md`, R-S1-01,
     OPEN; `steps/S10_two_transverse_photons.md:172–176`).
   - The longitudinal branch that S10 (`:249–251`) and S11 (`steps/S11_stray_longitudinal.md:15, 161–163, 173`)
     record as a departure, conditional on matter coupling to it, is its own row. A departure the ledger itself
     records, conditional on an open coupling, is an apparent conflict under the task's definition. Quote both
     sides and state the condition.
   - E10, E11 and E13–E15 use the same rule. `V3_STEP_PLAN.md:1183` states that the two polarizations are degenerate
     by symmetry in a homogeneous brane. State the conditions.
2. **The brane's `w` displacement.** Part 3 includes the plan-level sources that make the `w` displacement
   (`ξ_w = ℓh`, the `h`-branon) the charge field and the mediator of its static falloff: `V3_STEP_PLAN.md:453–461`
   and `:871–905`. Give their plan-level status. Source E02 to them. Add a short note separating the two uses of
   "transverse" in the sources: normal to the brane, and perpendicular to light's direction of travel.
3. **Polarization in S9b and S11c-d.** Part 3 includes:
   - the S9b spec's supplied polarization identification and its list of what is outside the step, cited at commit
     `ede8aa21` (`git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md`, lines 53–55 and
     106–109), because the working copy is under amendment;
   - `steps/S11c_d_profile_conditioned_scattering.md:58`.

   The rows on gravity and polarization (E22–E24) say that S9b takes polarization independence as an input, so S9b
   cannot reproduce it.
4. **The thermal count of radiating modes.** Part 2 has entries for the thermal (blackbody) count of radiating
   modes, laboratory and cosmological, from primary sources you open, with Part 4 rows.
5. **Part 1.** Give the reasons light has two polarization states: the quantum one (masslessness) and the classical
   one (in Maxwell's theory the longitudinal field is fixed by charges through Gauss's law and is not a wave). Label
   each.
6. **E21.** Quote what the cited version of record (DOI `10.1126/science.add0080`) reports. Preprint values appear
   only as the preprint's.
7. **E15.** Describe each number as the cited paper describes it. `0.20±0.08°` there is not an earlier version of the
   paper's own `β`; check what it is.

Update the classification counts and your final summary to match.

## Method and bounds
- As before: every repository claim carries its command and literal output in the appendix; every external fact
  carries a citation, and you say which sources you opened.
- Do not read `~/.claude/projects/`. Modify only `POLARIZATION_SURVEY.md`. Do not commit, push or touch git state.
- Write the file and exit. Your final message: the file's path, the number of experimental entries, the
  classification counts, and one line per item above saying what you changed.
