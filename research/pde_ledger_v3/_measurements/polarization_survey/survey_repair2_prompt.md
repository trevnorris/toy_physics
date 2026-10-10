# Repair 2: the polarization survey

A second pair of reviews found the problems below in `/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md`.
One reviewer cleared it. The other found four issues, and one of them traces to the task itself: its five class
names came with no criteria to separate them. The definitions below replace them. Repair the file in place so that
everything under "Must be true" holds. Keep your original task and method. Change nothing these items do not touch.

## The classes (these replace the task's class list)
Apply them in this order. A row takes the first class whose definition it meets. State each row's conditions in its
cell.
1. **Reproduced.** A v3 step computes the measured fact from inputs that do not themselves state it.
2. **In apparent conflict.** A model source states something opposed to the measured fact, or the ledger itself
   records a departure from it. An open condition, such as an unresolved coupling, makes the conflict conditional;
   state the condition. Against a reported signal that is not established, the conflict is conditional on the signal
   too.
3. **Testable.** A model source supplies or computes the object the measurement constrains (for example, the spectrum
   whose splitting a birefringence bound limits), but no step has compared the two.
4. **Required, no mechanism yet.** The fact concerns light within the model's current scope (propagation in the brane,
   and interaction with the brane's own structures and defects), and no model source supplies the object it
   constrains. A fact that follows only from a supplied input that states it also belongs here; name that input.
5. **Not addressed.** The fact depends on physics the model has not defined, such as the response of ordinary matter,
   material interfaces, or emission processes, and no model source bears on it.

## Must be true after the repair
1. Every row's class follows from these definitions. Update the class rule paragraph, the counts and your summary.
2. E32, E33 and E34 are classed by the same rule. E32's experimental side cites E33 and E34.
3. Part 3 includes these sources, each with its status as the source states it:
   - `research/pde_ledger_v3/steps/S11_stray_longitudinal.md:52–64`;
   - `research/pde_ledger_v3/steps/O2_steady_brane_balance.md:89–92, 110–111, 115–117`;
   - `research/pde_ledger_v3/steps/S11bB_interface_assembly.md:53–54`;
   - `git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md`, lines 117–119.

   E10–E15 state the parity/chirality condition these sources bear on, and how the cited β papers describe what β
   tests. Do not resolve it. The terminology note does not present the brane's `w` displacement only as the charge
   field.
4. Part 2 includes the photon-mass bounds it lacks: the limit the Particle Data Group adopts, and dispersion
   (time-of-flight) limits from fast radio bursts. Give each with its stated assumptions, from sources you open, and
   add Part 4 rows.

## Method and bounds
- As before: every repository claim carries its command and literal output in the appendix; every external fact
  carries a citation, and you say which sources you opened.
- Do not read `~/.claude/projects/`. Modify only `POLARIZATION_SURVEY.md`. Do not commit, push or touch git state.
- Write the file and exit. Your final message: the file's path, the number of experimental entries, the
  classification counts, and one line per item above saying what you changed.
