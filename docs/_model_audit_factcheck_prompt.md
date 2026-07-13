# FACT-CHECK this model audit — verify the DEFINED / IMPORTED / OPEN labels are honest

This is a toy 4D superfluid analog physics program (repo root: /var/projects/toy_physics). A recurring failure mode is **importing standard-EM structure while believing it was derived**. The attached audit is meant to state, honestly, exactly what the model has *defined natively* vs *imported* vs *left open*. **Your job is a FACT-CHECK, not a redesign or a physics brainstorm.** Do not propose new mechanisms. Verify the labels against the actual repo state and flag any that are wrong or anything missing.

## Read
- **`docs/model_definition_audit.md`** — THE audit under fact-check (status vocabulary: NATIVE-DEFINED / CALIBRATED / POSTULATED / IMPORTED / OPEN).
- `docs/conceptual_foundation.md` — the model's mechanisms (source of truth for what each force IS).
- `docs/em_medium_first_generative_plan.md` §8 — the native-`P` gate result (`NATIVE_P_NO_EMERGENT_GAUSS`, verified) = no native first-class U(1).
- `docs/em_charge_attribute_requirements.md`, `docs/em_gravity_native_ontology.md` — the EM reframes; note the electric sign / roll-vs-slip coupling being open.
- Skim `software/stage1_solver/` (decisions/directives/reports) for `pathA_29` (gravity localization), `pathA_36` (light 2 photons), `pathA_38` (electric Coulomb), `pathA_39` (magnetism) to check the "defined vs calibrated vs imported" claims are accurate.

## Check, concretely (cite file:line where you correct)
1. **Over-claims:** is anything labeled **NATIVE-DEFINED** or **CALIBRATED** actually **IMPORTED** or **OPEN**? (e.g., is gravity's localization really derived; is light's 2-photon mode really native given the postulated stiffness; is the electric-sign really as open as stated or is it more/less settled than the audit says?)
2. **Under-claims:** is anything labeled **IMPORTED / OPEN** actually more defined than the audit admits? (Don't let the audit be *too* pessimistic either — that's also dishonest.)
3. **The central finding (§C, §D):** is it accurate that the EM sign/law imports all trace to **one** missing definition — the **charge↔medium coupling** ("the lock" / roll-vs-slip / what the `±w` throat-body is as a dynamical object)? Or is that overstated (are there *independent* missing pieces), or understated?
4. **The `pathA_39` magnetism claim:** is it correct that the magnetic *sign* was imported (`j=sV`) rather than derived? Is the "swirl = circulation vs moving `±w` body" ambiguity real, and is the gravitomagnetism-vs-EM-magnetism distinction stated correctly?
5. **Omissions:** any force/phenomenon or open item that belongs in the table or §E and is missing entirely?
6. **Scope honesty:** does the audit correctly separate "form derived" from "constant calibrated" from "value imported"?

## Output
- **VERDICT:** `AUDIT_ACCURATE` or `AUDIT_NEEDS_CORRECTIONS`.
- A numbered list of any label corrections (row → current label → correct label → why, with file:line).
- Any missing item to add.
- One line on whether the §D central finding (one missing definition = charge↔medium coupling) is sound.
- Concrete over polite. If it's accurate, say so plainly.
