# BUILD — the S11b ledger card (TeX)

⛔⛔ **Write ONE card for S11b: the linear brane–bulk interface coupling law.**

## Deliverable

`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex`

Then add its `\input` to the light sector part file alongside the existing step cards, and ⭐ **build the
paper and report the page count and that it builds clean.**

## Sources of truth — ⭐ these two, and nothing else

- `research/pde_ledger_v3/steps/S11bA_interface_response.md`
- `research/pde_ledger_v3/steps/S11bB_interface_assembly.md`

⚠ **The second overturns a headline result of the first.** ⛔ Read both fully before writing a line.

⭐ Match the structure, voice and macro usage of the existing cards:
`paper/steps/S9_light_requires_shear.tex`, `S10_two_transverse_photons.tex`,
`S11_stray_longitudinal.tex`. ⛔ Do not invent a new card format.

## ⛔⛔ THE FRAMING RULE — this is the whole reason the card is written this way

**S11b is ONE step.** It was executed as two sub-steps (A then B) for tractability. ⛔ **The card must not
be structured as "part A, then part B", must not narrate the split, and must not read as one result
followed by its retraction.** ⭐ Write it as a single derivation of the interface coupling law.

⛔⛔ **The ledger keeps the SURVIVING solution plus departures.** Therefore:

- ⭐ **The surviving result is the interface law**: the bulk's face response, the permeable closure with its
  affinity, the assembled system, the admissibility and reciprocity conditions, the transverse null, and
  the breathing-mode stability boundary.
- ⛔⛔ **The frequency-dependent leak — "the velocity coupling dissipates iff `ωτ ≠ 0`", "slow disturbances
  see the leak, fast ones do not" — is DEAD. ⛔ It must NOT appear anywhere as a live result.** ⭐ It belongs
  in the card's **departures** section: derived, then **ruled out by thermodynamic admissibility**, with the
  reason stated in one or two sentences. ⚠ The algebra was correct; the physical claim did not survive an
  admissibility analysis that had not been performed when it was recorded.
- ⭐ A reader who stops halfway through this card must not be left holding a result the model forbids.

## What the card must contain

⭐ Take the values from the two step records; ⛔ do not re-derive and ⛔ do not restate anything they do not
say.

1. **What the object is** and why the step exists — S11 reduced three open questions to this one unbuilt law.
2. **The bulk's response to moving faces**, its regimes, and the added inertia. ⚠ Note the added inertia is
   `ρ_m/(2α)` **per face on the thickness coordinate** — ⛔ a pre-registered prediction of `ρ_m/α` was
   refuted.
3. **The permeable closure and its affinity**, with the affinity's normalization stated.
4. **The assembly** — the virtual-displacement rule, the balance-law route, and that the symmetry basis has
   **ten** independent quadratic invariants where the specification carried five. ⭐ Name the five new
   coefficients.
5. ⭐⭐ **Admissibility and reciprocity** — the conditions, and that reciprocity was reached by two
   independent routes.
6. **The transverse null** on a uniform background. ⛔ State explicitly that this does **not** settle whether
   confinement is unconditional.
7. **The breathing-mode stability boundary** `K₀`, and ⚠ **explicitly that the growing root is not an
   energy-conservation violation** — the stored energy simply has no minimum in that direction.
8. ⭐ **Departures**: the frequency-dependent leak (above), and the two refuted pre-registered predictions.
9. ⭐ **Known limits** — carry them from the S11b-B record. ⚠ **The first one is load-bearing**: the supplied
   objects cannot be falsified by the build, because both engines transcribed them; their verification comes
   from independent reviewer derivations. ⛔ Do not soften this, and ⛔ do not omit the slice scope of the
   root result.

## ⛔⛔ MACRO HAZARD — check this before writing

⚠ Some stage-field-style macros in `paper/macros.tex` are **SUPPRESSED in the default build**. ⛔ **Anything
a reader must see that sits in a suppressed field is invisible in the PDF.** ⭐ **Read `macros.tex` first**,
determine which fields render, and ⛔ **do not place the departures, the known limits, or the
does-not-settle caveat in a suppressed field.**

## ⛔ CONSTRAINTS

1. ⛔ **No overclaiming.** Every statement must be traceable to one of the two step records.
2. ⛔ **Do not present a verdict on the physics as following from the engines' `VERDICT: PASS`** — that
   means only that each script's internal checks did not contradict each other.
3. ⛔ **Do not edit the step records, the scripts, or the registry.**
4. ⛔ **Do not narrate process** — revisions, review rounds, defect counts and repair passes belong in the
   step records, ⛔ not in the ledger.

## Output

The card, the `\input` wiring, a clean build with page count, and a report **under 20 lines** stating: what
you placed in the departures section, which macro fields you used and which you avoided as suppressed, and
confirmation that the dead result appears nowhere as live.

⛔ Then stop and exit.
