# Independent review: targeted repair of the light-leakage scoping inventory (v3 → v4, Claude-authored)

## Artifact
`_scratch/light_leakage/light_leakage_scoping.md` (v4) in `/var/projects/toy_physics`. It is an inventory for a
later worst-case/best-case leakage bracket, written to the directive
`_scratch/light_leakage/light_leakage_scoping_directive.md`, and it contains no calculations. The previous
version is frozen at `_scratch/light_leakage/light_leakage_scoping_v3.md`. A fresh Claude author made v4 by
repairing two findings and nothing else. See the changes with:

```bash
diff _scratch/light_leakage/light_leakage_scoping_v3.md _scratch/light_leakage/light_leakage_scoping.md
```

## The two findings the repair addresses
**F1. Spontaneous Brillouin (and Raman) initiation.** v3 mapped it to a time-dependent background profile, the
background hold. But the records carry dynamical wave perturbations (`u`, `ζ_c`, `δW`, `θ`;
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:194–195`) on held material profiles (`:234–239`).
What must be true after the repair:
- Spontaneous and stimulated processes both point to the generic nonlinear program, whatever their initiation or
  gain class.
- The missing interaction, and the fluctuation statistics of the receiving branch, are each named as a gap.
- O7 keeps only independently time-varying material profiles.
- The source's thermal/vacuum distinction is quoted, not asserted for v3.
- The rows stay NOT ADDRESSED with no invented record status.

**F2. Friedland–Giannotti Eq. (8).** Its validity cell gave no domain. The source derives it for "an infinitely
small perturbation" that lifts the zero mode, using the unperturbed function (arXiv:0709.2164v1, §II).
What must be true after the repair:
- The validity cell states the small-lift domain (`m/k ≪ 1`), marked as an inference from the quoted derivation.
- It keeps "no v3 correspondence" and "no inherited bound".

## What to check
Only the changed material, and its consistency with the rest:
1. **Fidelity.** Does each changed cell meet the requirements above? Does it say only what the cited records and
   sources say? If you find that a requirement above is itself wrong against the records or sources, say so,
   with evidence.
2. **New defects.** Does any changed cell introduce a wrong route, status, mapping, destination, freeze, validity
   domain, yardstick, oracle correspondence or citation?
3. **Consistency.** Does any changed cell contradict a cell the repair left unchanged, for example a count, the
   survival matrix, or another row's pointer?
4. **Scope.** Is every change in the diff required by F1 or F2, or by consistency with them?
5. **No verdicts.** Does any changed cell state an outcome, a suppression requirement, or new physics?

## Required method
- Read the directive, then the cited records and sources, then the diff and the changed cells.
- Quote both sides for every finding.
- Each claim about a file or source needs the command or link and the literal text, or it is discarded. Write
  scratch files only under `/tmp`.
- Do not modify any file in the repository. Do not read other review outputs.

## Physics filter
Report a finding only if it changes what the later bracket would compute or claim from this inventory: a route,
a status, a mapping, a destination, a freeze, a validity domain, a yardstick, an oracle correspondence, or a
citation that would mislead the bracket. Do not report wording or formatting that changes none of these.

## Output
- **Verdict:** CLEAR or NEEDS REVISION.
- **Findings,** numbered. For each: the location, what is wrong, the evidence, and the fix.
- A short list of what you checked and found sound.

Keep it under 1,000 words, not counting command output.
