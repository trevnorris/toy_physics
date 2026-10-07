# Independent review: targeted repair of the light-leakage scoping inventory (v4 → v5, Claude-authored)

## Artifact
`_scratch/light_leakage/light_leakage_scoping.md` (v5) in `/var/projects/toy_physics`. It is an inventory for a
later worst-case/best-case leakage bracket, written to the directive
`_scratch/light_leakage/light_leakage_scoping_directive.md`, and it contains no calculations. The previous
version is frozen at `_scratch/light_leakage/light_leakage_scoping_v4.md`. A fresh Claude author made v5 by
repairing one finding and nothing else. See the changes with:

```bash
diff _scratch/light_leakage/light_leakage_scoping_v4.md _scratch/light_leakage/light_leakage_scoping.md
```

## The finding the repair addresses
In v4, rows 3–4 (Raman, Brillouin) and their taxonomy note (lines 26 and 41) justified routing these processes
to C3, the generic nonlinear program, with the reading that a coupling of transverse light to another wave
perturbation lies beyond the records' "first order in wave amplitude" truncation. The records put the
transverse↔thickness coupling at first order in wave amplitude, `O(εη)`
(`research/pde_ledger_v3/directives/S11c_decisions.md:123–126`). The inventory's own E1 route (line 57) is that
first-order interbranch operator.

What must be true after the repair:
- The rows 3–4 routing to C3 does not rest on any claim that couplings between wave branches as such lie beyond
  first order.
- Nothing in rows 3–4 or the taxonomy note contradicts E1's first-order interbranch operator.
- No v3 vertex is invented.

## What to check
Only the changed material, and its consistency with the rest:
1. **Fidelity.** Does each changed cell meet the requirements above? Does it say only what the cited records and
   sources say? If you find that a requirement above is itself wrong against the records, say so, with
   evidence.
2. **New defects.** Does any changed cell introduce a wrong route, status, mapping, destination, freeze, validity
   domain, yardstick, oracle correspondence or citation?
3. **Consistency.** Does any changed cell contradict a cell the repair left unchanged?
4. **Scope.** Is every change in the diff required by the finding, or by consistency with it?
5. **No verdicts.** Does any changed cell state an outcome, a suppression requirement, or new physics?

## Required method
- Read the directive, then the cited records, then the diff and the changed cells.
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

Keep it under 800 words, not counting command output.
