# Independent physics review — S11b SymPy interface-coupling-law engine (revised)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py`
— a SymPy engine (a REVISED version; review it fresh on its own merits). Running it flushes a tagged stdout
stream AND writes `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_exports.py` (the chain
record the next step imports). Both products are in scope.

## What to check — two burdens, discharge BOTH
**(A) Every emitted physical result is genuinely COMPUTED, not typed.** A payload that is a hand-authored CAS
object with no data dependence on the derivation is a defect invisible to code-reading. In particular the
**impedance** consumed by the assembly (face response, the regime table, added mass, port impedance,
dispersion, roots, the breathing slice) must be reached by computation from the bulk acoustic field, with
every regime/slice form descending from ONE such object; ⛔ no typed impedance literal, and ⛔ no separately
re-typed regime or `k=0`-slice form.

**(B) The export carries no tautological check.** For EVERY emitted check / residual / status / "identity"
tag, the residual must be able to FAIL.

## ⛔⛔ The decisive test for a check — ONE-SIDED CORRUPTION of the RESIDUAL, ⛔ not "the row moves"
For each check, corrupt **only** the live object it claims to police (the assembled equation / the bulk solve
/ the projection source), re-run, and look at the check's **residual or status** — the tautological element
itself. ⛔⛔ **"The emitted row's bytes changed" is NOT sufficient**: a tautological residual sitting in the
same payload as a moving ornament (e.g. a solved `omega²` beside an inert `x−x`) will move the row while the
residual stays `0`. Report any check whose **residual/status is inert** under a one-sided corruption of what
it polices, or whose two routes **share a constructor** (e.g. both call the same model-builder) so a form
change moves both together. Name the tag and the line.

## ⛔⛔ A FORM ablation is MANDATORY (rule 14)
Copy the engine to `/tmp` and ablate the COPY (⛔ never the working tree). Change the STRUCTURE of a
load-bearing object — flip a sign AND an off-diagonal, collapse two independent symbols, alter a governing
equation's form, corrupt the bulk-solve construction — re-run, and report the LITERAL diff. A COEFFICIENT
rescale tests arithmetic only; only a FORM change tests physics. For every load-bearing object ask
**WHICH LINE COMPUTED THIS?** and give the line, or report it uncomputed. Save every ablation script AND its
literal stdout to named absolute paths, and report them.

## Independent derivation
Write your OWN derivation script for the load-bearing objects (a coupling coefficient, the impedance, the
breathing root, the passivity form) BEFORE opening the engine; save the script and its literal stdout to
named absolute paths and report them. ⛔ A prose "I re-derived it" is discarded.

## Also flag regressions
This engine is a revision; a change may have introduced a NEW typed payload, a NEW tautology, or altered a
physics VALUE. Treat any such as a finding.

## What you are handed
- The engine and its product `S11b_exports.py` (paths above).
- The physics source of truth `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` (§§0–13; §6 energy
  accounting / convention checks; B7 dimensions) — what the engine must compute.
- The chain contracts (for wiring, ⛔ not physics answers): `F1`/`F9` in
  `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md`; `D1`/`D3` in `S9_REWRITE_PLAN.md`;
  `CLAUDE.md` (rule 2: a script PRINTS operands + residual, never asserts a conclusion; corollaries 1 & 3).
- ⚠ Your job is independent verification, ⛔ not confirming any expected result. Derive the physics yourself;
  add no expected value.

## Physics filter
Report a finding only if it catches a way the PHYSICS could be wrong, or a check that cannot fail, or a typed
load-bearing object, or a wiring defect corrupting the exported objects. ⛔ Not "wrong on a different input."

## Ablation sandbox
Copy to `/tmp`, ablate the COPY, ⛔ never modify the working tree. Pure Python — no CAS kernel, no seat, no
timeout constraint. Report every script + stdout path.

## Output
Findings most-severe first, each with WHICH-LINE evidence and its literal ablation/one-sided-corruption
output; then the FORM-ablation results with literal diffs and /tmp paths; then an explicit statement if
nothing survives the filter, naming which corruptions you ran and that the residuals BIT.
