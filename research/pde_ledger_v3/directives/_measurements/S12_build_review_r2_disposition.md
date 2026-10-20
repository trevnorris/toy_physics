# S12 build, review round 2 (scoped, after repair 2): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_build_review_r2_lookups.md`, beside this file
(generator `_scratch/s12/gen/build_r2_lookups.sh`). "D-n" cites a section of that file.

## The round

- **Why.** Round 1 confirmed two findings (`S12_build_review_r1_disposition.md`). All three findings in rounds 0 and 1
  traced to one gap in the decision list. As the third round, it was raised with the user, who chose: "Fix the list,
  then repair", with a scoped build review after; and "Leave it out for now" for order carried along a face. The list
  was revised and cleared as revision 4 at `b08673a2` (`S12_decision_list_review_p{4,5}_disposition.md`). The same
  builder then repaired the build against it.
- **Artifact.** The S12 build after repair 2, uncommitted at review on top of `b08673a2`; round 1 is preserved at
  `04a24c36`. It is committed unchanged with this file. All five baseline hashes matched after both legs, and the
  working tree was unchanged apart from this file and the lookups (D1).
- **Prompt and legs.** `_scratch/s12/build_review/s12_build_review_prompt_r2.md`, sha256 `4a8b9ac8…` (D1). It asked
  whether B3-1 to B3-5 hold for each engine, and for FORM controls that touch the face rule. Both legs received it
  identically and ran one after the other.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-build-r2-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-build-r2-grok/grok_stdout.txt`.
- **Verdicts.** Both legs: CLEAR, no findings (D2).

## What the legs checked

- **B3-1.** Each engine builds `n`, `u` and `J_χ` from separate live functions on each of the three regions. The
  order is live inside and zero outside, and no field has a face-concentrated part.
- **B3-2.**
  - The source's face parts, and what the total-material balance requires at each face, are reached by
    differentiation and coefficient extraction. Neither is typed.
  - The material requirement is printed. It is never solved or substituted.
- **B3-3.** Round 1's face-concentrated flux is gone from both engines and both transcripts (D3).
- **B3-4.**
  - Item 2 stays labelled and general in both offsets.
  - Dimensions are computed. The flux's units are solved from the order balance and printed as an inference.
  - The supplied balances are transcribed exactly, emission is objects only, and the Wolfram engine imports nothing.
  - Nothing is frozen.
- **B3-5.** The run record matches the guard records, and each guard stdout is byte-equal to its transcript (D4). Each
  leg's guarded rerun of copies reproduced both transcripts byte for byte.
- **Own derivations and controls (D5).**
  - Each leg derived the object in its own script before opening either engine, and found both engines equal to it.
    The Claude leg also checked the face terms by a second route, a weak form, with its own negative controls.
  - Every FORM ablation moved the printed output. They include:
    - fields that may not jump: `J` (SymPy) and `u` (Wolfram);
    - one shared density field (SymPy);
    - order live outside (Wolfram);
    - one offset for `G` (SymPy);
    - round 1's face flux re-added (Wolfram).
  - Every one-sided corruption of a supplied balance moved the printed output.

## Notes for the record and the comparison (no change to what is computed)

- SymPy's face objects carry a step in `w₊ − w₋`, which the supplied strict face order sets to one. Wolfram's do not.
  The cross-engine comparison applies that domain.
- SymPy forms each region's object by substituting that region's fields. Its region objects are therefore not a test
  of its support encoding (round 0 note).
- The flux's units are solved from the order balance, so homogeneity is definitional for `J_χ` and informative for
  the other terms (round 0 note).

## Outcome

No findings. Nothing outstanding changes what is computed or what may be claimed, so the S12 build is accepted as
reviewed in round 2, and committed with this file. Next come the cross-engine comparison, a separate step, and then
the S12 record and its register pass, each reviewed until clear.
