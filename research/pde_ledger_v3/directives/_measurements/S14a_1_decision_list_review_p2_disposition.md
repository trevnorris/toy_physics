# S14a part 1 decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_decision_list_review_p2_lookups.md`, beside this file
(generator `_scratch/s14a/gen/dl_p2_lookups.sh`). "K-n" cites a section of that file.

## The pass

- **Artifact.** Revision 2 of `directives/S14a_1_decision_list.md`: the one fold of pass 1, sha256 `592ad077…`. It is
  preserved unchanged in the commit that adds this file. It matched the baseline after the pass, and the working tree
  was unchanged apart from the lookups file itself (K1). Against the pass-1 text at `82fa9fa3`, the change is 40
  insertions and 13 deletions (K1).
- **Role.** As in pass 1: a pre-builder decision list whose content names the object, so it is reviewed until clear.
- **Prompt and legs.** `_scratch/s14a/dl_review_prompt_p2.md`, sha256 `89561235…` (K1). It is pass 1's prompt, with
  item 7 added: check each pass-1 repair against its "What must be true", and review the whole list. Both legs
  received it identically, in parallel. Codex (xhigh) wrote `/var/projects/_review_scratch/s14a1-dl-p2-codex/REVIEW.md`,
  and Grok (`grok-4.7`, effort high) wrote `/var/projects/_review_scratch/s14a1-dl-p2-grok/REVIEW.md`. Both exited 0
  (K1).
- **Verdicts.** Codex: CLEAR. Grok: FINDINGS, 1, which it rates narrow (K2).
- **What both legs found sound.** P1 to P5 are each repaired as their "What must be true" requires (K2). Each leg
  backs that with its own witnesses and literal stdout:
  - the face terms keep motion and tilt;
  - the native fields keep `w`;
  - the order lost is signed, with face parts;
  - the quotient is restricted to nonzero projected density;
  - the owners are as the accepted S12 record gives them.

  Codex's whole-list check found no other freeze, missing flux or overclaim. Its units audit puts the shell rates
  (face, core, outer, storage, order lost, `J`, difference) at one dimension, and the volume flux at another.

## Finding

**Q1. The printed combination of the face channel and the order lost does not state its order.** (Grok 1.)
**CONFIRMED, narrow.**
- **The list.** "the difference between the face channel and the order lost" (l.89, K3). No line says "the face
  channel minus the order lost" (count 0, K3).
- **The pass-1 line.** P3's "What must be true" required that "The printed difference uses that definition" (K3).
  It fixed the sign of the order lost, not the order of the subtraction.
- **The evidence.** Grok prints three readings of the combination (face minus lost, lost minus face, absolute value).
  For operands 2 and 5 they give `-3`, `3` and `3`. Both operands have the same dimension, so dimension does not
  choose (K3). Codex's witness read it as face channel minus order lost (K3), which is one reading, not the text.
- **Standing.** It changes the sign of one printed object. The operands are printed separately and are unaffected,
  and so are the balances, the two channels, `J`, storage and the outer flux. The two engines could still pick
  opposite orders, and the comparator would then show a convention mismatch, not physics.
- **What must be true.** The printed combination states its order of subtraction explicitly, in terms of the two
  printed operands.

## The leg disagreement

Codex found the list clear. Its witness used one order of the subtraction (K3), so it did not test the wording.
Grok's three readings (K3) bear out the finding.

## Adjudication

Q1 changes the sign of one printed object and nothing else. P1 to P5 are repaired, and both legs found no other
defect in revision 2. The next pass would be the third, so the user decides how to finish.
