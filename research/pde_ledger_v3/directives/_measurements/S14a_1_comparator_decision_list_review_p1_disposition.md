# S14a part 1 comparator decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_comparator_decision_list_review_p1_lookups.md`,
beside this file (generator `_scratch/s14a/gen/cmp_dl_p1_lookups.sh`). "D-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S14a_1_comparator_decision_list.md`, as first written, with its setup lookups
  `S14a_1_comparator_setup_lookups.md`. Both files matched the pass baseline after the pass, and the working tree was
  unchanged (D1).
- **Prompt and legs.** `_scratch/s14a/cmp/cmp_dl_review_prompt_p1.md` (D1). Both legs got the identical prompt, in
  parallel. The orchestrator wrote the list, so the legs are Codex + Grok (G1):
  - Codex (`gpt-6.1-sol`, effort xhigh, D1) wrote `/var/projects/_review_scratch/s14a1-cmpdl-p1-codex/REVIEW.md`. Its
    first launch was refused before any work ("Not inside a trusted directory"). It was relaunched with
    `--skip-git-repo-check` and the same prompt.
  - Grok (`grok-4.7`, effort high) wrote `/var/projects/_review_scratch/s14a1-cmpdl-p1-grok/REVIEW.md`.
- **Verdicts.** Codex: one finding. Grok: CLEAR (D1).

## Finding

**K1. A row binds whole tags, so an object one engine splits across tags and the other groups in one tag cannot be
fully compared.** (Codex 1, High.) **CONFIRMED.**
- **The list.**
  - Item 2: "Each row pairs one tag of each stream." "The table is injective." (D2)
  - Its component correspondence works within "the two payloads" of one row (D2).
- **The governing list.** Item 1 prints "every term that leaves at each face", and item 4 asks for "Dimensional
  homogeneity of the projected balances and of objects 2 and 3" (D3).
- **The gap.**
  - SymPy prints each face's closed-projection leaving term under its own tag, `face 0` and `face 1`. Wolfram prints
    both in one tag, `each closed-projection face term`, as a list (D4, D5).
  - SymPy prints the three in-plane flux dimensions under three tags, `flux 0`, `flux 1` and `flux 2`. Wolfram prints
    the three in one tag, `plane_flux` (D4, D5).
  - With one tag per side per row and an injective table, a grouped tag can sit on one row only. So at most one of the
    split tags can be compared with its counterpart (Codex's enumeration: 1 of 3 for the flux dimensions, D6). The
    others are left unjoined. A difference confined to one of them is never compared. The accounting shows the loss,
    but it does not make the comparison.
- **What must be true.**
  - A binding may join a tag or a component of a tag's payload on one side with a tag or a component on the other.
    Each end is named by its tag and its place in the payload, and is cited to the line that builds it.
  - Injectivity holds at those ends: no tag or component is bound twice, and a repeat is reported. Every component of
    a bound tag is bound, or printed unpaired with its reason.
  - An object that both engines print and the S14a part 1 decision list names is compared. It is left unjoined only
    where the cited lines build different objects, and its reason says so.
  - A control goes through split-against-grouped packaging. Each counterpart component is mutated on its own, and the
    control fails unless that component's comparison moves.

## Not findings

- **Grok's checks.** Grok checked items 1–7 against the streams and the engines and filed nothing. Its notes record
  facts the comparator meets on these streams, without a finding:
  - only one engine prints the projected velocity at a general radius;
  - the two engines print the dimension of an additive zero in different notations;
  - Wolfram's `{z, 0, 1}` integrals sit on a different printed object from the projected densities.
  The list's existing rules (unjoined with a reason; the same object read from both notations; a row only where the
  cited lines build the same object) cover each. K1's fold does not change them.

## Outcome

One finding, confirmed. The list as reviewed is preserved, with these dispositions, before the fold. The fold
repairs item 2, item 7 and the builder report for K1. It is the one fold G2 allows. Then the comparator build.
