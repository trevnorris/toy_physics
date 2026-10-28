# S15 reading decision list, review pass 3 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_decision_list_review_p3_lookups.md`, beside this
file (generator `_scratch/s15/gen/s15_reading_p3_lookups.sh`). "M-n" cites a section of that file. Every computed
result quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** `directives/S15_reading_decision_list.md`, revision 3, committed at `369ea6bb`. It repairs pass 2's
  finding N1 (`S15_reading_decision_list_review_p2_disposition.md`) in revision 2, which is preserved at `48d06300`.
  The list and the pass-2 disposition and lookups matched the baseline after the pass (M0).
- **Scope.** This is the third pass on the list's spec content. The user chose how to run it (AskUserQuestion,
  2026-10-28): "Scoped pass 3 (Recommended)". The option text read: "Codex and Grok check only the fix, and whether it
  changed anything else that affects what is computed or claimed." The prompt checks:
  - whether revision 3 makes N1's "What must be true" statements true;
  - whether the repaired text is well defined and comparable across engines, and consistent across item 2, the
    exactness paragraph, the record section and the claim;
  - whether it prescribes a derivation or hints at a value;
  - text the repair did not change, only where the repair changes it.
- **Prompt and legs.** `_scratch/s15/reading_dl/reading_dl_review_prompt_p3.md`, sha256 `be658886…` (M0). Both legs
  received it identically, in parallel:
  - Codex (`gpt-6.1-sol`, effort xhigh): `/var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-readdl-p3-grok/REVIEW.md`.
- **Verdicts.** Both CLEAR (M1).
- **The tree.** The status matched the pre-pass snapshot. No file in the repository root is newer than the prompt
  (M0).

## What the legs measured

- **The removed examples are not conditions the build needs** (M3, M4). Each leg built its own cases.
  - The expansion's response exists with `p′` not differentiable at a base that does not move: Codex
    `stationary_expansion_Eresp: 8*pi*Bp*R**3*d/3`; Grok `expansion_pprime_derivative_gap: 2` with
    `expansion_fixed_base_Eresp_over_volume: Bprime`.
  - The whole pressure has a separation derivative where `p` kinks at a density that is reached: Codex
    `whole_kink_separation_derivative: 8*pi*R**3*d/3`; Grok `kink_p_density_quotient_gap: 2`,
    `kink_composite_quotient_gap: 0`.
- **The new sentence that S15's conditions do not give every derivative holds** (M4). Grok's whole minus the
  expansion: `p` is differentiable at the total density and at the base (gaps `0`, `0`), the base moves, and the
  response quotients differ (`whole_minus_response_quotient_gap: -2`).
- **Each part carries its own condition** (M4). Grok's path example: geometric gap `2`, full gap `2*Abs(bprime + 1)`.
  The list's "each departure object" covers both parts separately.
- **Agreement on smooth shared data** (M3, M4). Two constructions agree, with residual `0`: Codex's whole minus the
  expansion and the named sums' parts; Grok's spherical mean, its separation derivative and its geometric departure.
- **The wording** (M2, M4). Grok counts `0` remaining instances of "wherever it exists", "wherever they exist" and
  the removed pressure sentences, and locates the new phrases at the repaired lines.

## Examined and not filed

- **Comparing conditions.** Codex: "Different sufficient conditions remain permissible engine disagreement. They must
  travel with the objects being compared; value agreement on shared coverage does not establish equal coverage or
  nonexistence elsewhere." (M3). Grok: "A comparator can compare those values … Unequal sufficient labels sit outside
  the value." Neither filed it. It is carried to the reading comparator's decision list: the comparator compares
  each value together with the condition it carries.

## Clearance

Both legs CLEAR. Nothing outstanding changes what is computed or what may be claimed. The S15 reading decision list
is cleared at revision 3 (`369ea6bb`), after three passes: pass 1 (F1–F4, folded once), pass 2 (N1, repaired) and pass
3 (scoped, CLEAR). Next is the build directive, from this list, with its leak and contract gate before launch.
