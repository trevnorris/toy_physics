# S16 decision list, review pass 3: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_decision_list_review_p3_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s16_dl_p3_lookups.sh`). "K-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** Revision 3 of `directives/S16_decision_list.md` (sha256 `29d3225e…`), with its setup lookups and the
  pass 2 dispositions. The list and lookups were uncommitted on `6ad2d816`. All three matched the baseline after the
  pass. The working tree was unchanged apart from the lookups file (K1).
- **Role.** The spec content of a pre-builder decision list, reviewed until clear (spec row). This is the third
  review of the spec content, a full one, which the user chose: "Repair, then full review (Recommended)".
- **Prompt and legs.** `_scratch/s16/dl/dl_review_prompt_p3.md` (K1). Both legs received it identically, in
  parallel:
  - Codex (`gpt-6.1-sol`, xhigh, workspace-write, working directory its scratch):
    `/var/projects/_review_scratch/s16-dl-p3-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s16-dl-p3-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: CLEAR. Grok: FINDINGS, 3 (K2).
- **What both legs found holding.** E2, D1, D2, D4, D5 and D6. Both found the material balance and the five option
  texts quoted faithfully. Codex also found E1 and E3 holding; Grok's findings 1 and 3 concern them.

## Findings

**F1. "At rest" does not say in which frame, and the two frames give two received momenta.** (Grok 1.)
**CONFIRMED.**
- **The list.**
  - The surfaces are spheres "about `X_c(t)` with live radii `r_in(t)` and `R(t)`" (`:87–89`).
  - Item 1: "Each surface term is printed in two parts, kept apart: the part present when that surface is at rest;
    the part due to that surface's motion" (`:105–107`).
  - Item 3 builds on that split (`:116–119`) (K3).
- **The source.** E1's statement: "The received momentum is the object the user's decision names, kept apart from
  the terms due only to the surfaces' motion" (K4).
- **The evidence.** Grok: the supplied balance can be written in coordinates moving with `X_c`. There it holds with
  flux `Π − (dX_c/dt) G`, and "at rest" removes only the change of the radii. Both readings reproduce `dP/dt`. They
  differ by `2*R(t)*g1(t)*Derivative(Xc(t), t)` in one dimension. On a shell they differ by
  `4*pi*(R(t)**3 - rin(t)**3)*U(t)*g1(t)/3`, the translational sweep (K5, F1). Codex's probe took the first reading:
  its printed sweep is the whole motion term (K5).
- **Effect.** Two engines can each follow the list and print different received momenta.
- **What must be true.** The list names the frame in which a surface is at rest. Either reading can then be recovered
  from what the build prints, and no reading is chosen beyond what the decisions say.

**F2. "All of these stay debts for S22's table" includes objects the build prints.** (Grok 2.) **CONFIRMED.**
- **The list** (`:18–26`). Five bullets name the printed count and moments, the printed boundary fluxes, `M`,
  compactness with a field of the surroundings, and the source normalization. Then: "All of these stay debts for
  S22's table" (K3).
- **The source.** The plan's sentence lists `M`, `a`, the multipoles, compactness, source normalization and
  boundary-flux suppression as "supplied inputs that belong in S22's debt table" (`:739–741`) (K4).
- **Effect.** A record that follows the sentence would file printed objects as supplied-input debts.
- **What must be true.** The list names as S22 debts only what this step does not deliver. It keeps those apart from
  the objects it prints.

**F3. The compactness requirement drops one of the source's two comparisons.** (Grok 3.) **CONFIRMED.**
- **The list** (`:150–152`): "its size, small against the scale over which a field of the surroundings varies. This
  is the source's "Compact support / small body" assumption (`4d_1pn_full.tex:697–700`)" (K3).
- **The source** (`4d_1pn_full.tex:697–700`): "small compared with the inter-body separation and with the scale over
  which the external field varies" (K4). The list contains "inter-body" 0 times (K3).
- **What must be true.** The requirement states both comparisons the cited assumption makes, or does not identify
  itself with that assumption.

## Agreement between the legs

Codex found nothing outstanding. Grok found F1–F3. They do not conflict on the facts:
- Codex's E1 probe used the coordinates of the supplied equation, one of F1's two readings.
- Codex quoted F2's and F3's lines as paired with their sources, but did not test them against the plan's
  supplied-input wording or the source's second comparison.

## Rounds

This is the third review of the spec content. F1 is a frame left unnamed in text the pass 2 repair kept. F2 is in
text the fold wrote, and F3 is in text the pass 2 repair wrote. The repairs are technical. They name no new
computed object beyond splitting a printed term, choose no law, and change no decision. Before any repair runs, the
user decides how to finish, as stated after pass 2.
