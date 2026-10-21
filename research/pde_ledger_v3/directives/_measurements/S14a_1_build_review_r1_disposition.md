# S14a part 1 build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s14a/gen/build_r1_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** The repaired build, under round 0's dispositions (`S14a_1_build_review_r0_disposition.md`, committed
  `0ee93e9b`). It is preserved unchanged in the commit that adds this file. The authors were split by engine:
  - the SymPy engine was repaired for B1 by its original author (resumed session `01a15e81…`);
  - the Wolfram engine was rewritten for B1–B3 by a fresh Codex author (session `01a15ecb…`), from the decision list.

  Both exited 0 (K1). The fresh author's log shows no `sed -n` or `cat` of the SymPy script (count 0, K1). Both final
  guarded runs exited 0 (K1). All five files matched the round-1 baseline after the round, and the working tree was
  unchanged apart from the lookups file itself (K1).
- **Prompt and legs.** `_scratch/s14a/build_review/s14a1_build_review_prompt_r1.md`, sha256 `fe211660…` (K1). It is
  round 0's prompt, with three changes:
  - item 10 added: check each round-0 repair, and review the build as a whole;
  - a guard log path under `_scratch/s11c`, one directory per leg;
  - the new file names.

  A fresh Claude agent (Opus) ran first, then Grok (`grok-4.7`, effort high), each with the identical prompt. Claude
  wrote `/var/projects/_review_scratch/s14a1-build-r1-claude/REVIEW.md`, and Grok wrote
  `/var/projects/_review_scratch/s14a1-build-r1-grok/REVIEW.md`.
- **Verdicts.** Claude: NEEDS REVISION, 1 finding, which it rates narrow. Grok: NEEDS REVISION, 1 finding (K2).
- **What both legs found sound.**
  - B1 is repaired in both engines, and B2 and B3 in the Wolfram engine (K2).
  - Claude evaluated both engines' printed objects on its two exact witnesses: 40 of 40 SymPy objects and 43 of 43
    Wolfram objects equal its own derivation.
  - All of each leg's ablations moved the output, including FORM ablations in each engine and one-sided corruptions.
  - Guarded reruns of both engines reproduce the transcripts byte for byte.
  - The Wolfram engine reads nothing from the SymPy side, and builds the projection by its own route.

## Findings

**C1. A SymPy residual is labelled as a two-route check of face pieces it cannot test.** (Claude 1.) **CONFIRMED.**
Narrow: it changes what may be claimed, not what is computed.
- **The engine.**
  - The bulk-direction face pieces (`bulk_parts`, the B1 repair) and the face-concentrated parts (`atoms`) enter both
    operands of the "projection integration" residual in the same way (`:177–178`, `:186–187`, `:199–205`, K3).
  - Its comparison line reads "Differentiation of projected integrals versus integration of the regional
    differential operand with full face parts" (K3).
  - The engine's own lines for those two pieces say "No second route; nothing to compare it with" (K3).
- **The evidence.** In Claude's copy with `bulk_parts` removed (which puts B1 back), and in its copy with `atoms`
  removed, the residual lines do not change (diff count 0 each). The S1 copy prints `0` on all four (K3). Its shell
  `J` and shell balance residual are wrong on Claude's witness (K3).
- **Why it matters.** The residual tests only the time and in-plane endpoint pieces. A record or comparator citing it
  as evidence for the face terms, or for the B1 repair, would claim more than the transcript shows.
- **What must be true.** The comparison line states which pieces the residual tests and which it shares between its
  operands, or the operands share no face-term construction. No physics object changes.

**G1. The Wolfram engine prints the order lost from the balance's left side, not from the projected order source.**
(Grok 1.) **CONFIRMED.**
- **The list.** Item 2: "the order lost over the shell: minus the projected order source integrated over the shell,
  including its face parts". Item 1 prints "the projected order source, including its face parts". The committed
  order balance has `n Γ_B` on its right side.
- **The engine.**
  - It stores two sources:
    - `required`, the `w`-integral of the interior divergence plus the face parts;
    - `suppliedSource`, the `w`-integral of `n Γ_B` plus the face parts.
  - Its projected order balance has `suppliedSource` on the right (`:163–165`, `:179`, K4).
  - The key `"source"` holds `required` (`:180`, K4). The order lost and the face channel minus the order lost are
    built from it (`:209`, `:215`, K4).
  - The `suppliedSource` form is printed separately, as "order lost using supplied bulk source" (`:210`, K4).
- **The evidence.**
  - `sgamma`, the engine's `Γ_B`, occurs zero times in the printed "order lost" and in "face channel minus order
    lost". It occurs once in "order lost using supplied bulk source" (K4).
  - The SymPy engine's "order lost" contains `GammaB` once (K4).
  - The two Wolfram forms differ by the printed source-law residual, which is a live expression (K4). They coincide on
    solutions of the order balance. That is why Claude's witness values for the order lost matched (K4).
- **Why it changes what is computed.** The object printed under the list's name is another object: equal to it on
  solutions, but not the same printed expression. A comparison of the two engines' order lost would then compare
  different objects under one name.
- **What must be true.** The object printed as the order lost is minus the shell integral of the projected order
  source: the `w`-integral of `n Γ_B` plus the face parts, the source on the right of the projected order balance. The
  face channel minus the order lost uses that object. The divergence form may be printed under its own name.

## The leg disagreement

Each leg found one item the other did not raise. On G1, Claude's witness evaluation could not see it, because the two
forms coincide on solutions. Grok's count and the source-law residual show it (K4). On C1, Grok described the SymPy
residuals as checks that "name two routes". Claude's ablations show this one does not test the face pieces (K3).

## Adjudication

C1 changes what may be claimed. G1 changes what is printed as two of the list's objects. Both repairs are local: a
comparison line in the SymPy engine, and which source the Wolfram engine's order lost is built from. Neither leg found
a defect in the other objects. On Claude's witnesses, every printed object matches an independent derivation. A third
review round would follow a repair, so the user decides how to finish.

## User decision (2026-10-21): "Fix both, scoped check"

Each engine's author repairs its own finding: the SymPy author C1, the Wolfram author G1. A fresh Claude agent and Grok
then check only those two edits, and that nothing else moved since `eaccd417`. The comparator follows.
