# S14a part 2 decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_2_decision_list_review_p2_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s14a2_dl_p2_lookups.sh`). "M-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** `directives/S14a_2_decision_list.md`, revision 2, committed at `945a1cd2` (sha256 `c7d8fd26…`), with
  `S14a_2_setup_lookups.md`. Both matched the baseline after the pass, and the working tree was unchanged apart from
  the lookups file itself (M1).
- **Role.** The physics content of a pre-builder decision list, reviewed until clear (spec row). The G2 pass and its
  one fold are complete (`S14a_2_decision_list_review_p1_disposition.md`; fold `945a1cd2`).
- **Prompt and legs.** `_scratch/s14a2/dl_review_prompt_p2.md`, sha256 `4df59427…` (M1). Item 0 checks the fold;
  items 1–7 review the list again. Both legs received it identically, in parallel: Codex (xhigh, workspace-write in
  its scratch), `/var/projects/_review_scratch/s14a2-dl-p2-codex/REVIEW.md`; Grok (`grok-4.7`, effort high),
  `/var/projects/_review_scratch/s14a2-dl-p2-grok/REVIEW.md`. Both exited 0 (M1).
- **Verdicts.** Codex: FINDINGS, 2. Grok: FINDINGS, 1 (M2).

## The fold, checked

- **Q3 is repaired.** Both legs. Codex's assembly audit closes the order budget with the regular return and drain
  totals, the net face part counted once, the face terms after the full concentrated part and the inner flux
  inward: residual `0` (M7). Grok's signed and orientation identities agree (M2).
- **Q1's pointwise condition is gone.** Both legs. Codex finds that its replacement identifies constant averages
  with the selected steady solution; that is R1 below.
- **Q2's operation is defined.** Grok finds it definite on the objects it acts on. Codex finds that its
  identification with the dynamical `ω → 0` limit, and its order relative to the far-radius limit, are outstanding;
  those are R1 and R2.

## Findings

**R1. The long-time average is not the `ω → 0` limit of the dynamical law.** (Codex 1.) **CONFIRMED.**
- **The list.** "The DC part of such an object `X` is its zero-frequency part, the long-time average", labelled
  "`ω → 0`" and "the user's steady-throughput limit" (l.74–93, M3).
- **The source.** The plan: the limit must be "**derived** as the `ω→0` limit of the dynamical law" (l.599–601,
  M8).
- **The evidence.**
  - Codex's family A solves both committed balances for every `ω`, with constant shell return and drain totals and a
    slab growing at constant speed. Its `J/(4π)` is `17/12 - cos(omega*t)`. The long-time average gives `17/12`. The
    family's `ω → 0` limit at fixed time gives `5/12`. The part the average drops contributes `-1` in that limit
    (M4).
  - Codex's family B, with a strictly positive density: the average of the whole volume-flux quotient is
    `-4*sqrt(15)/15`, and its `ω → 0` value is `-4/5`, in units of `4π` (M4).
- **Effect.** Both values belong to steady states of the committed balances. The list labels one of them as "the
  `ω → 0` limit" without deriving it as such.

**R2. Averaging at fixed radii before the far-radius limit can remove spreading and exchange that persist.**
(Codex 2.) **CONFIRMED.**
- **The list.** "The average is taken at fixed position and fixed radii. Any far-radius limit (item 3) is taken
  after it, on the DC objects" (l.94–95). "The DC part of each projected density's time rate is the slab's steady
  spreading, and it stays live" (l.92–93, M3).
- **The decision.** "Exchange live": keep "the exchange and storage between the throat and the sphere" (l.66–68).
- **The evidence.** Codex's spreading-front family solves both balances and the full shell budgets. Its global
  storage, net face source and `J` are constant in time: `1/6`, `7/12` and `11/6` in units of `4π`. Averaged at a
  fixed radius first, every finite shell gives storage `0`, face channel `0` and face source `0`, and the far limit
  keeps those zeros. Taking the far limit first gives storage `1/6` and face channel `-7/6`. The outer material
  flux is `-3` in one order and `-2` in the other (M5).
- **Effect.** The fixed order can remove, at every radius, the spreading and exchange that the decision keeps live.
  Existence of the averages and the far limits does not detect this.

**R3. "Nonzero over the average" is not one domain for the volume flux.** (Grok 1.) **CONFIRMED.**
- **The list.** The volume flux is printed "where the projected material density at the outer radius is nonzero
  over the average" (l.133–135, M3). Part 1's condition was `Σ(R,t) ≠ 0` (part 1 decision list l.93–98, M8).
- **The evidence.** Grok's solution has `Σ = c t`, nonzero for every `t > 0`, with no finite average. The average of
  the whole volume-flux quotient exists and equals `4*pi*beta/c`. Reading the gate as "the average of `Σ` is finite
  and nonzero" prints nothing; reading it as "`Σ ≠ 0` over the window and the quotient's average exists" prints
  `4*pi*beta/c`. On this solution the separate averages of the outer flux and the core channel diverge, while
  their balance with storage averages to `0` (M6).
- **Effect.** Item 3 prints a different object under each reading, and one of them omits the throughput on a
  steadily spreading solution.

## Agreement between the legs

Both legs found Q3 repaired and Q1's pointwise condition removed. On Q2, Grok found the operation definite on its
objects and raised only R3; Codex found R1 and R2. These do not conflict: Grok's finding concerns the domain of one
printed object, while Codex's concern the identification and order of the limit itself. Lookup sections M4 and M5 bear
out Codex's two witnesses as filed.

## Adjudication: the premise goes to the user

R1–R3 each change what is computed or what may be claimed. They are not folded here.

**The orchestrator's reading of passes 1 and 2 together.** Pass 1's Q2 and this pass's R1 and R2 have one root.
The committed law gives two balances for more unknown fields, with every form open, so it does not fix how a
solution evolves. A slow limit then needs a family of solutions and an order of limits, and the committed law
selects neither. Each definition written so far has picked one: a pointwise condition in revision 1, the long-time
average at fixed radii in revision 2. Each time the legs found solutions of the committed balances where another
choice gives another value (pass 1 K5; this pass M4, M5). This is a reading of the filed witnesses, not a computed
result of a build.

Under the standing rules on long loops and on many review rounds, the next step is a question to the user, not a
third revision: whether to narrow the claim, record the obligation as needing the slab's dynamics, or stop.
