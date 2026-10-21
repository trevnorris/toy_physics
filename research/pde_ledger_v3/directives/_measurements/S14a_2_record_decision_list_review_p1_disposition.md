# S14a part 2 record decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_2_record_decision_list_review_p1_lookups.md`, beside
this file (generator `_scratch/s14a/gen/s14a2_recdl_p1_lookups.sh`). "N-n" cites a section of that file. Every
computed result quoted here is a leg's saved stdout.

## The pass

- **Artifact.** `directives/S14a_2_record_decision_list.md`, orchestrator-written, sha256 `cf7a5445…`, with
  `S14a_2_record_setup_lookups.md` (`4a4e433e…`) and the filed witnesses `S14a_2_review_witnesses/` with their
  checksum file `S14a_2_review_witnesses.sha256` (`f090dde1…`). All three are preserved unchanged in the commit that
  adds this file. All matched the baseline after the pass, and the working tree was unchanged apart from the lookups
  file itself (N1).
- **Role.** A pre-builder decision list for a physics-bearing record (G2: one Codex + Grok pass, fold once). The
  record it governs is reviewed until clear by fresh Claude + Grok after Codex writes it.
- **Prompt and legs.** `_scratch/s14a2/recdl_review_prompt_p1.md`, sha256 `38bcd63f…` (N1). Both legs received it
  identically, in parallel. Codex (xhigh, workspace-write in its scratch) wrote
  `/var/projects/_review_scratch/s14a2-recdl-p1-codex/REVIEW.md`, and Grok (`grok-4.7`, effort high) wrote
  `/var/projects/_review_scratch/s14a2-recdl-p1-grok/REVIEW.md`. Both exited 0 (N1).
- **Verdicts.** Codex: FINDINGS, 1. Grok: FINDINGS, 2 (N2).
- **What both legs reproduced.** Each reran the seven filed witness scripts. All seven stdout files matched the
  filed ones, with empty stderr (N4).
- **What Grok checked and found usable as written:** the quoted user option; the plan's four S14a objects and the
  CONDITIONAL standing of S14, S15 and S20; the register's schema, rest-on test and both routes; that `R-S12-02`
  does not already state the new requirement; and the two carried citation notes.

## Findings

**S1. Item 2 gives every witness one common reach that one witness does not have.** (Codex 1, Grok 2.)
**CONFIRMED.**
- **The list.** "A statement reaches only as far as the witnesses reach: particular solutions of the committed
  balances on which one choice of family, average or order of limits gives one value and another choice another"
  (l.56–57, N3).
- **The source.** The pass-2 Grok witness's first section states its own scope: "Spherical projected balance, fixed
  faces, no face channel", with `Σ = c t` and `F = −c r/3 + β t/r²` checked against that one balance (N4). The
  pass-2 disposition's R3 uses it for "two readings" of the volume-flux gate (N4). That is neither a choice of
  family, average or order of limits, nor a checked solution of both committed balances. The other witnesses
  describe themselves in other terms ("A fixed admissible field configuration…", "Small admissible shell witness,
  not a construction of the proposed build", "Two exact solution probes, not an S14a construction or acceptance
  test", and so on; N4).
- **Effect.** The record could state, for the pass-2 Grok witness, a stronger verification than that witness
  carries.
- **Cause.** The list wrote one class from the orchestrator's reading of the witnesses, rather than pointing at
  what each witness and its disposition state.

**S2. Item 6's census misses three lines that still present part 1 as unaccepted, and its last bullet then
freezes two of them.** (Grok 1.) **CONFIRMED.**
- **The list.** "Every line that still describes them as written for, under or awaiting review states that
  acceptance instead, heading or not … (L7)", then "Nothing else in the part 1 record or the register's S14a part 1
  pass changes" (l.75–81, N3).
- **The sources.** L7's predicate prints four lines: record l.3, register l.8 and l.1483, `STATUS.md` l.11 (N5). A
  wider literal search of the same files also prints:
  - `STATUS.md:1`, "S14a part 1 record written for review";
  - `SUBSTRATE_REQUIREMENTS.md:1486`, "This record/pass is written, not yet accepted.";
  - `steps/S14a_1_projection.md:455`, "source-first record review remains next under the contract." (N5).
- **Effect.** The register's S14a part 1 pass would state acceptance in its heading and "not yet accepted" in its
  next sentence. The part 1 record would state acceptance at l.3 and pending review at l.455. Two authors could
  read the frozen lines differently.
- **Cause.** The list fixed the census by a phrase, not by the property the lines must meet.
- **Not affected.** Record l.359 and l.444 describe S12's earlier state and its repair. They are statements about
  S12, not about part 1 (N5).

## Agreement between the legs

The legs agree on S1, found independently. Only Grok found S2. The census in N5 bears it out.

## Fold

Both findings change what the record may claim or what it changes. G2 folds them once, in a separate commit after
this pass is preserved:
- **S1:** item 2 states that each witness is described by what its own source and its disposition say it constructs
  and checks, and by no more.
- **S2:** item 6 names the property: in the three files, every statement that presents part 1's record or register
  pass as not yet accepted, or its review as still to come, states the acceptance. The lookups list the lines found
  (N5) without claiming that the list is complete. Statements about S12's earlier state stay as they are.
