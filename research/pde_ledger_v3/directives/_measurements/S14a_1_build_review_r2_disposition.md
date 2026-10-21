# S14a part 1 build, scoped review round 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_build_review_r2_lookups.md`, beside this file
(generator `_scratch/s14a/gen/build_r2_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** The build after repair 2, under round 1's dispositions and the user's "Fix both, scoped check"
  (`S14a_1_build_review_r1_disposition.md`, committed `52c1f939`). Each engine's author repaired its own finding:
  - the SymPy author C1 (resumed `01a15e81…`);
  - the Wolfram author G1 (resumed `01a15ecb…`).

  Both exited 0 (K1), and both final guarded runs exited 0 (K1). All five files matched the round-2 baseline after the
  round. The working tree was unchanged apart from the lookups file itself (K1).
- **Prompt and legs.** `_scratch/s14a/build_review/s14a1_build_review_prompt_r2.md`, sha256 `54eae301…` (K1). It was
  scoped to six checks:
  - C1 and G1, each against its "What must be true";
  - nothing else moved since `eaccd417`;
  - no new defect;
  - guarded reruns;
  - a FORM ablation of the repaired Wolfram object.

  The round-1 transcripts were handed as copies taken from the annex. Their sha256s are the ones `eaccd417`'s run
  record gives (K1). A fresh Claude agent (Opus) ran first, then Grok (`grok-4.7`, effort high), each with the
  identical prompt. Claude wrote `/var/projects/_review_scratch/s14a1-build-r2-claude/REVIEW.md`, and Grok wrote
  `/var/projects/_review_scratch/s14a1-build-r2-grok/REVIEW.md`.
- **Verdicts.** Claude: CLEAR. Grok: CLEAR (K2).

## What the legs established

- **What moved.** The SymPy script changed 1 line and the Wolfram script 5 lines since `eaccd417` (K3). The SymPy
  transcript changed 4 lines. The Wolfram transcript changed 18 lines: 16 tag-only renames and the two order-lost
  payloads (K3).
- **C1 is repaired.**
  - The comparison line now says the residual compares the time and in-plane endpoint pieces. It also says the
    bulk-direction pieces and the face-concentrated parts are shared by both operands and not tested (K4).
  - Claude's ablations bear this out. Dropping `bulk_parts` or `atoms` changes no residual line, and all four still
    print `0`. Dropping the time or the in-plane endpoint pieces changes all four (K4).
  - Grok's symbolic version, with a tag on each piece, prints the same split (K4).
  - No physics object changed (K3).
- **G1 is repaired.**
  - The order entry's source is now `suppliedSource` (K5).
  - `sgamma`, the engine's `Γ_B`, now occurs in the printed "order lost" and in "face channel minus order lost" (once
    each, K5).
  - The divergence form is printed under its own name, "projected order divergence" (8 lines, K5).
  - Grok's parse of the transcript shows that the projected order balance's right side is the same expression as the
    projected order source (`SameQ` True). It also shows that the order lost is the same expression as minus that
    source's shell integral (K5).
  - Claude's FORM ablation W1 builds the order lost from the divergence again. Exactly the two order-lost lines move,
    and they return to their round-1 form. W2 drops the `n Γ_B` integral, and 9 lines move (K5). Grok's FORM ablation
    gives the same two-line result (K5).
- **Reproducibility.** Both legs' guarded reruns reproduce both transcripts byte for byte (K6).

## Known item, carried open

**N1. The run record's construction-order sentence is stale.** (Orchestrator, from K8; neither leg raised it.)
- **The sentence.** Line 5 of `directives/_measurements/S14a_1_build_runs.md`: "The Wolfram source was written before
  the SymPy source." It is present at `0ee93e9b`, at `eaccd417` and now (K8).
- **What is true.** The first build's author wrote it, about the first Wolfram engine, and there it was true. The
  current Wolfram engine was rewritten by a separate author after repair 1 of the SymPy engine had finished (K8). Its
  log shows no read of the SymPy script (round 1's lookups, K1).
- **Standing.** It changes no computed object. It is a provenance claim, so the S14a part 1 record states the actual
  authorship and order from these dispositions and does not repeat the sentence. The record's decision list also
  orders the run record's line corrected.

## Not findings

- **The run record's other edits.** Claude notes that the run record changed more than each engine's command, guard
  outcome, transcript and sha256: 29 changed lines, including the deleted "Construction:" bullet (K7). Both legs
  checked the new sentences against the artifacts and found them true, and the deletions only remove claims. The
  sentence "Neither engine reads the other engine or its output" is unchanged, and the Wolfram script has no `Get`,
  `Import` or `ReadList` (count 0, K7).
- **Claude's note on A5.** A one-sided interior control also moves the C1 residual, so the new comparison line
  understates what the residual tests. That claims less, not more.

## Adjudication

Nothing outstanding changes what is computed or what may be claimed from the engines' output. N1 is a provenance
sentence the record corrects. That meets the substantive exit (`CLAUDE.md` G4), under the user's "Fix both, scoped
check". The build is accepted as reviewed in this round: the five files byte-for-byte as in the round-2 baseline.

Three review rounds preceded acceptance: two full, one scoped. Three authors contributed:
- the first builder wrote both engines;
- the same author repaired the SymPy engine twice;
- a fresh author rewrote the Wolfram engine, then repaired it once.

Next is the cross-engine comparator.
