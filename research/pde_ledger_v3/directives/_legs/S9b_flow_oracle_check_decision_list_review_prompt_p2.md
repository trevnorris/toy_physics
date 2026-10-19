# Independent review: the S9b flow-only oracle-check decision list (orchestrator-written)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md`, uncommitted.
Its sha256 is in `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline_p2.sha256`.

It is the decision list a Codex builder will work from to write one check: S9b's flow-only light-bending condition
against the river/acoustic-flow prior art. It names the object the check computes. The builder sees only this list
and the sources it names.

This is a revision. The reviewed previous version is commit `504fa75c`;
`git diff 504fa75c -- research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md` shows the change. The
first pass's findings and their dispositions are in
`research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_decision_list_review_disposition.md`. Review the
whole revised list, not only the diff.

## What you are handed
- The artifact.
- Its sources, at the revisions it gives:
  - the S9b engine transcripts (git-annex content, already present);
  - the accepted S9b record `research/pde_ledger_v3/steps/S9b_light_bending.md` (`ff3bac63`), §3;
  - the scoping survey `research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md` (`36f59ff1`), §5;
  - `docs/medium_requirements_and_prior_art.md`, lines 240–270;
  - the user decision at `research/pde_ledger_v3/directives/_measurements/S9b_record_r2_review_disposition.md`
    (`4760e6e2`), lines 19–21.
- The two oracle manuscripts, if you can open them: arXiv gr-qc/0411060 (Hamilton and Lisle) and gr-qc/9712010
  (Visser).
- `CLAUDE.md`, sections M1–M3, E1, E2 and G2.

## What to check
0. **The first pass's findings.** Is each one resolved in the revised list? Did the revision introduce a new defect?
1. **Is this the right question at S9b's retained order, and not a convenient proxy?** Does it compare the objects
   the S9b engines print with the objects the oracles state? Does the residual it asks for tell us something about
   S9b's flow-only condition that the printed objects alone do not? Is any comparison missing that the oracles'
   stated conditions call for?
2. **Sources and tags.** Do the tags and fields the list names exist in the transcripts, in the form the list gives?
   Are the oracle citations, sections and equation numbers right? Are §5's conditions for each oracle all carried
   into the identification map?
3. **Live quantities (M3).** Does any decision hold fixed a quantity the sources keep live, or treat a supplied
   input as derived?
4. **Leakage and recipe (M2).** Does the list state or imply an expected value, sign or comparison outcome? Does it
   dictate how to compute something, where it should name the object instead? Can a builder reach a passing script
   by choosing an identification?
5. **Scope.** Are the out-of-scope items and their reasons sound? Should another identified oracle be in scope?

## Required method
- Your scratch directory is `/var/projects/_review_scratch/s9b-oracle-dl-p2-<your leg name, codex or grok>/`. Create
  it, and save every retrieval script with its literal stdout there (not under `/tmp`). Report the paths. A claim
  about a file's contents without its command and output will be discarded.
- Quote both sides for every finding.
- This is a review of decisions and their evidence. Retrieval scripts are enough; run no CAS. There is no script yet,
  so do not write or ablate one.
- Never modify the working tree.
- Do not read `~/.claude/projects/`.

## Filter
Report a finding only if it would change what the check computes, or what may be claimed from it.

## Output
- **Verdict:** SOUND or NEEDS REVISION.
- **Findings,** numbered. For each: the location, what is wrong, the evidence (script path and literal stdout), and
  the fix as a statement of what must be true.
- A short list of what you checked and found sound.

Keep it under 1,200 words, not counting command output.

## Bounds
Write your report and exit. Do not spawn agents.
