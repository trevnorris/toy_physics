# S15 item-4 Singular engine directive, review pass 2: disposition

## The gate record

- **Artifact.** `research/pde_ledger_v3/directives/S15_item4_singular_engine_directive.md`, revision 2, sha256
  `8c7e0035…`. It folds pass 1's findings G1–G3 (`S15_item4_singular_engine_directive_review_p1_disposition.md`).
  Revision 1 is preserved at `51fd98fa`. The fact check and leak gate for revision 2 is section R2 of
  `S15_item4_singular_engine_directive_lookups.md`.
- **Prompt.** `/var/projects/_review_scratch/s15-singular-dir-prompt/singular_engine_directive_review_prompt_p2.md`,
  sha256 `11b6c3aa…`, identical for both legs. It is pass 1's prompt with the artifact section updated, and it asks
  for a review of the whole revision, not only the fold.
- **Legs.** Codex (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, high), in parallel, 2026-11-09. Reports:
  `/var/projects/_review_scratch/s15-singular-dir-p2-codex/REVIEW.md` (sha256 `32946f56…`) and
  `/var/projects/_review_scratch/s15-singular-dir-p2-grok/REVIEW.md` (sha256 `88f710f7…`).
- **Verdicts.** Both CLEAR, with no findings.

## What the legs established

- **Object and scope.**
  - Grok checked that the quoted item-4 text is word for word the decision list's (`WORD_EQUAL True`).
  - Codex compared the directive's integrand test with the accepted Wolfram build's inventory. Both select 144 objects,
    the same 144, with none on only one side. The six two-drain region objects are among them.
- **Premises.** The pointed sections supply the balance, fields, region, domain and pressure conditions, the family,
  the separation derivative and the profile family. Note `4l`'s series is supplied as the identity it states.
- **Independence and leakage.** The construction is left to the builder. Nothing fixes a method or a truncation that
  both item-4 engines would share. The directive states no condition, value or expected outcome.
- **Claims and labels.**
  - Completeness is claimed only for returned conditions.
  - Every object in scope ends with its objects or a label, including after an interrupted attempt.
  - Uncomputed facts are printed as stated premises beside the steps that rely on them.
- **Operation.** The paths, the guard pool, the manifest contents, the README section and the stopping point are
  well defined.

## Disposition

Nothing outstanding changes what is computed or what may be claimed (G4). Revision 2 is **accepted** as the builder
directive for the S15 item-4 Singular engine.

## Next

A fresh Codex builder runs the directive. Its build is Codex-written, so the build review's legs are a fresh Claude
agent and Grok, reviewed until clear. The Grok launch for that review must let Grok take the guard's lock (round-5
build review, recorded). Then the item-4 comparator, and the S15 close.
