# S16 build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s16_build_r1_lookups.sh`). "N-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout or an engine's transcript line.

## The round

- **Artifact.** The S16 build after repair 1, with the same five files as round 0. The round 0 version is preserved
  at `b60e38fb`.
  - The original Codex session (`01a16176…`) made the repair, exiting 0 at 268,074 tokens. It was resumed at its
    default model, `gpt-6.1-sol`, effort high, with full access. The build itself ran at `gpt-6-astra` (N1).
  - The repair changed both engines and the run record (N1). Each guard run that wrote the final transcripts exited 0
    (N1).
  - All five files matched the baseline after the round. The working tree was unchanged apart from the lookups file
    (N1).
- **Prompt and legs.** `_scratch/s16/build_review/s16_build_review_prompt_r1.md`, sha256 `4f00ce5e…` (N1). It is the
  round 0 prompt, with a first item that tests the repairs. The legs were a fresh Claude agent (Opus) and Grok
  (`grok-4.7`, effort high). Both got the identical prompt and ran in parallel, with one Wolfram kernel at a time.
  - Claude: `/var/projects/_review_scratch/s16-build-r1-claude/REVIEW.md`.
  - Grok: `/var/projects/_review_scratch/s16-build-r1-grok/REVIEW.md`, exit 0.
- **Verdicts.** Claude: CLEAR, no findings, four notes. Grok: NEEDS REVISION, one finding (N2).

## The repairs, as the legs checked them

- **B1 holds at the rate.** Each engine prints the reference term apart from the flux (N3):
  - SymPy: `dipole.0.bulk.reference.centerMotion = Volume(-Sigma(x1, x2, x3, t)*Derivative(Xc1(t), t))`, beside
    `dipole.0.bulk.weight = Volume(F1(x1, x2, x3, t))`.
  - Wolfram: the same pair under `dipole.1.bulkReference.centerMotion` and `dipole.1.bulkWeight`.

  Both legs undid the split on a copy. The term then returns to the weight, and their checks fail (N5):
  - Claude's U1 gives 27 nonzero lines in each engine, the first `dipole.0.bulk.weight | … | residual=258879*pi/80`.
  - Grok's `RESIDUAL b1 undone weight minus int F1` is nonzero, `936*pi` at its point.
- **B1's second part and B2 hold.**
  - Claude evaluated every printed object of each engine against its own ground truth, on two seeds of live fields
    (N5). All 418 had residual 0. The pieces summed to their derivative in 44 of 44 cases. The test covered every
    piece of all 44 derivatives, including those whose integrand carries a velocity.
  - Undoing the volume-part split (U2) puts `−Ẍ_c N` back into `localEvolution`:
    `Drate0.termRatePart.bulk.reference.centerMotion.localEvolution | … | residual=258879*pi/64` (N5).
  - Grok matched the pieces to an independent Reynolds derivative on two witnesses, including squared and mixed
    velocity products: `RESIDUAL A printed pieces minus reynolds = 0`, and the same for witness B (N5).
- **B3 holds.** Both legs found the run record's route and domain text matching the code and transcripts (N7).
- **The other repair changes.** Undoing them (Claude U3 and U4) changes no value (N5). They are refactors.

## Grok's finding

**G1. The `localEvolution` part of a motion term's derivative carries that term's undifferentiated velocity.**
(Grok 1.) **The fact is CONFIRMED. It is NOT CONFIRMED as a build defect, and it is carried to the record as a limit
on claims (C1).**
- **The fact.**
  - In each transcript there are 24 `localEvolution` lines under a term of the motion. Nineteen of them contain a
    velocity of `X_c` or of a radius.
  - The 20 `localEvolution` lines under other terms contain none (N4).
  - Claude's note (d) states the same fact.
- **The two readings.**
  - **Grok** reads "the part present when that surface is at rest" as the value with those velocities set to zero.
    It asks that `localEvolution` hold no velocity factor, and that the undifferentiated velocity move into the motion
    piece its parent term is named for.
  - **Claude** judged this disposition's round 0 reading consistent with the list. Read literally, "no jet" cannot
    be met by any surface term, because the surface position `X_c + r n̂` is in every integrand pulled back to the
    surface. Under that reading every term that comes only from the motion stays identifiable under its parent's
    name.
- **Adjudication.**
  - The list's split into a rest part and a motion part, in two pieces, is a split of each surface term (item 1,
    N6). Both legs found the rates satisfying it.
  - Item 2 asks that the terms that come only from the motion be kept apart, in two pieces (N6). The derivative of a
    term of the motion is printed under that term's name, and that name is already one of the two pieces
    (`…centerMotion` or `…radiusMotion`). So every part of it, `localEvolution` included, is kept apart and filed in
    its piece.
  - The list names no parts within a derivative. The engines' split by which time dependence is differentiated is
    finer than the list asks.
  - Grok's partition is the engines' partition with each motion term's `localEvolution` added to that term's own
    piece. It can be read off the printed parts by that sum, but the engines' split cannot be recovered from Grok's.
    So the change Grok asks for would remove a printed split and give the record nothing it cannot already compute.
  - Grok's own evaluation found the totals and the pieces correct (N5).
  - So nothing outstanding changes what is computed. The limit on claims is C1.

## Carried to the record and the comparator (not findings)

- **C1. What `localEvolution` is.** `localEvolution` is the derivative at fixed geometry. Under a term of the motion
  it carries that term's undifferentiated velocity. So it is free of the motion only under the rest, flux, face,
  source and weight terms (N4). The record must not add `localEvolution` parts across parent terms and call the sum
  free of the motion (Grok 1; Claude note d).
- **C2. What the engines' own partition residuals test.** Each engine's 44 partition residuals test that the pieces
  sum to the derivative. They do not test which piece holds a term, because they share the splitter that makes the
  pieces (Claude note a; Grok). The attribution rests on the legs' independent splits (N5).
- **C3. A shared partition.** Both engines split the derivatives with the same jet-slot construction, which the same
  repair author wrote (Claude note b). A cross-engine comparison cannot see a convention they share. The legs'
  independent splits agree with both engines. The comparator compares the pieces one by one, as round 0 noted.
- **C4. Where homogeneity is printed.** The added zero reference terms leave some top-level dimension residuals empty:
  `Missing["zeroTerm"]` appears on 42 lines of the Wolfram transcript, against 8 in round 0 (N8). Homogeneity is still
  printed per object, one unit per sum, for example `count.rate.dimensions` and `dipole.*.rate.dimensions` (N8). The
  record cites those lines.
- **C5. The repair's model.** The repair ran at `gpt-6.1-sol`, the build at `gpt-6-astra`, in the same session (N1).

## Disposition: cleared

Nothing outstanding changes what the build computes, or what the record may claim within C1–C4. The S16 build is
accepted at the commit that adds this file.

## The build's review history

| Round | Claude | Grok | Outcome |
|---|---|---|---|
| 0 | NEEDS REVISION, 3 | CLEAR | B1–B3 confirmed; preserved `b60e38fb`; repair 1 by the original session |
| 1 | CLEAR | NEEDS REVISION, 1 | G1: the fact confirmed, not a defect; C1–C5 carried; accepted |

## Next

The cross-engine comparator, with its own decision list. It compares the pieces one by one (C2, C3). The step
record and its register pass follow.
