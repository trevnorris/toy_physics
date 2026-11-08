# S15 reading build review, round 4: disposition

Round 4 reviews repair 4 of the S15 reading build: the SymPy and Wolfram engines, the Wolfram runner, the two final
transcripts and the run record, as the builder left them when its 4o session ended (exit 0, 2026-11-08 05:34Z).
Codex wrote all of them, so the legs are a fresh Claude agent and Grok (G1). Commands and their literal output are in
`S15_reading_build_review_r4_lookups.md` (L1–L9).

## The gate record

- **Artifact.** The eight files in `S15_reading_build_review_r4_baseline.sha256`, hashed before the first leg launched.
  They are unchanged now (L2). The two transcripts are not committed; their reviewed copies are at
  `scripts/out/preserved/S15_separation_departure.r4.out` and `mathematica/out/preserved/S15_separation_departure.r4.out`,
  with the same sha256 (L3). Their run manifests are `S15_separation_departure_sympy_run_manifest.json` and
  `S15_separation_departure_wolfram_run_manifest.json`.
- **Prompt.** `research/pde_ledger_v3/directives/S15_reading_build_review_prompt_r4.md`, sha256 `0b2dc58d…` (L2). Both legs
  read this file. It is round 3's prompt with a new head (the repair, notes 4l, 4n and 4o) and item 0 rewritten for this
  repair. The notes `resume_s15_reading_repair_4.md` and `4d`, which the prompt names at their old location, are now in
  `research/pde_ledger_v3/directives/` with the same sha256.
- **Legs.** Grok (`grok-4.7`, effort high), 2026-11-08 05:36Z–07:22Z, report sha256 `003b20bf…`. Fresh Claude agent
  (Opus), launched after Grok ended, with Grok's directory unreadable to it during the run; report sha256 `b4cb5179…`
  (L1). Neither leg saw the other's output. Both verdicts: NEEDS REVISION, on the same single finding.

## Z1 (both legs; confirmed): the six two-drain region objects get item-4 conditions computed in the final Wolfram run, while the label and the run record say otherwise

**The legs.** For each `region/component_{1,2,3}/{momentum,rate}_labels_{drain1, drain2}`, the final Wolfram run started
an item-4 kernel. That kernel computed the sufficient parameter branches, verified them by substitution, and printed
them as `condition_N` with `Q_under_condition_N`, before it printed the third-kind label. Sixty conditions in all.
The label says "the exponent-sensitive attempt was not made in this final run. The retained sufficient conditions are
printed above." The run record says the six "were not attempted in these final segments; each retained its previously
printed sufficient conditions and corresponding Q objects." Both legs found the printed conditions true. The Claude leg
found them identical to round 3's by value, recomputed and not carried over.

**Verified.**
- The runner runs an item-4 kernel for every object, the six included, and passes the label text to it (L5).
- The engine prints `condition_N` and `Q_under_condition_N` (`.wl:1364–1365`) before it reads the label
  (`.wl:1370–1373`). The exponent-sensitive attempt (`.wl:1377`) comes after the label's `Return` (L6).
- The run manifest records six item-4 kernels for the six, on source `366298d8…`, return code 0, 15–71 s each (L7).
- The transcript prints 4, 16, 4, 16, 4 and 16 conditions for the six, each followed by one `label_kind` and one
  `status` line (L8).
- The run-record sentence is as quoted (L9).

**Cause.** The cause is on the orchestrator's side: my notes conflict. Note 4 says "Every condition printed now that is
true stays printed." Note 4n says "These six objects are not attempted again in this build," and asks for a label
that says "the object was not attempted in that run." (L4). The builder kept the true conditions printed by computing
them again in the final run, and then described them as retained and the six as not attempted. That description is
not what the run did. The defect is in what is claimed, not in any value: both legs found the 60 conditions true.

**Disposition: repair.** The sixty conditions are true, and each object's kernel ran in about a minute. Note 4's rule
keeps a true condition printed, and the user's standing rule is not to delete computed results on a judgement of
their worth. So the conditions stay printed, and the claims are corrected. What must be true, for each of the six, in
the final Wolfram transcript and in the run record:

- they say that the final run computed the printed sufficient conditions, and `Q` under each;
- they say that the exponent-sensitive attempt was not made, and that the complete conditions were not obtained;
- the label is the third kind, set by the user's decision;
- nothing says that the object was not attempted in that run, or that its conditions come from an earlier run.

For these six objects, this replaces 4n's requirement that the label say the object was not attempted in that run.
Nothing else about the six changes. This is a technical repair under the standing go-ahead for repair rounds. The user
decision behind it was to stop the long attempts on the six, and it still holds: no exponent-sensitive attempt is made
on any of them.

## Checked and found sound by the legs

Each item names the leg or legs that checked it. Scripts and literal stdout are in each leg's scratch directory, under
`/var/projects/_review_scratch/s15-reading-r4-{claude,grok}/`.

- **Wolfram item 4 is complete.**
  - Claude: on 12 representative objects, Wolfram's condition sets equal the leg's own exact sets in both directions.
    The 12 cover carried, centre and radius at both surfaces, region momentum and rate, and the mixed carried group.
  - Grok: on four drain-1 objects, the printed conditions are the leg's own exponent relations.
- **`Q` under each condition** equals the leg's own value.
  - Claude: four objects, 376 conditions, ratio 1.
  - Grok: carried, radius, centre and region momentum, by exact routes. One region-momentum `Q` is confirmed
    numerically only, because the symbolic reduction left an `AppellF1` unreduced.
- **Series.** The `4l` series is not truncated in place of any object (both legs).
- **Dependence.** On copies, Wolfram's item-4 conditions depend on its item-3 closed form.
  - Claude: dropping the special-exponent branches cuts carried c1 from 12 conditions to 2; the corrupted family
    `r^(2a)` moves every exponent condition.
  - Grok: replacing the geometric derivative with 0 removes the exponent condition.
  - SymPy FORM ablations move the printed output in both legs.
- **Labels.**
  - Wolfram carries no 6-hour and no out-of-memory label (both legs).
  - The three objects that ran out of memory at a shared 8 GiB cap were rerun alone at 56 GiB, with no other guarded
    unit running. Each returned in under 6 hours (Grok: `journalctl --user` and the guard records).
  - SymPy prints the fourth-kind label on all 336 item-4 objects, with no conditions (both legs).
- **Nothing else changed.** Items 1, 2, 3, 5 and 6, and Wolfram's item 4 where `d` is absent, are unchanged from
  round 3 by value (both legs).
- **Reproducibility.** Claude: a guarded SymPy rerun matches by value; a Wolfram rerun of three objects matches line
  for line after renaming.
- **Transcripts.** Both transcripts end complete. Every Wolfram kernel ran source `366298d8…`. Their sha256 match the
  run record and the manifests (both legs).
- **Run record.** It agrees with the guard records and the attempt history, except for the sentence in Z1 (Claude).

## Recorded, no repair required

- **Partial blocks in the Wolfram transcript.** The three shared-cap attempts that ran out of memory left their
  partial blocks inline and unmarked. The Claude leg measured each as an exact prefix of its object's final block, so
  no value differs. Grok found the same pattern. A reader of the transcript sees those objects' first conditions
  twice, with the same values. This changes no computed value or claim.
- **A zero by construction.** `Q_restriction_residual` is zero by construction (`.wl:555`). It appears 288 times among
  `exact_attempt` intermediate lines, and nothing claims it (Claude).

## Carried to the S15 close (not part of this repair)

- **Item 4 has one engine.** By the user's decision in 4o, item 4 has only the Wolfram engine in this build. The S15
  record must say so.
- **More-memory list.** The six two-drain objects go to `DEFERRED_HEAVY_RUNS.md`, from the run record.
- **Scratch citations.** The builder's run record, attempt history and Wolfram run manifest cite paths under
  `_scratch/` (guard records and work files). By the user's rule of 2026-11-02, nothing committed cites `_scratch`.
  What the S15 record relies on moves out before that record is accepted.

## Next

1. This commit preserves the reviewed baseline. It is not acceptance.
2. Write repair note 4p. It is orchestrator-written, so its legs are Codex and Grok.
3. The builder resumes under 4p.
4. Round 5 of the build review follows.
