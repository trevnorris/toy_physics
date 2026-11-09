# S15 reading build review, round 5: disposition

Round 5 reviews repair 4 of the S15 reading build after note `resume_s15_reading_repair_4p.md`, as the builder left it
when its 4p session ended (exit 0, 2026-11-08 21:37Z). The repair changed the Wolfram runner, the six two-drain region
objects' item-4 label lines in the final Wolfram transcript, the run record and both run manifests. The engines did
not change (L6). Codex wrote all of it, so the legs are a fresh Claude agent and Grok (G1). Commands and their literal
output are in `S15_reading_build_review_r5_lookups.md` (L1–L7, L4b).

## The gate record

- **Artifact.** The ten repository files in `S15_reading_build_review_r5_baseline.sha256`, hashed before the first leg
  launched; all ten are unchanged now (L2). The two transcripts are not committed. Their run manifests are
  `S15_separation_departure_sympy_run_manifest.json` and `S15_separation_departure_wolfram_run_manifest.json`; the
  final Wolfram transcript's sha256 is `7d804465…`, the SymPy transcript's `10160f89…` (unchanged from round 4).
- **Prompt.** `/var/projects/_review_scratch/s15-reading-r5-prompt/S15_reading_build_review_prompt_r5.md`, sha256
  `65268367…` (L2), read by both legs. It is round 4's prompt with a new head (repair note 4p, the round-4 baseline
  `0d9675eb` and its preserved transcripts) and item 0 rewritten for this repair. Items 1–8 were limited to changes;
  items 9–12 applied in full to the changed code.
- **Legs.** Grok (`grok-4.7`, effort high), 2026-11-08 21:47Z–22:46Z, report sha256 `dd790545…`. Fresh Claude agent
  (Opus), launched after Grok ended, with Grok's directory unreadable to it during its run; report sha256 `0d51f3ef…`
  (L1). Neither leg saw the other's output.
- **Verdicts.** Both CLEAR, with no findings (L1).

## What the legs established

Each leg derived the six objects' item-3 closed forms itself before opening the engines, and checked its derivation
against a second route. Scripts and literal stdout are under `/var/projects/_review_scratch/s15-reading-r5-{claude,grok}/`.

- **The six objects' lines.** Each block prints its sufficient conditions and `Q` under each, then the third-kind
  label, which now says that the run computed the printed sufficient conditions and `Q` under each, made no
  exponent-sensitive attempt, and did not obtain the complete conditions. Neither leg found any of the ruled-out
  wording anywhere in the final transcript.
- **The 60 conditions are true, and `Q` under each is right**, by each leg's own computation (Claude: three routes, with
  negative controls that come out nonzero). The Claude leg's probes found exponent conditions that are not printed,
  for which the geometric departure vanishes with `Q` nonzero. This is what the label says: the complete conditions
  were not obtained.
- **How the final Wolfram transcript was made.** Both legs audited the join byte by byte: the round-4 transcript with
  the six objects' item-4 blocks replaced by the six parts of guard run `repair5_wolfram_labels_01` (Wolfram source
  `366298d8…`, runner `a5b6c68c…`). Inside the blocks only the label lines differ, plus 300 `WORKING_ON` lines that
  are equal once the kernel's generated names are normalized. Everything outside the blocks is byte-identical to
  round 4. The run record and both manifests say so.
- **Nothing else changed.** No object in either transcript changed by value outside the six objects' label lines.
- **Reproducibility.** Claude reran the six kernels from copies under the guard and got all six blocks byte for byte.
  Grok reran one and matched it on every condition and `Q`.
- **Controls.** On copies:
  - Removing the label for one of the six sends its kernel into the exponent-sensitive attempt (Claude).
  - Routing it to another label kind changes `label_kind` (both legs).
  - A FORM change to drain 1's density or flow in the engine changes the computed conditions (both legs).

## Recorded, no repair required

- **Stale-looking sha256 values in the Wolfram manifest (Claude).** `execution`, `same_engine_baseline_audit` and
  `final_transcript_production.preserved_transcript_sha256` give `5779d92e…`, the round-4 transcript's sha256, which
  the preserved round-4 copy still has (L3). Those sections describe the 4o run and an audit made on that file, so the
  value is accurate history. The top-level `transcript_sha256` and the join verification give the current file's
  `7d804465…` (L3). Both legs found every byte outside the six objects' label lines unchanged, so the audit's content
  still holds. This changes no computed value or claim.
- **The attempt-history inventory (Claude)** does not list the final runs' sufficient-condition kernels on the six.
  The run record's new section and the Wolfram manifest record them. `4n` asks the run record to carry that history,
  and it does.
- **Grok ran without the memory guard (L4b).** Its sandbox could not open the guard's lock file, so its derivation,
  rerun and ablations ran under an 8–16 GiB address cap instead. No other guarded job was running, and the Claude leg
  ran after Grok ended. Grok's evidence is its scripts and literal stdout, and the guard does not bear on their
  values. This is a defect in how I launch Grok. A Grok leg that must run CAS jobs through the guard needs a launch
  that lets it take the guard's lock.

## Disposition

Nothing outstanding changes what is computed or what may be claimed (G4). The S15 reading build, after repair 4, is
**accepted**.

## Carried to the S15 close (unchanged from round 4)

- **Item 4 has one engine.** By the user's decision in `4o`, item 4 has only the Wolfram engine in this build. The S15
  record must say so.
- **More-memory list.** The six two-drain region objects go to `DEFERRED_HEAVY_RUNS.md`, from the run record.
- **Scratch citations.** The builder's run record, attempt history and run manifests cite paths under `_scratch/`, and
  so do two quotes in this round's lookups (L2's command, L4b). By the user's rule of 2026-11-02, nothing committed
  cites `_scratch`. What the S15 record relies on moves out before that record is accepted.
- **Local transcripts.** The final transcripts and the round-4 preserved copies are deleted when S15 closes.
