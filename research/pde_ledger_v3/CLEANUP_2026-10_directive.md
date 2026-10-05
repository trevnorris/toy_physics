# Cleanup and consolidation directive: S11c era (2026-09-09 → 2026-10-04)

**Author:** Claude (orchestrator), 2026-10-04. **Executor:** Codex. **Reviewer:** Claude, at every STOP below.

## Goal

Four weeks of S11c work added 4,434 files and about 1,200 commits, but only one step record (the 50-line
`steps/S11c_PARTIAL_CLOSEOUT.md`). The user cannot tell what was learned. This directive covers three things:

1. Consolidate everything into a clear ledger: what was established, what was not, and what is open and which later
   section owns it.
2. Remove what the ledger doesn't need.
3. Rebuild the branch history as a short, readable sequence.

The full history is preserved in the annotated tag `archive/pre-cleanup-2026-10-04` (→ `2f56b303`), already on
GitHub and GIN. Every old commit ID keeps resolving through it.

**Success looks like:** the user can read one short document and know where the v3 light sector stands. The
ledger's step records, plan, STATUS and TeX say the same thing. Each kept claim points at the file that backs it.

## Ground rules (they apply to every phase)

1. **No new science.** Do not run CAS engines, numerical workers or new calculations. Do not launch reviews. Make no
   new physics claims. If two records conflict, or a result looks wrong, list it under *Conflicts and questions*.
   ⛔ Do not resolve it yourself.
2. **Cite everything.** Every statement in FINDINGS or in a ledger record cites its source file (path plus line or
   section). Quote numbers exactly. Carry every qualification the source attaches: model point, frozen quantities,
   single- vs dual-engine, and review status (who reviewed it, and the literal verdict).
   - ⛔ Never upgrade a status. "Unresolved" is not a result. "Not detected" is not "no leak". "Scoped clear" is not
     "cleared".
3. **Classify every finding** as exactly one of:
   - ESTABLISHED (with its qualifications);
   - CONDITIONAL (true given stated inputs);
   - UNRESOLVED;
   - WITHDRAWN/FAILED;
   - OPEN, deferred to a named later step (S12, S22, Q2, R10, …).
4. **Plain language first.** The user is a programmer, not a physicist. Every document opens with a short
   plain-language summary. Jargon goes in the body, with the summary saying what it means.
5. **Commit hygiene.**
   - Make a few commits per phase, each with a message that says what changed.
   - ⛔ No launch, admission, readiness, receipt, hook, guard or "record that X started" commits or files.
   - ⛔ No new JSON process records.
   - ⛔ Create no documents beyond the ones named in this directive.
6. **Off-limits:**
   - ⛔ Do not edit `CLAUDE.md`, `AGENTS.md` (Claude supplies the new text, see Phase 3b), `.claude/` or the `pde_ledger_v2` tree.
   - ⛔ Do not edit the hash-chained `scripts/*_exports.py`.
   - ⛔ Do not move or delete the archive tag.
   - ⛔ Never push or force-push. Claude or the user pushes after review.
7. **Annex, ignore rules and attributes.**
   - Use `datalad save`.
   - ⛔ Never `git add -f` anything. A file that `.gitignore` excludes stays untracked: no logs, no packet
     tarballs, nothing under `_scratch/`.
   - ⛔ Never annex a `*_exports.py`.
   - Kept annexed `.out` stay annexed.
   - ⛔ Never edit `.gitattributes` or `.gitignore`, except for the Phase 3d restore.
8. **No deletion before Phase 3.** Deleting means removing a file from the tree. The archive tag keeps it.
9. **Run in the foreground.** No time limits, watchers, schedulers or completion hooks.
10. **At every STOP, write a short report (at most one page)** covering:
    - what was done;
    - commit IDs;
    - which files the reviewer should read;
    - the *Conflicts and questions* list.

    Then stop and wait.

## Phase 1: Inventory and findings → STOP

Working folder: `research/pde_ledger_v3/cleanup_2026_10/`.

**1a. Inventory.**
- Write a small committed script that lists every file added (A) or modified (M) in
  `dada3b7d..archive/pre-cleanup-2026-10-04`. That's 4,434 A and 24 M; `dada3b7d` is the last commit before S11c-d
  work began.
- Write the list to `INVENTORY.tsv` with these columns: path, A/M, bytes, workstream, kind.
- *kind* is one of: spec/amendment/contract, decision list or build directive, result report, review report or
  disposition, script that produces a reported result, worker/continuation/tooling, process record, output
  (`.out`/`.json` result), Lean proof/contract, Lean CAS bridge, doc/note, ledger front matter.
- Every file gets exactly one workstream and one kind. Put the counts by workstream × kind at the top of
  `FINDINGS.md`.

**Seed workstreams.** Refine as needed, and make sure every file is assigned:

| # | Workstream | Anchors |
|---|---|---|
| 1 | S11c-d physics spec and its amendments and contracts | spec cleared `399a8516`; pole contract `7c98b8ee`; total-transverse-loss amendment `d648e1f4`; scattering-form amendment |
| 2 | S11c-d symbolic build | directive `64b85989`; program brief; symbolic scattering/current engine; ω=1 fixed point |
| 3 | Upstream repairs and their reviews | inertia `a74da30a`; mechanical load `c643112a`; thickness coordinate (find it); c2 face trace `a05b05e3`; sheet `f9e28f5f`; Wolfram repair audit `946acc84`/`5acfdf30`; 10-01 joint review `6a95a0e9` (open c2 `[0,2]` mixed term) |
| 4 | Numerical radiating track and the ω=3 benchmark | `136b73c5` |
| 5 | Near-unity track | uniform check `098d910b`; grazing; defect weak-form, packet, receiving and centre-drive work, 10-01 → 10-04 |
| 6 | Light no-leak clean-condition packet | A/B v1–v5, paused at `af9d4d55` |
| 7 | Lean S10/S11 and the formalization policy | `c2bdb663`, `c215b0d6`; the S10 CAS byte bridge under `lean/s10/S10Audit/CAS/` |
| 8 | Exploratory throat/EM notes | archived in `2e62f660` |
| 9 | Muonium gravity corrections | `00c296a7` |
| 10 | Execution infrastructure | guards, runners, resource/review policies, AGENTS.md history |
| 11 | Ledger front matter and other notes | STATUS, V3_STEP_PLAN, CLAUDE.md, skills, prior-art and remote-compute notes |

**1b. `FINDINGS.md`.** This is the clear picture. Read the reports, dispositions, reviews, specs and step records,
not the process JSON.
- **A. Plain-language summary (at most one page).** Where the v3 light sector stands now. What S11c set out to
  answer, what it found, and why it closed as partial. Name the next step (S12).
- **B. Per workstream, in at most about half a page each:**
  - the question it asked;
  - its findings, each classified per rule 3 and cited;
  - what failed or was withdrawn;
  - what is open and which later step owns it;
  - review status;
  - the evidence files that any kept claim depends on.
- **C. Conflicts and questions:** places where records disagree with each other, or with
  `steps/S11c_PARTIAL_CLOSEOUT.md`, `STATUS.md` or `V3_STEP_PLAN.md`.
- **D. Proposed ledger changes:** which records to create, rewrite or correct in Phase 2, and why. Prefer folding
  into existing canonical records over new files: one canonical record per step.
- **E. Proposed keep rules:** which classes of file the ledger needs. No per-file list yet.

**STOP.** Claude reviews FINDINGS. The user approves Part D.

## Phase 2: Ledger records → STOP

Apply the approved Part D. At minimum:
- **The S11c-d record and the S11c closeout.** One canonical S11c-d record, and the closeout reconciled with it.
  Include the uniform near-unity result, the ω=3 benchmark exactly as its report states it, the repair review
  status, and the open c2 `[0,2]` term.
- **The S11c-b and c2 step records.** Correct them for the repairs made after they closed, with each repair's
  review status.
- **`V3_STEP_PLAN.md`, S11/S11c sections.** Rewrite them to the current state. ⛔ Rewrite superseded prose; do
  not stack banners or annotations on it.
- **`STATUS.md`.** Rewrite it as a short front door of about 80 lines or fewer, down from 1,289: where we are, what
  is next (S12), and open debts with their owners. Old clauses are deleted, not appended to. They live in the tag.
- **Paper.** The S11c content in `paper/parts/part01_light.tex` must match the records. Build the PDF. Check that
  no reader-critical content sits in a macro field the default build suppresses (`paper/macros.tex`).

**STOP.** Claude reviews the records, with Grok as the second leg.

## Phase 3: Keep list, AGENTS.md and pruning → STOP

**3a. Keep list, by tier.** What we keep follows from what the approved ledger records rely on, not from what
exists. The archive tag is the durable home of everything else. **When unsure, use Tier 2.**
- **Tier 1: claims a record states as ESTABLISHED or CONDITIONAL that later work may use.** Examples: the repaired
  S11c-b/c2 engines and exports, the uniform near-unity result, and the Lean proofs. Keep full provenance: the
  result report, the producing script, its inputs and outputs, and their transitive imports and hash pins. A Tier 1
  result must be re-runnable from the kept tree.
- **Tier 2: UNRESOLVED, WITHDRAWN, paused or exploratory work.** Keep only the final report or disposition that a
  record cites, plus any literal review report the record cites. ⛔ Do not keep its transitive chain. Cite anything
  else as `archive/pre-cleanup-2026-10-04:<path>`. A hash pin inside a pruned or Tier 2 file does not need to
  resolve.
- **Process records** are pruned unless a Tier 1 script reads them when it runs. These are launch, admission,
  readiness, receipt, hook, guard, resume, journal and checkpoint files, and review launch prompts and records.
- **Directives:** keep each step's cleared governing spec and its cleared build directive. Prune drafts, superseded
  versions, per-round review prompts and preservation copies; they are in the tag.
- **Scratch references.** If a Tier 1 file cites something under `_scratch/` or `/tmp/`, either move that artifact
  into the tree (annexed if large) or re-point the citation to the archive tag. List each case. Tier 2 citations are
  re-pointed or left as they are.
- **Lean.** Keep what `lake build` needs for the S10 proof, its coverage contract and its mutation controls.
  Remove the byte-level CAS bridge (`S10Audit/CAS/PY*.lean`, `WL*.lean`, `*Rerun*` and their bindings) only if
  the build still passes without it. The reason is CLAUDE.md L1 and L-LEAN. Report what was removed.
- Write `KEEP.tsv` and `PRUNE.tsv`. Give each kept file its tier and the record claim that needs it. Include counts
  by workstream × kind, plus a list of ambiguous files.

**3b. AGENTS.md.** Claude supplies the replacement text at this STOP. Commit it unchanged, in its own commit.

**3c. Prune**, on a local branch `cleanup/pruned` cut from the Phase 2 tip, so the step is reversible.

**3d. Restore `.gitattributes` and untrack ignored files**, on the same branch:
- Restore `.gitattributes` byte-for-byte to its `dada3b7d` version: the 8-line annex policy. The 214 added lines go,
  including:
  - 164 per-file `-diff` entries;
  - 24 entries for files under the ignored `_scratch/`;
  - 6 entries that hid diffs of the hash-chained `scripts/*_exports.py`.

  The annex policy (`out/*.out` annexed, everything else plain git) must stay as it was.
- Untrack every file that matches `.gitignore` and was added after `dada3b7d`, using `git rm --cached`. That's 33:
  22 `*.log` and 11 `*.tar.gz` review packets. They stay on disk and in the archive tag.
- Leave the 20 older tracked-but-ignored files from July/August alone, including the 8 in the v2 tree. List them in
  the report.
- Afterwards, `git ls-files -ci --exclude-standard` must list only those 20.

Then verify:
- every path cited in a kept record exists;
- every hash pin in a kept **Tier 1** file resolves to a kept file with a matching hash;
- `python3 -m py_compile` passes on every kept script;
- any existing export-chain integrity check passes (⛔ no engine reruns);
- `lake build` passes;
- the paper builds.

Report each check with its command and its literal output.

**STOP.** Claude reviews the keep and prune lists and the checks.

## Phase 4: History rebuild (local only) → STOP

- Create branch `ledger-v3-rebuild-clean` from `dada3b7d`.
- Add a small number of commits, grouped by workstream, that reproduce the approved pruned tree exactly:
  `git diff cleanup/pruned ledger-v3-rebuild-clean` must be empty. Changes to `CLAUDE.md` and `AGENTS.md` each go
  in their own commit.
- Annex pointers must be unchanged, and the `git-annex` branch must not be touched.
- **Fresh-clone check**, done under `_scratch/` (not `/tmp`):
  1. Clone the local repo.
  2. Add the `gin` remote.
  3. Run `git annex get --from gin` on every kept annexed `.out`.
  4. Run `lake build`.
  5. Run `py_compile` on the kept scripts.
  6. Build the paper.

  Report each step with its command and its literal output.

**STOP.** Claude verifies the result and the user approves. Then Claude or the user force-pushes to both remotes.

## Phase 5: Scratch → STOP

List every top-level item under `_scratch/` with its size, and say whether any kept file references it. Propose
what to delete. Delete only after the user approves; this is local disk only.

## Appendix: measured facts (2026-10-04)

These were measured as follows:
- `git diff --name-status dada3b7d 2f56b303 | cut -c1 | sort | uniq -c` → 4,434 A and 24 M;
- `git rev-list --count 5ab9caf1..HEAD` → 1,206 commits at the time;
- `git diff --name-only --diff-filter=A 63a06590 HEAD`, with each file's extension extracted, sorted and counted →
  1,969 json, 1,057 md, 682 py, 254 lean, 173 txt and 152 out among the 4,376 files added since `63a06590`;
- `git annex find --not --in gin | wc -l` → 0, and `git annex find --in gin | wc -l` → 152;
- `du -sh _scratch` → 181G;
- `wc -l STATUS.md` → 1,289;
- `wc -l .gitattributes` → 222, against 8 at `dada3b7d`; `grep -c ' -diff$' .gitattributes` → 164;
  `grep -c '_scratch' .gitattributes` → 24;
- `git ls-files -ci --exclude-standard` → 53 tracked files that match `.gitignore`. Of these, 33 were added after
  2026-09-08, by commit date of addition. `git check-ignore -v --no-index` attributes 22 of them to
  `.gitignore:10:*.log` and 11 to `.gitignore:128:*.tar.gz`, about 1.9 MB in total.

`steps/S11c_PARTIAL_CLOSEOUT.md` is the only file under `steps/` added in this period.
