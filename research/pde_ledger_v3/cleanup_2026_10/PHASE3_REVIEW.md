# Phase 3 review (Claude) and fix instructions

**Reviewed:** `cleanup/pruned` at `bec92067`: commits `0dcc53dd`, `a811df3b`, `bd1119bb`, `9b98bd79` and `bec92067`.
**Reviewer:** Claude (orchestrator), as the directive's Phase 3 STOP specifies.
**Verdict: NEEDS REVISION.** Most of Phase 3 is accepted. Two workarounds and one wrong claim must be fixed before
Phase 4.

## Accepted (verified by Claude)

- **G1:** the N6 premise caveats are in the STATUS open table (line 32).
- **G2:** `DEFERRED_HEAVY_RUNS.md:38–42` has the c2 self-energy entry and points at the c2 record.
- **AGENTS.md:** byte-identical to `AGENTS_replacement.md`, in its own commit `a811df3b`.
  Command: `cmp AGENTS.md research/pde_ledger_v3/cleanup_2026_10/AGENTS_replacement.md`.
- **.gitattributes:** byte-identical to `dada3b7d`. `git ls-files -ci --exclude-standard | wc -l` → `20`.
- **Protected paths and tag:** CLAUDE.md, `.claude/` and the v2 tree are unchanged since `394767fd`. The tag resolves
  to `2f56b303`.
- **Citation edits:** in `bd1119bb`'s 59 modified files, the only non-repoint changes are removed link markup and the
  Lean README's description of the removed bridge.
  Check: a word-diff of each file, excluding added text that contains `archive/pre-cleanup-2026-10-04`.
- **Lean:** keeping `PY`, `WL`, `Support`, `BasisCompletion` and `Bindings` is correct. The reviewed S10 contract
  names the Bindings D3 link (`lean/s10/COVERAGE.md:119–122`), and CLAUDE.md L1 allows a targeted certificate for a
  named fidelity link. Removing the other 50 modules is accepted.
- **Tier 2 choices** in FINDINGS G, "Ambiguous cases", are accepted.

## Fixes

**H1. Remove the `.blob` files (`9b98bd79`).**
- `.gitignore:119` has excluded `*.pickle` since `a3262f51` (2026-06-12, repo hygiene). Renaming pickles to `.blob`
  so that they get tracked has the same effect as `git add -f`. AGENTS.md forbids it: "A file that `.gitignore`
  excludes stays untracked."
- Fix: remove the eight `.blob` files from the tree.

**H2. Do not store caches as `scripts/out/*.out`.**
- The seven `scripts/out/S11c_uniform_retained_<sha>.out` files (about 173 MB) are SQLite/pickle caches. They use a
  `.out` name so that the transcript annex rule applies. Six of them exist only in the local annex
  (`git annex whereis` → 1 copy), so the Phase 4 fresh clone could not get them from GIN.
- One of them is the committed transcript `scripts/out/S11c_d_mixing_scattering_sympy_audit.out` under a new name
  (`R100` in `bd1119bb`).
- Fix:
  - Remove the six caches.
  - Restore the transcript to its original path if a kept record cites it. Otherwise prune it under the tier rules.
  - Do not rename kept files.

**H3. Replace the uniform "replay museum" with the cited output.**
- FINDINGS G says the uniform worker uses absolute paths and refuses to reuse an existing run directory. A replay
  needs an isolated checkout rebuilt from a relocation map. So the copies under `_measurements/retained_uniform/`
  (107 files) do not make the result re-runnable from the kept tree.
- The inputs were only ever in `_scratch/`, and the archive tag does not contain them either. The honest version of
  Tier 1 here is:
  - keep the worker script, the result report, and the **output files the report cites**;
  - running the check again means regenerating its intermediate caches from the kept b/c2 engines.
- Fix:
  - Keep only the output files that the uniform result report cites, as plain JSON, in one directory named after the
    result (for example `_measurements/S11c_d_near_unity_uniform_output/`). Re-point the report's citations to it.
  - Remove the rest of `_measurements/retained_uniform/`, along with the relocation-map rows and the
    replay-reconstruction text in FINDINGS G.
  - Add one STATUS open item: "Replaying the uniform check needs its intermediate caches regenerated from the kept
    b/c2 engines and the S10/S11 export chain repaired (H4). The original caches are not in git. They stay on local
    disk under `_scratch/` (Phase 5)." Owner: S11c-d if reused.
  - In the Phase 5 list, mark the scratch run directories that hold those inputs as **keep on disk**.
- **User option B, used only if the user asks:** annex the caches under an honestly named path. That needs a
  `.gitattributes` line the user approves and a GIN push before the Phase 4 fresh clone.

**H4. The S10/S11 export mismatches are not pre-existing.**
- At `dada3b7d` both export chains matched:
  - `git show dada3b7d:research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py | sha256sum` →
    `7de1764c…`, equal to the pin;
  - the S11 producer → `352fc502…`, also equal to its pin.
- The producers were changed inside the cleanup window, and the exports were not regenerated:
  - `56595cf7` (2026-09-11) added the S10 `TRANSVERSE_RANK_DROP` audit family, plus `argparse` and the removal of
    `unlink`;
  - `035bb654` (2026-09-16) changed the S11 producer.

  Command: `git log dada3b7d..HEAD -- <producer>`.
- Fix:
  - Correct FINDINGS G: replace "pre-existing" and "just as before cleanup" with the facts above.
  - Remove `_measurements/retained_uniform/frozen-export-inputs/`. The original bytes are at `dada3b7d`, which will
    be the root of the rebuilt branch.
  - Add one STATUS open item: "S10/S11 exports pin their pre-2026-09-11 producers. Regenerate them from the current
    producers under review, or revert the producer edits. The S10 edit may be what the S10 Lean D3 link relies on,
    so a revert is not automatic." Owner: S10/S11 maintenance. Not run in this cleanup.

**H5. Do not edit the frozen method and keep a copy.**
- `frozen-method/SCATTERING_FORM_AMENDMENT.original` exists only because the hash-pinned
  `directives/S11c_d_SCATTERING_FORM_AMENDMENT.md` was edited.
- Fix:
  - Restore the directive to its original bytes and delete the copy.
  - Exempt hash-pinned frozen files from the link check. List them in FINDINGS G. Their old links resolve through
    the archive tag.

## After the fixes

- Regenerate KEEP and PRUNE, and rerun `inventory.py --check-selection`.
- Rerun py_compile, `lake build` and the paper build. Report each with its command and literal output.
- `git status --ignored` must show the eight pickles and the six caches as on disk and untracked, not deleted.
- Report the new total of tracked files against `dada3b7d` (`git ls-files | wc -l`: 9,472 at `dada3b7d`; 10,089 now).

Then STOP. Claude checks H1–H5 before Phase 4.

## Round 2: check of the H1–H5 fixes (`1848323b`)

**H1–H5: verified fixed** by Claude:
- **H1:** `git ls-files | grep -c '\.blob$'` → `0`.
- **H2:** no `uniform_retained` file is tracked. `scripts/out/S11c_d_mixing_scattering_sympy_audit.out` is back at
  its original path, and `git annex whereis` reports 2 copies.
- **H3:** `_measurements/retained_uniform/` is gone. The 31 files in `_measurements/S11c_d_near_unity_uniform_output/`
  are byte-identical (`cmp`) to `_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01/complete`.
  The nine scratch run directories in FINDINGS G exist on disk with 0 tracked files. STATUS line 31 carries the
  replay debt.
- **H4:** FINDINGS G names `56595cf7` and `035bb654` and withdraws "pre-existing". STATUS line 32 carries the
  regeneration debt.
- **H5:** the amendment is byte-identical to its `0dcc53dd` version. FINDINGS G lists it as the one frozen
  exemption.
- **Invariants:**
  - `git ls-files | wc -l` → `10001`;
  - `git ls-files -ci --exclude-standard | wc -l` → `20`;
  - `.gitattributes` and `AGENTS.md` still match.

**Phase 3 is accepted. Proceed to Phase 4** of the directive:
- branch `ledger-v3-rebuild-clean` from `dada3b7d`;
- workstream-grouped commits, with `git diff cleanup/pruned ledger-v3-rebuild-clean` empty;
- CLAUDE.md and AGENTS.md changes each in their own commit;
- annex pointers unchanged, and the `git-annex` branch untouched;
- the fresh-clone check under `_scratch/`.

The H4 export mismatch is a known OPEN item. Report its check as such, not as a failure of the rebuild. Then STOP.
