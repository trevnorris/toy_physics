# Codex — apply batch I.1 paper-alignment resolutions

You answered 7 paper-vs-script questions in `redteam/resolutions/batch_I1_paper_alignment.md`. The user (Trevor) agreed with all 7 recommendations. Your job now is to **double-check each recommendation** (re-read the relevant files to make sure your earlier rationale still holds), then **apply the edits** across the relevant scripts and paper files.

## Authorization scope

For this session you are authorized to edit:

- `paper/stages/stage_*.tex`
- `paper/parts/*.tex` (only if a paragraph genuinely needs to move between part-level files)
- `paper/appendices/stage_appendix_part*.tex` (only if a row needs to be re-worded to match a paper-card scope change)
- `scripts/moving_throat_pde_stage*.py`
- `mathematica/moving_throat_pde_stage*.wl`

You are **NOT** authorized to edit:

- `notes/**/*.md` — leave notes alone. If a note needs to change, append a `Blocked: paper_misalignment` note in the questions file under that Q's recommendation block and skip the notes edit.
- `redteam/**/*` — leave the audit pipeline files alone (including this prompt and the questions file). The only exception: you may append a `## Apply notes` paragraph at the bottom of `redteam/resolutions/batch_I1_paper_alignment.md` summarizing what you actually did vs. what you recommended.
- Anything under `.claude/`, `notes/MEMORY.md`, `.git/`, build artifacts.

## Procedure (per question)

For each Q1..Q7 in `redteam/resolutions/batch_I1_paper_alignment.md`:

1. **Re-read the cited files** in the question's "Files to read for context" block, plus your earlier rationale.
2. **Decide if the recommendation still holds.** If yes, proceed. If a closer read reveals the recommendation was wrong, mark the Q with a `## Apply: revised` block explaining what changed and apply the corrected fix instead. If you cannot confidently apply, mark `## Apply: blocked` with a one-line reason and skip — Trevor will resolve manually.
3. **Apply the edits** across the relevant scripts and paper files. Cite file:line for each edit you make.
4. **Run the affected scripts** to confirm they exit 0:
   - SymPy: `python3 scripts/moving_throat_pde_stage{NNN}_*_sympy_audit.py`
   - Mathematica: `math -script mathematica/moving_throat_pde_stage{NNN}_*_mathematica_audit.wl`
   - Only one Mathematica process at a time across the system; the project has a strict single-seat rule.
5. **If the script fails after your edit, iterate**: read the failure, diagnose, fix, re-run. Up to ~5 iterations per stage before marking blocked.
6. **Append a per-Q result block** to `redteam/resolutions/batch_I1_paper_alignment.md` under each Q's `## Recommendation` block:

   ```yaml
   ## Apply: <done|revised|blocked>
   files_changed:
     - <path>: <one-line summary>
   sanity_check: <"sympy + mathematica exit 0" | "sympy only (no .wl change)" | "no run needed (paper only)" | error description>
   notes: <one-line, only if relevant>
   ```

## Cross-stage considerations

These questions are not independent — please plan the order:

- **Q1 + Q2 are both on stage 001** (same scripts, same paper card). Apply both, then run stage 001 sympy + mathematica once.
- **Q3 (stage 006)** is independent.
- **Q4 (stage 007)** is independent but requires writing new `H(w)` profile code in both engines. Stage 008 (`scripts/moving_throat_pde_stage008_*.py`, `mathematica/moving_throat_pde_stage008_*.wl`) already handles `H(w)` concretely per your earlier rationale — reuse its conventions where possible (`H(w) = exp(-w^2/rho^2)` with a new positive `rho` parameter) so the cross-stage convention stays consistent.
- **Q5 + Q6 are both on stage 010** (same scripts, same paper card). Apply both. For Q5, paper gets a `δP_n` paragraph (with `\stagefield{Output}` updated and new equation labels — use `eq:stage010-dP2`, `eq:stage010-dP4` analogous to the existing `eq:stage010-du2`/`eq:stage010-du4`); scripts get the missing `δu_2`/`δu_4` assertions. For Q6, paper grows paragraphs describing the 7 clusters with equation labels (use stable labels — they'll be cited downstream in 011, 012, 013–014).
- **Q7 (stage 011)** is the most complex: trim the 5 clusters from stage 011's scripts, and verify whether stages 022/023/024 already publish the relevant signatures (per your earlier rationale, they do). If a destination stage's paper card already covers the signature, no paper move is needed — just trim the source script. If a destination paper card is missing a particular signature, add a brief equation/paragraph there with the appropriate label. If you discover a destination script does NOT currently verify a moved cluster, leave the source verification in place and append a note explaining why (we'll triage that separately).

## Important global rules

- **Single Mathematica process at a time.** Never run two `math -script` invocations concurrently. If you need to validate multiple `.wl` files, do them serially.
- **Do not run `$RT exec-sympy` or `$RT exec-mathematica`** — they race on `redteam/MANIFEST.yaml`. Use direct `python3` and `math -script` invocations.
- **No commentary `python3 -c` scripts** for "verifying" things. Read and reason.
- **The paper builds with LaTeX** — if you touch `paper/*.tex`, prefer minimal edits that respect existing `\stagefield{...}` macros, `\StageFile{...}` directives, and the `eq:stageNNN-...` label convention. Don't introduce new commands.
- **Equation labels you add should follow the convention** `eq:stage{NNN}-{short-name}` so downstream stages can `\eqref{}` them cleanly.
- **For Q7 cross-stage moves**, the destination stage's paper card and/or script may legitimately already cover the material — in that case the action is just a trim at the source, not a move.

## When done

1. Make sure every Q has either an `## Apply: done` / `## Apply: revised` / `## Apply: blocked` block.
2. Append a `## Apply notes` paragraph at the very bottom of `redteam/resolutions/batch_I1_paper_alignment.md` with:
   - Overall summary: "X applied, Y revised, Z blocked"
   - Any files you touched that weren't in the original question file:line list (cross-stage moves, label additions, etc.)
   - Any concerns Trevor should review before the next batch
3. Exit. Do NOT commit anything to git — Trevor + Claude handle the commit at the end after review.

## Project context (carried from prior session)

- This is the moving-throat PDE project: a toy superfluid analog. Internal/private framing freely uses unification language; the paper itself is strict toy-analog. Don't change the paper's framing in your edits.
- Mathematica `.wl` lives in `mathematica/`. SymPy `.py` lives in `scripts/`. Their saved outputs go in `mathematica/output/` and `scripts/output/` respectively.
- The red-team pipeline is the source of these resolutions. Your job here is to land them; the verifier pass happens separately later.
