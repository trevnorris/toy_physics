# Codex — apply batch I.2 paper-alignment resolutions

You answered 6 paper-vs-script questions in `redteam/resolutions/batch_I2_paper_alignment.md`. The user (Trevor) cross-checked your recommendations against the destination stages' paper cards and scripts, **agreed with Q3, Q4, Q5, Q6 as you wrote them, and revised Q1 and Q2 from (c) acknowledgement to (b) trim**. The revised rationale is in the questions file. Re-read the file before starting.

Your job now is to **double-check each recommendation against the current paper + scripts** (re-read the relevant files to make sure the rationale still holds), then **apply the edits** across the relevant scripts and (for Q6) potentially the scripts on the receiving end of the trim if necessary.

## Authorization scope

For this session you are authorized to edit:

- `scripts/moving_throat_pde_stage*.py`
- `mathematica/moving_throat_pde_stage*.wl`

You are **NOT** authorized to edit in this session:

- `paper/**/*.tex` — all batch I.2 questions resolved (b) or (a-add-assertion); no paper edits needed. (Stage 010 paper already exports δP_2/δP_4 and Y20 lane ratios from the I.1 v2 batch; stage 017 paper already exports the wall-only obstruction and lane signature.) If reading the paper reveals a paper edit is genuinely needed (e.g., a label that needs renaming), mark the Q with `## Apply: blocked` and skip — Trevor will resolve manually.
- `notes/**/*.md` — leave notes alone.
- `redteam/**/*` — leave the audit pipeline files alone (including this prompt and the questions file). Exception: you MAY append a `## Apply: <done|revised|blocked>` block under each Q's `## Recommendation` block in `redteam/resolutions/batch_I2_paper_alignment.md` summarizing what you actually did.
- Anything under `.claude/`, `notes/MEMORY.md`, `.git/`, build artifacts.

## Procedure (per question)

For each Q1..Q6:

1. **Re-read the cited files** in the question's "Files to read for context" block AND the destination stage's paper + script (to confirm the destination genuinely verifies the moved content before you trim the source). The destinations per question are:
   - Q1 (stage 013 trim) → destinations: stage 010 owns δP_2/δP_4 (verify in `paper/stages/stage_010.tex` eq:stage010-dP2/dP4 and stage 010 script asserts them); stage 014 owns the mechanism sieve (verify in `scripts/moving_throat_pde_stage014_*_sympy_audit.py:77-89,129-138`)
   - Q2 (stage 014 trim) → destinations: stage 010 owns δP_2/δP_4 and Compat transport (verify in `paper/stages/stage_010.tex` and `scripts/moving_throat_pde_stage010_*.py`); stage 013 owns Xi_load (verify in `paper/stages/stage_013.tex:35-37`)
   - Q3 (stage 015 trim) → destination: stage 017 owns Y20 lane signature + b=3a + wall-only obstruction (verify in `paper/stages/stage_017.tex` and `scripts/moving_throat_pde_stage017_*.py:19-110`). Stage 010 also owns Y20 lane ratios.
   - Q4 (stage 018 trim) → destinations: stage 019 owns one-pole closure + compatibility (verify in `scripts/moving_throat_pde_stage019_*.py:35-56`); stage 020 owns gate determinant + Xi_1 + wall-slope solve (verify in `scripts/moving_throat_pde_stage020_*.py:33-70`)
   - Q5 (stage 020 Y20 trim) → destination: stage 010 and stage 017 both own Y20 lane ratios.
   - Q6 (stage 021 add composed assertion) → no destination check; this is a script-side addition only.
2. **CRITICAL GUARDRAIL FOR Q1-Q5 (b-trim recommendations):** Before deleting content from the source script, GREP the destination script to confirm the destination actually verifies the equivalent assertion. If the destination has the paper claim but the destination's script does NOT verify it, do NOT delete — instead, mark the Q with `## Apply: blocked` and write `reason: destination <stage> script does not verify <claim>; trim would silently remove coverage`. Trevor will decide whether to move the assertion to the destination script or keep at source.
3. **Decide if the recommendation still holds.** If yes, proceed. If a closer read reveals the recommendation was wrong, mark the Q with a `## Apply: revised` block explaining what changed and apply the corrected fix. If you cannot confidently apply, mark `## Apply: blocked`.
4. **Apply the edits** across the relevant scripts. Cite file:line for each edit.
5. **Run the affected scripts** to confirm they exit 0:
   - SymPy: `python3 scripts/moving_throat_pde_stage{NNN}_*_sympy_audit.py`
   - Mathematica: `math -script mathematica/moving_throat_pde_stage{NNN}_*_mathematica_audit.wl`
   - Only one Mathematica process at a time across the system; the project has a strict single-seat rule.
6. **If the script fails after your edit, iterate**: read the failure, diagnose, fix, re-run. Up to ~5 iterations per stage.
7. **Append a per-Q result block** under each Q's `## Recommendation` block:

   ```yaml
   ## Apply: <done|revised|blocked>
   files_changed:
     - <path>: <one-line summary>
   sanity_check: <"sympy + mathematica exit 0" | "sympy only (no .wl change)" | "no run needed" | error description>
   destination_verified: <"yes — <file:line> asserts <claim>" | "blocked — <reason>">
   notes: <one-line, only if relevant>
   ```

## Per-question scope detail

### Q1 — Stage 013 trim
- **Trim from sympy** (`scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`): the δP_2/δP_4 definitions and their ∂/∂G_w' assertions (lines 106-119 approx, plus the assertions around 122-124 and 144-145), the mechanism sieve (lines 126-151).
- **Keep** the Taylor primitive map (z_0, z_2, z_4, n_0, mu1 moments, Xi closed-form), Xi_load `∂Xi/∂Pprime = 2/P` check.
- **Trim from mathematica** (`mathematica/moving_throat_pde_stage013_*.wl:126-166`): same δP_2/δP_4 + sieve blocks.
- **Keep n_2, n_4 derivations** in both engines — they support the Taylor map narrative even though no paper card formally owns them yet.

### Q2 — Stage 014 trim
- **Trim from sympy** (`scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`): Xi_load definition + assertion (line 67, line 126), δP_2/δP_4 definitions + assertions (lines 74-75, 127-128), `Compat` printed-only check (find and remove).
- **Keep**: K_1/H_even sieve (lines 77-89, 129-138, 144-146), compensation denominator factor checks (lines 131-133), sign-flip mutation (line 133).
- **Trim from mathematica** (`mathematica/moving_throat_pde_stage014_*.wl:113,117-126,138-140`): same Xi/δP/Compat blocks.

### Q3 — Stage 015 trim
- **Trim from sympy** (`scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:103-208`): wall-only K1/Hev specializations, Gaussian overlap concretizations, even-gate determinant, perturbed solves, Y20 ratios, grouped trace identities.
- **Keep**: K_eta exact-quadratic-recovery block (A1-A8 approx), promoted-action setup, concrete boundary-discharge probe, helpers `real_y20_square_ratio` and `grouped_trace_anomaly` can be REMOVED if no other block uses them after trimming.
- **Trim from mathematica** (`mathematica/moving_throat_pde_stage015_*.wl:104-196`): same trim — M4-M9 wall-only + Y20 + grouped.

### Q4 — Stage 018 trim
- **Trim from sympy** (`scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:25-88`): one-pole numerator identity, K_Σ closure conditions + compatibility, 2x2 even-gate determinant = 1/27, wall-stiffness/wall-inertia slope solve, Xi_1 residual.
- **Keep**: Gaussian-profile concrete sanity checks for M_Σ, K_Σ branch integrals (the trailing A18/A19 / M8 block — these are what the paper actually exports).
- **Trim from mathematica** (`mathematica/moving_throat_pde_stage018_*.wl:24-122`): same trim.

### Q5 — Stage 020 trim
- **Trim from sympy** (`scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:45-50`): the three Y20 overlap lane assertions and the helper function if no longer used.
- **Keep**: all the wall-slope solve and Xi_1 algebra (lines 30-70 ex).
- **Trim from mathematica** (`mathematica/moving_throat_pde_stage020_*.wl:34-67`): same Y20 block + `GauntIntegral` helper if unused elsewhere.

### Q6 — Stage 021 add composed assertion
- **Add to sympy** (`scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py`): after `Dcorr` is defined around line 239, add an `expect_zero` assertion that substitutes `Gamma_port → a^5/(27*c_s**5)` and compares to the paper's closed form `-I * N0 * a**5/(27*c_s**5) * omega**5`. Use the variable names already in scope; reuse Section IV's `a, c_s` symbols.
- **Add to mathematica** (`mathematica/moving_throat_pde_stage021_*.wl`): mirror the same assertion using the local symbol names (`radius`, `cS`, `n0`, `omega`).
- This is the only Q that ADDS rather than trims.

## Important global rules

- **Single Mathematica process at a time.** Never run two `math -script` invocations concurrently. Serialize across all 11 stages if needed.
- **Do not run `$RT exec-sympy` or `$RT exec-mathematica`** — they race on `redteam/MANIFEST.yaml`. Use direct `python3` and `math -script` invocations.
- **No commentary `python3 -c` scripts** for "verifying" things. Read and reason.
- **The paper builds with LaTeX** — you are NOT editing paper this session, but if reading the paper reveals a label change is needed (rare), mark blocked rather than editing.
- **Equation labels in scripts** — if a script comment references `eq:stage013-xi` etc., keep the reference intact; just trim the code being removed, not the comment anchors to other stages' labels.
- **Helper function removal** — when trimming a block that used a top-level helper like `real_y20_square_ratio` or `grouped_trace_anomaly`, check if any remaining code uses it before deleting. Better to leave an unused helper than to break a downstream block.

## When done

1. Make sure every Q has either an `## Apply: done` / `## Apply: revised` / `## Apply: blocked` block.
2. Append a `## Apply notes` paragraph at the very bottom of `redteam/resolutions/batch_I2_paper_alignment.md` with:
   - Overall summary: "X applied, Y revised, Z blocked"
   - Per-stage script exit status (sympy + mathematica)
   - Any blocked items Trevor should resolve manually
3. Exit. Do NOT commit anything to git — Trevor + Claude handle the commit at the end after review.

## Project context (carried from prior session)

- This is the moving-throat PDE project: a toy superfluid analog. Internal/private framing freely uses unification language; the paper itself is strict toy-analog.
- Mathematica `.wl` lives in `mathematica/`. SymPy `.py` lives in `scripts/`. Their saved outputs go in `mathematica/output/` and `scripts/output/` respectively.
- The red-team pipeline is the source of these resolutions. Your job here is to land them; the verifier pass happens separately later.
- Pitfall list from past Mathematica sessions is in `.claude/skills/redteam-audit/prompts/codex.md` — preemptively patch `expectZero` with the `ConditionalExpression` strip, avoid `Part[]`-on-pattern-parameter inside `Do[Module]` blocks, use `1/pi1 == 0` for pole checks not `=!= Infinity`. These should not bite this session (no new `Solve` or `Module[Do]` blocks expected), but worth knowing.
