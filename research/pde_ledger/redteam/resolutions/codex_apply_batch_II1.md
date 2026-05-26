You are applying the user-approved resolutions for batch II.1 v2 paper_misalignment items. The user has explicitly authorized all three directions (Q1=a, Q2=b, Q3=a) from the recommendations you wrote earlier in `redteam/resolutions/batch_II1_paper_alignment.md`.

# Per-question scope and authorization

## Q1 — Stage 029 docstring/banner relabel (Stage 12 → Stage 029)

**Authorized to edit:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` (docstring lines 3 and 5 only)
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` (banner string at line 33 only)

**NOT authorized to edit:**
- The notes file (`notes/stages/moving_throat_pde_stage029_dynamic_loading.md`) — leave legacy "Stage 11"/"Stage 12" prose alone.
- Any other line in the scripts.
- Inline comments inside the scripts that reference physical results from "Stage 11"/"Stage 12" descriptively — those are content references, not file-identification labels.

**Specific edits:**
- SymPy docstring line 3: change `moving_throat_pde_stage12_dynamic_loading_sympy_audit.py` → `moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`.
- SymPy docstring line 5: change "SymPy audit for Stage 12 of the moving-throat PDE program." → "SymPy audit for Stage 029 of the moving-throat PDE program."
- Mathematica banner at line 33: change `banner["STAGE 012 — DYNAMIC LOADING"]` → `banner["STAGE 029 — DYNAMIC LOADING"]`.

After editing, run both scripts via `python3 scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` and `math -script mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` (RESPECT THE MATHEMATICA SINGLE-SEAT RULE — one math -script at a time across the whole system; only run Mathematica after no other Mathematica is running). Confirm exit 0 for both. Record results.

---

## Q2 — Stage 029 alpha_crit trim (script-only)

**Authorized to edit:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` (lines 189-194 region — the `alpha_crit` assertion block)
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` (lines 159-166 region — the `det(alpha_crit) == 0` assertion block)
- The script's docstring/comments may need minor edits to remove the scope bullet referencing `alpha_crit` so they don't dangle. Be surgical — don't restructure the script.

**NOT authorized to edit:**
- The paper card `paper/stages/stage_029.tex` — no paper change for Q2 (the user chose trim-scripts, not add-to-card).
- The notes file (leave alone).
- Stage 031 (or any other stage) — out of scope for this Q.
- Any unrelated assertion in stage 029.

**Destination-verification guardrail (MANDATORY)** — before deleting, grep stage 031 to confirm it owns the `alpha_crit` derivation:
1. Run `grep -n "alpha_crit\|alphaCrit" scripts/moving_throat_pde_stage031_*.py mathematica/moving_throat_pde_stage031_*.wl paper/stages/stage_031.tex`.
2. Confirm: stage 031's paper card boxes `alpha_crit = AB/(B*kappa_0^2 + A*kappa_1^2)` at `paper/stages/stage_031.tex:43` (Output line 65); stage 031 sympy script defines `alpha_crit = A*(A+DK)/T0` at line 87, asserts `lam_-(alpha_crit)` at line 94, and uses it through line 116; stage 031 Mathematica defines `alphaCrit` at line 59, asserts `lam_-(alphaCrit)` at line 61, uses it through line 71.
3. Record the verification in the applied block: `destination_verified: yes — <file:line>` for each of (.tex, .py, .wl).
4. If you can't confirm equivalent assertion in destination, BLOCK rather than delete.

**Specific edits:**
- SymPy `.py` lines 189-194: delete the `alpha_crit` definition + the `expect_zero` / `assert` block. If there's an adjacent narrative comment, delete or edit it to remove the scope claim.
- Mathematica `.wl` lines 159-166: delete the `alphaCrit` definition + the `expectZero` / determinant assertion block.
- Docstring/scope bullet: remove any line in either file that lists "exact softening threshold alpha_crit" as a stage 029 deliverable (stage 031 owns it now in script terms too).

After editing, run both scripts. Confirm exit 0 for both. Record results.

---

## Q3 — Stage 035 paper + notes coefficient fix

**Authorized to edit:**
- `paper/stages/stage_035.tex` (line 71 only — the bracket polynomial inside `\left(...\right)`)
- `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md` (line 86 only — the same polynomial)

**NOT authorized to edit:**
- The scripts (`scripts/moving_throat_pde_stage035_..._sympy_audit.py`, `mathematica/moving_throat_pde_stage035_..._mathematica_audit.wl`) — they are correct.
- Any other part of either the paper TeX or the notes file.

**Specific edits:**
- `paper/stages/stage_035.tex:71`: change `81\delta^3+72\delta^2+206\delta^2\xi+297\delta\xi^2+138\xi^3` → `81\delta^3+72\delta^2+189\delta^2\xi+297\delta\xi^2+121\xi^3`.
- `notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:86`: change `81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3` → `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`.

After editing, no script re-run is required (scripts already verify the correct form). Just confirm the edits landed and the surrounding LaTeX still compiles trivially (no need to actually rebuild PDF — just verify no obvious syntax break).

---

# Output format

Append to the existing `redteam/resolutions/batch_II1_paper_alignment.md` a new `## Apply log` section listing for each Q:

```
### Q<N> applied
- direction: <a|b|c|skip>
- files modified: <list with paths and line ranges>
- destination_verified: <yes — file:line | n/a>
- post-edit checks: <script exit codes, or "n/a — paper/notes only">
- notes: <anything surprising>
```

# Critical rules (recap)

1. **Per-question scope is strict.** Do NOT edit files outside the per-Q authorization list.
2. **Mathematica single-seat.** Only one `math -script` invocation at a time across the whole system.
3. **No JSON output.** YAML or markdown only.
4. **No fake commentary scripts.** Read and reason, then make edits.
5. **Destination-verification guardrail on Q2.** If you can't confirm equivalent assertion in stage 031 scripts/paper, BLOCK Q2 — do not delete.
6. **Confirm scripts pass after Q1 and Q2.** For Q3, no re-run needed.
7. **Report blocked items explicitly** — don't silently skip.

# Working directory

You are at `/var/projects/toy_physics/research/pde_ledger`.
