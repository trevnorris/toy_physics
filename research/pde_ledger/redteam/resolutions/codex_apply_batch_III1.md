You are applying the user-approved resolutions for batch III.1 v2 user-gate items. The user has explicitly authorized all four directions (Q1=a, Q2=a, Q3=b, Q4=a) from the recommendations you wrote earlier in `redteam/resolutions/batch_III1_paper_alignment.md`.

# Per-question scope and authorization

## Q1 — Stage 043 Mathematica D_phi sign fix

**Authorized to edit:**
- `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl` (the `dPhi` / `dPhiExpected` definition lines, around line 52-53)

**NOT authorized to edit:**
- The SymPy script (already uses the correct kappa-first convention)
- The notes file (already uses the correct convention)
- The paper card (doesn't mention `D_phi` sign)

**Specific edits:**
- Change `dPhi = Det[{{y0, y1}, {kappa0, kappa1}}]` → `dPhi = Det[{{kappa0, kappa1}, {y0, y1}}]` (swap rows so `(kappa_0, kappa_1)` is in the first row, matching the notes' `D_phi := kappa_0 y_1 - kappa_1 y_0`)
- Update `dPhiExpected` to remove the leading minus sign (no longer needed since the determinant convention now matches notes/SymPy directly)
- Verify all assertions still pass

After editing, run `math -script mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl` to confirm exit 0 (RESPECT Mathematica single-seat rule — orchestrator is not running other math -script).

---

## Q2 — Stage 045 import Stage-044 residual for F3

**Authorized to edit:**
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py` (the F_tr collapse block around lines 172-189)
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl` (the analogous block around lines 125-136)

**NOT authorized to edit:**
- Stage 044 scripts (read-only — that's the upstream source)
- Paper TeX or notes files

**Specific approach:**
1. Read `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90` to get the canonical `F_cont` / `R_target - F_cont(xi_phys)` residual.
2. In stage 045's SymPy script, import or restate the `F_cont` expression (as a comment block citing Stage 044's file:line for traceability). Substitute tracking conditions `R_phi → R_U` (and `lambda_0 → 2/9`), then assert the result equals the notes' `F_tr` form (the displayed expression at `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:218-220`).
3. The new SymPy assertion should look like:
   ```python
   # Import Stage-044 F_cont residual; substitute tracking (R_phi -> R_U) + D/N (lam0 -> 2/9).
   # See: scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146
   F_cont_stage044 = <restate the canonical residual here using stage 045's symbols>
   F_tr_from_stage044 = sp.simplify(F_cont_stage044.subs([(R_phi, R_U), (lam0, sp.Rational(2, 9))]))
   expect_zero(
       "F_tr collapse from Stage-044 residual",
       sp.simplify(F_tr_from_stage044 - F_tr_expected),
   )
   ```
4. Mirror in Mathematica with `FullSimplify` + same substitution chain.
5. Keep or remove the existing `F_track`/`F_tr_expected` algebraic-substitution check at your discretion — if it's a useful subsidiary anchor, keep it; if it's now redundant with the Stage-044 check, remove it (cite the rationale in the Applied block).

After editing, run both stage 045 scripts. Confirm exit 0 for both. Single-seat Mathematica.

---

## Q3 — Stage 045 label relabel (scripts + notes)

**Authorized to edit:**
- `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py` (docstring line 3, banner line 31, final-print line 191)
- `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl` (banner line 26)
- `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md` (the "Stage 28" reference at line 232)

**NOT authorized to edit:**
- Any other line in either script (descriptive references to "Stage 28" in content prose are NOT file-identification labels; only metadata labels are in scope)
- Paper TeX (already correct)
- Any other stage's files

**Specific edits:**
- SymPy line 3: `Stage 28 SymPy audit.` → `Stage 045 SymPy audit.`
- SymPy line 31: `STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT` → `STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT`
- SymPy line 191: `All Stage-28 symbolic checks passed.` → `All Stage-045 symbolic checks passed.`
- Mathematica line 26: `STAGE 028 — COHERENT LOCAL TRACKING` → `STAGE 045 — COHERENT LOCAL TRACKING`
- Notes line 232: `## 7. Best current theorem statement after Stage 28` → `## 7. Best current theorem statement after Stage 045`

After editing, run both stage 045 scripts. Confirm exit 0.

---

## Q4 — Stage 046 notes coefficient fix (paper-only)

**Authorized to edit:**
- `notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md` (lines 87-91, 131-133, 136-139)

**NOT authorized to edit:**
- Any script (scripts are correct)
- The paper card (doesn't quote the coefficients)
- Any other file

**Specific edits in notes:**
- `P_R` (lines 87-91): change `230 R delta^3` → `162 R delta^3` and `230 R delta xi^2` → `162 R delta xi^2`
- `P_1` (lines 131-133): change `248 R delta^2 xi` → `180 R delta^2 xi` and `230 delta^3` → `162 delta^3`
- `P_2` (lines 136-139): change `237 R^2 xi^4` → `220 R^2 xi^4`

No script re-run needed. Just verify edits landed cleanly and surrounding markdown reads correctly.

---

# Output format

Append to `redteam/resolutions/batch_III1_paper_alignment.md` a new `## Apply log` section listing for each Q:

```
### Q<N> applied
- direction: <a|b|c|skip>
- files modified: <list with paths and line ranges>
- destination_verified: <yes — file:line | n/a>
- post-edit checks: <script exit codes, or "n/a — paper/notes only">
- notes: <anything surprising>
```

# Critical rules

1. **Per-question scope is strict.** Do NOT edit files outside the per-Q authorization list.
2. **Mathematica single-seat.** One `math -script` invocation at a time.
3. **No JSON output.** YAML or markdown only.
4. **No fake commentary scripts.** Read and reason directly.
5. **For Q2 — Stage 044 source cross-verification (MANDATORY)**: before importing the residual, grep stage 044's actual SymPy and Mathematica scripts to confirm the `F_cont` expression matches what you cite. Record this in the Applied block: `destination_verified: yes — file:line`.
6. **Report blocked items explicitly** — don't silently skip.

# Working directory

`/var/projects/toy_physics/research/pde_ledger`
