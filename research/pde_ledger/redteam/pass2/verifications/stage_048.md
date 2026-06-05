---
unit_id: 048
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T08:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 048

## Per-finding outcomes

### F1 — stale_output (self-label + stale committed transcripts)

**Classification:** resolved

**What changed:**
Single label-only edit in the SymPy source. Per the captured diff
(`exec_logs/stage_048_diff.patch`), the ONLY changed line is the docstring header at
`scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py:3`:
`Stage 31 SymPy audit.` → `Stage 48 SymPy audit.`. Confirmed in the post-fix file (line 3
now reads `Stage 48 SymPy audit.`). The `.wl` was not touched (absent from the diff; mtime
unchanged at 2026-06-03 15:59:11). Both committed `.txt` outputs were refreshed by the
orchestrator re-run.

**Assessment:**
Correct and complete. The diff is strictly label-only:
- Exactly one hunk, one logical line changed; only the stage-number token (`31`→`48`) in the
  docstring header. Surrounding docstring text (`SymPy audit.`, the numbered "Checks:" list)
  is byte-identical.
- No equation, value, variable name, symbol declaration, or assertion was modified. The
  closing pass-line (`All Stage-048 symbolic checks passed.`, py:127 / sympy out L52) and the
  already-canonical banner (`STAGE 48`, py:31 / sympy out L8) were correctly left untouched,
  matching the directive's explicit DO-NOT-TOUCH list.
- The adjacent in-source labels the auditor flagged for deferral (2-digit banner zero-pad;
  `.txt` banner stale `STAGE 031`) are correctly out of scope — the banner stale-`.txt` is
  cured by the re-run; the zero-pad is deferred to the numbering pass and was not (and should
  not have been) altered here.

Refreshed outputs confirm the fix:
- sympy `.txt` line 8: `STAGE 48 — SUPPORT COMPENSATION THEOREM AUDIT` (was `STAGE 31`).
- mathematica `.txt` line 8: `STAGE 048 — SUPPORT COMPENSATION THEOREM` (was `STAGE 031`).
- mathematica `.txt` final line 47: `Stage 048 Mathematica audit passed.` (already canonical).

Every pre-existing residual/PASS survives the refresh (no math regression): sympy shows all
`= 0` residuals (dG_tr/dxi, F_tr(0)-1, (1-xi)F_tr coeff, M_crit-G_tr, S(0)-1, dS/dzeta,
inverse maps S(zeta_req)/S(zeta_crit), pole/branch margins, dxi_phys/dzeta) and the closing
"All Stage-048 symbolic checks passed." Mathematica shows all 13 `PASS:` lines and every
`formula = 0` residual. Both engines agree on every printed closed form.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 48 — SUPPORT COMPENSATION THEOREM AUDIT` (banner now canonical)
- `inverse map S(zeta_req)-S_req = 0`, `inverse map S(zeta_crit)-S_crit = 0` (load-bearing
  inverse-map checks still pass)
- `All Stage-048 symbolic checks passed.` / `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- `STAGE 048 — SUPPORT COMPENSATION THEOREM` (banner canonical)
- `PASS: inverse map S(zeta_req)-S_req`, `PASS: branch margin formula`, `PASS: dxi_phys/dzeta formula`
- `Stage 048 Mathematica audit passed.` / `# exit_code: 0`

**Output freshness:** confirmed. Both `.txt` mtimes are 2026-06-05 08:12:44 — newer than the
sympy `.py` (07:54:26) and the untouched `.wl` (2026-06-03 15:59:11). The log header dates
(08:12:33 sympy / 08:12:37 mathematica) confirm a post-fix re-run.

## Material-change assessment

`material_change`: false. The sole edit is a docstring stage-number token; no derived result,
equation, value, or assertion changed. The `.wl` was untouched, both engines still exit 0 with
identical residuals/PASS set, and the engines still agree. No downstream unit can depend on a
docstring label. No `upstream_stale` propagation warranted on math grounds.

## Side observations (non-blocking)

- The sympy `.txt` banner is 2-digit (`STAGE 48`) while the mathematica `.txt` is zero-padded
  (`STAGE 048`); the auditor already flagged the zero-pad as deferred to the dedicated
  numbering pass. Not a blocker here.
- sympy out L31 prints the unguarded `-oo*sign(eps - 1)` before L32 resolves it to `oo` under
  `0<eps<1` — expected/known behavior noted by the auditor, not introduced by this fix.

## Verdict justification

The single low-severity `stale_output` finding is fully resolved. The captured diff is strictly
label-only (one docstring token, `31`→`48`); no equation, value, or assertion changed, and the
`.wl` plus the canonical banner/pass-line were correctly left untouched per the directive. Both
refreshed transcripts now show the canonical `STAGE 48`/`STAGE 048` banners, every prior
PASS/`= 0` residual is intact, both engines exit 0, and both `.txt` mtimes post-date their
scripts. No regressions, no tautology introduced. `verified`, `material_change: false`.
