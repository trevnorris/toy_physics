---
unit_id: 104
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 104

## Per-finding outcomes

The original auditor report (`redteam/reports/stage_104.md`) has `verdict: clean` with `findings_count: 0`, so there is no `redteam/directives/stage_104.md` file (confirmed: directory listing shows `stage_102.md`, `stage_105.md` but no `stage_104.md`). No per-finding rework was required.

The only edit that landed for this unit is the Cluster C banner sweep applied to both engines. Inspected via `git diff`:

- `scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py` line 18: `banner("STAGE 87 — ...")` → `banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")`.
- `mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl` line 26: `banner["STAGE 087 — ..."]` → `banner["STAGE 104 — EXACT OUTGOING l=2 DTN FINGERPRINT"]`.

Both edits are single-line cosmetic banner updates; no derivation, assertion, or symbolic line was touched. The two diffs are exactly the one-line banner change each.

## Exec log assessment

**SymPy:** exit=0 (inferred — log ends with the `RESULT:` block and the script raises `AssertionError` on any non-zero residual; all nine `expect_zero` checks print `= 0`). Notable lines (from `scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt`):

- L3: `STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT` (banner reflects the edit).
- L22–27: `static DtN slot = 0`, `Y z^2 coefficient = 0`, `Y z^4 coefficient = 0`, `Y imag z^5 coefficient = 0`, `Y z^6 coefficient = 0`, `Y imag z^7 coefficient = 0`.
- L35–37: `omega^2 coefficient = 0`, `omega^4 coefficient = 0`, `imag omega^5 coefficient = 0`.

**Mathematica:** exit=0 (script ends with `Exit[0]`; the log's final line `Stage 104 Mathematica audit passed.` confirms no `Exit[1]` from `fail[]`). Notable lines:

- L3: `STAGE 104 — EXACT OUTGOING l=2 DTN FINGERPRINT`.
- L9–19: nine `PASS:` lines (`static DtN slot`, `Y z^2`/`z^4`/`imag z^5`/`z^6`/`imag z^7`, `omega^2`/`omega^4`/`imag omega^5`).
- L20: `Y_2^out(omega) ... = 1 + (aThroat^2*omega^2)/(9*cSound^2) + (4*aThroat^4*omega^4)/(81*cSound^4) + ((I/27)*aThroat^5*omega^5)/cSound^5` — matches the notes' canonical ω-frame form.

**Output freshness:** confirmed.

- sympy script mtime: 2026-05-27 15:08; sympy output mtime: 2026-05-27 15:18 (output newer than script).
- mathematica script mtime: 2026-05-27 15:08; mathematica output mtime: 2026-05-27 15:24 (output newer than script).

Both outputs were regenerated after the banner edit.

## Material-change assessment

`material_change`: false.

The banner sweep is a cosmetic display string and changes no derived quantity, assertion residual, or symbolic expression. Downstream units that consume any numeric output from stage 104 (e.g., the canonical odd coefficient `a^5/(27 c_s^5)` used in stage 105 onward) see byte-identical values; the only diff is in the banner header printed at the top of the log.

## Side observations (non-blocking)

- The SymPy banner spells the DtN suffix as `DtN` while the Mathematica banner spells it `DTN`. This was already present pre-edit (the prior STAGE 87 banners had the same casing split) and is cosmetic only; the auditor did not flag it and it is out of scope here.
- The Mathematica script uses symbol names `aThroat`/`cSound` while SymPy uses `a`/`c_s`; this is a long-standing engine-side naming choice (no impact on the numeric/symbolic identity) and is unchanged by this edit.

## Verdict justification

The auditor closed unit 104 with `verdict: clean` and zero findings. The only post-audit edit was a Cluster C banner sweep (`STAGE 87`/`STAGE 087` → `STAGE 104`) in both engines, confirmed in the working-tree diff. Re-generated exec logs (post-edit mtimes) show the updated banner and every assertion (9 SymPy + 9 Mathematica) reports a zero residual / `PASS`. No derivation, assertion, or numerical value changed, so there is no downstream impact. Verified.
