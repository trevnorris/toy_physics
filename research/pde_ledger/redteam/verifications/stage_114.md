---
unit_id: 114
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

# Verification — unit 114

## Per-finding outcomes

The original auditor report (`redteam/reports/stage_114.md`) returned `verdict: clean` with `findings_count: 0`. No directive was generated under `redteam/directives/stage_114.md` because there were no findings to address. The only Codex action was a Cluster C banner sweep edit applied independently of any per-stage directive:

- `scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py` line 18: banner string updated from `STAGE 97 ...` to `STAGE 114 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL`. Verified by re-reading the file.
- `mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl` line 26: banner string updated from `STAGE 097 ...` to `STAGE 114 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL`. Verified by re-reading the file.

Both edits are cosmetic header-string changes only; no algebraic content, no assertion text, no symbol definitions, and no targets were touched. They resolve the stage-numbering inconsistency the auditor flagged as a relabeling artifact (not a math finding) in the "Attacks tried that failed" section of the audit report.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L3: banner now reads `STAGE 114 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL` (banner sweep landed)
- L11: `Schur form identity = 0`
- L12: `low-frequency normalized outlet identity = 0`

Both assertions (A1, A2) simplify to literal `0` and pass `expect_zero`. No error or warning text.

**Mathematica:** exit=0. Notable lines:
- L3: banner now reads `STAGE 114 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL`
- L7: `PASS: Schur form identity`
- L9: `PASS: low-frequency normalized outlet identity`
- L17: `Stage 114 Mathematica audit passed.`

Both assertions pass; the trailing "passed" sentinel is emitted.

**Output freshness:** SymPy script mtime 2026-05-27 15:08:30 → output mtime 2026-05-27 15:18:08 (fresh, post-edit). Mathematica script mtime 2026-05-27 15:08:32 → output mtime 2026-05-27 15:24:45 (fresh, post-edit). Both outputs were regenerated after the banner-sweep edit.

## Material-change assessment

`material_change`: false.

The banner edits are pure print-string changes inside `banner(...)` / `banner[...]` calls. No symbol definitions, target expressions, matrix entries, assumptions, or assertions were altered. The derived results (`δΛ(D)`, `ρ_c`, `σ_c`, `κ_c`, `γ_c`) are byte-identical to the pre-edit run, as confirmed by comparing the assertion identities and the printed identifications against the auditor's report tables. No downstream unit can possibly be affected.

## Side observations (non-blocking)

None worth flagging beyond what the auditor already noted. The auditor's own self-test section already covered the stage-numbering artifact (item (c) under "Attacks tried that failed"), and Cluster C addressed precisely that.

## Verdict justification

Audit returned 0 findings with `verdict: clean`; no per-stage directive was needed. The only edit applied was a Cluster C banner sweep updating the stage-number prints in both engines (`.py` L18 and `.wl` L26) from the legacy "STAGE 97/097" label to "STAGE 114". Both edits are cosmetic. Both exec logs exit 0 with all assertions passing and outputs regenerated after the edits. No algebraic content moved, so `material_change` is false. Verification is `verified`.
