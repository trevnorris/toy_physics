---
unit_id: 053
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 053

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Single label-only edit captured in `redteam/pass2/exec_logs/stage_053_diff.patch`:
`scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:3` docstring
`Stage 36 SymPy audit: ...` → `Stage 53 SymPy audit: ...`. Confirmed live in the
file (py:3 now reads `Stage 53 SymPy audit: exact overlap-boost window for the
lowest support lane.`). Both committed transcripts were regenerated
(`scripts/output/...sympy_audit.txt` and
`mathematica/output/...mathematica_audit.txt`, mtime 2026-06-05 12:22, newer than
both scripts).

**Assessment:**
Correct and minimal. The edit is strip-the-number identical to HEAD — only the
integer "36"→"53" in the docstring changed; no symbol, assertion, or numeric
expression touched (diff is exactly one `-`/`+` line pair on the docstring).

- (a) Directive's required change applied: yes — py:3 docstring `Stage 36`→`Stage 53`,
  2-digit format preserved (reads "Stage 53", not "Stage 053").
- (b) Label-only / strip-the-number identical to HEAD: yes — confirmed by the
  diff; the only delta is the docstring digit.
- (c) Already-correct 2-digit banner py:25 `banner("STAGE 53 — ...")` was correctly
  LEFT UNPADDED — it is unchanged (still "STAGE 53"), exactly as the directive
  instructed ("Do NOT pad the already-correct STAGE 53 banner").
- (d) No math/assertion/value changed: confirmed — all seven SymPy assertions and
  all seven Mathematica assertions still present and PASS in the refreshed
  transcripts; final symbolic forms (`Omega_alpha`, `I_W`, `pi^2/4`, endpoints,
  linear coeff) unchanged.
- (e) material_change false.

The auditor's secondary observation (the SymPy output banner self-label reading
"STAGE 53" vs the .wl's "STAGE 053") is a deliberately-preserved 2-digit format,
not a defect — the directive explicitly chose to keep py:25 unpadded. No rework.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 53 — EXACT OVERLAP-BOOST WINDOW` (banner; unpadded as intended)
- `Omega_max - pi/2 = 0`, `A_I,max - pi^2/4 = 0` (ceiling asserts pass)
- `Omega_alpha closed form = 0`, `linear coefficient - (4-pi)/(2pi) = 0` (boxed
  formula + small-alpha coeff pass)

**Mathematica:** exit=0. Notable lines:
- `STAGE 053 — EXACT OVERLAP-BOOST WINDOW` (canonical banner)
- `PASS: Omega_alpha closed form`, `PASS: linear coefficient - (4-Pi)/(2Pi)`
- `Stage 053 Mathematica audit passed.` The `Limit::alimv` warnings are the known
  benign limit-variable-assumption notices; limits still evaluate to `1` and
  `Pi/2`.

**Output freshness:** confirmed. Both `.txt` outputs have mtime 2026-06-05 12:22,
newer than the .py (12:22) and the .wl (2026-06-03 15:59). The refreshed SymPy
transcript banner reads `STAGE 53` (matching py:25), the Mathematica transcript
reads `STAGE 053` (matching wl:32). The stale 036/36 banners flagged in the audit
are gone.

## Material-change assessment

`material_change`: false. The only edit is a docstring digit; no derived result,
symbol, or assertion changed. No downstream unit can depend on this.

## Side observations (non-blocking)

None. The SymPy banner remaining at 2-digit "STAGE 53" while the .wl uses "STAGE
053" is an intentional, directive-sanctioned format difference, not an
inconsistency to flag.

## Verdict justification

The sole low-severity stale_output finding is fully resolved: Codex applied the
exact number-only docstring correction (Stage 36→Stage 53), preserved the 2-digit
format, correctly left the already-correct STAGE 53 banner unpadded, and the
orchestrator's independent re-run refreshed both transcripts (both exit 0, all
PASS lines present, banners now 53/053). The diff shows no collateral edits and no
math/assertion/value change. verdict: verified, material_change: false.
