---
unit_id: 057
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T12:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 057

The original auditor report raised two low-severity `stale_output` items (F1 stale
transcripts, F2 stale SymPy self-labels). The directive consolidated these into a
single applied finding (directive `findings_count: 1`): the F2 self-label fix plus
the F1 output refresh. I verify them together below.

## Per-finding outcomes

### F1 (directive) — stale_output (SymPy self-labels + output refresh)

**Classification:** resolved

**What changed:**
The captured diff (`exec_logs/stage_057_diff.patch`) touches exactly two source
lines in `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`:
- py:3 docstring `Moving-Throat PDE — Stage 40 SymPy audit.` → `Stage 57 SymPy audit.`
- py:134 `print("\nStage 40 audit passed.")` → `print("\nStage 57 audit passed.")`
Plus the regenerated `scripts/output/...sympy_audit.txt`. I confirmed the current
file state: py:3 reads `Stage 57`, py:134 reads `Stage 57 audit passed.`.

**Assessment:**
Correct and exactly matches the directive's required change.
- (a) Required change applied: both self-labels are now `Stage 57`, 2-digit format
  preserved (not padded to 057), as the directive specified.
- (b) Label-only / strip-the-number identical to HEAD: the diff body is purely the
  string `40`→`57` on two comment/print lines; no code, symbol, assertion, or
  numeric expression is touched. Stripping the stage number from both versions
  yields identical source.
- (c) The already-correct cross-refs at py:83-84 (`Stage 056` Pe-monotonicity
  carry-forward comment) were left untouched — verified in the current file.
- (d) The `STAGE 57` banner (py:33) was correctly left un-padded per the directive's
  explicit "Do NOT pad" instruction; it remains `STAGE 57`.
- (e) No math/assertion/value changed: the diff shows only the two label lines, and
  the refreshed SymPy `.txt` result lines (x, A_K, zeta_phys, both partials, the two
  identity=0 lines, zeta_max, kappa_max, Omega_req, y_req, kappa_req, both
  solve-vs-closed-form identities) are byte-identical to what the auditor recorded,
  with every `expect_zero` reporting residual 0 and both sweeps PASS.

The F1-from-report output-freshness sub-item is also satisfied: both `.txt` mtimes
(2026-06-05 12:22) are newer than their scripts (.py 11:56, .wl 2026-06-03), and the
banners now read `STAGE 57` (sympy txt:3) / `STAGE 057` (math txt) with closing lines
`Stage 57 audit passed.` (sympy txt:26) / `Stage 057 Mathematica audit passed.`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 57 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP` (now correct number)
- `partial_kappa identity = 0`, `partial_y identity = 0`, `kappa_max identity = 0`,
  `kappa_req identity = 0`, `y_req identity = 0` — all residual-0 identities hold.
- `partial_kappa zeta < 0 ... PASS`, `partial_Pe zeta > 0 ... PASS`
- closing `Stage 57 audit passed.`

**Mathematica:** exit=0. Notable lines:
- `STAGE 057 — ...` (banner already canonical pre-fix, untouched).
- 14 `PASS:` lines incl. the independent physical-stiffness `A_K(physical)` route,
  both partial identities, both sweeps, `zeta_max - limit`, `kappa_max identity`,
  the `zeta_max(kappa_max) - zeta_req` round-trip, and `kappa_req`/`y_req` identities.
- Two `Limit::alimv` warnings (benign, pre-existing; assumptions on the limit
  variable ignored — does not affect the residual-0 limit result).
- closing `Stage 057 Mathematica audit passed.`

**Output freshness:** confirmed. `.txt` mtimes (12:22:12) post-date both scripts.

## Material-change assessment

`material_change`: false. The only source edits are two stage-number labels in a
docstring and a closing print; no derived quantity, assertion threshold, symbol, or
numeric value changed. Downstream units cannot depend on a label string. The refreshed
transcripts carry identical result lines to the auditor's hand-verified record.

## Side observations (non-blocking)

None affecting verification. The two `partial` corroboration items the auditor noted
(`partial_Pe>0` 6-point numeric sweep as a Stage-056 carry-forward; eta-monotonicity
via the external IFT `dy/deta>0` sign) remain as honestly-annotated scaffolding — out
of scope for this label-only fix and not introduced/altered here.

## Verdict justification

The single applied finding is fully resolved: the directive's required change (py:3 +
py:134 `Stage 40`→`Stage 57`, 2-digit format kept) was applied exactly, with no
collateral edits — the captured diff is purely those two label lines. The already-correct
`Stage 056` cross-ref (py:83-84) and the un-padded `STAGE 57` banner were left as-is per
the directive. Both engines re-ran to exit 0 with every assertion PASS and identical
result values, and the committed transcripts were refreshed (mtimes newer than scripts).
This is a clean label-only / strip-the-number change with `material_change: false`.
Verdict: **verified**.
