---
unit_id: 191
batch: V.3
verifier_model: Opus 4.8 (1M context)
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 191

## Per-finding outcomes

The audit report (`reports/stage_191.md`) returned `verdict: clean` with `findings_count: 0`.
There is no `directives/stage_191.md` (correctly absent — nothing to fix), and the captured
`exec_logs/stage_191_diff.patch` is empty (0 bytes), confirming Codex made no source-code edit.
With zero findings there are no per-finding rows to classify; verification reduces to confirming
the clean disposition holds and the committed outputs are clean.

## Disposition confirmation (no-change clean)

**Outputs clean.** Both committed `.txt` outputs carry the canonical banner `STAGE 191`:
- SymPy: `STAGE 191 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM` (line 3); the
  closing ledger banner `STAGE 191 LEDGER` (line 333). Every residual reads `= 0` (response/prefactor
  compilers, pole defect, grouped inverses, projector matrices `[[0,0,0],...]`, isotropic-branch
  coords, `Delta_branch on isotropic one-pole normalized branch = [0;…;0]`, all Packet B round trips,
  orbit-lock collapses).
- Mathematica: `STAGE 191 -- MINIMAL PDE DATA PACKET MATHEMATICA AUDIT` (line 3); `STAGE 191
  MATHEMATICA LEDGER` (line 113). 40 `PASS:` lines, every residual `{0,…}` / `{{0,0,0},…}`, and the
  single boolean closure check `formal closure zero-set iff both packet zero-sets vanish = True`.
  No `FAIL` / `False` / `$Failed` / `Unequal` lines present.

**Disposition holds — genuinely independent `.wl`.** The audit's independent-derivation section
documents that for every load-bearing deliverable the `.wl` reaches the result by a distinct
operation: native `Series`/`Coefficient` plus an `D[]`-derivative cross-route for the Taylor moments
(vs `.py`'s closed-form-vs-`sp.series`), a `Solve`-based weighted-basis decomposition and
metric-built `Outer[Times,v,v.G]/(v.G.v)` projectors (vs `.py`'s formula + hardcoded literal
projectors), and a `Log`/`Solve` log-map inversion for Packet B (vs `.py`'s hand-written
`exp(...)` inverse). Several deliverables carry a third internal cross-route the `.py` lacks. The
only shared algebra is definitional packet ordering / substitution dictionaries, which any second
engine must restate verbatim. The B9 / printed-py V "home-stretch theorem" is a by-design propositional
tautology mirroring the notes' "iff", carries no physics weight, and is correctly not flagged.

**Value reconciliation.** 11 deliverable values checked against the notes carrier (the `.tex` card
lists deliverable names only), 0 misaligned. The numeric constant `54 G c_s^5/(5 a^5 c^5)` matches
notes verbatim on both engines; the isotropic-one-pole target substitution (`D4 = -3 D2^2/D0`,
`N0 = D0 P0_target`) is the genuine `Delta_pole=0`/`Delta_norm=0` surface, so the vanishing is an
earned consequence, not a rigged definition.

## Exec log assessment

**SymPy:** exit=0. No `fail`/`error`/`traceback`/`exception` lines; 20 `= 0` residual lines, all
exact symbolic zeros.

**Mathematica:** exit=0. 40 `PASS:` lines, e.g.:
- `PASS: response Series/Coefficient - derivative route` (`{0, 0}`)
- `PASS: Delta_branch native route - SymPy packet` (`{0, 0, 0, 0, 0, 0, 0, 0}`)
- `PASS: Delta_orbit solve zero-set - quotient lock` (`{0, 0, 0}`)
- `PASS: formal closure zero-set iff both packet zero-sets vanish` (`True`)
No failures.

**Output freshness:** confirmed. Both `.txt` outputs are dated 2026-06-09 14:05:09, newer than
both scripts (`.py` 2026-06-01 11:30:24, `.wl` 2026-06-01 11:35:04). Re-generated post-audit.

## Material-change assessment

`material_change`: false. No source-code change (empty diff patch); no derived result altered. No
downstream unit is affected.

## Side observations (non-blocking)

Card line 11 ("Mathematica audit: none yet.") lags the now-present, clean `.wl`. The audit correctly
classified this as a paper-side `.tex` doc-sync item out of red-team (scripts-only) scope, and the
USER has DEFERRED it to the paper-cleanup tracker. Non-blocking; no directive, no edit (Codex must
not touch `paper/`). Does not affect this verdict.

## Verdict justification

`verified`. The stage was audited clean with zero findings, no directive was generated and no source
edit was made (empty diff patch). Both committed outputs carry the canonical `STAGE 191` banner and
read all-PASS / `= 0` / `True`; the SymPy and Mathematica exec logs exit 0 with no failures; outputs
are fresh (re-generated 2026-06-09 14:05, newer than both scripts). The independent-derivation
disposition holds (distinct `Series`/`D[]`/`Solve`/metric-projector/log-map routes) and value
reconciliation is complete (11/11 match). The card-text "Mathematica audit: none yet." lag is a
USER-DEFERRED paper-side item, non-blocking. Nothing looks off.
