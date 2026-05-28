---
unit_id: 165
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 165

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:41,44`. The no-op short-circuit
`dlog_LW = sp.simplify(sp.diff(sp.log(LW_law), sp.log(a)) if False else 1)  # documentary only`
was replaced with the genuine logarithmic derivative `dlog_LW = sp.simplify(a * sp.diff(sp.log(LW_law), a))`,
and a new assertion `expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)` was added at line 44 after the two
existing `print` lines (which are unchanged). Diff matches the directive's required-change block verbatim.

**Assessment:**
Correct and non-tautological. `LW_law = πa√(1+r²)/(2√3)` is linear in `a`, so `a·d/da log(LW_law) = a·(1/a) = 1`,
and `dlog_LW - 1` simplifies to exactly 0 — confirmed by SymPy log line 7 `d ln L_W - d ln a at fixed r_* = 0`.
The assertion is sensitive: if the `a`-exponent of `LW_law` were anything other than 1 the residual would be nonzero
and `expect_zero` would raise. SymPy now exercises the headline deliverable D1 to the same degree as the Mathematica
script (log lines 6-7, `PASS`). The directive's note about not differentiating w.r.t. `sp.log(a)` was heeded
(differentiation is w.r.t. `a`). No collateral edits in this block.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Both engines. SymPy `:82-89`: a four-line comment block was added explaining that `channel_g`/`channel_r` are the LHS
of `eq_g`/`eq_r` at the solved `(dv,dT)` and vanish by construction; the two `expect_zero(...)` assertions were
demoted to `print("fixed-g channel (solver consistency) =", channel_g)` / `print("fixed-r channel (solver consistency) =", channel_r)`.
Mathematica `:68-74`: parallel three-line comment block added; the two `expectZero[...]` calls were replaced with
`Print["fixed-g channel (solver consistency) = ", fmt[channelG]]` / `Print["fixed-r channel (solver consistency) = ", fmt[channelR]]`.
The `channel_g/channel_r` (and `channelG/channelR`) definitions, `eq_r/eq_g/sol`, and all substantive ratio/product/n=5
`expect_zero`/`expectZero` calls are untouched. Diff matches the directive verbatim.

**Assessment:**
This is the intended demotion, not a regression. The directive's "minimal, safe variant" was selected because the
Stage-249 δ⊥ coefficients are not available from this stage's notes, so fabricating a substantive δ⊥ construction would
be unsafe; downgrading to a clearly-labelled solver-consistency print is the correct, conservative resolution. The exec
logs confirm the intent: SymPy log lines 23-24 and Mathematica log lines 36-37 now read
`fixed-g channel (solver consistency) = 0` / `fixed-r channel (solver consistency) = 0` as plain prints with NO `PASS:`
prefix. The transcript no longer over-credits a tautology as a verified PASS. The substantive PASS lines (D1, ratio,
product, n=5) remain present and unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L7 `d ln L_W - d ln a at fixed r_* = 0` (new F1 assertion residual — passes)
- L16 `ratio law - (2 d ln c_s - d ln a) = 0`, L17 `product law - (dZ + 3 dcsw - drho - 4 da) = 0` (substantive, unchanged)
- L22 `n=5 product law - (dZ + 5 drho - 4 da) = 0` (substantive, unchanged)
- L23-24 `fixed-g channel (solver consistency) = 0` / `fixed-r channel (solver consistency) = 0` (F2 demotion — plain prints, no PASS)
No traceback / AssertionError; the script ran to completion through the carry-forward block, so exit 0.

**Mathematica:** exit=0. Notable lines:
- L6-7 `d ln L_W - d ln a at fixed r_* = 0` then `PASS: d ln L_W - d ln a at fixed r_*` (D1, genuine, pre-existing)
- L20-23 ratio/product `PASS` lines (substantive, unchanged)
- L30-31 `n=5 product law ... = 0` then `PASS` (substantive, unchanged)
- L36-37 `fixed-g channel (solver consistency) = 0` / `fixed-r channel (solver consistency) = 0` (F2 demotion — plain prints, no PASS)
- L47 `Stage 165 Mathematica audit passed.`
The `fail[]` helper exits 1 and prints `FAIL:`; none appears, the success line printed, and the script ends with `Exit[0]` → exit 0.

**Output freshness:** confirmed. Script mtimes 2026-05-28 15:54:48 (`.py`) / 15:54:55 (`.wl`); output mtimes
2026-05-28 16:10:12 (sympy `.txt`) / 16:11:19 (mathematica `.txt`). Both `.txt` outputs are ~15 min newer than their
scripts, so they were regenerated after Codex's edits. (The MANIFEST's recorded mtimes at lines 5517/5521/5525/5529 are
stale pre-fix values and were not relied upon; the on-disk stat values and log content are the authoritative evidence.)

## Material-change assessment

`material_change`: false.

Neither edit changes any derived result. F1 adds a verification of a law (`δlnL_W = δln a`) that was already used as
`subs(dLW, da)` and already verified by the Mathematica engine; the residual is identically 0 for the unchanged
`LW_law`, so no downstream constant or drift law moves. F2 only changes how two by-construction-zero quantities are
reported (PASS → plain print); their values (both 0) and all downstream-carried formulas (`dv_DN`, `dT_DN`, ratio,
product, n=5 forms in log lines 14-22 / 18-30) are unchanged and still match the auditor's hand-derivation and the
notes' boxed forms. No downstream unit depends on a value that changed.

## Side observations (non-blocking)

- The copy-paste banner artifact `"STAGE 148 — EXACT LOWER-BRANCH DRIFT LAWS"` (sympy L30, wl L26) noted by the auditor
  is still present in both scripts. The auditor explicitly did not raise it as a finding, so it is out of scope here;
  flagging only for the record. It is cosmetic and does not affect any assertion.
- The Mathematica script remains a line-by-line transliteration of the SymPy script (the auditor's noted-but-not-raised
  `mathematica_transliteration` observation). As the auditor reasoned, this carries little physical content for a unique
  2x2-linear-solve procedure; not a verification blocker.

## Verdict justification

Both findings are `resolved`. F1's no-op was replaced with a genuine, sensitivity-preserving logarithmic-derivative
`expect_zero` whose `= 0` residual is confirmed in the refreshed SymPy log, bringing SymPy coverage of the headline D1
deliverable up to parity with Mathematica. F2's two tautological assertions were correctly demoted to clearly-labelled
solver-consistency prints in both engines (no `PASS:` prefix in either log), which is the intended fix rather than a
regression. Both diffs match the directive's required-change blocks with no collateral edits, both refreshed exec logs
show exit 0 with all substantive PASS lines intact, and the outputs are confirmed fresher than the scripts. No
regressions, no material change.
