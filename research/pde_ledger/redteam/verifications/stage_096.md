---
unit_id: 096
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T14:40:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 096

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
SymPy script `moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py` lines 88-91 now read
`eps_2 = sp.Integer(0); eps_4 = sp.Integer(0); print("eps_2 =", eps_2); print("eps_4 =", eps_4)` — the two
`expect_zero("eps_2", …)` / `expect_zero("eps_4", …)` calls are deleted and replaced by `print` statements, exactly
as the directive prescribed. The previously tautological `expect_zero("zeta_req - c_pole/c_geom", …)` line is
also gone (no longer present; original line 117 deleted; the script jumps from `expect_zero("zeta_req - 1/3", …)`
at line 116 straight to the FINAL LEDGER block at line 118). Mathematica mirror `…_mathematica_audit.wl` lines
57-58 retain `eps2 = 0; eps4 = 0;` but lines 74-75 are now `Print["eps2 = ", fmt[eps2]]; Print["eps4 = ", fmt[eps4]];`
in place of the two `expectZero[…]` calls, again exactly per directive.

**Assessment:**
Correct. The deleted assertions were trivially-true by construction (`Integer(0) == 0`, `(c_pole/c_geom) - (c_pole/c_geom) == 0`).
The `print` replacements preserve the carried numerical values in the transcript without dressing them as assertions.
Output transcripts confirm: SymPy `.txt` lines 35-36 show `eps_2 = 0`, `eps_4 = 0` (bare values, not residuals);
Mathematica `.txt` lines 48-49 show `eps2 = 0`, `eps4 = 0` (no `PASS:` prefix). No stray `PASS: eps_2` / `PASS: eps2`
or `zeta_req - c_pole/c_geom = 0` lines remain anywhere in either transcript (grepped). No collateral edits.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy script lines 127-131 now contain the four-line `HYPOTHESIS CARRIED` print block (blank line + header + three
condition lines), inserted immediately after the FINAL LEDGER block ending at line 125. Mathematica `.wl` lines
85-89 mirror it (`Print[""]; Print["HYPOTHESIS CARRIED"]; Print["These results are conditional…"]; Print["and the…"]; Print["not an unconditional…"];`)
placed between `Print["Stage 096 Mathematica audit passed."];` and `Exit[0];`, matching the directive verbatim.

**Assessment:**
Correct. These are pure prints, not assertions; both engines exit 0 unchanged. SymPy `.txt` lines 57-60 and
Mathematica `.txt` lines 63-66 both display the four-line block verbatim. The conditioning language echoes the
paper card's downstream-use clause exactly as the directive specified.

### F3 — stale_output (banner)

**Classification:** resolved

**What changed:**
SymPy script line 53 now reads `banner("STAGE 096 — GEOMETRY-LANE CHECK VERDICT")`; Mathematica `.wl` line 26 now
reads `banner["STAGE 096 — GEOMETRY-LANE CHECK VERDICT"];`. Additionally, the directive's note documents an
orchestrator-direct docstring fix: SymPy script line 9 reads "Stage 092 obstruction formula" (no longer the stale
"Stage 75" reference). The "Stage 77" reference noted as acceptable in the original report is no longer present
either — the docstring now only references Stage 092 upstream.

**Assessment:**
Correct. Grep for `STAGE 079`, `STAGE 075`, `Stage 75`, `Stage 77` in both scripts returns no hits. Output
transcripts now both lead with `STAGE 096 — GEOMETRY-LANE CHECK VERDICT` at line 3 (SymPy) and line 3 (Mathematica).
Pure cosmetic correction, no math touched.

## Exec log assessment

**SymPy:** exit=0 (script-side `expect_zero` raises on residual ≠ 0; transcript ends cleanly with the HYPOTHESIS
CARRIED block at line 60 and no `AssertionError`). Notable lines: `<Y00|Y20> = 0` (line 5), `c_pole - 1/4 = 0`
(line 42), `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] = 0` (line 44), `zeta_req - 1/3 = 0` (line 46).
5 non-trivial `expect_zero` residual lines remain after F1's deletions; orthogonality block (lines 5-34) contributes
15 residual = 0 lines, each printed twice (script first `print`s the value then `expect_zero` re-prints it), for
30 total ortho lines — matches script structure.

**Mathematica:** exit=0. 20 `PASS:` lines as expected (15 in Section I orthogonality + 5 in Section II arithmetic);
no `PASS: eps2` or `PASS: eps4` lines remain (F1 cleanup confirmed). Notable lines: `PASS: <Y00|(-Delta)Y22s>`
(line 38), `Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)] = 0` / `PASS:` (lines 54-55), `Stage 096
Mathematica audit passed.` (line 61).

**Output freshness:** SymPy script mtime 2026-05-27 11:14:56; SymPy output mtime 2026-05-27 14:28:40 (newer,
post-fix run). Mathematica script mtime 2026-05-27 11:15:02; Mathematica output mtime 2026-05-27 14:30:25 (newer,
post-fix run). Both outputs fresh.

## Material-change assessment

`material_change`: false.

No derived numerical result changed. `c_pole = 1/4`, `c_geom = 3/4`, `Yhat_Q^cons = 3/4 + (1/4)/(1-ω²/Ω_Q²)`,
`rho_alpha = 4/3`, `zeta_req = 1/3` all still hold and are still asserted via the (substitution-style but
non-tautological) residual checks. The edits removed tautological asserts, added pure prints, and corrected a
banner label — none of these propagate to downstream units. Downstream stages quoting Stage 096's deliverables
see identical values.

## Side observations (non-blocking)

- The SymPy script still prints each orthogonality residual twice per mode (once via `print(...)` at lines 78-80,
  once via `expect_zero(...)` re-print at lines 82-84). This is a transcript-formatting redundancy carried in from
  before the audit; not part of F1/F2/F3 and not affected by the fix. Flagging only because the duplicated lines
  inflate the transcript a bit.
- The directive's "Applied: F1+F2+F3" block is correctly appended at the end of `directives/stage_096.md` with
  files_changed / summary / deviation fields populated — orchestrator state advanceable.

## Verdict justification

All three findings (F1 tautological-asserts deletion + value-print substitution; F2 HYPOTHESIS CARRIED
annotation block; F3 STAGE 079 → STAGE 096 banner sweep, plus the bonus Stage 75 → Stage 092 docstring fix)
are applied verbatim to both engines. The two `eps_2`/`eps_4` tautologies are gone, the `zeta_req - c_pole/c_geom`
restatement is gone, the carried values still appear in transcripts as pure prints, both engines still exit 0
with all 20 Mathematica `PASS:` lines and all 5 surviving SymPy `expect_zero` residual = 0 lines, and the
HYPOTHESIS CARRIED block appears in both transcripts. Outputs are fresh (mtimes post-script-edit). No material
change to downstream-relevant numbers. Verdict: `verified`.
