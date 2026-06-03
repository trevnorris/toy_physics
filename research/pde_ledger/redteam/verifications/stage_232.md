---
unit_id: 232
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 232

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script; user-resolved to prefactor 100)

**Classification:** resolved

**What changed:**
The single authorized notes edit was applied in
`notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md`:
- line 153: `\Xi_\chi = 168\,\Theta_w^{(\chi)}\Lambda_\ell^2` → `\Xi_\chi = 100\,\Theta_w^{(\chi)}\Lambda_\ell^2`
- line 157: `\Xi_J = 168\,\Theta_w^{(J)}\Lambda_\ell^2` → `\Xi_J = 100\,\Theta_w^{(J)}\Lambda_\ell^2`

The `git diff HEAD` shows exactly these two lines changed (a 4-line diff: 2 removed, 2 added). The quoted decimals
`\approx 5.5548332017764099\times 10^5` and `\approx 1.2663707072528143\times 10^5` were left intact, as directed.

**Assessment:**
Correct and minimal. The edit aligns the notes' written prefactor with the value the notes' own decimals require
(`100·4.06863235008162·1365.2827 ≈ 555483.32`, matching the quoted `5.5548332017764099e5`; `168·…` would give `≈933212`,
which contradicts it). The SymPy script (`scripts/…_sympy_audit.py`) was NOT touched — `git diff HEAD` on the `.py` is
empty, and lines 149–150 still read `100 * Theta_w_* * Lambda_ell**2`, consistent with the resolution. No `paper/`
file was touched, and no other notes line was altered. The F1 resolution is a documentation correction; it changes no
derived value (the verified pipeline already used 100), hence no regeneration of downstream literals.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
A new Mathematica audit was created (untracked):
`mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl`,
covering claims M1–M7. Exec log exits 0 with 38 PASS / 0 FAIL.

**Assessment — independence (the load-bearing check):**
The `.wl` is a GENUINELY INDEPENDENT route, not a port of the SymPy `Delta_closed` algebra.
- M3/M5 obtain `Delta(Pe)` by NATIVE symbolic integration: line 91 `deltaIntegralExpr = Integrate[kernel[xVar]·sourceStable[xVar,pVar], {xVar,0,1}, …]`. Lines 97–99 assert the integral actually evaluated (`!FreeQ[…, Integrate]`); the log line `M3 Delta integral expression head = Times` confirms Mathematica produced its own closed form. `deltaByIntegral` (line 101) and the M5 `FindRoot` (line 134) both use THIS native integral — they never reference the SymPy `Jc`/`Js` decomposition (which is the `.py`'s hand-derived `Delta_closed`, lines 117–128). This is structurally distinct algebra.
- The M3 endpoint cross-checks compare the analytic endpoint formulas against the LIMITS of the native integral: `delta0FromLimit = Limit[deltaIntegralExpr, pVar→0]` (line 105), `deltaInfFromLimit = Limit[…, pVar→∞]` (line 106), plus an independent `Integrate[kernel,{x,0,1}]` uniform-source check (line 107). These verify the formula AGAINST the integral, exactly the cross-check the directive required to catch an algebra error in `Delta_closed`. Log residuals are 0 to ~56–79 digits.
- M5 uses `FindRoot` (line 134), not the SymPy bisection. M2 uses `FindRoot` for the Robin root (line 63), not bisection.
- Prefactor `cPrefactor = 100` (line 121, M4), as required.
The seeds (`n80["1.5"]`, `n80["11155"]`, `n80["2505"]`) are FindRoot guesses, not the asserted answers; asserts compare against independent decimal literals. Engine cross-check: every Mathematica value agrees with the SymPy log to the stated tolerance (e.g. `zeta_max` 2.4675297457…, `Pe_*^(chi)` 11155.726586…, ceiling gaps 9.78e-8 / 1.94e-6 match both logs).

**Non-vacuous FAIL path:** Yes. `fail[]` (lines 8–12) prints `FAIL:` and calls `Exit[1]`. All four `expect*` helpers route to `fail` on the negative branch. The residuals printed in the log are genuine nonzero small numbers (not `0`-vs-`0` of a literal against itself), confirming the comparisons exercise real computed-vs-target deltas.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `Verified Delta(Pe) closed form against direct quadrature and endpoint bounds.`; `Pe_*^(chi) = 11155.7265863205…`; `zeta_max - zeta_phys^(chi) = 9.78…e-8`; `All Stage 232 checks passed.` No failures/tracebacks (grep count 0).

**Mathematica:** exit=0. 38 PASS, 0 FAIL. Notable: `M3 Delta integral expression head = Times` (native integral evaluated); `PASS: M3 Delta_0 formula vs Pe->0 integral limit`; `PASS: M5 chi Pe - Xi Delta(Pe)`; `All Stage 232 Mathematica checks passed.`

PASS-count manifest reconciliation (no silent parser-skip): 36 `expect*` CALL sites (the `grep` raw 40 includes the 4 `:=` definition headers at lines 22/28/34/40, which are not calls); the 2 calls inside `solveBranch` (lines 144–145) fire twice → 36 + 2 = 38 runtime PASS. Exactly matches. Per-claim: M1=2, M2=5, M3=5, M4=2, M5=6, M6=2, M7=16 → 38.

**Output freshness:** confirmed. `scripts/output/…232…sympy_audit.txt` mtime 2026-06-02 22:01:47 and `mathematica/output/…232…mathematica_audit.txt` mtime 2026-06-02 22:01:47 are both newer than the `.py` (2026-05-11) and the `.wl` (2026-06-02 21:53:33). The committed mathematica `.txt` tail matches the exec log (`All Stage 232 Mathematica checks passed.`).

## Material-change assessment

`material_change`: false.

F1 is a notes-only typo correction to a value the verified pipeline already used (100); no derived number changed. F2 adds a second engine that reproduces the existing SymPy results — it confirms, it does not alter. The auditor's robustness note checks out: `Omega_Pe → π/2` as `Pe → ∞`, so `zeta_phys → A_K·π²/4 = zeta_max` and the support side saturates near the ceiling regardless of the prefactor; the qualitative verdict (support/source non-bottlenecked) is unaffected. No downstream unit consumes a numeric constant from this stage that changed. No staleness flag warranted on substantive grounds.

## Side observations (non-blocking)

- The orchestrator-captured `redteam/exec_logs/stage_232_diff.patch` is 0 bytes (empty). I verified the actual edits directly via `git diff HEAD` and `git status` instead: the F1 notes 2-line change is present, the SymPy `.py` is unchanged, and the `.wl` is the only new file for this unit. This is a tooling artifact in the capture step, not a problem with the fix; flagging so the orchestrator can check the diff-capture path.
- The SymPy script retains its own `Delta_closed` hand-derived form (lines 117–128); that is fine — the `.wl`'s native `Integrate` is precisely the independent corroboration of that closed form.

## Verdict justification

Both findings are resolved. F1 is the exact authorized notes edit (168→100 on lines 153/157, decimals intact, script and paper untouched). F2 is a genuinely independent Mathematica route: `Delta(Pe)` is obtained by native `Integrate` of `kernel·source` (head `Times`, not unevaluated), the endpoint formulas are cross-checked against the integral's own limits, and the fixed-point/Robin roots use `FindRoot` rather than the SymPy bisection — all at prefactor 100. Both engines exit 0, 38/38 Mathematica PASS reconcile exactly with the call sites (no parser-skip), the FAIL path is a non-vacuous `Exit[1]`, outputs are fresh, and the two engines agree to the stated tolerances. No material change; downstream is not stale.
