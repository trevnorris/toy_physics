---
unit_id: 068
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T21:45:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: true
---

# Verification — unit 068

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
SymPy (`scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:44-110`):
- Lines 53-60 introduce symbols `A_in, C_sym`, form `A_trans = C_sym * A_in`, square it, substitute `A_in**2 -> W_wall` and `C_sym**2 -> C2`, and assert `W_res_derived - C2*W_wall == 0`.
- Lines 64-68 introduce `Cres = sp.symbols("C_res", positive=True)`, set `Pres_derived = 1/Cres**2`, and assert `Pres_derived - 1/Cres**2 == 0` and `(1/Cres**2)*Cres**2 - 1 == 0`.
- Lines 81-110 replace the postulated `Wfail_match = Pe_req/Delta_inf` etc. with `sp.solve(sp.Eq(W_match_sym*Delta_sym, Pe_req), W_match_sym)` for the matched branch and `sp.solve(sp.Eq(C2*W_prof_sym*Delta_sym, Pe_req), W_prof_sym)` for the profile branch. Adds the new non-trivial cross-relations `Wfail_res * C2 - Wfail_match == 0` and `Wsuff_res * C2 - Wsuff_match == 0`, then keeps the resonance-substitution assertions `Wfail_res(C2->1/Pres) - Pres*Wfail_match == 0` and the Wsuff counterpart.

Mathematica (`mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:28-66`):
- Mirrors the SymPy structure via independent `Solve[...]` calls: `Solve[Wres == C2*Wwall, Wres]` (line 36), `Solve[Pres*Cres^2 == 1, Pres]` (line 41), `Solve[Wmatch*DeltaSym == PeReq, Wmatch]` (line 46), `Solve[C2*Wprof*DeltaSym == PeReq, Wprof]` (line 51). Same assertion set as SymPy.

**Assessment:**
The diff matches the directive's required-change block essentially verbatim in both engines. The new `Wfail_res * C2 - Wfail_match` and `Wsuff_res * C2 - Wsuff_match` assertions are genuinely non-tautological: `Wfail_match` is solved from `W_match*Delta = Pe_req` (no `C2` involvement), while `Wfail_res` is independently solved from `C2*W_prof*Delta = Pe_req`. The two solutions come from two different posed equations, and their `C2`-weighted equality is a real cross-relation between the matched and profile branches (i.e. M3/M4 of the claim manifest).

Caveats noted but not blocking:
- The `W_res - C2*W_wall` derivation in SymPy is partly cosmetic — `(C_sym * A_in)**2` simplifies algebraically to `C_sym**2 * A_in**2`, then two substitutions rename to the audit symbols. SymPy *does* perform the symbolic squaring step, so it is not strictly tautological, but it is the weakest of the four assertions.
- The `P_res - 1/C_res^2` assertion compares `Pres_derived` (defined one line earlier as `1/Cres**2`) to `1/Cres**2`. By itself this is the identity `x - x = 0`. The directive explicitly specifies this exact code, so the verifier marks F1 as resolved per directive compliance; however, the practical non-triviality of *this particular* assertion is low. The cross-relation assertions (M3/M4) carry the load.
- The Mathematica `Solve[Wres == C2*Wwall, Wres]` likewise returns `Wres -> C2*Wwall` and the subsequent assertion is then trivially zero; same caveat as the SymPy mirror.

Overall: F1 is resolved because (a) Codex applied the directive's verbatim replacement, (b) the new M3/M4 cross-relations are non-tautological, (c) the structural separation between matched-branch and profile-family solves is now in place, and (d) the docstring claims are now exercised by independent `solve`/`Solve` calls rather than postulated symbols.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
SymPy (`scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:112-134`):
- Way A: `success_band_widthA = sp.simplify((Wsuff_res - Wsuff_match).subs(C2, 1/Pres))` and the failure analog.
- Way B: solves `sp.Eq(Wsuff_match + gap_sym, Pres * Wsuff_match)` for `gap_sym` (and the failure analog).
- Asserts `success_band_widthA - success_band_widthB == 0` and `failure_band_widthA - failure_band_widthB == 0`.

Mathematica (`mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:68-84`):
- Way A: `FullSimplify[(WsuffRes - WsuffMatch) /. C2 -> 1/Pres, ...]`.
- Way B: `gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap]`.
- Same two A-vs-B assertions.

**Assessment:**
This is genuinely non-trivial. Way A pulls from the F1-derived `Wsuff_res = Pe_req/(C2*Delta_0)` (which came from the profile-family Solve) and `Wsuff_match = Pe_req/Delta_0` (from the matched-branch Solve), takes their raw difference, then substitutes `C2 -> 1/Pres` to get `Pe_req*(Pres-1)/Delta_0`. Way B solves a different equation (`x + gap = Pres*x`) for `gap`, giving `(Pres-1)*Wsuff_match = (Pres-1)*Pe_req/Delta_0`. The agreement is not built into either definition — it tests the relation between the F1-derived profile threshold (at the resonance point) and the gap-equation algebra. If F1's profile-threshold solve had returned the wrong functional form (e.g. `Pe_req/(C2^2*Delta_0)`), Way A would no longer equal Way B and the assertion would fail.

The literal anti-pattern from the original report — "`Pres*X - X` followed by `X/Pres - (Pres-1)`" — is no longer present in either script.

Stylistic note (non-blocking): the Mathematica file contains a dangling line `gapSym;` (line 74) that does nothing — the substantive `gap /. First@Solve[..., gap]` uses the symbol `gap`, not `gapSym`. The directive itself contained this same stray line, so Codex faithfully copied it. The assertion semantics are unaffected.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The Mathematica file now contains six distinct `Solve` calls (lines 36, 41, 46, 51, 75, 76) — none of which appear in the SymPy script as `Solve[...]`-shaped operations. The SymPy script uses `sp.solve(sp.Eq(...))` and `.subs(...)` rather than the same algebraic operations. The `expect_zero`/`expectZero` helper signatures still parallel each other, but the underlying algebraic path differs.

The verbatim "FINAL LEDGER" prose block (originally a word-for-word copy of the SymPy script's interpretation text) has been replaced with three lines:
```
W_res derived as |C|^2 W_wall.
P_res derived as 1/C_res^2.
Band widths cross-checked two ways.
```
This satisfies the directive's structural requirement that the .wl no longer mirror the .py's interpretation text.

**Assessment:**
F3 is resolved as a side-effect of F1+F2 plus the explicit prose rewrite. Concretely: (a) ≥4 distinct `Solve[...]` calls present (✓), (b) FINAL LEDGER prose no longer matches the SymPy script verbatim (✓), (c) Mathematica still calls `Exit[0]` at the bottom (✓ line 91).

One residual transliteration concern: the assertion labels in `expectZero` strings (e.g. `"Wfail_res * C2 - Wfail_match"`, `"P_res - 1/C_res^2"`) still match the SymPy `expect_zero` labels verbatim. The directive did not require renaming these, and they are diagnostic strings rather than derivation steps, so this is a cosmetic match rather than a transliteration of derivation logic. Not blocking.

## Exec log assessment

**SymPy:** exit=n/a. Both exec_logs/stage_068_sympy.log and stage_068_mathematica.log are absent from `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/` — only `stage_068_diff.patch` is present. The orchestrator appears not to have re-run the scripts after Codex's edits, or the log capture was skipped for this batch.

However, the output `.txt` files were regenerated:
- `scripts/output/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.txt` (mtime 20:00, vs. source mtime 19:57)
- `mathematica/output/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.txt` (mtime 20:00, vs. source mtime 19:58)

Both output files show all assertions printing residual `= 0`, and the Mathematica output shows `PASS:` for every `expectZero` call. Specifically the SymPy output (lines 5-21) contains:
```
W_res - C2 * W_wall = 0
P_res - 1/C_res^2 = 0
P_res*C_res^2 - 1 = 0
Wfail_res * C2 - Wfail_match = 0
Wsuff_res * C2 - Wsuff_match = 0
Wfail_res(C2->1/Pres) - Pres*Wfail_match = 0
Wsuff_res(C2->1/Pres) - Pres*Wsuff_match = 0
success band A vs B = 0
failure band A vs B = 0
```
And the Mathematica output (lines 5-32) shows the same nine residuals, each followed by a `PASS:` line.

**Mathematica:** exit=n/a (log missing). Output file shows all `PASS:` lines and reaches the FINAL LEDGER banner, implying `Exit[0]` was hit successfully.

**Output freshness:** Confirmed. Both output `.txt` files are newer than their source scripts (20:00 vs. 19:57/19:58).

The missing exec logs are a process gap; the regenerated output files provide circumstantial evidence that the scripts ran to completion without raising, but the verifier flags this as an audit-trail caveat. Given the output content (every assertion = 0, every PASS recorded, FINAL LEDGER reached), the verifier accepts the substantive outcome.

## Material-change assessment

`material_change`: true.

Stage 068's substantive content changed: the audit script now derives `W_res = C^2 W_wall`, `P_res = 1/C_res^2`, and the matched/profile threshold pair from `solve`-based premises rather than postulating them. Downstream units that consume the `P_res = 1/C_res^2` relation or the resonance-amplification scaling could in principle be re-audited against a tightened upstream. However, since the *symbolic content* of the derived expressions is identical to the previously postulated forms (the same `W_res = C^2*W_wall` and `Wfail_res = Pe_req/(C2*Delta_inf)` come out of the new solves), no numerical or functional result actually changed — only the derivation route. The orchestrator's policy of marking units > 068 as `upstream_stale: true` is appropriate, but a targeted re-audit on downstream units is unlikely to surface new defects from this edit alone.

## Side observations (non-blocking)

- The SymPy banner still reads `"STAGE 51 — RESONANCE-CORRECTED THRESHOLDS"` (line 37) and the Mathematica banner reads `"STAGE 051 — RESONANCE-CORRECTED THRESHOLDS"` (line 26). Both should plausibly say `STAGE 068`. The original auditor report did not flag this and the directive did not ask for renaming, so it is out of scope here. Worth a future cleanup pass.
- The Mathematica file has a stray no-op statement `gapSym;` on line 74 that the directive itself contained. Cosmetic only.
- The `expectZero` label strings in the Mathematica script still match the SymPy label strings verbatim. This is non-blocking under F3 since labels are diagnostic, not derivation logic.
- The SymPy `W_res` derivation (lines 53-60) and the Mathematica `WresRule` derivation (line 36) are both close to identity-style operations under the hood. The non-triviality of stage 068 now rests primarily on M3/M4 (`Wfail_res * C2 - Wfail_match == 0` and friends) and on the band-widths A-vs-B cross-check, which are the genuine load-bearing assertions.

## Verdict justification

All three findings are `resolved`: Codex applied each directive's required-change block essentially verbatim in both engines; the new cross-relations (`Wfail_res * C2 - Wfail_match`, `Wsuff_res * C2 - Wsuff_match`, band A vs B) are non-tautological because they connect independently-solved expressions from differently-posed equations; the Mathematica script now uses six distinct `Solve` calls and a Mathematica-idiom FINAL LEDGER, breaking the transliteration. Exec logs were not captured for this batch, but the regenerated output files show every assertion evaluating to zero with explicit `PASS:` markers and both scripts reaching their final banners. The previously-flagged tautology pattern (`Pres*X - X` followed by `X/Pres - (Pres-1)`) is gone. Verified.
