---
unit_id: 026
batch: II.1
verifier_model: claude-opus-4-7
verify_date: 2026-05-22T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 026

## Per-finding outcomes

### F1 — tautological_check (overlap collapse)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:129-146` — the three duplicate `integrate(u0*f0, ...)` instantiations were collapsed into a single `overlap_u0_f0` definition and a single `overlap_u0_u0` definition; the four `I_(...)` labels are now aliases to those two underlying overlaps; the four bogus assertions were replaced with two non-tautological ones (`overlap_u0_f0 - kappa` and `overlap_u0_u0 - 1`).
- `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:97-114` — the same restructuring applied; assertions are `overlap_u0_f0 - kappa` and `overlap_u0_u0 - 1`.

**Assessment:**
The edits match the directive verbatim. The output (`scripts/output/.../sympy_audit.txt` lines 40-47 and `mathematica/output/.../mathematica_audit.txt` lines 44-53) shows exactly two assertions instead of four, with the four `I_(...)` lines still printed for traceability. The two remaining assertions are non-tautological: `overlap_u0_f0` is computed by `integrate(u0*f0, s, 0, L)` whereas `kappa` is computed via a separate path `integrate(u0*f_n, s, 0, L)` then substituted `n -> 0`, so the equality is a real consistency check between the two integration paths. `overlap_u0_u0 - 1` independently checks the constant-mode normalisation (this is the integral of 1/L over [0,L], which equals 1; the assertion would fail if SymPy/Mathematica miscomputed the bound).

### F2 — tautological_check (`_expected` rebuilds)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py` — the 16-line block defining `Delta_expected, Q_expected, P_expected, B0_expected, P0_expected` and the five corresponding `expect_zero(...)` calls were deleted (visible in `stage_026_diff.patch` lines 187-201 removed). Section III.2 now ends with the print statements; the function returns directly afterwards.
- `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:93` — the corresponding `deltaExpected, qExpected, pExpected, b0Expected, p0Expected` symbols were removed from the `Module[{...}]` declaration. The block defining and asserting them (formerly lines 121-131) is gone.

**Assessment:**
Matches the directive exactly. The saved outputs no longer contain any `Delta - Delta_expected = 0`, `Q - Q_expected = 0`, etc. lines (search of both `.txt` files confirms zero hits for `_expected`). Closed-form symbolic results for Delta/Q/P/B0/Z0/N0/D0/P0 are still printed (sympy output lines 52-63; mathematica output lines 58-69), so the diagnostic content is preserved. The tautological assertions are gone.

### F3 — tautological_check (K_req back-substitution)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:195-202` — `K_req_expected = sp.simplify(B0 + Q/Delta + mhat**2 * P**2 / (target * Delta**2))` was removed; the assertion is now `expect_zero("residual @ K_req", residual.subs(K, K_req))`.
- `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:148-162` — `kReqExpected` removed from the `Module[{...}]` symbol list (line 148) and from its definition/usage; the assertion is now `expectZero["residual @ K_req", residual /. k -> kReq]`.

**Assessment:**
Matches the directive. The sympy output line 122 shows `residual @ K_req = 0`; the mathematica output lines 80-81 show `residual @ K_req = 0` and `PASS: residual @ K_req`. The check is non-tautological in the directive's sense: it asks "does substituting the solver's output back into the original residual yield zero?" rather than "does the solver's output match a hand-rearrangement of the residual?". As the auditor noted in their self-test, the back-substitution still passes by construction (any correct solver output will), but it is no longer a literal algebraic restatement; it would fail if `sp.solve`/`Solve` returned a wrong branch or a typo'd expression.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:53-91` (`overlapLaw[]`) — the single `Integrate[u0*fN, {s, 0, l}]` call was replaced by (i) an indefinite-integral path `indef = Integrate[u0*fN, s]` followed by boundary evaluation `(indef /. s -> l) - (indef /. s -> 0)`, and (ii) the analytic short form `kappaN = Sqrt[2]/((n + 1/2)*Pi)` typed in directly with comment justification. The assertion now compares the two paths: `kappa_n (analytic) - (fundamental thm)`.

**Assessment:**
Matches the directive's required change verbatim. The output (mathematica `.txt` lines 28-31) shows both paths yielding `(2*Sqrt[2])/(Pi + 2*n*Pi)` and their difference being 0. The `.wl` file no longer contains `Integrate[u0*fN, {s, 0, l}]` with explicit limits (only the indefinite form). Hand-checking the analytic short form: `u0*fN = Sqrt[2]/l * Sin[(n+1/2) Pi s/l]`; antiderivative is `-Sqrt[2]/((n+1/2) Pi) * Cos[(n+1/2) Pi s/l]`; at s=l the cosine vanishes (cos of half-integer multiple of pi), at s=0 it is 1; so the boundary contribution is `Sqrt[2]/((n+1/2) Pi)`, matching the typed-in `kappaN`. The independence is modest (both paths still use Mathematica's `Integrate`, just with different limits and a typed-in literal as the second leg), but it does what the directive asked: it exercises a different code path from the SymPy `integrate(u0*f_n, (s, 0, L))` call.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator did not deposit `redteam/exec_logs/stage_026_sympy.log` — only `stage_026_diff.patch` is present. Per the verifier instructions ("the saved outputs are already fresh — read them"), I read the saved `.txt` instead. Notable lines from `scripts/output/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.txt`:
- L20: `int u0^2 - 1 = 0`
- L28: `kappa_n - expected = 0`
- L46-47: `overlap_u0_f0 - kappa = 0` and `overlap_u0_u0 - 1 = 0`
- L122: `residual @ K_req = 0`

No `Delta - Delta_expected` or `K_req - expected` strings appear anywhere in the output.

**Mathematica:** exit=n/a. Same situation — no `stage_026_mathematica.log` was captured; I read the saved `.txt`. Notable lines from `mathematica/output/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.txt`:
- L21, L23: `PASS: int u0^2 - 1`, `PASS: int f0^2 - 1`
- L31: `PASS: kappa_n (analytic) - (fundamental thm)`
- L51, L53: `PASS: overlap_u0_f0 - kappa`, `PASS: overlap_u0_u0 - 1`
- L81: `PASS: residual @ K_req`
- L84: `Stage 9 Mathematica audit passed.`

No `PASS: Delta - Delta_expected` etc. appears, and no `PASS: K_req - expected` appears. The script's `Exit[0]` would not have been reached if any `expectZero` had failed (failure path is `fail[name, res]` → `Exit[1]`), so the "audit passed" banner is a positive signal even without a captured log.

**Output freshness:** Confirmed.
- `scripts/.../sympy_audit.py` mtime: May 21 16:59
- `scripts/output/.../sympy_audit.txt` mtime: May 21 17:01 (newer — fresh)
- `mathematica/.../mathematica_audit.wl` mtime: May 21 17:00
- `mathematica/output/.../mathematica_audit.txt` mtime: May 21 17:02 (newer — fresh)

## Material-change assessment

`material_change`: false.

The fixes are restricted to:
1. Removing tautological assertions (F2 and the bogus three-quarters of F1).
2. Replacing one tautological assertion with a back-substitution check (F3) with the same `K_req` symbolic value as before.
3. Restructuring the Mathematica overlap derivation into two equivalent paths (F4) producing the same `kappa_n` and `kappa` values.

No printed symbolic result was changed: `kappa`, `kappa_n`, `Delta`, `Q`, `P`, `B0`, `Z0`, `N0`, `D0`, `P0`, `K_req`, and `K_geom` are identical between pre- and post-fix outputs (verified by spot-checking the closed forms in the sympy `.txt`: `Delta = Omega_U**2*Omega_W**2 - 8*lambda_R**2/pi**2`, `B0 = 8*lambda_B**2/(pi**2*varpi**2)`, `K_req`'s top-level structure, etc. — all match the expressions the auditor's engine-cross-check table already noted). Downstream units that consume any of these symbolic outputs are unaffected.

## Side observations (non-blocking)

1. In `wl:97-98`, `overlapU0F0` and `overlapU0U0` are assigned without being declared in the enclosing `Module[{...}]` symbol list (line 93). They leak into the `Global`` context. The script still runs and the assertions still pass (the saved output confirms this), but a cleaner edit would have added them to the Module declaration. Not a finding — purely cosmetic.

2. The F4 fix achieves the directive's stated goal but with a relatively shallow form of independence: both paths still use Mathematica's `Integrate`, and the "analytic short form" is typed in literally rather than derived. A stricter Fourier-projection independence (as the original auditor suggested as the gold-standard option) was not attempted. The directive's `Required change` section explicitly accepted the indefinite-integral + boundary-evaluation approach, so this is consistent with the directive, but a future reviewer may want to revisit if cross-engine independence is being stressed harder.

3. The orchestrator did not capture exec logs for this stage (only the diff patch). The verifier prompt explicitly tells the verifier to fall back to reading the saved outputs in that case, which I did. Flagging for the orchestrator's awareness in case log capture should be hardened.

## Verdict justification

All four findings were applied non-tautologically and the saved outputs corroborate the edits: the four tautological lines (three byte-identical overlap assertions, five `_expected` rebuilds, one `K_req - expected` rearrangement) are gone, replaced by two genuine overlap checks, an empty section (closed forms printed only), a back-substitution check, and an indefinite-vs-typed-form comparison for `kappa_n`. The `Stage 9 Mathematica audit passed.` banner and final `residual @ K_req = 0` confirm both engines completed without `expectZero` failures. No regressions visible in the diff (only deletions and label-renames in the asserted region; no closed-form expressions were touched). Verdict: `verified`.
