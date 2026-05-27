---
unit_id: 070
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T13:25:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 070

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex inserted a new `Module[...]` block in `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:75-120`, between the existing `expectZero["Xi - W_wall", ...]` at line 73 and the closing `banner["STAGE 053 THEOREM LEDGER"]` at line 122. The block matches the directive's prescribed text verbatim. It introduces:
- `fProf[xi_] := Sech[xi]` with `fp` and `fpp` derivatives (lines 80-82)
- Two `NIntegrate[..., {xi, -Infinity, Infinity}, WorkingPrecision -> 30]` calls computing `IfNum` and `IgNum` (lines 83-84)
- Numeric substitutions `ruleNum = {a -> 1, L -> 1, ell -> 1/10, ...}` (line 96)
- Three `Abs[... - ...] < tol` assertions wrapped in `If[...pass, ...fail]` for `kappa`, `W_wall`, and `Xi = W_wall` (lines 105-119)

**Assessment:**
The fix is correct and substantively resolves the finding. The new block uses `NIntegrate` — a genuinely independent numerical code path (quadrature) that does not appear in the SymPy script. The numerical assertion `Abs[kappaNum - kappaCmp] < 10^-10` is non-tautological because:
- `kappaNum = KxNum*L^2/TxNum` is computed by assembling the thin-shell forms with the quadrature moments and dividing out
- `kappaCmp = 4*(m*cSw*L/hbar)^2 + (IgNum/IfNum)*(L/ell)^2` is the closed-form expression with the same quadrature moments substituted

The exec log shows `kappa_num = 144.00000000000000000000000000000039923598` and `kappa_closed = 144.00000000000000000000000000000039923599` — they agree to ~10^-29, which would not occur if either the thin-shell substitution or the closed-form `kappa` decomposition were wrong. Similarly `W_wall_num = W_wall_closed = 400.` and `Xi_num = 400.`. All three new `PASS` lines appear in the exec log (lines 37, 40, 42). The Mathematica script exits 0.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:51-64`): inserted directly after `J1 = sp.simplify(If / Hw)` (line 50). The new block defines a fresh local `xi`, builds `chi_phi = d/dxi sech(xi)`, computes `If_sym = sp.integrate(chi_phi**2, (xi, -sp.oo, sp.oo))`, then asserts the structural identity `I1_constH/J1_stage48 - 4*pi*a**2*ell == 0` via `expect_zero`. A comment block anchors the Stage-48 vs Stage-47 normalization conventions explicitly.
- Mathematica (`mathematica/...wl:87-94`): three `Print` statements inside the F1 `Module` document the Stage-48 and Stage-47 normalization conventions, followed by an `expectZero` for `(4*Pi*a^2*ell*IfMoment/Hw)/(IfMoment/Hw) - 4*Pi*a^2*ell`.

**Assessment:**
The new SymPy assertion is substantively distinct from the original `Xi - W_wall == 0` because it directly targets the structural ratio `I_1/J_1 = 4 pi a^2 ell` that was the implicit dependency. The check now uses a concrete profile-derived moment (`If_sym = 2/3` from a real integral of `(sech*tanh)^2`), confirming that the constant-`H_w` approximation in `I_1` and the Stage-48 normalization in `J_1` differ only by the shell-measure factor. The exec log shows:
- `I1 / J1 - 4 pi a^2 ell (sech profile) = 0` (line 16)
- `sech-profile moment I_f = 2/3` (line 17, matching analytic expectation)

This is non-tautological: had `J_1` been defined with the shell measure included (e.g. `J1 = Nphiphi/Hw`), the assertion would fail because the ratio would not equal `4 pi a^2 ell`. The original `Xi - W_wall == 0` could not distinguish that error because it cancels symmetrically.

The Mathematica mirror's `expectZero` uses symbolic `IfMoment` (the abstract symbol), so it is a pure algebraic identity check — strictly it confirms the cancellation algebra. The substantive numerical content comes via the F1 numeric block, which uses the same Sech profile and a non-trivial `ell -> 1/10`, exercising the structural ratio implicitly through the numeric `kappa` and `W_wall` agreement. The `Print` statements explicitly state the convention. Acceptable — the mirror is documented.

The directive's `applied: true` is set in the frontmatter (with an `orchestrator_notes` block explaining the manual update due to the directive's unusual "do NOT edit directive prose" instruction). I have confirmed the edits are present in both files (via the diff patch and direct file reads) and the exec logs show passing results, matching the orchestrator's manual confirmation.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `redteam/exec_logs/stage_070_sympy.log`:
- `kappa - expected = 0`
- `I1 / J1 - 4 pi a^2 ell (sech profile) = 0`
- `sech-profile moment I_f = 2/3  (expected 2/3)`
- `W_wall - expected = 0`
- `Xi - W_wall = 0`

All four original `expect_zero` calls plus the new F2 anchor pass cleanly. No warnings.

**Mathematica:** exit=0. Notable lines from `redteam/exec_logs/stage_070_mathematica.log`:
- `PASS: kappa - expected`, `PASS: W_wall - expected`, `PASS: Xi - W_wall` (original symbolic checks intact)
- `I_f (sech profile) = 0.6666666666...` (matches analytic 2/3 to ~10^-29)
- `PASS: I_1 / J_1 - 4 pi a^2 ell (independent of profile, symbolic check)`
- `kappa_num = 144.000...` and `kappa_closed = 144.000...` agree to ~10^-29
- `W_wall_num = W_wall_closed = 400.` (exact)
- `Xi_num = 400.`, matching `W_wall_num`
- `PASS: kappa numeric profile check`, `PASS: W_wall numeric profile check`, `PASS: Xi = W_wall numeric profile check`

**Output freshness:** The exec logs are dated `2026-05-26T13:13:53-06:00` (sympy) and `2026-05-26T13:13:58-06:00` (mathematica). Both scripts were last modified at `May 26 13:03`, ~10 minutes before the exec runs — logs are fresh relative to scripts. Note: the per-script `output/*.txt` files under `scripts/output/` and `mathematica/output/` are dated `May 22 20:09` (pre-fix) and were NOT regenerated by the orchestrator's exec calls; those exec calls produced the logs at `redteam/exec_logs/` only. This is consistent with the orchestrator's policy of routing audit-cycle runs through `exec_logs/` rather than overwriting committed output snapshots. The logs themselves are the authoritative post-fix transcripts and they show passing assertions.

## Material-change assessment

`material_change`: false.

No load-bearing derived result was altered. The original three symbolic identities (`kappa - expected`, `W_wall - expected`, `Xi - W_wall`) all still pass with identical values. The carry-forward quantities `T_X`, `K_X`, `kappa`, `W_wall = Xi` are unchanged. The new checks are additional corroboration only; they do not modify any downstream-visible expression. Downstream units (071+) do not need re-auditing on account of stage 070.

## Side observations (non-blocking)

1. **Misstated analytic value in the I_g comment.** The Mathematica `Print` at line 86 of the post-fix `.wl` writes `"   (analytic 8/15 = ", N[8/15, 30], ")"`, but the actual `I_g = NIntegrate[(d^2/dxi^2 sech(xi))^2, ...]` is `14/15 ≈ 0.9333`, not `8/15 ≈ 0.5333`. The exec log makes this visible:
   - `I_g (sech profile) = 0.9333333333... (analytic 8/15 = 0.5333333333...)`
   
   Independent verification via direct SymPy integration: `sp.integrate(sp.diff(1/sp.cosh(xi), xi, 2)**2, (xi, -oo, oo)) = 14/15`. The auditor's original report at lines 122 and 229 (the basis for the directive's prescribed text) stated `8/15` as the analytic value, which was wrong. Codex faithfully implemented the directive as written.
   
   This does NOT invalidate the F1 numeric assertions because both `kappa_num` and `kappa_closed` use the same `IgNum` from `NIntegrate`, so they agree (to ~10^-29) regardless of what the closed-form value of `I_g` "should be" — the check tests algebraic consistency of the assembled vs closed-form `kappa` decomposition, not the value of `I_g` against its analytic closed form. Worth fixing the printed reference value in a future cleanup pass; not blocking verification.

2. **Frontmatter `applied: true` was set manually.** The directive's frontmatter shows `applied: true` plus `orchestrator_notes` explaining the manual update because the directive contained both `Do NOT run python or mathematica. Only edit files.` and `Do NOT introduce new features...` — Codex applied both findings correctly to the script files (visible in the diff patch) but did not append the per-finding `## Applied: F1` / `## Applied: F2` blocks. This is process-only; the actual code edits are correct and present. Per the user's note, the orchestrator confirmed via `codex_logs/070_iter1.txt`.

3. **Stage label inconsistency in transcripts.** The exec logs print "STAGE 53" and "STAGE 053" banners in places where the unit is now stage 070. This is a pre-existing labeling artifact (the file headers and banners use the old stage-53 nomenclature even though the file is named stage070). Pre-existing, not introduced by this edit. Out of scope.

## Verdict justification

Both findings are resolved. F1's required Mathematica numeric profile block was inserted verbatim per the directive, uses `NIntegrate` (an independent quadrature code path from SymPy's symbolic integration), and the new numeric assertions agree to ~10^-29 — a non-tautological consistency check between the assembled thin-shell forms and the closed-form expressions evaluated on a concrete `Sech` profile. F2's normalization-convention anchor exercises the previously implicit `I_1/J_1 = 4 pi a^2 ell` structural identity via a profile-derived moment in SymPy and via documented symbolic substitution in Mathematica; the new assertion would catch a normalization error that the original `Xi - W_wall == 0` could not. Both engines exit 0. The original symbolic checks remain intact. No regressions. The only side observation (the auditor's typo-grade `8/15` vs actual `14/15` for the sech I_g) does not affect the substance of the assertions.
