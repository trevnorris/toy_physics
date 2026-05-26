---
unit_id: 057
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 057

## Per-finding outcomes

### F1 — paper_misalignment (user-resolved direction (a) — script edit)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:83-91` — new block computing `dPe = sp.simplify(sp.diff(zeta_phys, Pe))` and sweeping `Pe in {1/10, 1/2, 1, 2, 5, 10}` at `kappa=1, y=pi/4` with `if val <= 0: raise AssertionError` and a final `print("partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS")`.
- `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:96-107` — analogous `Module[{pevals, signOk}, ...]` with `pevals = {1/10, 1/2, 1, 2, 5, 10}`, `D[zetaPhys, pe] /. {pe -> #, kappa -> 1, y -> Pi/4}`, `AllTrue[..., ... > 0]`, dispatching `pass[...]` or `fail[...]`.

**Assessment:**
Both engines now anchor the Pe-monotonicity sign locally. Each block carries the required carry-forward comment naming Stage 056 (sympy lines 83-85 cite "Stage 056 (notes §4: dOmega_Pe/dPe > 0 ... via Cov_Pe(chi_0, s) > 0)"; wl lines 96-98 mirror the same comment). The sweep covers six Pe values spanning small (`1/10`) to moderately large (`10`), with `kappa=1, y=pi/4` consistent across engines. Sympy uses strict `if val <= 0` (rejects zero) and Mathematica uses strict `> 0`; both are non-tautological because they evaluate the actual `D[zeta_phys, Pe]` symbolic form at concrete points. Transcripts show `partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS` (sympy line 15 of .txt) and `PASS: partial_Pe zeta > 0 on constructive branch (numerical sweep)` (mathematica line 21 of .txt). F1 direction (a) is satisfied.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- Mathematica `kappaMax` derivation (wl:110-111): replaced literal `kappaMax = FullSimplify[Pi^4/(4(4 zetaReq - Pi^2)), ...]` with `kappaMaxSol = Solve[zetaReq == (Pi^2/4)(kappa + Pi^2/4)/kappa, kappa, Reals]; kappaMax = FullSimplify[kappa /. First[kappaMaxSol], ...]`. The existing `kappa_max identity` check at wl:119 now compares Solve-derived expr to the literal closed form — non-tautological.
- Mathematica `kappaReq` derivation (wl:124-125): replaced literal `kappaReq = FullSimplify[(omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2), ...]` with `kappaReqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2), kappa, Reals]; kappaReq = FullSimplify[kappa /. First[kappaReqSol], ...]`. The `kappa_req identity` check at wl:134-137 is now solve-vs-closed-form.
- Mathematica `y_req` derivation (wl:138-140): replaced substitution-cancellation `expectZero["y_req defining equation", zetaReq - FullSimplify[...]]` with `yReqSqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + ysq), ysq, Reals]; yReqSqSolved = FullSimplify[ysq /. First[yReqSqSol], ...]; expectZero["y_req identity", yReqSq - yReqSqSolved]`.
- SymPy `y_req` derivation (py:125-132): replaced `expect_zero("y_req defining equation", zeta_req - sp.simplify(...))` with `y_req_sq_solved = sp.solve(sp.Eq(zeta_req, Omega_Pe**2 * (kappa + sp.pi**2/4)/(kappa + y**2)), y**2)[0]; expect_zero("y_req identity", y_req_sq - y_req_sq_solved)`.

**Assessment:**
All three Mathematica closed-form quantities (kappaMax, kappaReq, yReqSq) are now derived via `Solve` and compared against the literal closed forms. The SymPy `y_req` check is now a real solve-vs-literal comparison (the SymPy `kappa_req` derivation already used `sp.solve` pre-edit at py:112-117 and remains unchanged — that pre-existing check was non-tautological, so no edit was needed there). Transcripts confirm: sympy line 24 prints `y_req identity = 0`; mathematica lines 30-31, 39-42 print `PASS: kappa_max identity`, `PASS: kappa_req identity`, `PASS: y_req identity`. The renames from "defining equation" to "identity" appear in both transcripts. Edit boundaries respect the directive — no collateral changes to surrounding definitions, prints, or unrelated assertions. The leading `kappa_req defining equation` check (Mathematica wl:130-133) was retained from the original script as an unchanged consistency check and still passes (txt line 38).

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy py:75-81: added the `partial_kappa` sign sweep block (`for y_val in (sp.pi/8, sp.pi/6, sp.pi/4, sp.pi/3, sp.Rational(7,16)*sp.pi)` ... `if val >= 0: raise AssertionError` ... `print("partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS")`).
- Mathematica wl:36 — `$Assumptions` line now includes `y < Pi/2`.
- Mathematica wl:85-94: added the `Module[{yvals, signOk}, ...]` block sweeping the same y-values with strict `< 0` and `pass[...]`/`fail[...]` dispatch.

**Assessment:**
Both new sweeps exercise the actual symbolic `dkappa = sp.diff(zeta_phys, kappa)` / `D[zetaPhys, kappa]` expression at concrete y-values; these are non-tautological (a wrong derivative would not reliably yield negative values at all five sample points). The y-values cover small (pi/8), moderate (pi/4), and near-upper-boundary (7 pi/16 ≈ 0.9 * pi/2) cases. Mathematica `$Assumptions` now declares `y < Pi/2` per the directive. Transcripts show `partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS` (sympy line 14) and `PASS: partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)` (mathematica line 20).

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- Mathematica wl:43-56: inserted an independent physical-operator derivation block `Module[{KW, Kphi0, aKPhys}, KW = KX + Pi^2 TX/(4 LL^2); Kphi0 = KX + TX y^2/LL^2; aKPhys = FullSimplify[KW/Kphi0, ...]; aKKappaFromPhys = FullSimplify[aKPhys /. KX -> kappa TX/LL^2, ...]; expectZero["A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)", aKKappaFromPhys - (kappa + Pi^2/4)/(kappa + y^2)]]`. The existing `aKX`/`xSub`/`aKKappa` substitution chain at wl:58-60 remains intact as a consistency check.

**Assessment:**
Block is inserted between the `omegaPe` definition (wl:38-41) and the existing `aKX` chain (wl:58), exactly as the directive specified. Uses the symbols `KX`, `TX`, `LL` that do not appear in the SymPy script — confirms independence. The check compares the kappa-substituted physical-operator ratio against the canonical `(kappa+Pi^2/4)/(kappa+y^2)` form; passes only if the physical-operator algebra independently yields the same closed form (non-tautological). Mathematica transcript line 6 shows `PASS: A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)`. Original substitution chain is preserved and still passes (txt line 11: `PASS: A_K - (kappa+Pi^2/4)/(kappa+y^2)`).

## Exec log assessment

**SymPy:** exit=n/a (per-stage exec log `stage_057_sympy.log` not present; only `stage_057_diff.patch` is in `redteam/exec_logs/`). Fallback transcript at `scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt` (mtime 2026-05-26 03:06) shows the script completed to its final `Stage 40 audit passed.` line with all `expect_zero` residuals = 0 and both new PASS lines for the F1 (Pe sweep) and F3 (kappa sign sweep). Notable lines:
- `A_K - (kappa+pi^2/4)/(kappa+y^2) = 0` (line 7)
- `partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS` (line 14)
- `partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS` (line 15)
- `y_req identity = 0` (line 24)

**Mathematica:** exit=n/a (per-stage exec log `stage_057_mathematica.log` not present). Fallback transcript at `mathematica/output/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.txt` (mtime 2026-05-26 03:06) ends at `Stage 057 Mathematica audit passed.` with all `expectZero` PASS dispatches. Notable lines:
- `PASS: A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)` (line 6) — new F4 line
- `PASS: partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)` (line 20) — new F3 line
- `PASS: partial_Pe zeta > 0 on constructive branch (numerical sweep)` (line 21) — new F1 line
- `PASS: kappa_max identity`, `PASS: kappa_req identity`, `PASS: y_req identity` (lines 31, 40, 42) — F2 lines now backed by `Solve`

Two `Limit::alimv` warnings appear at lines 25 and 27 around the `zeta_max - limit(Pe->inf, y->0)` check — these are pre-existing carry-forward warnings unchanged by the F1–F4 edits, and the corresponding `expectZero` still passes (line 29). Not blocking.

**Output freshness:** both `.txt` outputs are mtime 2026-05-26 03:06; source files are mtime 2026-05-26 03:04. Outputs are 2 minutes newer than sources — confirms regeneration after Codex edits. Per-stage `stage_057_sympy.log` / `stage_057_mathematica.log` are absent from `redteam/exec_logs/`; only `stage_057_diff.patch` is captured there. The diff patch matches the changes I observe in the live source files, and the saved `.txt` transcripts contain all post-fix PASS lines, so the fallback is sufficient.

## Material-change assessment

`material_change`: false.

The four edits add new local verifications (Pe-monotonicity sweep, partial_kappa sign sweep, A_K physical-operator derivation) and tighten three threshold checks from tautologies to genuine Solve-vs-closed-form comparisons. None of these change a derived closed form that downstream stages consume — every closed form (Omega_Pe, A_K, zeta_phys, zeta_max, kappa_max, kappa_req, y_req^2, Omega_req^2) is bit-identical pre- and post-edit (the script transcripts print the same expressions for them). The added Mathematica `$Assumptions` clause `y < Pi/2` strengthens the domain but does not change any algebraic result the engines previously simplified to. Downstream units depending on the closed forms (e.g., any stage that consumes the placement-map outputs) are unaffected.

## Side observations (non-blocking)

- The Mathematica `kappa_req` Solve result (txt line 36) is wrapped in `ConditionalExpression[..., (omegaPe^2 < zetaReq < ...)]`. The `expectZero` helper at wl:20-30 explicitly strips `ConditionalExpression[e_, _] :> e` so the residual still reduces to 0 and the check passes (PASS at line 39). No action required — the wrapper-stripping idiom is the project's standard per `feedback_mathematica_script_idioms`.
- Per-stage exec logs are not landing in `redteam/exec_logs/` for stage 057; only the diff patch is captured there. The orchestrator likely needs to update its log-capture path or the verifier needs to keep relying on the `scripts/output/` and `mathematica/output/` `.txt` transcripts as ground truth. Not blocking for unit 057, but worth flagging if downstream verifiers expect the per-stage logs to exist.
- The banner inside both scripts still reads "STAGE 40" / "STAGE 040" (the legacy numbering); the audit, directive, and verification all use the new `stage_057` identifier. Cosmetic only — not part of any finding.

## Verdict justification

All four findings have `## Applied` blocks recording the requested edits; the live source files at `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py` and `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl` contain the exact code blocks described in the directive (and shown in `stage_057_diff.patch`); the regenerated `.txt` transcripts (mtimes newer than source) show every required new PASS line and no failures. The Pe-monotonicity sweep added under direction (a) includes the carry-forward comment naming Stage 056 in both engines (sympy py:83-85, mathematica wl:96-98). F2's tautology fixes use `Solve`/`sp.solve` non-tautologically. F3's sign sweeps exercise the actual derivative symbolically before sign evaluation. F4's physical-operator derivation introduces `KX`, `TX`, `LL` exclusively on the Mathematica side. No regressions are visible in the diff: all pre-existing PASS lines remain (line-count diff in transcripts is purely additive). Verdict: verified.
