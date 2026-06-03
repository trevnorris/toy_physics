---
unit_id: 238
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T08:45:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 238

## Per-finding outcomes

### F1 — tautological_check (support-blindness, reframed iter 2)

**Classification:** resolved

**What changed:** `scripts/...stage238..._sympy_audit.py:155-202`. Codex added the live `M_tr` support channel and three discriminating guards before the retained `∂=0` observable checks:
- (a) Negative control: `M_tr = M_mix*(1 + zeta*(1-eps)/(1-zeta*eps))` (line 155); `dMtr_dzeta` asserted `== M_mix*(1-eps)/(1-zeta*eps)**2` (158-161) and `dMtr_dMmix` asserted `== 1 + zeta*(1-eps)/(1-zeta*eps)` (162-165); plus an explicit `if dMtr_dzeta == 0 or dMtr_dMmix == 0: raise` non-vacuity guard (166-167).
- (b) Leak detector: `Rtr_leak = Rtr*M_tr/M_mix` with `if simplify(diff(Rtr_leak, zeta)) == 0: raise` (170-172).
- (c) Exclusion: loop over `{R_tr, T^2, eps_eta, q_tr, q_nt_factored, q_eta}` raising if `zeta`/`M_mix` in `free_symbols` (176-187), then the retained six `∂_zeta=0` / six `∂_{M_mix}=0` asserts (189-202).

**Assessment:** Correct and non-vacuous. All three guards present and discriminating. The negative control is genuinely nonzero (closed forms match notes §4 and are enforced `!=0`). The leak detector proves a support leak into an observable WOULD trip the `∂_zeta` check (its witness in the Mathematica log is a nonzero rational). The exclusion + retained `∂=0` checks now carry weight because the symbols are proven live and in-scope. Codex correctly did NOT route `M_tr` into the observables (that was the blocked fabrication path). The `q_nt` used is `q_nt_factored = log(T2/T2_ref) - (B/C)q_tr` — still support-free, consistent with the original script's choice. Faithful to the reframe.

### F2 — insufficient_verification (rigid first-order q_nt)

**Classification:** resolved

**What changed:** `...py:142-150`. `q_eta_first` is now recomputed from a perturbed `eps_eta*exp(h*dlnepseta)` via `diff(log(...), h)|_{h=0}` (142-143) instead of the trivial `= dlnepseta`. The duplicate `q_nt_factorized_rigid = dlnT2_calc` / `dlnT2_expected` check is removed; line 150 now asserts `q_nt_first_rigid - dlnT2_calc == 0`, where `q_nt_first_rigid = q_nt_first_direct.subs(dlnRtr_calc, 0)` (147) is anchored to the `_calc` perturbation chain (`q_nt_first_direct` is built from `dlnRtr_calc` and `dlnRtarget_calc`).

**Assessment:** Substantive, not a relabeled tautology. The assertion now exercises `q_nt` against the `diff`-derived `dlnT2_calc` (not `dlnT2_expected`), so a perturbation error would surface.

### F3 — tautological_check (tracking gate)

**Classification:** resolved

**What changed:** `...py:209-219`. The factorization assert now references `dlnRtr_calc` (the `diff`-derived line-102 quantity), not `dlnRtr_expected` (211). A second gate-locus assert was added: `dlnRtr_calc.subs(dlndelta, -(1+deltaU)/(1+chi0)*dlnchi)` simplifies to 0 (216-219).

**Assessment:** Substantive. Both asserts depend on the perturbation derivation; the gate-locus substitution exercises the gate's physical content (drift vanishes when `tracking_condition = 0`).

### F4 — missing_verification_script (new independent Mathematica)

**Classification:** resolved

**What changed:** New file `mathematica/...stage238..._mathematica_audit.wl` covering M1–M7, with `expectZero`/`expectNonZero`/`expectSupportFree` helpers (the first strips `ConditionalExpression[0,...]`), exit 0.

**Assessment:** Genuinely independent route, not a transliteration:
- M4 first-order compilers: `.wl` pre-expands the observable logs ADDITIVELY by hand (`logRtr = Log[1+chi0/(1+deltaU)] - Log[1+chi0]`, etc., lines 132-137) then perturbs and `D[...,h]/.h->0` — the directive's requested contrast against the `.py`'s `diff(log(simplify(rational_obs)), h)` path.
- M1 identity: `FullSimplify[Together[...]]` directly on `Rtarget T2 - Lambda0(1-epsEta)` (vs `.py`'s `simplify(factor(...))`).
- M3 factorization/rigid: log-difference / `PowerExpand` additive forms (vs `.py`'s `exp(...) - ratio` comparison).
- M2 mirrors the reframed non-vacuous form: negative controls (`expectNonZero[D[Mtr,zeta]]`, `D[Mtr,Mmix]`), leak detector (`expectNonZero[D[RtrLeak,zeta]]`, witness nonzero in log), `FreeQ` exclusions, retained `D[obs,*]==0` claims.
- M6 anchored to derived `dlnRtrCalc` + gate-locus substitution.
- Var names (`Bstar/Cstar`, `epsEta`, `dlnRtrCalc`), assertion choreography, and decomposition all differ from the `.py`. M5 rigid `qNtFirstRigid` is phrased as `-dlnRtargetCalc - epsEta/(1-epsEta)dlnepsEta` (the `dlnRtrCalc->0` slice of the direct form), still anchored to the derived `dlnRtargetCalc`/`dlnT2Calc` chain — legitimate, non-vacuous.

## Exec log assessment

**SymPy:** exit=0. 31 `[ok]` lines incl. the new "live support channel depends on zeta and M_mix", "support leak detector for R_tr", "structural support exclusion...", "tracking gate factorization", "tracking gate drift vanishes on gate locus". Matches the saved `.txt` line-for-line.

**Mathematica:** exit=0. All M1–M7 covered with `[ok]`; negative-control/leak witnesses print genuinely nonzero closed forms (e.g. zeta witness `-(((-1+eps)*Mmix)/(-1+eps*zeta)^2)`, leak witness `-(((1+chi0+deltaU)*(-1+eps))/((1+chi0)*(1+deltaU)*(-1+eps*zeta)^2))`), structural-exclusion lines print `True`, all residuals `= 0`.

**Output freshness:** confirmed. `.py` mtime 08:22:38, `.wl` mtime 08:24:46, both `output/*.txt` mtime 08:31:15 (newer than the scripts). SymPy `.txt` matches the exec log exactly.

## Material-change assessment

`material_change`: false. No derived constant, sign, or closed form changed — the underlying math (M1–M7 identities, the three first-order compilers, the transfer-shape/dressing gates) was already correct per the auditor; the fixes only made the support-blindness (F1), rigid first-order `q_nt` (F2), and tracking-gate (F3) checks actually test it, and added an independent Mathematica engine (F4). Downstream units do not consume any newly-altered value. No specific downstream re-audit concern.

## Side observations (non-blocking)

- The verify prompt's primary input paths carried a spurious `toy_projects/` segment; the real tree is `/var/projects/toy_physics/research/pde_ledger/...`. No impact on the verdict.
- The F1 exclusion/`∂=0` checks use `q_nt_factored` (not the raw `q_nt` containing `Rtarget`); both are support-free, and this matches the original script's design choice — noted only for completeness.

## Verdict justification

All four findings are `resolved`. F1 (the load-bearing reframe) carries the three required discriminating guards — live negative control (enforced nonzero), a leak detector proving the `∂` check would catch contamination, and structural exclusion — so the retained `∂=0` asserts are no longer vacuous, and `M_tr` was correctly NOT fabricated into the observables. F4 is a genuinely independent Mathematica route (additive log pre-expansion, `FullSimplify`/`PowerExpand` log-difference forms, distinct var names/choreography) whose M2/M6 mirror the same non-vacuous form. F2 and F3 are substantive re-anchors to the `diff`-derived `_calc` chain, not relabeled tautologies. Both engines exit 0, FAIL paths are live (`AssertionError` / `Exit[1]`), pass counts match the manifest, outputs are fresh, and no `notes/`/`paper/` files were touched.
