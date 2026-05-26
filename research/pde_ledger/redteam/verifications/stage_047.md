---
unit_id: 047
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T02:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 047 (batch III.1 v2)

This verification supersedes the prior 2026-05-22 verification of unit 047. The 2026-05-26 re-audit caught that the prior F1/F2 fixes (split num/den + channel-saturation guard, and `mSupp = mMix*supportLoadFactor` + `sEnhance = mTr/mMix`) had themselves introduced tautologies. The directive verified here removes those tautologies.

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:41-53` — Removed `rho0_num/rho0_den/rho0`, `sigma0_num/sigma0_den/sigma0`, the two `assert sp.simplify(...) == 0` no-op channel-saturation guards, and the two `expect_zero("rho_0 - chi_0", ...)` / `expect_zero("sigma_0 - chi_0", ...)` calls. Replaced with a comment block explaining the Stage-28 upstream provenance of the saturation rule, retained `chi0 = sp.simplify(gamma * c_etaU / KU)`, and renamed the banner to "1. Coherent interference ratio" (singular). Only `print("chi_0 =", chi0)` remains in the trace.
- `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:37-47` — Mirror edit: removed `rho0Num/rho0Den/rho0`, `sigma0Num/sigma0Den/sigma0`, the two `If[TrueQ[FullSimplify[...]]]` no-op guards, and the two `expectZero["rho_0 - chi_0", ...]` / `expectZero["sigma_0 - chi_0", ...]` calls. Replaced with the parallel comment block, retained `chi0 = FullSimplify[gamma*cEtaU/kU, ...]` and `Print["chi_0 = ", fmt[chi0]]`.

**Assessment:**
The change matches directive option (b) verbatim. The diff (`stage_047_diff.patch:5-49, 92-139`) confirms only the lines named in the directive were touched. The captured exec logs show `chi_0 = c_etaU*gamma/KU` in sympy (log:14) and `chi_0 = (cEtaU*gamma)/kU` in mathematica (log:10) with no `rho_0 =`, `sigma_0 =`, `PASS: rho_0 - chi_0`, or `PASS: sigma_0 - chi_0` lines anywhere, satisfying verification condition (b) of the directive. Downstream usage of `chi0` is unchanged (sympy:80, 89-90; .wl:68, 75, 82, 92), so removing rho0/sigma0 has no propagation hazard. No tautological assertion remains in this block.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:77-90` — Replaced the `supportLoadFactor = zeta*(1 - eps)/(1 - zeta*eps); mSupp = FullSimplify[mMix*supportLoadFactor, ...]` block plus `sEnhance = FullSimplify[mTr/mMix, ...]; sClosedForm = ...; expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm]` with an independent closed-form `mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)), ...]` matching notes §4, kept `mTr = FullSimplify[mMix + mSupp, ...]`, and redefined `sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), ...]` from the closed-form S of Eq. app-stage047-S. The tautological `expectZero["S from ratio agrees with closed-form S", ...]` is deleted.

**Assessment:**
The change matches the directive's "After" block exactly. The retained `expectZero["M_tr - M_mix S", mTr - mMix*sEnhance]` (.wl:99) is now substantive: both `mMix` (line 75) and `mSupp` (line 82) are written as independent closed forms in the dimensionless ratios `chi0, epsEta, eps, zeta`, so `mMix + mSupp` is *not* algebraically `mMix * (1 + zeta(1-eps)/(1-zeta*eps))` by string construction — the equality is now a real algebraic identity that `FullSimplify` must establish. The exec log line 28 confirms `M_tr - M_mix S = 0` and line 29 shows `PASS: M_tr - M_mix S`. The exec log also lacks the line `PASS: S from ratio agrees with closed-form S`, satisfying directive verification condition (c). All other PASS lines required by directive condition (d) are present: `PASS: M_tr - M_mix S` (log:29), `PASS: product law` (log:33), `PASS: support-loaded R_target reconstruction` (log:35), `PASS: dR_target_loaded/dzeta` (log:37), `PASS: R_target_loaded(zeta) - R_target_loaded(0)` (log:39), `PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2` (log:41), `PASS: S(zeta=0)-1` (log:43). The printed `M_supp` (log:24) is algebraically equivalent to the sympy `M_supp` (sympy log:39) after distribution.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `chi_0 = c_etaU*gamma/KU` (no rho_0/sigma_0 lines — F1 resolved)
- `M_tr - M_mix*S = 0` (still substantive in sympy)
- `product law = 0`, `support-loaded R_target reconstruction = 0`, `dR_target_loaded/dzeta = 0`, `R_target_loaded(zeta) - R_target_loaded(0) = 0`, `dS/dzeta - (1-eps)/(1-zeta eps)^2 = 0`, `S(zeta=0)-1 = 0`
- `spoiled dR_target/dzeta at a probe point = -44.866...` (negative-control probe A14 still firing non-vacuously)
- `All Stage-30 symbolic checks passed.` (banner still says "Stage-30" — see side observations)

**Mathematica:** exit=0. Notable lines:
- `chi_0 = (cEtaU*gamma)/kU` (no rho_0/sigma_0 lines — F1 resolved)
- `M_supp = (-72*(cEtaU*gamma + kU)^2*lamW^2*(ell^2*kU + Pi^2*tU)*zeta) / ((cEtaU^2 - kEtaEff*kU)*(9*kU*kWEff*Pi^2*(ell^2*kU + Pi^2*tU) - 8*gamma^2*lamW^2*(11*ell^2*kU + 9*Pi^2*tU)*zeta))` — expanded normal form matches the sympy `M_supp` after distribution
- `PASS: M_tr - M_mix S` (F2 resolved — now substantive)
- All other required `PASS:` lines present, no `PASS: S from ratio agrees with closed-form S`
- `Stage 047 Mathematica audit passed.`

**Output freshness:**
The exec logs at `redteam/exec_logs/stage_047_sympy.log` and `stage_047_mathematica.log` were captured 2026-05-26 01:52, post-edit (scripts mtime 01:46). However, the persistent saved transcripts `scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt` and `mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt` retain their 2026-05-22 12:57 mtime — they were NOT re-saved post-fix. The captured exec logs do reflect the post-fix runs, so verification stands on the logs, but the stale `output/*.txt` files still contain the pre-fix `rho_0 =`, `sigma_0 =`, `PASS: rho_0 - chi_0`, `PASS: sigma_0 - chi_0`, and `PASS: S from ratio agrees with closed-form S` lines. See side observations.

## Material-change assessment

`material_change`: false.

Both edits are script-only and remove or replace assertions while leaving the symbolic definitions used downstream (`chi0`, `Mmix`, `Msupp`, `Mtr`, `S`/`sEnhance`, `Rtarget`) algebraically unchanged:
- The sympy script's `chi0` is identical pre/post (only the `rho0`/`sigma0` aliases and their assertions were removed).
- The .wl script's `mSupp` is rewritten from `mMix*supportLoadFactor` to the closed form `8*zeta*zW*(1+chi0)^2/(Pi^2*(1-epsEta)*(1-zeta*eps))` — these are algebraically equal under `FullSimplify`.
- `sEnhance` is rewritten from `mTr/mMix` to `1 + zeta*(1-eps)/(1-zeta*eps)` — algebraically equal.

No derived result that downstream units (048-050) consume is changed. The orchestrator's default `upstream_stale: true` flag on units > 047 can be a no-op cross-check.

## Side observations (non-blocking)

1. **Persistent output transcripts not refreshed.** `scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt` and `mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt` still contain pre-fix content (mtime 2026-05-22 12:57). The exec logs in `redteam/exec_logs/` were captured post-fix and are authoritative for this verification, but the persistent `.txt` outputs that future auditors inspect for "output freshness" are now stale. The orchestrator should regenerate these via `$RT exec-sympy 047` and `$RT exec-mathematica 047` (single-seat Mathematica policy in effect) so the saved transcripts agree with the current scripts.
2. **Banner labels still say "Stage 30" / "Stage 030".** The sympy script's banner (`STAGE 30 — COHERENT KERNEL MAP AUDIT`) and closing line (`All Stage-30 symbolic checks passed.`) and the .wl `STAGE 030 — COHERENT KERNEL MAP` banner are stale relative to the unit's renumber to 047. Not part of any finding; flagging for future cleanup.
3. **Channel-saturation provenance pointer.** The replacement comment in both files cites "Stage 28" as where the matching condition is established. The auditor's original report referenced "Stage-28 derivation" in its required-change text, so this is consistent with the directive; not verified against the actual Stage 028 unit by this verifier (out of scope).

## Verdict justification

Both findings are `resolved`. F1's no-op guards and tautological assertions in sympy:41-74 and .wl:37-69 are removed and replaced with the directive's prescribed comment+chi0-only block; downstream usage is unaffected. F2's `mSupp = mMix*supportLoadFactor` and `sEnhance = mTr/mMix` construction in .wl:97-115 is replaced with independent closed forms for `mSupp` and `sEnhance` from notes §4 / Eq. app-stage047-S; the kept `expectZero["M_tr - M_mix S", ...]` is now a substantive algebraic identity that `FullSimplify` must establish. Both engines exit 0; the captured exec logs show all required `PASS:` lines present and the tautological lines gone. No regressions in the diff, no collateral edits beyond the directive's named line ranges. `verified`.
