---
unit_id: 047
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T13:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 047

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py:41-67`: The single-line constructions of `rho0` and `sigma0` were replaced with named `rho0_num/rho0_den` and `sigma0_num/sigma0_den` subexpressions followed by `sp.simplify(num/den)`. Two new `assert` statements check the channel-saturation rule `num*KU == den*(gamma*c_etaU)`.
- `mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:37-63`: Mirror change — separate `rho0Num/rho0Den/sigma0Num/sigma0Den` definitions, FullSimplify on the ratio, plus two `If[...Exit[1]]` guards that print "FAIL: rho0 channel-saturation rule violated" / "FAIL: sigma0 ..." if the saturation rule fails. The existing `expectZero["rho_0 - chi_0", ...]` / `expectZero["sigma_0 - chi_0", ...]` calls at lines 68-69 are preserved as the directive required.

**Assessment:**
Edit matches the directive verbatim, including comments and structure. Codex implemented the directive's minimal-change variant (split numerator/denominator and add channel-saturation guards) rather than the auxiliary-alias fallback. The directive itself explicitly anticipated that the new check would test whether num/den pre-satisfy the saturation rule rather than constructing rho_0 from a physically independent kernel, and the auditor accepted that as a defensible minimal change. The post-fix sympy/mathematica outputs still print `rho_0 = c_etaU*gamma/KU` (post-simplify) and the assertions all pass. No collateral edits outside the cited line range; the rest of the script (assertion list at lines 73-74 / wl 68-69) is unchanged. Mtimes confirm both outputs were regenerated post-fix (12:57 > script mtime 12:55).

Side caveat (non-blocking, recorded under "side observations"): the new guard `rho0_num*KU - rho0_den*gamma*c_etaU` is itself a symbol-cancellation identity (`gamma*lamW*c_etaU*KU - KU*lamW*gamma*c_etaU` = 0 by associativity), so it is non-tautological only in the sense that any future regression in `rho0_num` or `rho0_den` would be caught before the final assertion. This matches the directive's stated semantics and is not a regression.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.wl:97-115`: `mSupp` is now constructed as `FullSimplify[mMix * supportLoadFactor]` with `supportLoadFactor = zeta*(1-eps)/(1-zeta*eps)`, instead of being copied from the SymPy script's closed form `8*zeta*zW*(1+chi0)^2/(Pi^2*(1-epsEta)*(1-zeta*eps))`. `sEnhance` is now derived as `FullSimplify[mTr/mMix]` rather than written directly. A new cross-check `expectZero["S from ratio agrees with closed-form S", sEnhance - sClosedForm]` (line 115) verifies the ratio-derived `sEnhance` equals the closed-form `1 + zeta*(1-eps)/(1-zeta*eps)`.

**Assessment:**
Edit matches the directive verbatim. The Mathematica path to `mSupp` and `sEnhance` is now distinct from the SymPy path (sympy builds `Msupp` from the closed form and `S` from the closed form independently; Mathematica builds `mSupp` from the support-loading rule applied to `mMix` and `sEnhance` from the ratio `mTr/mMix`). The new `PASS: S from ratio agrees with closed-form S` line is present in the captured `.txt` output. The final closed-form expressions still agree with the SymPy outputs after FullSimplify (verified by inspection of `M_supp = ...` and `S(zeta;eps) = ...` lines in both `.txt` files — they are algebraically equivalent rational functions of the same kernel symbols, just rewritten under different simplification choices). All existing `PASS: M_tr - M_mix S`, `PASS: product law`, `PASS: support-loaded R_target reconstruction`, `PASS: dR_target_loaded/dzeta`, `PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2`, `PASS: S(zeta=0)-1` lines remain present. No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (canonical transcript at `scripts/output/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.txt`; the rendered prompt's `redteam/exec_logs/stage_047_sympy.log` path does not exist — orchestrator wrote the post-fix transcript to the canonical output instead). Notable lines:

```
rho_0 = c_etaU*gamma/KU
sigma_0 = c_etaU*gamma/KU
rho_0 - chi_0 = 0
sigma_0 - chi_0 = 0
M_tr - M_mix*S = 0
spoiled dR_target/dzeta at a probe point = -44.8664510213650197722046979701
All Stage-30 symbolic checks passed.
```

The new internal `assert` statements pass silently (no AssertionError text in the transcript), and the final "All Stage-30 symbolic checks passed." line indicates the script ran to completion past the assertions and the spoiler probe.

**Mathematica:** exit=0 (canonical transcript at `mathematica/output/moving_throat_pde_stage047_coherent_kernel_map_mathematica_audit.txt`). Notable lines:

```
PASS: rho_0 - chi_0
PASS: sigma_0 - chi_0
PASS: S from ratio agrees with closed-form S
PASS: M_tr - M_mix S
PASS: product law
PASS: support-loaded R_target reconstruction
PASS: dS/dzeta - (1-eps)/(1-zeta eps)^2
Stage 047 Mathematica audit passed.
```

No `FAIL:` lines, no `$Failed`, no traceback. The new `PASS: S from ratio agrees with closed-form S` line (F2 cross-check) is present at line 25 of the transcript. No `FAIL: rho0 channel-saturation rule violated` / `FAIL: sigma0 ...` lines (F1 guards passed silently).

**Output freshness:** confirmed — both `.txt` outputs (12:57) are newer than both edited scripts (12:55). Diff at `redteam/exec_logs/stage_047_diff.patch` (12:58) post-dates the script edits.

## Material-change assessment

`material_change`: false.

The edits in F1 split intermediate expressions but preserve the value of `rho0`, `sigma0`, and `chi0` exactly (they reduce to the same simplified `c_etaU*gamma/KU`). The edits in F2 reorder how `mSupp`, `mTr`, and `sEnhance` are constructed in the Mathematica engine but produce algebraically equivalent results to the SymPy engine — `M_tr - M_mix S = 0`, `product law = 0`, `support-loaded R_target reconstruction = 0`, and `R_target` closed forms continue to agree across engines. No derived result that a downstream unit might depend on has changed.

## Side observations (non-blocking)

- The F1 channel-saturation guard `rho0_num*KU - rho0_den*gamma*c_etaU` reduces to `0` by associativity/commutativity of `*` alone, without invoking simplify cancellation. So while it is structurally an "independent" check on the saturation rule, in practice the only way it could fail is if a future edit corrupted one of the four symbols. The directive accepted this as the minimal-change variant; flagging here only because a hypothetical iteration-2 polish could replace the guard with an independent channel-saturation rule defined via a `Subs` / replacement rule. Not blocking verification.
- The .py script's docstring and banner still call the unit "Stage 30" rather than "Stage 047". The .wl banner reads "STAGE 030 — COHERENT KERNEL MAP". Naming drift noted by the auditor in the inventory but not raised as a finding; verifier does not raise it as a new finding.
- The sympy script's spoiler probe (`eta_bad`, lines 144-177) has no Mathematica counterpart. The auditor noted this under F2(c) but the directive did not require porting it; verifier does not block on it.

## Verdict justification

Both findings F1 and F2 are resolved: Codex applied the directive verbatim across both engines, the post-fix transcripts run to exit 0 with all expected PASS lines (including the new `PASS: S from ratio agrees with closed-form S` for F2), output mtimes are fresh, and the diff shows no collateral edits outside the cited line ranges. The values of all downstream-relevant quantities (`rho0`, `chi0`, `M_tr`, `R_target`, `S`) are preserved bit-for-bit in simplified form; the edits restructured how they are constructed without changing what they evaluate to. Verdict: `verified`.
