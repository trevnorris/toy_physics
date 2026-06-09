---
unit_id: 180
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 180

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
The `.wl` was fully re-authored (diff touches only the `.wl`). The four deliverables are now reached by routes the `.py` does NOT use:

- **D1 multi-port collapse** (`.wl:45-52`): `teff2Base = t1^2 + t2^2`; `xiEff = scaledFirstVariation[teff2Base, {{t1, t1*tau1}, {t2, t2*tau2}}]`, where `scaledFirstVariation[expr, incs] = Σ_i (∂expr/∂x_i · δx_i)/expr` (`.wl:36-39`). This is an analytic first-variation of the bare base sum with port-load increments — NOT the `.py`'s `exp(2·eps·lam·tau)`-perturbed `diff(log(Teff2),eps).subs(eps,0)/lam` (`.py:41-42`). The exponential perturbation tag is gone.
- **D2 one-port continuum** (`.wl:60-74`): builds an equation system `continuumEquations` (`k0Sym·muEta==kEta`, `betaSym·muEta·kW·(1-epsW)^2==muW·kEta·zW·(1+rho)^2`, `transferSq·k0Sym==betaSym`) and `Solve`s for `transferSq`, then compares against both the `muW/KW` and `Ω_W²` forms. The `.py` instead hand-builds `beta0/K0` and `.subs`-es (`.py:61-71`). Branch-equation solving is a genuinely different route.
- **D3 selected-branch** (`.wl:83-91`): `selectedEquations` system + `Solve` for `selectedTransferSq`, vs the `.py`'s hand `Lambda`/`Rtarget_def`/`.subs` (`.py:82-89`). Independent.
- **D4 slope laws** (`.wl:103-117`): `scaledFirstVariation` partial-derivative route on the *static* closed forms `directTransferShape`/`selectedTransferShape`, vs the `.py`'s `exp(e·lam·…)`-tagged `T2_pert` then `diff(log,e).subs(e,0)/lam` (`.py:100-118`). Different operator, no perturbation tag.

**Assessment:**
The 173-lesson check passes. The old verbatim mirror (`D[Log[hand-coded teff2/t2Pert], …]/.→0` on the same objects) is removed in every section. The strongest independence gain is D2/D3, which now route through `Solve` of an equation system rather than `/.`-substitution of the same hand-form — a mis-transcribed closed form on the SymPy side would NOT be silently reproduced. D1/D4 use the analytic first-variation `Σ ∂expr·δx / expr` of static forms, a distinct symbolic operation from the `.py`'s exp-tag-then-log-differentiate. All six `expectZero` assertions remain non-tautological (constructed side vs separately-written closed form). Re-author is genuinely independent, not a cosmetic relabel or a port. Banner reads `STAGE 180`. Exec log exits 0 with all six PASS lines.

### F2 — stale_output

**Classification:** resolved

**What changed:**
No source edit (correct — directive said none needed). Both committed `.txt` were regenerated post-fix: mtimes are `2026-06-09 00:21:36` (sympy) and `00:21:36` (mathematica), newer than both scripts (`.py` 2026-06-03 15:59; `.wl` 2026-06-08 22:25). Line 8 (the banner) of both outputs now reads `STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE` (confirmed in both exec logs).

**Assessment:**
Stale-banner condition cleared. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE` (banner refreshed)
- `multi-port effective-shape identity = 0`
- `T^2 = ZW(1+rho)^2 / [OmegaW^2 (1-epsW)^2] = 0`; `selected-branch slope law = 0`

**Mathematica:** exit=0. Notable lines:
- `STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE`
- six `= 0` residuals each followed by `PASS:` (multi-port, beta0/K0→muW/KW, Ω_W², selected-branch T², direct slope, selected-branch slope)
- `Stage 180 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` mtimes (2026-06-09 00:21:36) are newer than the `.wl` (2026-06-08 22:25:38) and the `.py` (2026-06-03 15:59:11).

## Material-change assessment

`material_change`: false. The edit is method-only on the Mathematica side: all five deliverable values are unchanged and still verify to residual 0 by both engines. No derived numeric/symbolic result changed; the `.py` (the value-bearing engine) is untouched. Downstream units do not depend on a changed result. (Orchestrator may still flag units >180 `upstream_stale` per policy, but there is no specific value concern here.)

## Side observations (non-blocking)

- D1's `kA` factor in `nPort1=kA·t1²`, `nPort2=kA·t2²` cancels in `rho1=nPort1/(nPort1+nPort2)`, so the `kA` carrier is cosmetic but harmless; the identity still genuinely exercises `tau1`/`tau2` via `scaledFirstVariation`, so the check is real.
- `expectZeroUnder` for the Ω_W² form correctly supplies `omegaW2==kW/muW` as an assumption rather than `/.`-substituting, which is itself a cleaner independent treatment than the `.py`.

## Verdict justification

The USER-AUTHORIZED re-author (F1) is genuinely independent: the verbatim `D[Log[hand-coded …]]/.→0` mirror is gone in all four deliverables, replaced by equation-system `Solve` routes (D2/D3) and analytic first-variation partials of static closed forms (D1/D4) — routes the SymPy script does not use, so a shared hand-form mis-transcription could now fail. All six residuals are 0 with PASS, banner is STAGE 180, `math -script` exits 0, the `.py` is unchanged (diff is `.wl`-only), and both outputs are freshly regenerated (F2). Every deliverable is preserved and method-only. Verdict: verified.
