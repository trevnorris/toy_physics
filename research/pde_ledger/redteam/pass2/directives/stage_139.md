---
unit_id: 139
batch: IV.5
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T16:35:07Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 139

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`,
`summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what F1 requires. Do NOT touch
paper.tex, notes/, or any prose document. The SymPy reference script (`...stage139_..._sympy_audit.py`)
must stay UNCHANGED — this finding is about the Mathematica engine only.

After editing, RUN the Mathematica script (`math -script <wl>`) and iterate until it exits 0 with all
in-file checks passing. (The orchestrator independently re-runs both engines afterward.) Run under
`timeout 600`; a timeout is a failure — reformulate, do not raise the cap.

## F1 — mathematica_transliteration (USER-AUTHORIZED FULL RE-AUTHOR, 2026-06-08)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl` (whole file)

**Issue:**
The current `.wl` is a numerical mirror of the `.py`, not an independent second-engine audit. Concretely:
- it imports `piStar = SetPrecision[1.50882951349316, 30]` (wl:31) and `sQStar = SetPrecision[0.658075937605429, 30]` (wl:32) as hardcoded literals from Stage 134 — the same literals the `.py` imports (py:8-9);
- its "independent" `S_q` reconstruction (wl:82-84) evaluates the SAME hardcoded Stage-134 kernel closed form `Pi*(kappa*Tanh[kappa] + Pi*(Exp[-Pi]*Sech[kappa]-1))/((1-Exp[-Pi])*(kappa^2-Pi^2))` that the `.py` hardcodes (py:67-69) — so it is not independent between engines;
- the file self-describes the branch line as "a sanctioned mirror of the SymPy route" (wl:42).
Everything the `.wl` computes, the `.py` computes identically; Mathematica adds no independent verification. Per the dual-engine rule (a `.wl` is required wherever Mathematica CAN independently verify; test = "is it possible," not "is it necessary"), this is a finding. An independent route IS feasible (every other IV.5 dual-engine stage establishes Pi_* by root-finding and S_q by an independent integral): 142/144/146 obtain Pi_* via `FindRoot` on its defining equation, and 146/147 obtain `S_q(Pi)` by direct quadrature of the source moment `Integrate[Sigma*K_q, {x,0,1}]`.

**Required change (requirement + acceptance — YOU design the exact Mathematica route):**
Re-author the `.wl` so its load-bearing inputs are derived in-engine, not imported as literals:

1. **Derive `Pi_*` in the `.wl`** by root-finding its Family-1 compensation defining equation
   `g(Pi) == g_minus`, where `g(Pi) = 2*Pi*(2*Pi*Exp[Pi]+Pi_const)/((4*Pi^2+Pi_const^2)*(Exp[Pi]-1))`
   is the cos-channel source moment (note: `Pi_const` here means the mathematical constant π — keep
   the symbol bookkeeping straight; this is the same `gFormula` used in stages 142/144/146/147) and
   `g_minus = rF - Sqrt[1+rF^2]/2` is the lower compensated branch. Use a **cleared-denominator
   bracketed `FindRoot`** (the 144/146 pattern: bracket `{p, 1.4, 1.6}`, `WorkingPrecision -> 80`),
   NOT a hardcoded literal. Assert the cleared-denominator residual is `~0` at the found root (a
   genuine root self-check). Watch the sign: the cleared-denominator residual numerator carries
   `(Exp[p] - 1)`, NOT `(1 - Exp[p])` (the §131-F3/144 trap).
2. **Derive `S_q(Pi_*)` independently** — reconstruct it by direct quadrature of the mixed-channel
   source moment `Integrate[Sigma*K_q, {x, 0, 1}]` at `Pi -> Pi_*`, with
   `Sigma = Pi*Exp[-Pi*x]/(1 - Exp[-Pi])`, `K_q = Cosh[kappa*(1-x)]/Cosh[kappa]`, `kappa = Pi_const/2`
   (the 146/147 route), evaluated to high precision. (Equivalently you may `DSolve` the mixed D/N BVP
   for the channel profile and read the boundary slope, per stages 133/134 — your choice; the
   requirement is only that `S_q(Pi_*)` is NOT taken from the hardcoded Stage-134 kernel closed form
   as its source of truth.) The hardcoded kernel closed form may survive ONLY as the RHS of an
   `expectApprox`/`expectZero` cross-check, never as the computed value feeding the gains.
3. The gain algebra is irreducible and stays: `M_s = Pi_*/(1 - R_q*S_q)`, `M_q = -R_q*M_s`, on both
   the natural branch (`R_q^nat = (1-rF)^2/(1+rF^2)`) and the compensated branch (`R_q^comp = 1/4`).
   Now it is fed by the independently-derived `Pi_*` and `S_q(Pi_*)`.
4. **Remove the "sanctioned mirror" comment (wl:42)** and any other self-description as a port; the
   file should read as an independent route. `rF` stays the Stage-121 closed form (it is independently
   anchored against `1.77799353547498` already).

**Acceptance criteria (the verifier will check ALL):**
- The `.wl` no longer contains `SetPrecision[1.50882951349316, ...]` or `SetPrecision[0.658075937605429, ...]` as the SOURCE of `Pi_*` / `S_q(Pi_*)`; both are derived in-engine (FindRoot / quadrature or DSolve).
- No occurrence of the phrase "sanctioned mirror" (or "mirror of the SymPy"/"port of") remains.
- The derived `Pi_*` matches `1.50882951349316` and the derived `S_q(Pi_*)` matches `0.658075937605429` to `< 1e-12` (via `expectApprox`), so the independent route corroborates the Stage-134 values.
- All existing printed deliverables still match the notes literals to `1e-12`: `r_F1 = 1.77799353547498`, `R_q^nat = 0.145454452260421`, `g_-^F1 = 0.758035078944662826919680890414`, `M_s^nat,* = 1.66854252965624`, `M_q^nat,* = -0.242696939724365`, `M_s^comp,* = 1.80594111095636`, `M_q^comp,* = -0.451485277739090`, and `R_q^comp - 1/4` to `1e-25`.
- The `g_-^F1` branch-discrimination anchor (distinguishes lower 0.758 from upper ~2.79) and the `R_q^nat != 1/4` content are preserved.
- `math -script <wl>` exits 0; the committed `.py` is byte-identical to HEAD.

**Verification command:**
`redteam exec-mathematica 139` (and `exec-sympy 139` to confirm the `.py` is unchanged + still passes); confirm the new FindRoot/quadrature blocks appear, no hardcoded Pi_*/S_q literal feeds the gains, no "mirror" phrase remains, and all `expectApprox` lines PASS.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit to derive `Pi_*` by cleared-denominator `FindRoot` and `S_q(Pi_*)` by source-moment quadrature before computing the gains.
- deviation: none
