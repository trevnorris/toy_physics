# Directive — Path-A derivation #1: return sources `S_η^(ψ,A)` + the static self-consistent throat balance

**Program:** Path A (promote `S_Σ[R]` to a parent dynamical throat field; close the matter/gauge→wall return
loop). This is the ANALYTICAL derivation that must ground the Path-A pre-registration. **Target-blind analysis
only — NO numerical solve, NO comparison to any target, NO `R_norm`/`P0_target`/GR-constant value enters.**

**Role split:** Codex DERIVES + writes the verification (you choose the engine; symbolic variational calculus —
Mathematica `.wls` is well-suited and the 2nd engine, SymPy acceptable; cross-check in both where cheap, in the
dual-engine spirit). Claude REVIEWS + runs a transliteration-fidelity audit. Iterate your scripts until they
exit 0; the verifier reviews substance afterward. Scripts you RUN get `timeout 600`; do not wrap your own
session in a timeout.

## Context you must read first (ground every step in these — cite file:line)

- Parent action + promotion: `notes/moving_throat_pde_program_compact.md:196-218, 408-431, 480-528` (the
  `S_total = S_ψ[ψ,A;Σ] + S_EM[A] + S_Σ[R]` promotion; `S_Σ` compact form L510-528; `K_η` stiffness L523-528).
- Matter Lagrangian: `compact.md:556-588` (`L_ψ`, the single interaction term `−V_conf(X;Σ)ρ` at L564;
  `Σ = r − R(Ω,w,t)`, `R = R0(w)+η`). Maxwell sector + current: `compact.md:590-630, 632-685`.
- Wall PDE + sources: `compact.md:1259-1305` (`S_η^(2)` action, modal operator), `compact.md:1383-1455`
  (linearized skeleton, `δV_conf = −(V_wall'(Σ0/ℓc)/ℓc)η` L1424-1429, the boxed wall eq with
  `S_η^(ψ)+S_η^(A)+f_ext` L1441-1451).
- Schematic matter return: `research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md:344-356`
  (`S_η^(ψ) ~ −(V_wall'(Σ0/ℓc)/ℓc)δρ + …`).
- Reduced wall↔gauge algebra (where the sources feed): `research/pde_audit/notes/stage_v2_09_maxwell_mixed_
  kernel_derivation.md` (full); `compact.md:3946-4010` (`Δ, Q, P, D0=K_*−Q/Δ, N0=(P/Δ)², P0=N0/D0`).
- The crux resolution this derivation must honor: `software/stage1_solver/decisions/08_pathA_crux_resolution.md`
  (Fork B narrow: reciprocity ⇒ the return sources add **NO new numerator magnitude**; Path A's lever on
  `R_norm` is the **self-consistent background / `D0`**, via the **derived `R0(w)`** replacing M1c's hand-frozen
  smoothstep).

## Deliverables

A derivation note (markdown, committed-quality, at `software/stage1_solver/derivations/pathA_01_return_source_
and_balance.md` — create the dir) + the verification script(s) under `software/stage1_solver/mathematica/`
(and/or a SymPy module) that EXECUTE the symbolic checks. Bulk run artifacts may go under the gitignored
`runs/`; the note + script are the committed reproducer.

### D1 — Matter return source `S_η^(ψ)` (confirm exact reciprocity, complete it to the order the action fixes)
Vary the single interaction term `−V_conf(X;Σ)ρ` (`compact.md:564`) with `δΣ = −δR = −η`. Show **symbolically**:
- the forward wall→matter coupling is `δV_conf = −(V_wall'(Σ0/ℓc)/ℓc)η` acting on `ρ0` (the M1b/BdG `C_η`);
- the return matter→wall source is the SAME kernel against `δρ`: `S_η^(ψ) = −(V_wall'(Σ0/ℓc)/ℓc)δρ + [the
  higher terms the action actually fixes]` — derive the "+…" explicitly (do NOT leave it schematic; state
  exactly which terms the quadratic expansion of `−V_conf·ρ` around `(ρ0,R0)` produces, incl. any
  `−½ V_conf''·ρ0·η`-type self-term that belongs to `K_η`, vs the genuine matter-return cross term).
- Establish the **adjoint/self-adjoint structure**: forward and return derive from one symmetric bilinear in
  `(η, δρ)` ⇒ no new free coupling magnitude is introduced by the return. This is the key reciprocity claim.

### D2 — Gauge return source `S_η^(A)` (the genuinely new piece — "named, not formulated")
Derive `S_η^(A)` from the η-variation of the parts of `S_ψ + S_EM` that depend on the throat geometry through
the gauge sector. Determine **which** geometry dependence sources it:
- the matter current `J_ψ^M[ψ,A]` depends on the confinement (hence on `R`) through `ρ, ψ` → an A-coupled
  return via `δ(−A_M J_ψ^M)/δη` (matter-mediated, reciprocal of the forward `η→δJ_ψ→δA`); AND/OR
- any geometry dependence of the Maxwell sector itself. **Check explicitly whether the localization weight
  `Z(w)` is geometry-independent** (it appears as `Z(w)` only — `compact.md:594, 682, 1434`); if so, state that
  there is **no direct `η→Z` gauge↔wall coupling** in the declared action (consistent with decision-04 F2), and
  that `S_η^(A)` is therefore matter-mediated (reciprocal of the forward gauge route) — again **no new free
  magnitude**. If you find a genuine direct geometric gauge↔wall coupling (e.g. via moving boundaries of the
  open finite-radius conduit `R(0)=a, R(L)>0`, `compact.md:310-325`), formulate it and say so explicitly.
- Verify the result respects the mixed gauge invariants `E_w, C_a` (`stage_v2_09…:19-47`).

### D3 — The static self-consistent throat balance (where the physics lives)
Write the **stationary, isotropic (ℓ=0) wall-balance equation** of the PROMOTED closed system: the static part
of the wall PDE with the static return source, i.e.
`−∂_w(T_w∂_w R0) + [S_Σ restoring force](R0) = S_η^(ψ,A)|_static[ρ0(R0), A0(R0)]`,
whose solution `R0(w)` is the **derived self-consistent equilibrium** (NOT a free choice). Make explicit:
- exactly what the promotion of `S_Σ[R]` adds (the `μ_Σ R_t²`, `T_{w,Σ}`, `T_{Ω,Σ}`, `U_Σ` terms, `compact.md:
  510-528`) and how its quadratic reduction reproduces the SAME effective `μ_η/T_w/T_Ω/K_η` (so these stay
  posited free_choice — promotion does NOT secretly re-derive them);
- the difference from M1c: M1c froze `R0(w)` (the smoothstep, `decisions/07`); here `R0(w)` solves the balance.
- how the derived `R0(w)` then feeds the SAME downstream chain (BdG `B_n`, wall `K,M`, gauge `Z0=Q/Δ`, `N0`) →
  a self-consistent `D0 = K − B0 − Z0`.

### D4 — The no-new-numerator-knob confirmation (honor decision-08)
Show **explicitly** that closing the loop does NOT introduce a new term in the `N0` numerator `P = Ω_U²G_W +
R·G_U`: the return sources change only (i) the self-consistent background `(ρ0, A0, R0)` and (ii) the
denominator/Schur structure `D0`. State plainly that Path A can move `R_norm` ONLY through `D0`/background, and
that the GR target corresponds to the softening edge `D0 → 0` (`K → B0+Z0`; v2_09 Schur boundary
`stage_v2_09…:245-263, 497-507`). Do NOT compute any residual; just locate the lever.

## Acceptance criteria (the verification script must check, target-blind)
1. **Reciprocity:** forward kernel (`δV_conf` on `ρ0`) and return kernel (`S_η^(ψ)` on `δρ`) are the SAME
   symmetric bilinear — verified symbolically (e.g. the cross-Hessian `∂²S/∂η∂ρ` is symmetric). No free
   magnitude beyond what the forward chain already uses.
2. **Gauge invariance:** `S_η^(A)` is built from gauge-invariant data (`E_w, C_a` / the covariant current) —
   verified.
3. **Dimensional consistency** of every source term vs the wall-PDE LHS `μ_η∂_t²η` etc.
4. **No-numerator-knob:** symbolic confirmation that the return sources do not appear in `P = Ω_U²G_W + R·G_U`;
   they enter only background + `D0`.
5. **Reduction consistency:** the quadratic reduction of promoted `S_Σ[R]` reproduces the effective
   `μ_η/T_w/T_Ω/K_η` form (`compact.md:1269-1305, 523-528`) — promotion adds dynamics, not a re-derivation of
   the frozen constitutive functions.
6. **Honest open-items list:** anything the declared action does NOT fix (e.g. a genuine `S_Σ` constitutive
   posit, any boundary-term ambiguity) is listed as posited/open, not silently closed.

## Recurring sins to avoid (this project's calibrated failure modes)
- NO can't-fail / tautological "checks" presented as derivations (reciprocity must be a genuine symbolic
  identity, not asserted).
- NO hardcoded values; NO target leakage (no `10.8`, `54/5`, `R_norm`, `P0_target` anywhere).
- Do not dress a partial source as complete: if the "+…" terms or `S_η^(A)` are only partially fixed by the
  action, say exactly how far the derivation goes and what remains posited.
- Label machinery vs physics: this is a target-blind derivation of STRUCTURE; no physical claim about whether
  the target is met.
