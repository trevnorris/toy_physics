# Decision 11 — Path-A GATE-A freeze sheet (S_Σ constitutive family + DOF + field→D0 map)

**Date:** 2026-06-18
**Status:** **FROZEN — `frozen: YES`** (user-authorized 2026-06-18, pre-reg §M Stage-2 freeze gate). Was:
Claude+Codex determination COMPLETE, pending freeze.
**Freeze record (G4 / decision-07 analog):**
- `frozen: YES`; `signed_off_by: user`; `freeze_date: 2026-06-18`.
- `candidate_freeze_hash: ed358569393fed5fc29c0c13286a07cd438db467da6c1bc663a09bb04b1691c9` (SHA-256 over the §8
  canonical spec; byte-reproducible — two stamp runs gave the identical hash, independently re-verified).
- `source_revision (git_head at freeze): 8bd82b9`; frozen packet:
  `software/stage1_solver/frozen/pathA_gate_a/ed358569393fed5fc29c0c13286a07cd438db467da6c1bc663a09bb04b1691c9/`
  (`freeze_sheet.json` + `freeze_hash.txt`). Frozen family registered as the
  `homogeneous_isotropic_hooke_v1` provider in `src/stage1_solver/patha_static_balance.py` (transliteration-fidelity
  audited FAITHFUL to §1; 11 solver tests pass; firewall untouched; sheet carries NO target residual / `R_norm` /
  `D0` / pass-flag). Freeze stamped BEFORE any calibration (freeze-before-solve discipline).
This is the §M **item-1 gating artifact** that converts the Path-A pre-registration DRAFT into a freeze-ready
record (`docs/pathA_preregistration.md` §E/§K/§M; decision-07 analog).
**Mechanism:** Claude+Codex consult ([[claude-codex-resolve-math]]). Prompt + full log:
`_scratch/pathA_gateA_sigma_family_consult_prompt.md` / `…_consult.log`. Conceptual forks put to + authorized by
the USER (freeze-prep, 2026-06-18): naturalness = **report value, no hard exclusion**; `T_Ω` = **tied/isotropic**;
and — the decisive call — **NO inserted anharmonicity DOF**: the primary model is the strict harmonic wall
(`g=0`, 1 material scale `τ`), and a miss is treated as a **missing-physics signal handled κ_PV-style**, never as
a cue to add a rescue parameter.
**Target-blind:** this sheet SPECIFIES the family + DOF + extraction map + calibration protocol. It does NOT run
calibration, compute `R_norm`/`D0`, or export. No target residual / pass-flag / target value appears here.

## Consult provenance + the methodology call (why no `g`)
- **Codex** opened with a strict single-scale Hookean wall (1 DOF `τ`, constant moduli, `R_*=a`, no quartic,
  `T_Ω=τ/a²`). Good calls retained: pin `R_*`, drop the sinusoidal-in-`w` placeholder modulations (fitting
  artifacts), tie `T_Ω` (isotropic), the eigenproblem K-reduction (§5).
- **Claude** observed that with constant moduli `U_{Σ,RR}=τ/a²`, `T_{w,Σ,R}=0` ⇒ the l=2 stiffness operator
  `−τη''+7τ/a²η` contains **no `R0`** ⇒ `K=τ·κ̂` (κ̂ fixed geometric) is **independent of the self-consistent
  background**, so decision-08 **mechanism #2** ("the self-consistent background changes `K`") acts only weakly
  (through `B0/Z0` and `R0(τ)`), not through the dominant denominator term. The fix CONSIDERED was a single
  anharmonicity `g` to make `K=K[R0]` background-dependent.
- **USER call (the integrity decision):** `g` would be an **inserted DOF reached for to make the machinery behave**
  (revive mechanism #2), not a contribution any derivation produced — exactly the fudge the program forbids. By
  the κ_PV precedent (`n=5 → β=3` required, two known κ's summed to 1.5; the missing 1.5 was NOT a free knob but a
  **derived** third contribution `κ_PV`), a discrepancy must diagnose **missing physics whose contribution is
  determined**, never a tuned parameter. **Decision: drop `g` (g=0), keep the strict harmonic wall, and treat a
  miss as a missing-mechanism signal (§4b), not a DOF to add.** Honest cost (disclosed, not hidden): with `g=0`,
  `K` is background-independent — itself a reportable finding ("linear elasticity alone cannot move `R_norm`
  through the self-consistent background").

## 1. Frozen constitutive family (target-blind)
Promoted wall action `S_Σ[R]=∫ dt dw dΩ L_Σ`,
`L_Σ = ½ μ_Σ R_t² − ½ T_{w,Σ} R_w² − ½ T_{Ω,Σ}|∇_Ω R|² − U_Σ`. The 4-term form is the symmetry-allowed,
two-derivative effective action (derivation-backed); the coefficients are the posited material law fixed here.
**Family id:** `homogeneous_isotropic_hooke_v1` (strict harmonic; calibrated parameter: `τ` only).

| function | FROZEN form | domain | units (G=c=c_s=1, action∼L²) | class | admissibility |
|---|---|---|---|---|---|
| `μ_Σ(R,w)` | `τ` | `w∈[0,L]`, `R>0` | `[μ_Σ]=1` | F (1-scale tie) | `τ>0` |
| `T_{w,Σ}(R,w)` | `τ` | `w∈[0,L]`, `R>0` | `[T_{w,Σ}]=1` | F (1-scale tie) | `τ>0` |
| `T_{Ω,Σ}(R,w)` | `τ/a²` | `w∈[0,L]`, `R>0` | `[T_{Ω,Σ}]=L⁻²` | F (isotropic tie) | `τ>0` |
| `U_Σ(R,w)` | `½(τ/a²)(R−a)²` | `R>0` | `[U_Σ]=1` (dimensionless; density in action∼L²) | F | `τ>0` (bounded below, single well at `R=a`) |

Constitutive derivatives (analytic; `[U_{Σ,R}]=L⁻¹`, `[U_{Σ,RR}]=L⁻²`):
`U_{Σ,R}=(τ/a²)(R−a)`; `U_{Σ,RR}=τ/a²` (constant); `T_{w,Σ,R}=T_{w,Σ,RR}=0`; `μ_{Σ,R}=0`; `T_{Ω,Σ,R}=0`.

**Physical rationale (target-blind):** a passive, **homogeneous, isotropic, linearly-elastic (Hookean) wall
membrane** — one elastic scale `τ` (inertia + longitudinal + transverse moduli) and a harmonic confinement well
with preferred radius `R_*=a`. This is the leading symmetry-allowed effective wall; nothing is inserted beyond it.
The inhomogeneous `R0(w)` profile is **DERIVED** from the static balance + boundary conditions + the matter/gauge
return source. (Anharmonicity is NOT posited; if a miss later shows it is physically required, its coefficient is
to be DERIVED, not fitted — §4b.)

## 2. Boundary / structural class (frozen; NOT calibrated DOF — disclosed per §K)
- Mouth: Dirichlet `R0(0)=a`. Exit: natural zero-traction `T_{w,Σ}R0'(L)=0` (`Y_L=0` open limit). No hard cap;
  require `R0(w)>0`.
- Perturbation (modal) BCs: `η(0)=0` (anchored mouth), `T_{w,Σ}η'(L)=0` (natural exit).
- No counterpressure / `ρ_ref` subtraction (absolute `ρ0` source, decision-10). No moving-boundary EM term
  (excluded; cannot be added post-hoc, §K).
- Geometry pins (branch, decision-07): `a=1`, `L/a=37/20` ⇒ `L=1.85`.
- **Source/port conventions that scale `R_norm` (pinned, per decision-07 lines 33–34): `m̂0=1`
  (point-particle natural source-map limit), `S_port=1` (N0 in gravitationally-normalized port convention).**
- Units: `G=c=c_s=1`; benchmark stays the GR quadrupole `54 G c_s⁵/(5 a⁵ c⁵)`.

## 3. Complete DOF count (honest; §E/§K)
**Calibrated continuous DOF: `{τ}` = 1.** `[τ]=1`.
**Pinned (NOT calibrated):** `g=0` (no anharmonicity); `R_*=a`; the ties `μ_Σ=τ`, `T_{w,Σ}=τ`, `T_{Ω,Σ}=τ/a²`;
w-homogeneity of all moduli; `m̂0=1`, `S_port=1`.
**Declared discrete/auxiliary choices that must be COUNTED (the "hidden DOF" §K demands):**
- modal channel = `l=2`; mode branch = **lowest positive** `L₂` eigenmode; mode normalization (`∫χ²=1`) +
  orientation (`∫χ>0`) fixed before extraction;
- **root branch = the stable Schur side `D0>0`** when selecting the `τ` that solves `R_norm=0` (v2_09 stability
  gate `D0>0`; a branch-selection switch);
- optimizer/root-finder + tolerances + mesh/grid-convergence ladder, frozen before calibration;
- **no candidate-family search and no rescue DOF** (this single harmonic family is the pre-committed model; a miss
  is handled by §4b, which adds NO free parameter).
- NOTE: 1 DOF vs 1 anchor ⇒ calibration is UNIQUE → there is **no naturalness-penalty tie-breaker** and **no `w_g`
  weight** to count.

## 4. Calibration objective + optimizer (1 DOF)
- **Anchor (open):** `R_norm = 0`, i.e. the GR quadrupole `P0 = N0/D0 = 54/5` under `m̂0²S_port=1` (Peters 1964 /
  `benchmarks.yaml`). `R_norm` is the calibration anchor, **NOT a held-out prediction** (decision-09; pre-reg §H).
- **Procedure:** 1 DOF vs 1 anchor ⇒ solve `R_norm(τ)=0` by deterministic root-finding on the converged closed
  background (each evaluation = a full closed solve + §5 extraction), selecting the **stable-side root `D0>0`**.
  Unique; no tie-breaker, no regularization weight.
- **Naturalness reporting (user choice):** report `τ*` and the naturalness diagnostic (`|ln τ*|`; the implied
  cancellation ratio `K/(B0+Z0)` and digit count). **Report value, no hard exclusion**: a `τ*` far from O(1)
  (ultrasoft wall) or a many-digit `K≈B0+Z0` near-pole is **reported as fine-tuned, NOT a clean physical pass.** A
  soft pre-declared "natural band" for `τ` accompanies the report as context only (non-exclusionary).

## 4b. Miss-response protocol (κ_PV discipline — NO rescue DOF)
A miss = the GR anchor is reachable only at an unnatural/fine-tuned `τ*` (the decision-08 honest prior:
`K=4.06 ≫ B0+Z0≈0.0047`, an ~`1.6e7` cancellation). **A miss is NOT rescued by inserting a constitutive DOF**
(no free `g`, no extra wall coefficient tuned to the target). It is treated as a signal that a physical
contribution is unaccounted for — diagnosed and **DERIVED**, κ_PV-style (the missing `1.5` in `β` was a derived
third term, not a free knob). Candidate re-examinations, each yielding a **determined** (not fitted) contribution
that is then re-frozen target-blind before any re-test:
- a missing term/coupling in `S_Σ` the symmetries permit and the microscopic physics demands;
- a **derived** anharmonicity — the actual quartic from matching `U_Σ` to the microscopic gauged-GPE free energy
  near the wall (a determined coefficient; legitimate precisely because it is derived, unlike a posited/ fitted
  `g`);
- a missing contribution to `K`, `B0`, or `Z0` (the direct κ_PV analog: the bundle, not a knob);
- re-examination of the throat / boundary / source physics.
The map κ_PV↔Path-A: `β=3` ↔ `D0→0` (`K→B0+Z0`); "known κ's sum to 1.5" ↔ `K≫B0+Z0`; "don't fudge the gap" ↔
"don't soften `K` with a free parameter"; "find + derive the missing contribution" ↔ re-examine the wall bundle.

## 5. Closed-background field→coefficient extraction map (§M item 2)
How `K, M, B_n, Z_n, N_n, D_n, P_n, R_pole, R_norm` are read off the converged closed solution `{ψ0, A0, R0}` (the
not-yet-built extraction module; SPECIFIED here, BUILT post-freeze as calibration machinery; never touches the
firewall `research/pde_audit/simulation/` or `physical_export_permitted`).

**(a) Static wall stiffness `K`.** On the converged `R0(w)`, form the l=2 wall stiffness operator (second variation
of `S_Σ` about `R0`, frozen BCs):
```
K_η(w) = U_{Σ,RR}(R0,w) − ∂_w[T_{w,Σ,R}(R0,w) R0'] + ½ T_{w,Σ,RR}(R0,w) (R0')²
       = τ/a²                                        (g=0; U_{Σ,RR} constant; the T_{w,*} terms vanish)
L₂ η = −∂_w(T_{w,Σ}(R0,w) ∂_w η) + [K_η(w) + 6 T_{Ω,Σ}(R0,w)] η ,   η(0)=0, T_{w,Σ}η'(L)=0
     = −τ η'' + 7τ/a² η
```
Discretization-ready: assemble symmetric stiffness `A₂` for `L₂`, the norm matrix `W` and mass matrix `M_μ` (with
`μ_Σ`); solve for the **lowest positive** mode, fix normalization/orientation:
```
A₂ χ = K W χ ,   χᵀ W χ = 1 ,   ∫χ dw > 0       ⇒   K = χᵀ A₂ χ = τ·κ̂ ,   M = χᵀ M_μ χ
```
**Honest note (disclosed):** for the harmonic family `K = τ·κ̂` with κ̂ a fixed geometric eigenvalue — `K` does
NOT depend on the solved `R0(w)`. The self-consistent background enters `D0` only through `{B0, Z0, N0}` (overlaps
with the solved `{ψ0,A0}`) and through `R0(τ)` shifting those overlaps. This is the reportable finding of §0/§6,
not a defect to patch.
**(b) Back-reaction + transfer moments** (reuse the SAME `χ` in every overlap; kernels frozen target-blind).
Field→coupling map: BdG `c_j = λ_B I_{η,φ_j}`; mixed-port `g_U = λ_U I_{η,u}`, `g_W = λ_W I_{η,w}`,
`R_mix = λ_R I_{u,w}` (overlaps of `χ` with the solved `{ψ0,A0}` mode profiles). Moments (ω⁰ AND the held-out
higher moments needed for `R_pole`, `P2`, `P4`):
```
B0 = Σ_j c_j²/ϖ_j²,  B2 = Σ_j c_j²/ϖ_j⁴,  B4 = Σ_j c_j²/ϖ_j⁶
per mixed port:  Δ = Ω_U²Ω_W² − R_mix²,  Q = g_U²Ω_W² + 2 g_U g_W R_mix + g_W²Ω_U²,  P = Ω_U² g_W + R_mix g_U
Z0 = Σ Q/Δ ;  Z2, Z4 ;  N0 = Σ P²/Δ² ;  N2, N4         (Z_n,N_n per the v2_21/v2_22a moment expansions)
D0 = K − B0 − Z0 ;  D2 = −(M + B2 + Z2) ;  D4 = −(B4 + Z4)
P0 = N0/D0 ;  P2 = (D0 N2 − 2 D2 N0)/D0² ;  P4 = (D0² N4 − 2 D0(D2 N2 + D4 N0) + 3 D2² N0)/D0³
R_norm = m̂0² S_port P0 − 54 G c_s⁵/(5 a⁵ c⁵) ;  R_pole = D0(B4+Z4) − 3(M+B2+Z2)²
```
This is the existing V2-21/V2-22a extraction algebra
(`research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py`,
`…stage_v2_22a_profile_to_coefficient_adapter.py` — READ-ONLY references; formulas verified faithful 2026-06-18,
incl. `D2/D4/P2/P4/R_pole`). `D0 = K − B0 − Z0` reconciles with pre-reg §F `D0 = K_* − Q/Δ` since
`K_* = K − B0(BdG)`.
**Solved-field-dependent inputs:** `{ψ0,A0,R0}, ρ0`, `χ`, the BdG/mixed mode spectra + overlaps, hence
`{M, B_n, Z_n, N_n}` and via them `{D_n, P_n}` (note `K=τκ̂` is NOT background-dependent for this family).
**Frozen inputs:** the constitutive family + ties + `g=0`, boundary class, `m̂0=1`, `S_port=1`, channel/
mode-selection + normalization, the port/kernel conventions (`λ_B, λ_U, λ_W, λ_R, Ω_U, Ω_W, ϖ_j, …`), the
extraction formulas, constants/units, optimizer/tolerances.

## 6. Mechanism statement + honest prior (decision-08)
Path A moves `R_norm` **only** through the denominator `D0 = K − B0 − Z0` / the self-consistent background — never
through a new numerator (reciprocity ⇒ `N0` forward-determined). For the harmonic family the operative levers are
**the `τ` scale (`K=τκ̂`, mechanism #1) and the near-pole (#3)**; the background acts on `D0` only through
`{B0,Z0,N0}` (mechanism #2, partial). The test: does the self-consistent throat equilibrium reach the Schur
softening edge `D0→0` (`K→B0+Z0`) at a **natural** `τ`? **Honest prior:** `K=4.06`, `B0+Z0≈0.0047` ⇒ reaching the
anchor needs `K→~0.0047`, an ~`1.6e7` cancellation, so the likely outcome is a MISS at an absurd `τ*`. Under the
chosen policy that miss is **not** rescued — it launches the κ_PV re-examination (§4b). A target-blind natural
pass would be a genuine breakthrough; a clean documented miss is the strongest evidence yet against the model and
the well-posed start of "what physics are we missing," not a failure to paper over.

## 7. §J validation status (item 4)
- **DONE (chunks 1a–1c, committed):** static-balance operator + dual-engine MMS (1a); closed self-consistent
  Newton over `{ψ,A,R0,μ}` (1b); self-consistent-balance validation + ≥3-level closed grid-convergence +
  conservation + genuinely-independent flux-discretization gate (1c). These ran on the PLACEHOLDER families.
- **Remaining for the frozen family (post-freeze calibration machinery):** re-run the closed solve with the frozen
  `homogeneous_isotropic_hooke_v1` (`g=0`) forms; build + validate the §5 extraction module; the `R_norm(τ)=0`
  root-find; calibration-covariance propagation into the held-out residuals; margin-to-Schur (`D0`) error bars.
  The closed solve is still required (it supplies `R0(τ)` and the `{B0,Z0,N0}` overlaps) even though `K=τκ̂`.
  "A result without §J does not count" (pre-reg §J).

## 8. Target-blindness + freeze-hash basis (`candidate_freeze_hash_def`, G4 analog)
Canonical spec hashed at freeze, over: `parent_action_status = path_A_promoted_S_Σ`; branch identity;
`family = homogeneous_isotropic_hooke_v1` + the 4 forms + ties + **`g=0`** + domains/units; `R_*=a`; geometry
`{a, L, boundary_class, Y_L=0, no_hard_cap}`; **source/port `{m̂0=1, S_port=1}`**; gauge convention; the §4
calibration objective + root-finder + tolerances + stable-side branch selection; the §5 extraction formulas;
mesh/convergence ladder; channel/mode-selection + normalization. **Hashed:** the calibration *objective* including
the external GR-benchmark anchor constant `54Gc_s⁵/(5a⁵c⁵)` (the thing calibrated TO — per pre-reg §K and
decision-07). **Excluded:** any target *residual*, pass/fail field, computed `D0`/`R_norm` value, or
held-out-prediction target value.

## 9. Conceptual flag + the gate
Promoting `S_Σ[R]` to a dynamical field is THE conceptual change (already user-GO'd at the §M Stage-1 conceptual
gate). This sheet is the §M item-1 gating artifact; item 2 (extraction map, §5) and item 3 (engine = extend the
torch WP1, decision-10) are satisfied; §7 covers item 4. **NEXT = USER `frozen: YES`** (Stage-2): set
`frozen: YES` + compute `candidate_freeze_hash` over §8 BEFORE any calibration. Until the user freezes, the build
stays target-blind: NO calibration / `R_norm` / `D0` / export / GATE-A-against-anchor.
