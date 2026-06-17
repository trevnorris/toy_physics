# Directive — Spike-2: open-exit DtN + forward δJ_ψ transfer + Z0/N0 via the ω⁵ self-energy

You are the **CODER** (Codex codes / Claude reviews). Build Spike-2 of 3 (decision 05): the open-radiation
boundary + the forward wall→gauge transfer + the **basis-invariant outgoing transfer `N0`** (and `Z0`) that
makes the primary `R_norm` a clean test. Spike-1 (the genuine ℓ=2 VSH Maxwell operator) is DONE + committed
(`edb6830`); reuse its verified weak operator. NO V2-22B packet / NO `R_norm` yet (Spike-3). Pre-freeze,
target-blind, smoke background (so the numbers are MACHINERY, not physical — label loudly).

Follow the accepted design: `software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md` (D2, D3, D4,
D6) and `_scratch/mixed_maxwell_spike_design.log` (those sections). Implement D2/D3/D4; design the numerical
route. Reuse: the Spike-1 operator `mathematica/mt15_03_spike1_vsh_maxwell_operator.wls` (the weak `AAction`
VSH Maxwell operator — the MMS-verified one), the M1b BdG sector `mathematica/mt15_02_bdg_wall_derivation.wls`
(the BdG response + wall mode), and the step-8c phasor current.

## What to build

### A. Open-exit DtN / impedance BC (D2)
Replace Spike-1's closed test BC with the genuine OUTGOING radiation condition at `w=L` (and the radial outer
boundary): `Z C_a = Y_out,l(ω) E_a` on the physical tangential lanes (NOT `A=0`, NOT a hard cap). The
outgoing admittance/DtN branch `Y_out,l` is a **target-blind modeling choice where the canonical record is
silent — declare it and FREEZE it before computing any transfer number.** Mouth + radial-outer per D2
(no injected Maxwell drive at the mouth; one scalar gauge anchor). The operator becomes complex/non-Hermitian
at ω≠0 (outgoing loss) — that is correct.

### B. Forward wall→gauge source current `j_eta` (D3, no Path A)
Build the chain `η → δV_conf → δψ → δJ_ψ`: take the wall mode (from the M1b wall sector), apply the
`δV_conf = −4 V_radial r⁴/R0⁵ · η` drive to the M1b BdG operator to get the matter response `δψ`, then form
the linearized current with the **step-8c phasor Fréchet formula** (NOT `_matter_number_current`):
`δρ = 2(ψ_R0 δψ_R + ψ_I0 δψ_I)`, `δj_i = (ℏ/m)(δψ_R ∂_iψ_I0 + ψ_R0 ∂_iδψ_I − δψ_I ∂_iψ_R0 − ψ_I0 ∂_iδψ_R)
− (q/m)(δA_i ρ0 + A_i0 δρ)`, `δJ^0 = q δρ`, `δJ^i = q δj_i`. Confirm `S_η^(A)` (the return source) is NOT used.

### C. Z0/N0 via the retarded Maxwell self-energy (D4) — the crux
Solve the Spike-1 Maxwell operator (with the open DtN BC) for `δA` given source `j_eta`, at small target-blind
ω samples. Form `Σ_A(ω) = μ0 ⟨j_eta, G_A(ω) j_eta⟩`:
- `Z0 = Re Σ_A^cons(0)` (+ `Z2, Z4` from the low-ω expansion of the conservative part).
- `N0 = − lim_{ω→0} Im Σ_A^ret(ω) / [ Γ_port · ω⁵ ]` (+ `N2, N4` from higher orders), where `Γ_port` is the
  canonical outgoing-port normalization.

### D. ⚠ NORMALIZATION PINNING — DERIVE `Γ_port`, do NOT fit it to chi_Q=1 (acceptance criterion #1)
The design wrote `Γ_port = a⁵/27c_s⁵`; the compact's canonical `σ_Q^can = 9/(8 Ω_Q⁵) = 4a⁵/27c_s⁵` with
`Ω_Q = 3c_s/2a`. **Resolve this factor by DERIVING `Γ_port` from the physical outgoing-power (Poynting) flux
definition (D2) + the compact's canonical σ_Q^can — NOT by tuning it so chi_Q comes out 1.** chi_Q is the
TARGET; fitting the normalization to chi_Q=1 would rig `R_norm`. Instead: derive `Γ_port` independently, then
compute `chi_Q` on a reference canonical configuration as an INDEPENDENT CONSISTENCY CHECK — on the canonical
outgoing branch it should come out 1 *because the physics matches*, not because you set it. Report the derived
`Γ_port` (state whether it is `a⁵/27c_s⁵` or `4a⁵/27c_s⁵` and WHY, with the Poynting/σ_Q^can derivation), and
report chi_Q as a check with this caveat made explicit.

### E. Validation (D6 — these are the real gates)
- **Basis-invariance:** `N0` (and `Z0`) invariant under random finite-dimensional basis rotations of the
  internal gauge block (decision 04: U/W are a basis choice — the transfer must not depend on it).
- **V2-09 finite-block regression:** with a finite block `K_int=[[A,−R],[−R,B]]`, `g=(g_U,g_W)`, `ℓ=(0,1)`,
  reproduce `Z = Q/Δ` and `N0 = (Ω_U²g_W+R g_U)²/Δ²` (decision 05 D4 consistency check).
- **Convergence:** `Z0/N0` stable under mesh refinement AND ω-window refinement (report orders).
- **Outgoing flux positive** under the DtN BC (no hard-cap behavior); pure-gauge source → zero physical
  transfer.
- **Transliteration-fidelity targets:** the `δJ_ψ` current vs the 8c phasor formula; the DtN BC vs D2; the
  self-energy/ω⁵ extraction + `Γ_port` vs D4 + the compact.

## Scope OUT / honesty
NO V2-22B packet, NO `R_norm`, NO V2-22B extension (Spike-3). Do NOT invent `S_η^(A)`. The background is still
the M1b smoke (non-self-consistent) → `Z0/N0` are MACHINERY/method-demonstration, NOT physical (loud label;
physical values need M1c). Do NOT touch firewalled `research/pde_audit/simulation/`. No GATE-A.

## Working rules
≤2 `math`/`wolframscript` kernels (1 here); `timeout 600` per script RUN; iterate to exit 0; never wrap your
session in a shell `timeout`. Files under `software/stage1_solver/mathematica/` (e.g.
`mt15_04_spike2_transfer_n0.wls` + report); artifacts → gitignored `runs/`. Deterministic substantive outputs.
No `git add`/commit. No network/GPU. Target-blind.

## Report back
The DtN BC + frozen `Y_out,l`; the forward `j_eta` chain (+ that it matches the 8c phasor current, no `S_η^(A)`);
the derived `Γ_port` (a⁵/27c_s⁵ vs 4a⁵/27c_s⁵, with the Poynting/σ_Q^can derivation) and chi_Q as an
independent check (NOT a fitted knob); `Z0/N0` (+Z2/Z4,N2/N4) with the loud machinery label; basis-invariance
result; V2-09 regression match; convergence orders. Then Claude runs the fidelity audit (current/DtN/ω⁵/Γ_port)
+ adversarial (is Γ_port derived-not-fitted? basis-invariance + V2-09 genuine? any can't-fail?) + arbiter.
