# Decision 05 — mixed-Maxwell port spike: derivation design (REVIEWED + ACCEPTED) + build plan

**Date:** 2026-06-17
**Mechanism:** Claude+Codex DESIGN consult (read-only) — `_scratch/mixed_maxwell_spike_design_prompt.md` /
`…_design.log` (full D1-D6 derivation). User chose "design-first, then build". Claude reviewed the design
critically (skeptical read, not rubber-stamp) and ACCEPTS it as the build plan, with the flags below.
**Goal:** derive the mixed-Maxwell outgoing transfer (`N0`, `Z0`) from the parent action so the primary
`R_norm` becomes a CLEAN falsification test (decision 04: reachable WITHOUT Path A).

## The accepted design (condensed)

- **D1 — ℓ=2 vector-harmonic Maxwell operator.** Decompose δA in vector spherical harmonics: `δA_0=φ Y`,
  `δA_w=a_w Y`, brane 3-vector `= a_r Y e_r + a_E Ψ_lm + a_B Φ_lm` (`Ψ=∇_Ω Y` poloidal/electric,
  `Φ=e_r×∇_Ω Y` toroidal/magnetic). Assemble from the quadratic PARENT weak form
  `S_A^(2)=(1/2μ0)∫r²dr dw Z[|E|²−|B|²−|C|²−ξ⁻¹|G|²] − ∫A·δJ` (angular weights `L=l(l+1)=6` on the Ψ/Φ
  lanes), NOT from scalarized components. This restores the vector-Laplacian cross-terms the torch smoke
  dropped (`−2a_r/r²`, `2L a_E/r²`, `2a_r/r²`). Lane split: `a_B` dynamical+gauge-invariant; poloidal
  combination dynamical after constraint projection; pure gradients = gauge (deflate); `φ`+longitudinal
  enforce Gauss/Lorenz; `A_w` physical only via `E_w, C_a`.
- **D2 — open BCs + inner product.** `H=Z` gauge; `boundary_class=open_impedance` (no hard cap); r=0
  VSH-parity regularity; **open exit `w=L` = outgoing impedance/DtN** `Z C_a = Y_out,l(ω) E_a` (NOT A=0),
  DtN frozen before target eval; mouth = no injected Maxwell drive + one scalar gauge anchor (record silent →
  declared target-blind choice). Physical energy norm + Poynting outgoing flux on the `G=0` quotient;
  pure-gauge modes must have zero physical flux/norm (else reject).
- **D3 — forward transfer (no Path A).** `η→δV_conf→δψ (BdG, from M1b)→δJ_ψ→δA (Maxwell Green)`, using the
  **step-8c phasor Fréchet current** (`δρ=2(ψ_R0δψ_R+ψ_I0δψ_I)`, `δj_i=(ℏ/m)(…)−(q/m)(δA_iρ0+A_i0δρ)`) —
  NOT `_matter_number_current`. `S_η^(A)` confirmed NOT needed for the forward number.
- **D4 — basis-invariant N0/Z0 (the crux).** `Σ_A^cons(ω)=μ0⟨j_eta,G_A^cons j_eta⟩`, `Z0=Re Σ_A^cons(0)`;
  `Σ_A^ret(ω)=Σ_A^cons(ω) − i N0 (a⁵/27c_s⁵)ω⁵ + O(ω⁷)`, `N0=−lim_{ω→0} Im Σ_A^ret/[(a⁵/27c_s⁵)ω⁵]`. I.e.
  **N0 = the ω⁵ retarded-self-energy coefficient = the canonical quadrupole radiation-reaction number** (ties
  to `chi_Q`, the layer-2 "single new retarded number"). PROVEN consistent with V2-09: finite block
  `K_int=[[A,−R],[−R,B]]`, `g=(g_U,g_W)`, `ℓ=(0,1)` ⇒ `Z=Q/Δ`, `N0=(Ω_U²g_W+R g_U)²/Δ0²`. Avoids the U/W
  basis entirely (decision 04: R≠0 ⟹ U/W not eigenmodes).
- **D5 — V2 packet mapping.** Option (b): a minimal **direct-derived-coefficient** path. V2-21 ALREADY
  supports direct `{K,M,Bn,Zn,Nn}`; add a V2-22B `derived_maxwell_transfer` branch carrying `Z0/Z2/Z4,
  N0/N2/N4` + flux normalization + residuals + hashes, generating a V2-21 `direct_coefficients` lane →
  V2-21 computes `D0=K−B0−Z0`, `R_norm=mhat0²S_port N0/D0 − 54Gc_s⁵/(5a⁵c⁵)`. Cleaner than faking a U/W port.
- **D6 — acceptance criteria** (the real tests): MMS on all 5 lanes detecting the dropped vector-Laplacian
  terms; pure-gauge→zero F/flux/deflated; conservative operator self-adjoint under `Z r²` inner product;
  open BC positive outgoing flux + no hard-cap; Green + gauge residuals converge under refinement; low-ω
  `Z_n/N_n` extraction stable under ω-window refinement; **basis-invariance under random basis rotations**;
  **V2-09 finite-block regression reproduces `Q/Δ`, `P²/Δ²`**; term map cites the parent eqs.

## Claude review verdict: ACCEPT as the build plan

Sound, standard-correct multipole electrodynamics; the N0↔ω⁵↔chi_Q connection + the V2-09 consistency proof
make it genuinely consistent with the established framework, not a bolt-on. **Three carried flags:**
1. **ω⁵ NORMALIZATION must be pinned exactly** — design writes `a⁵/27c_s⁵`; compact `σ_Q^can=4a⁵/27c_s⁵`
   (compact §… Ω_Q=3c_s/2a). Same quantity or a factor-4 convention gap? This DIRECTLY scales `R_norm`
   (primary observable). The build MUST pin it against `chi_Q=1` / the V2-09 port coefficient. #1 acceptance
   item.
2. **D5 extends the V2-22B validator** (layer-2 audit infra) — must be ADDITIVE, preserve every existing gate
   + the recursive target-blindness/forbidden-field enforcement, and get its own fidelity check. The
   `mixed_ports` path stays intact.
3. **Design is a PLAN; the build verifies it** — MMS + the transliteration-fidelity audit on the VSH
   reduction are load-bearing.

## Build plan — DECOMPOSE into 3 reviewable chunks (each gets fidelity + adversarial + arbiter)

- **Spike-1 (foundation, highest fidelity risk):** the ℓ=2 VSH Maxwell operator (5 lanes) + gauge
  constraint/deflation + conservative self-adjointness under `Z r²` + MMS (must detect the dropped
  vector-Laplacian terms) + pure-gauge-deflation check. No transfer yet.
- **Spike-2 (the N0 physics):** open-exit DtN BC + the forward `δJ_ψ` chain (reuse 8c phasor) + conservative
  & retarded Green solves + `Z0/Z2/Z4, N0/N2/N4` extraction via the ω⁵ self-energy, with the **#1
  normalization pinning** + basis-invariance + V2-09 regression.
- **Spike-3 (integration → clean R_norm):** the minimal reviewed V2-22B direct-derived extension + emit the
  packet + run the V2 chain → first **CLEAN** `R_norm` (still on the smoke background → a clean structural
  TEST, not physical until M1c).

## Remains open AFTER the spike (do not overclaim)
Self-consistent background selection (M1c); the Path-A return source `S_η^(A)` for wall self-dynamics /
reciprocity; any microscopic derivation of localized Maxwell beyond its parent-sector status. PENDING USER GO
on Spike-1.
