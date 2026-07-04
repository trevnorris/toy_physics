# pathA_39 — Magnetism = the moving 4D throat-body interaction: the O(V) charge-coupled force + the scalar-admixture crux

**Status:** ⭐ v1 — §2 RESOLVED (the Claude↔Codex→GLM setup gauntlet is complete; resolution =
`_scratch/pathA_39_s2_resolution_v5.md`, Codex `SOUND_AS_A_RESOLUTION`). **Type:** brane/throat EM sector — the fourth and last
force-sector. **Downstream of** `pathA_38` (electric charge = the static 4D throat-body interaction → `1/r²` Coulomb) and
`pathA_36` (light = the brane shear wave, `c_γ²=μ_R/ρ_br`, and `B_eff=ρ_B0²/χ_c>0`). **Next:** Codex design-review of THIS
directive → dual-engine execution (Stage 0+1 first) → tri-review. **Author:** orchestrator (scaffolding).

> **⚠️ READ THE HONEST FRAME FIRST (§4).** This gate does NOT chase clean isolated Maxwell — that is EXCLUDED up front on the
> provenance-fixed branch (`M_h>0` from pathA_38 forces `h` to propagate; `B_eff>0` from pathA_36). The gate COMPUTES *which
> characterized departure* magnetism is: the charge-coupled magnetic force (conditional on a declared source), the two-magnetism
> separation, and — the crux — whether the charge couples to the density pole (a detectable 5th-force we do NOT observe) or not.
> A characterized departure is a first-class result for a unified toy model, NOT a defeat. Do not rescue, do not oversell.

---

## §0. The crux in one paragraph (reframed)
Magnetism is **NOT** "does the curl of the shear field `u` close Maxwell as a gauge field" (the v2 framing; GLM showed that is a
predetermined fail because pathA_36 already computed the longitudinal stiffness `B_eff>0`). Magnetism **IS the velocity-dependent
(`O(V)`/`O(V₁V₂)`) part of the pathA_38 4D throat-body interaction** — what two *moving* throats do to each other — which is
charge-coupled by construction (the ±w body, not the mass-flow). THE CRUX: (i) does that `O(V)` interaction give a brane-localized
charge-coupled magnetic force with the right sign and falloff; (ii) does it stay separated from **gravitomagnetism** (the weak 3D
flow-swirl of the mass-inflow); and — the load-bearing question — (iii) **does the charge source have nonzero residue on the
longitudinal/density (`c_s`) pole?** A nonzero charge residue = a detectable 5th-force-like scalar we do NOT observe = a real,
falsifiable departure; zero residue = clean coexistence of transverse EM with the gravity/density channel.

## §1. The leverage — the electric-sector route, applied again
- **`pathA_38`** earned the STATIC electric force as the interaction of two throats' 4D bodies, mediated by the brane-localized
  embedding Goldstone `h` (`E=-∇h`, `1/r²`, sign `±w`, ODD parity). **Magnetism = the `O(V)` part of that same interaction.** This
  keeps it charge-coupled and avoids the gauge-closure mold. Banked and IMPORTED: `M_h>0`, `c_E` (from pathA_38's dynamic Green
  `exp(iRω/c_E)/(4πR)`), `q_h`, `f_h(w)`.
- **`pathA_36`** established the transverse shear photons (`c_γ²=μ_R/ρ_br`, 2 pols, `f_u(w)`) AND the real longitudinal density
  sector `B_eff=ρ_B0²/χ_c>0`. Both IMPORTED.
- **Two magnetisms:** gravitomagnetism = 3D flow-swirl of the mass-inflow (EVEN, mass, `c_s`, weak) vs EM magnetism = the 4D
  throat-body swirl (ODD/±w, charge, strong). Consistency checks: moment ∥ spin ∥ angular momentum; EM sense flips with ±w while
  gravitomagnetic sense follows mass (electron vs positron). **This split is settled at the SOURCE level (parity of the declared
  source); the OPERATOR level (Stage 3) is UNPROVEN — for a localized moving throat, `O(V)` operator mixing is generic, so
  `FAIL_GRAVITOMAGNETIC_CONTAMINATION` is the expected outcome until Stage 3 computes otherwise.** The strength hierarchy (EM ≫
  gravitomagnetism, from the swirl concentrating toward the finite mouth radius) is a physical picture AWAITING derivation, not a
  computed result.

## §2. Resolved setup (full detail = `_scratch/pathA_39_s2_resolution_v5.md`)
- **The moving-throat source law is NOT free.** pathA_38 earned the STATIC projection, not the moving-worldtube boundary variation.
  The gate runs on a **parameterized ansatz** `S_s(t,x,w;X,V)=S₀^s(x−X,w)+V_i S_{i,s}^{(1)}(x−X,w)+O(V²)`, projected onto native
  modes — **no Biot–Savart / `V×E` / pre-written `B` inside.** The first-principles `S_i^(1)`, `q_A`, `q_L`, `δÔ_V`, `w`-compactness,
  and hierarchy magnitude are **new-analytic-or-sim-deferred**; the magnetic-force pass is **CONDITIONAL on the declared source**.
- **`M_h>0`, `c_E` are IMPORTED from pathA_38 (not free)** → `h` is a propagating charge-coupled scalar → the exact-Maxwell branch
  is KNOWN UNREACHABLE, and `FAIL_EXTRA_H_BRANON` is EXPECTED.
- **Order-parameter hygiene:** the vector object the moving source needs is `u_T` (pathA_36 shear), NOT a vector `χ_B` — a scalar
  Z₂ wall HOSTS `h`, `u_T`, and the velocity projection, but does NOT explain by symmetry why `q_A/q_h` is Maxwell or why residues
  vanish. Do NOT smuggle a director order parameter through "orientation-lock" language; state whether pathA_38's `(η,n,τ)`
  auxiliaries are physical OP components or wall-normal bookkeeping.

## §3. The staged executable gate (dual-engine, empirically able-to-fail)

**Shared input:** the parameterized `S_s`; `M_h>0, c_E, B_eff, c_γ²=μ_R/ρ_br, q_h, f_h, f_u` imported from pathA_38/36; NO target
readback anywhere.

**Stage 0 — scalar-admixture PARAMETER SCAN (not decisive alone).** Build
`L_scalar = [pathA_36 longitudinal density block, B_eff] + [h block (M_h/2)ḣ²−(K_h/2)k²h², M_h & c_E IMPORTED>0] + [allowed
h–u_L/density mixing] + sources (ρ_q q_h, j_L q_L, mass source)`. Compute `G_scalar(ω,k)`, the charge amplitude `A_qq`, and the
scalar poles `(speed, residue_to_charge, residue_to_mass)`.
- Stage 0 is a **PARAMETER SCAN over `(q_L, mixing)`** (`M_h/c_E` fixed inputs, NOT scanned) → the danger-zone map. Its
  `SCALAR_COEXISTENCE_CLEAN` verdict requires `q_L=0`, which is DERIVED in Stage 1 → if `q_L=0` is set here without a Stage-1
  derivation, flag **`PASS_BY_DECLARATION`** (not a real pass). **Cheapest DECISIVE test = Stage 0+1 together.**
- **FAIL (reachable):** `FAIL_OBSERVABLE_SCALAR_ADMIXTURE`, `FAIL_EXTRA_H_BRANON` (expected, `M_h>0`), `FAIL_CHARGE_COUPLED_CS_SCALAR`.
- **Controls:** inject `q_L=ε` → must fire; closed steady current `∇·j=0` → no longitudinal response (the WIRE limit — not the
  point-throat limit); import `B_eff>0` (don't tune away); `M_h>0,q_h≠0` → extra-scalar flag unless zero residue proven after
  diagonalization; verify Stage-0 `c_E` matches pathA_38's dynamic Green speed (else not using the banked result).
- **Physical prior:** a point-like moving throat has `∇·j≠0` → `q_L≠0` expected → the clean region is likely UNREACHABLE without a
  transverse-projection mechanism (new-analytic-or-sim). The `h`–`u_L` mixing is not forbidden by `w→−w` (h odd; `u_L` in-brane), so
  "clean" needs BOTH `q_L=0` AND zero charge residue on the `h`-dominated eigenpole.

**Stage 1 — moving-source projection.** Project the declared `S_{i,s}^{(1)}`: `q_A^T=⟨f_uT|Π_T S_i⟩`, `q_L=⟨f_uL|Π_L S_i⟩`,
`q_bulk`, `q_even_to_A`, anisotropy. **Pass:** `q_A^T≠0` (charged); flips under `s→−s`; `V=0`→no source; neutral ±composite→zero
monopole EM source; uncharged even drain→zero charge-coupled source; `q_bulk=0` at monopole order (else report tail). **FAIL:**
`FAIL_SWIRL_NOT_SOURCED`, `FAIL_WRONG_AXIS_OR_ANISOTROPIC_SWIRL`, `FAIL_GRAVITOMAGNETIC_CONTAMINATION`, `FAIL_BULK_WAKE_LEAKAGE`,
`FAIL_SOURCE_NOT_COMPACT` (may stay sim-deferred; declares itself conditional).

**Stage 2 — primary magnetic force (one genuine quadratic-action derivation in two languages; no second-method overclaim).**
Integrate out the `u_T` and `u_L` fields sourced by `S_i^(1)` to obtain the `O(V₁V₂)` exchange kernel. **Pass content:**
brane-localized current–current kernel (`~1/R` potential, `1/R²` point-force); `s` sets sign; obtained from the action/source
projection, not inserted. **Corrected Stage-2 landing:** transverse `u_T` exchange gives like-current attraction for `μ_R>0`, and
longitudinal `u_L` exchange has the same attractive sign for `B_eff>0`; the old `q_L²/q_A^T²` cancellation/crossover table is retired.
The scalar-current channel is therefore an unavoidable attractive admixture, not a tunable sign cancellation. **True Stage-2 FAILs:**
`FAIL_WRONG_FALLOFF`, `FAIL_TARGET_READBACK` (perturb the propagator functional form/source → the derived force MUST change; if it
stays Biot–Savart, it read the target back). **`FAIL_NO_LORENTZ_FORCE` is a DIAGNOSTIC, not a fail:** the
medium has a preferred rest frame; the `O(V)` force is `∝V_medium×B_medium`, and a Lorentz-force form requires `c_E=c_γ` (a Stage-4
question). A preferred-frame current–current force with the right sign/falloff is a valid departure.

**Stage 3 — operator parity under motion (compute or explicitly defer).** `δÔ_i=∂Ô/∂V_i|₀`; test `[δÔ_i,P_w]`, `⟨odd|δÔ_i|even⟩`,
`⟨even|δÔ_i|odd⟩`. If it needs the nonlinear throat solve, label the pass `PASS_CONDITIONAL_ON_NO_OPERATOR_PARITY_MIXING`. Do NOT
promote source-level parity to operator-level; `O(V)` mixing is generic → the two-magnetism separation is conditional on this stage.

**Stage 4 — field-coupling classification (diagnostic; replaces v2's exact-Maxwell entry ticket).**
`FIELD_EXACT_MAXWELL_LINEARIZED` (`M_h=0,B_L=0`, first-class generator) — **KNOWN UNREACHABLE (`M_h>0`); listed only to mark it
excluded** / `FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY` (transverse EM + `B_L=B_eff>0`, allowed ONLY if Stage-0 charge residue
vanishes) / `FIELD_SCALAR_VECTOR_DEPARTURE` (the EXPECTED landing: propagating `h`-branon) / `FIELD_CHARGE_COUPLED_DENSITY_DEPARTURE`.
Note: `c_E=c_γ` (needed for a Lorentz/Maxwell form) is the retired `η_T=C_hu²/(ρ_u K_h)=1` condition in disguise; its
kernel-structure question reopens only if anyone tries the exact-Maxwell sub-branch. Maxwell would need BOTH `c_E=c_γ` AND `c_γ=c_s`
(the `λγ` cone-lock) — neither is tested here; both are deferrals.

**Integrity (mandatory):** dual-engine SymPy+Mathematica `ENGINE_AGREE`; dimensional firewall + ablations (omit `V`; omit compact
`ρ_b`; use `c_s` where `c_γ` required; mix mass/charge-current dims); source + parity controls; pass-by-construction ablations
(derive, don't read back); delocalized/noncompact + ghost/wrong-sign controls; `FAIL_TARGET_READBACK`. Every FAIL branch empirically
reachable — force a bad input, the classifier must EMIT the `FAIL_*`, not raise ([[feedback-negative-verdict-short-circuit]]).
**Retired from v2:** mandatory `B_L=0` / `M_h=0` as whole-gate pass conditions; `A_μ=(h,u)` as an assumed dictionary; "Fork S" as an
automatic name for any non-Maxwell coupling.

## §4. Honest scope + expected landing
- **In-scope algebra now:** Stage 0 (scalar block) + Stage 1 (source projection, conditional on the declared `S_i^(1)`) + Stage 2
  (the `O(V₁V₂)` force from the declared source) + Stage 4 (classification), all reusing pathA_38/36/29 machinery. Stage 3
  (operator parity) computed if tractable, else explicitly deferred.
- **Sim-deferred / new-analytic:** the first-principles moving-worldtube source `S_i^(1)` (hence `q_A`, `q_L` values/signs), the
  operator perturbation `δÔ_V`, source `w`-compactness, and the strength-hierarchy magnitude.
- **Expected landing (honest):** exact isolated Maxwell EXCLUDED (`M_h>0`, `B_eff>0`); the likely result is a **characterized
  departure** — a propagating charge-coupled `h`-branon + the scalar-admixture-residue verdict + a source-level (not operator-level)
  two-magnetism split. A charge-coupled `c_s` scalar (Stage-0 fire) would be a detectable 5th-force-like departure. **This is
  first-class for a unified toy model; a clean coexistence would be a stronger-than-expected win. Do not rescue, do not oversell.**
- **This is the FOURTH sector.** Its result does NOT by itself do the consistency knit (`λγ=c_γ/c_s=1` + NG gauntlet + pde_ledger
  assembly), the separate closing step ([[project-brane-existence-defect-structure]]).

## §5. References (read first)
- `_scratch/pathA_39_s2_resolution_v5.md` (the full resolved setup) + the review trail (`pathA_39_s2_convergence.md` v2,
  `..._v3/v4` , `..._glm_review.md`, `..._v4_glm_review.md`).
- `directives/pathA_38_throat_body_electric_localization.md` + `reports/pathA_38_*` (the static charge `E`; `M_h/c_E`, `q_h`, `f_h`).
- `reports/pathA_36_c5_phase_potential.md` + `_results.yaml` (light `c_γ`, `f_u`, `B_eff=ρ_B0²/χ_c`).
- `reports/pathA_29_results.yaml` (localization machinery + the `c_s` residual).
- `docs/conceptual_foundation.md` ⭐ FOUR-SECTOR CHAIN + §3/§4.
- Memories: `[[project-brane-existence-defect-structure]]`, `[[feedback-native-em-mechanisms]]`,
  `[[feedback-negative-verdict-short-circuit]]` (able-to-fail / pass-by-construction), `[[project-calibrated-pde-goal]]`,
  `[[project-mhd-reconnection-parked]]` (the downstream 4D-reconnection payoff).

## §6. Changelog
- **v1 (2026-07-04)** — §2 RESOLVED. Reframed from "does `∇×u` close Maxwell (gauge closure)" to "magnetism = the `O(V)` part of the
  pathA_38 4D throat-body interaction," after a physical-picture conversation (user) + the full Claude↔Codex→GLM setup gauntlet
  (resolution v5, Codex `SOUND_AS_A_RESOLUTION`; GLM `SOUND_WITH_CONCERNS` folded). Staged gate (Stage 0 scalar-admixture screen →
  Stage 1 source projection → Stage 2 force diagnostic → Stage 3 operator parity → Stage 4 field classification). Exact Maxwell
  excluded up front (`M_h>0`); expected landing = a characterized departure. Cheapest first computation = Stage 0+1.
- v0 DRAFT (2026-07-03) — carry-over spec (the gauge-closure framing, now superseded); nailed the crux + leverage + the §2 setup
  sub-questions that the v1 gauntlet resolved.
