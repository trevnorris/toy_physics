# Directive — M1c-run: the PHYSICAL frozen falsification run → physical R_norm

You are the **CODER** (Codex codes / Claude reviews). This is the PHYSICAL run (GATE A is GO, decision-07
structured freeze). Produce the first PHYSICAL `R_norm` on the frozen, target-blind, pre-registered
effective-closure branch. The machinery all exists (M1c-prep `mt15_06` + `m1c_background_export.py` + the V2
chain). **Airtight discipline: FREEZE FIRST (write + hash the freeze sheet), THEN solve; no value mutated
after any residual is seen; do NOT tune anything to `R_norm`.**

## The frozen branch (decision-07, GATE-A GO) — set EXACTLY these, NOT the smoke placeholders
- Geometry: `a = R_mouth = 1`, `L = 37/20 = 1.85`, **`R_exit = 3/2`**, **`R0(w) = 1 + (1/2)(3x²−2x³)`,
  `x = w/L`** (cubic smoothstep a→1.5); boundary class = mouth Dirichlet / open-impedance(Robin) exit / NO
  hard cap; `ε_r=1/20`, `ℓ=a/20` only if a radial support layer is used.
- Wall constitutive (flat): `μ_η = T_w = T_Ω = K_η = 1` (positivity `μ_η>0,T_w>0,K_η+6T_Ω=7≥0` ✓).
- Constants (natural units): `G = c = c_s = 1`, `mhat0 = 1`, `S_port = 1`, `theta_tail = 1` INACTIVE
  (`tail_sector_active = false`). Gauge `H = Z`. EOS quintic `h=(5K/4)ρ⁴` (do not simplify).
- `chi_Q` / `N_Q`: **EXTRACT-and-compare, do NOT freeze `chi_Q=1` as an input** (it's the target). Freeze only
  the relation `N_Q = chi_Q⁻¹`.
- **NOTE the current torch `BranchSmokeConfig` is smoke (`r_mouth=1.2, r_exit=0.9`, …) — DO NOT use it; build
  the frozen M1c config from the freeze sheet.**

## Sequence (strict)
1. **Freeze sheet FIRST.** Emit the GATE-A freeze sheet as canonical JSON (the full frozen set above +
   `parent_action_status=effective_closure`, branch identity, gauge, extraction formulas, solver
   tols/mesh/backend/source-revision, validation protocol) → a **content-only freeze hash** (NO live
   timestamp; NO target residuals / pass-fail fields). This hash names the run.
2. **Torch WP1 solve** with the frozen config (reuse `run_branch_continuation`; pick a resolution that
   converges within `timeout 600` and supports a grid-convergence / Richardson error estimate). Confirm the
   coupled stationary residual is SMALL (self-consistent). Export the background (the M1c-prep exporter path,
   content-only hashing).
3. **Derived chain** (the M1c-prep pipeline `mt15_06` machinery, now on the frozen background): BdG+wall +
   Spike-2 transfer (gauge-active) → the direct-derived V2-22B packet, carrying the derived
   `K,M,B,Z,N` + provenance + `m1c_physical_frozen_run=true`,
   `claim_scope=pre_registered_effective_closure_branch`. **Do NOT touch `physical_export_permitted` / the
   firewalled `research/pde_audit/simulation/`.**
4. **Emit the FROZEN packet to a TRACKED location** `software/stage1_solver/frozen/m1c/<freeze_hash>/`
   (escape the `software/**/runs/` ignore — this is the falsification ARTIFACT; it must be committable +
   byte-reproducible).
5. **V2 chain** → compute the PHYSICAL `R_norm = mhat0²·S_port·N0/D0 − 54Gc_s⁵/(5a⁵c⁵)`, `D0=K−B0−Z0`, plus
   `R_pole, P2, P4`, and EXTRACT `chi_Q` then compare to 1.
6. **§J error budget.** Compose the numerical floors from the chain (solver/Newton, discretization/Richardson,
   boundary, conservation) + the **M1c-prep interp-method sensitivity** (N0 ~1.6% from interp order) + the
   **background-resolution** error at this resolution. State the budget LOUDLY scoped: it is NUMERICAL only —
   it EXCLUDES the free_choice/posit uncertainty (the result is CONDITIONAL on the frozen target-blind posits,
   NOT a numerical error) and model-form. Be honest if the resolution is modest (a tighter run is GPU/cluster
   later).

## Acceptance criteria
- Freeze sheet written + hashed BEFORE the solve; no value mutated after any residual.
- Torch WP1 residual small (self-consistent) + a resolution/convergence estimate for §J.
- Frozen packet validates through the V2 chain (open/stability/forbidden-scan), in the TRACKED `frozen/` dir,
  byte-reproducible (two reruns → identical freeze hash + packet hash + `R_norm`).
- `R_norm` (+ `R_pole, P2, P4, chi_Q`) computed and reported HONESTLY whatever they are — pass = within the §J
  budget of the GR target; miss = honest falsification of THIS branch+posits. NO tuning.
- `chi_Q` extracted+compared (not frozen as input). Loud CONDITIONAL-on-posits framing; loud resolution scope.

## Working rules
≤2 `math`/`wolframscript` kernels; `timeout 600` per script RUN (incl. the torch solve — pick a resolution
that fits); never wrap your session in a shell `timeout`. Code under `src/stage1_solver/` + `mathematica/`;
the FROZEN artifact under tracked `software/stage1_solver/frozen/m1c/<hash>/`; transient run output →
gitignored `runs/`. No `git add`/commit. No network/GPU. Do not touch `physical_export_permitted` or
`research/pde_audit/simulation/`.

## Report back
The freeze hash (written before the solve); the torch WP1 residual + resolution/convergence; the derived
coefficients; the PHYSICAL `R_norm`/`R_pole`/`P2`/`P4`/`chi_Q` (honest, whatever they are) with the §J budget
and the conditional-on-posits + resolution scope; the tracked frozen-packet path + byte-reproducibility
(two-rerun identical hashes). Then Claude reviews: fidelity (frozen config matches decision-07; chain
faithful) + adversarial (freeze-before-solve honored, target-blind/not-tuned, conditional framing honest,
budget honest, no firewall/guard touched, no can't-fail) + arbiter (re-run → R_norm + hashes reproduce).
