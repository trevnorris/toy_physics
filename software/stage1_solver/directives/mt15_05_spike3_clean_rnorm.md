# Directive — Spike-3: V2-22B direct-derived-coefficient extension → first CLEAN R_norm

You are the **CODER** (Codex codes / Claude reviews). Final spike of 3 (decision 05 D5): wire the DERIVED
mixed-Maxwell transfer (Spike-2) + the DERIVED wall/BdG sectors (M1b) into the V2 residual chain via a
minimal direct-derived-coefficient path, producing the first **CLEAN `R_norm`** (no posited ports in the
`N0` path). Pre-freeze, target-blind; the background is still M1b smoke → the `R_norm` NUMBER is MACHINERY
(a clean test STRUCTURE, not a physical result; M1c required). NO GATE A, NO frozen run, NO `S_η^(A)`.

Follow decision 05 D5: `software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md` (D5) +
`_scratch/mixed_maxwell_spike_design.log` (D5). Inputs: M1b derived `K, M, B0/B2/B4`
(`mathematica/mt15_02_bdg_wall_derivation.wls`) and Spike-2 derived `Z0/Z2/Z4, N0/N2/N4`
(`mathematica/mt15_04_spike2_transfer_n0.wls`). Use the DERIVED values (re-run/import — do NOT re-posit).

## Part 1 — minimal, ADDITIVE V2-22B direct-derived extension (⚠ layer-2 audit infra)
The V2 chain is under `research/pde_audit/scripts/` (the validator/extractor — NOT the firewalled
`research/pde_audit/simulation/`, which you must not touch). V2-21 already supports direct coefficients
`{K,M,Bn,Zn,Nn}`; V2-22B currently REQUIRES posited `mixed_ports` (+ ≥1 A_w port). Extend V2-22B to ALSO
accept a **direct-derived-coefficient** packet:
- Add an OPTIONAL schema branch (e.g. `derived_maxwell_transfer` carrying `Z0/Z2/Z4, N0/N2/N4` + provenance:
  `status=derived_green_function_transfer`, gauge convention, `flux_normalization` (Γ_port), operator/gauge
  residuals, source hashes). When present, the validator uses it to generate a V2-21 `direct_coefficients`
  lane and **relaxes the mixed_ports requirement** for that path.
- **STRICTLY ADDITIVE:** the existing posited-`mixed_ports` path must keep working UNCHANGED. Every existing
  gate stays: open-impedance/`R_exit>0`/no-hard-cap, `K>0/M>0`, monotone grid, profile norms, and the
  **recursive forbidden-target-field scan** (`R_pole/R_norm/P0_target/gamma_*/...` still rejected — the new
  `Z0/N0` are INPUTS, not forbidden outputs; do NOT weaken the scan).
- **REGRESSION (mandatory):** after the change, the existing committed fixtures still behave —
  `…/fixtures/stage_v2_22b_sample_solver_output_valid.json` validates, the hard-cap control is rejected — and
  the prior M1a/M1b/Spike packets on the mixed_ports path still validate. Run these and report.
- If the V2 scripts have a derivation note / tests, update them minimally to cover the new path (additive).

## Part 2 — Mathematica producer → CLEAN R_norm
Emit a V2-22B **direct-derived** packet integrating the DERIVED sectors: `K, M` + `B0/B2/B4` (M1b) and
`Z0/Z2/Z4, N0/N2/N4` (Spike-2), plus constants `G, c_s, c, a, mhat0, S_port` (target-blind placeholders,
NOT tuned). Run the real V2-22C/V2-21 chain on it → compute `R_norm = mhat0²·S_port·N0/D0 −
54Gc_s⁵/(5a⁵c⁵)`, `D0 = K − B0 − Z0`. Report `R_norm` (and `R_pole, P2, P4`) with the LOUD label: **CLEAN
structure (no posited ports in the N0 path), but a MACHINERY number on the smoke background — NOT a physical
falsification result until M1c.** Do NOT tune anything to make `R_norm` small; whatever it is, report it.

## Acceptance criteria
- The V2-22B extension is additive; existing fixtures + prior packets still pass; the forbidden-field scan
  still rejects a target-leak (teeth-test it on the new path too: a direct-derived packet with an injected
  `R_norm` field must be rejected).
- The direct-derived packet validates and flows through V2-21 → a finite `R_norm` (+ R_pole/P2/P4).
- The `Z0/N0` used are the Spike-2 DERIVED values and `K/M/B` the M1b DERIVED values (provenance recorded);
  nothing re-posited in the N0 path.
- Loud machinery label; no target-tuning; target-blind constants.
- Term/provenance map: which derived inputs came from M1b vs Spike-2; the extension's additive nature.

## Scope OUT
M1c / self-consistent background; GATE A; the frozen physical run; `S_η^(A)`/Path A; any change to the
firewalled `research/pde_audit/simulation/`; reporting `R_norm` as physical.

## Working rules
≤2 `math`/`wolframscript` kernels; `timeout 600` per script RUN; iterate to exit 0; never wrap your session in
a shell `timeout`. WL + producer files under `software/stage1_solver/mathematica/` (e.g.
`mt15_05_spike3_clean_rnorm.wls` + report); V2 extension edits under `research/pde_audit/scripts/` (+ its
note/tests). Run artifacts → gitignored `runs/`. No `git add`/commit. No network/GPU. Target-blind.

## Report back
The V2-22B extension (what was added, that it's additive, the fixture/prior-packet regression results, the
forbidden-field teeth-test on the new path); the direct-derived packet (derived inputs + provenance); the
computed `R_norm`/`R_pole`/`P2`/`P4` with the loud machinery label; the term/provenance map. Then Claude
reviews: fidelity (the extension preserves the chain; the packet uses derived values) + adversarial (additive?
gates intact? forbidden-scan intact? R_norm not tuned? machinery label honest?) + arbiter (re-run, fixtures
still pass, R_norm reproduces).
