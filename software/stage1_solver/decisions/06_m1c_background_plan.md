# Decision 06 — M1c self-consistent background + frozen physical run: design (REVIEWED + ACCEPTED) + plan

**Date:** 2026-06-17
**Mechanism:** Claude+Codex design consult (read-only) — `_scratch/m1c_background_design_prompt.md` /
`…_design.log` (M1-M6). User chose "design the M1c background first." Claude reviewed + ACCEPTS, with flags.
**Goal:** turn the spike's CLEAN-but-machinery `R_norm` into a PHYSICAL falsification number: a self-consistent
WP1 background (not the smoke ρ0~243) + GATE A (freeze free_choice target-blind) + a frozen, byte-reproducible,
COMMITTED packet → physical `R_norm`.

## Accepted design (condensed)

- **M1 — background:** produce the self-consistent WP1 bundle (`ψ_R0,ψ_I0,A_00,A_r0,A_w0,ρ0,R0(w)`) on the
  `(r,w)` FV grid satisfying the coupled stationary residual (small, not ~243); first acceptance gate =
  imported background residual small + grid-convergence evidence.
- **M2 — engine choice: REUSE the validated torch WP1 (Recommendation A).** The torch coupled solver is
  already validated (MMS order~2, Newton residual~1e-8, JVP 9.1e-12; §7 steps 1-4); only its PARAMS are smoke.
  Cross-engine interface: torch exports a canonical JSON background (grid/fields/ρ0/R0/Z/config/residuals/hash)
  → MMA imports + interpolates/projects, replacing the hardcoded smoke formulas in mt15_02 + mt15_04; torch
  replays the residual, MMA runs consistency checks. (B fresh MMA-FEM = rebuild a large validated solver, NO;
  C dual-engine background = reserve as an optional adversarial cross-check.) **NOTE: the torch WP1 must be
  RE-SOLVED with the GATE-A frozen values — "reuse" = reuse the validated MACHINERY, not the smoke solution.**
- **M3 — derived-operator consistency:** replace M1b's + Spike-2's separate hardcoded smoke grids/formulas with
  ONE shared imported background; conservative projection for ρ0/weights, point interp for smooth fields,
  record interp order + one refinement.
- **M4 — GATE A + sequence:** freeze (target-blind, before any residual) the wall free_choice
  (`μ_η,T_w,T_Ω,K_η`), geometry (`a,L,R_mouth,R_exit,R0(w)` family, boundary class; `L/a=37/20` is free_choice
  not a target), source/extraction constants (`G,c_s,c,a,mhat0,S_port,theta_tail`), + solver
  consts/tols/mesh/source-rev → all into the freeze hash. Strict sequence: **GATE A freeze → torch WP1 solve →
  export/import residual checks → M1b BdG/wall on imported bg → Spike-2 transfer → Spike-3 packet → V2-22C →
  R_norm/error budget**; no value touched after any residual is seen. **GATE B: do NOT mutate the torch
  `physical_export_permitted` guard (firewalled simulation/); the MMA packet instead self-declares
  `parent_action_status=effective_closure`, `m1c_physical_frozen_run=true`,
  `claim_scope=pre_registered_effective_closure_branch`** (Path B = scoped Stage-1 test; Path A separate).
- **M5 — frozen-packet discipline:** byte-reproducible (content-only hashing, REMOVE the `DateString[Now]`
  the WL scripts currently use) + COMMITTED to a TRACKED path `software/stage1_solver/frozen/m1c/<freeze_hash>/`
  (escape the `software/**/runs/` ignore); acceptance = two reruns → identical packet/observable/freeze hashes.
- **M6 — open after M1c:** Path A / `S_η^(A)` (wall self-dynamics); R_pole/P2/P4/chi_Q/WP3 secondary
  (not re-ranked after R_norm).

## Claude review verdict: ACCEPT. Carried flags:
1. **#1 RISK — nonzero-spatial-gauge revalidation.** The smoke background had `A_r0=A_w0=0`, so the phasor
   current's gauge terms (`−(q/m)(δA_iρ0 + A_i0 δρ)`) + the VSH-Maxwell gauge couplings were DORMANT. On the
   real self-consistent background `A_i0≠0` → those terms ACTIVATE for the first time. M1c MUST re-run the
   Spike-1/2 fidelity + operator checks (current Fréchet, pure-gauge, basis-invariance, V2-09, Green residuals,
   N0>0) WITH `A≠0`, watching for spatial-gauge SIGN errors. This is the analog of the spike's #1 flag.
2. **Torch re-solve, not reuse the smoke solution** (see M2 note) — emphasize in the build.
3. **Resolution / error budget:** the torch WP1 resolution (CPU direct-LU caps ~100k DOF) sets the §J error
   budget; a first M1c at modest resolution with an honest budget is fine — the point is residual-small
   self-consistency, not maximal resolution.
4. **Large-artifact management:** the exported background fields can be large; commit the PACKET + a
   hash-referenced/modest-resolution background, not necessarily full field arrays at high res.

## Plan — DECOMPOSE (mirrors the established build-machinery-then-freeze discipline)
- **M1c-prep (target-blind, NO GATE A):** build the torch background export + MMA import/projection + swap the
  shared background into M1b/Spike-2 + the **A≠0 revalidation** of the derived chain — all on smoke/placeholder
  values, proving the cross-engine pipeline + that the gauge terms are faithful when active. NO frozen values,
  NO physical R_norm.
- **🚪 GATE A:** freeze the free_choice values target-blind. **Split of roles (user-clarified 2026-06-17):**
  the user's GATE-A role is the **GO authorization for the physical run + the target-blind guarantee**, NOT
  choosing the physics values. The actual free_choice **VALUES** (what the prereg §E / ANSATZ_LEDGER already
  pins vs what must be set, and the functional forms) are a **MATH/documentation determination → resolved by a
  Claude+Codex consult + agreement** ([[claude_codex_resolve_math]]; value-classification is math, not a user
  decision). Codex + Claude agree + record it (a decision record); escalate to the user ONLY if it changes the
  conceptual nature. Values frozen WITHOUT reference to making R_norm come out 0.
- **M1c-run (physical):** torch re-solve with frozen values → export → MMA derived chain → frozen packet
  (tracked frozen/ dir, byte-reproducible) → V2 chain → PHYSICAL `R_norm` + §J error budget. Pass = target-blind
  GR-quadrupole match; miss = honest falsification of this pre-registered effective-closure branch.

PENDING USER decision on how to proceed (M1c-prep build vs GATE-A-values discussion first) + the GATE-A freeze.
