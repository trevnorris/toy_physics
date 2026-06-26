# Moving-Throat PDE — Completion Ladder (what "finish the PDE" actually requires)

**Why this doc exists.** "Complete the moving-throat PDE" is **not one step** — it is a ~6-gate ladder, each
gate a pathA_29-scale effort (directive → review gauntlet → dual-engine → tri-review → user gate). This is the
master checklist so no session forgets the full distance. Maps the roadmap Phases (`..._roadmap.md`) + scaffold
§11 (`..._phase1_linearized_scaffold.md`) + the cross-ℓ gate from the reconciliation note onto a single
progress-trackable list. **Calibrations are allowed (toy model) but every claim must be earned, not asserted;
the goal is a legitimate brane/bulk/soliton structure supporting EM+GR.**

Status legend: ☐ not started · ◐ in progress · ☑ done.

---

## The ladder

### Gate 1 — Frozen-wall D/N unit test  ☑ DONE = `DN_UNITTEST_BC_DEPENDENT`  (`pathA_30`)
- **Tests:** does the distributed geometry lift `Σ=r−R(Ω,w,t)` + `V_conf(Σ)`, frozen (`η=0`) and linearized,
  reduce the longitudinal matter channel to the const-coeff Helmholtz resonator with mouth DtN
  `Z₀₀(ω)=−(ω/c_s)tan(ωL₀/c_s)`, pole ladder `(n+½)πc_s/L₀`?
- **Able-to-fail core:** the *reduction* (wrong lift → variable-coeff operator, renormalized speed, shifted
  domain, wrong BC). Robin-family counterfactual + BC-provenance default to `BC_DEPENDENT`.
- **Why first:** upstream of every ℓ-sector; fails cheapest if the geometry lift is wrong.
- **Difficulty:** low. **Decisive?** no (necessary-not-sufficient plumbing).
- **RESULT (2026-06-25):** the reduction is genuine — const-coeff Helmholtz, bulk `c_s=√(5Kρ⁴/m)`, domain
  `L₀`, DtN `−(ω/c_s)tan`, half-shifted ladder, Robin counterfactual (one symbolic-`α` determinant recovering
  D/N at `α→0`, D/D at `α→∞`), dual-engine agreed (SymPy + independent Mathematica transfer-matrix). The
  Bogoliubov `k⁴` term is honestly retained + deferred under `kξ≪1`. **D/N is a banked calibration input**
  (`bc_derivation_emitted=False`) — earning it → `PASS` is an optional later upgrade. Tri-review CLEAN (arbiter
  re-run both engines + fidelity + adversarial mutation-tested: a wrong answer fails, PASS is reachable,
  engine agreement independent). Artifacts: `software/stage1_solver/{tools/pathA_30_*, reports/pathA_30_*}`,
  directive `directives/pathA_30_dn_unit_test.md`.
- **Known NITs (non-blocking, fold via Codex at next touch):** (1) in `compute_verdict`, `engine_agreement=False`
  short-circuits to `ENV_BLOCKED` *before* the FAIL checks → a genuine engine *disagreement* would mislabel as
  `ENV_BLOCKED` instead of `FAIL_ENGINE_DISAGREE` (conservative — errs toward re-run, never a false PASS).
  (2) `dim_check` uses hand-coded `{M,L,T}` tuples (verified non-tautological); directive §5 now blesses this.

### Gate 2 — Scalar breathing on (`η₀₀`)  ☑ DONE = `BREATHING_CALIBRATED`  (`pathA_31`)
- **Tests:** switch on only the monopole channel; verify the distributed Hellmann–Feynman force reproduces
  the **old `(a,L)` collective closure** in the lowest-mode truncation (scaffold §11.2).
- **Able-to-fail core:** the geometry-*on* coupling must regenerate the already-frozen `(a,L)` dynamics; a
  mismatch means the distributed lift is not the right generalization.
- **Difficulty:** medium. **Decisive?** no (consistency with the existing reduced closure).
- **RESULT (2026-06-25):** the ℓ=0 reduction genuinely reproduces the legacy `(a,L)` Lagrangian — profiles
  `α_a,α_L` derived (`L₀`-harmonic lifting), `M_AB/K_AB` from operator projection (not `∂²E_geom`), the
  distributed Hellmann–Feynman force matches the legacy `F_a,F_L` via two genuinely-independent routes, the
  Hessian pattern matches, dual-engine agreed (max Δ 5.7e-13). The 2-mode truncation is clean
  (generalized-eigenproblem `V₂`-overlap `o_1=0.993` at the physical `β_L0=1.85`) across the whole order-unity
  wall-stiffness window (`K_η/T_w≲2.6`); fails only for sharp walls. **CALIBRATED** because `μ_η,T_w,K_η` are
  calibration inputs. **Caveat:** the overlap measure passes any mouth-BC profile — the anti-wrong-profile
  teeth are the recomputed HF mismatch; clean truncation is shown *for an assumed order-unity wall stiffness*.
  Artifacts: `software/stage1_solver/{directives,tools,reports}/pathA_31_*`.
- **Hard-won (process):** v1 was REJECTED by tri-review (HF `x−x`, hardcoded counterfactual, chosen-to-pass
  truncation threshold); remediation also caught a Parseval-muddled measure + double-`μ_η`; a final surgical
  fix made the STRUCTURE gate computed-not-hardcoded (with an able-to-fail probe). The Gate-1 "fold via Codex
  at next touch" NIT (`engine_agreement` ordering) still applies to the stage1_solver verdict idiom.

### Gate 3 — Grouped-`P2` (ℓ=2) sector  ☐
- **Tests:** switch on `η_2m` one channel at a time; build the grouped response matrix `Z^eff_{AB}(ω)`; test
  the **isotropy gate** `a₂=b₂=0` on the reference branch (scaffold §11.3, handoff §11).
- **Able-to-fail core:** the three grouped lanes `{20,21,22}` must collapse to common coefficients on the
  isotropic branch; anisotropy that won't vanish = the branch is wrong. **This is where the quadrupole sector
  first appears.**
- **Difficulty:** hard. **Decisive?** partial (isotropy is a real gate).

### Gate 4 — Full coupled extraction → quadrupole normalization  ☐  ⭐ THE PAYOFF
- **Tests:** on the isotropic passive/outgoing branch, extract the canonical pair `(K̄₀, K̄₂)`, compute the
  outgoing prefactor `P₀=N₀/D₀` and `γ_quad^eff`, and hit the **single sharp target**
  `m̂₀²·P₀ = 54Gc_s⁵/(5a⁵c⁵)` ⟺ `γ_quad^eff = 2G/(5c⁵)` (handoff §12.3; scaffold §11.4, §8.3).
- **Able-to-fail core:** the branch data `(K_A,M_A,B_{A,n},Z_{A,n},N_{A,n})` of the *actual solution* must
  produce the right `P₀` **form**.
- **CALIBRATION BOUNDARY:** `54/5` is a calibration anchor (Gate-4 = `GENUINE_BLOCKED`). The PDE delivers the
  **form/branch** of `P₀` + passes isotropy; the overall normalization of `G` stays a labeled calibration
  knob. **We are not deriving Newton's constant.**
- **Difficulty:** hardest. **Decisive?** yes (this is the ledger's single OPEN item, "actual branch
  realization").

### Threaded through Gates 2–4 — PN match-back  ☐  (decisive falsification)
- **Tests:** at each gate, the extracted low-frequency data must match the **already-audited PN ladder**
  (2PN/3PN/2.5PN/4PN in `research/4d_*pn*` + `research/1pn_orbital_dynamics`). Roadmap Phase 4. **Do NOT
  re-derive the PN ladder** (memory `project-pn-gravity-ladder`).
- **Why it matters:** this is where the PDE can actually **die** — the cheapest decisive falsifier on the
  whole push.

### Gate 5 — Scalar/dipole side conditions + cross-ℓ unification  ☐
- **Tests:** from the *same* PDE, (a) confirm scalar soft-mode suppression / projection-locking, (b) confirm
  the dipole outgoing branch is cubic + finite-size demoted (roadmap Phase 5), and (c) **unify pathA_29's
  ℓ=0/1 brane↔bulk return residual** (`R0=−M0, R1=−D1`, bounded residual `∝ε₀`) with the ℓ=2 channel — the
  **new cross-ℓ consistency gate** from the reconciliation note.
- **Able-to-fail core:** one `S_return` must *simultaneously* leave the bounded ℓ=0/1 residual **and** deliver
  the ℓ=2 `P₀`. A return that "cleanly cancels at all ℓ" is **suspicious**, not a win
  (`feedback-falsification-is-the-goal`).
- **Difficulty:** hard. **Decisive?** yes (it's the GR-departure prediction's survival).

### EM half — Maxwell mixed channels  ☐  (switches on ~Gates 2–3)
- **Tests:** keep `A_w / J^w / F_{μw}` alive (do **not** zero them as ontology; scaffold §1.3); the localized
  Maxwell sector must ride the same lift and reproduce the brane-effective EM reduction without breaking
  Magnus (magnetism). This is the **EM** half of "supports EM+GR" — not just gravity.
- **Difficulty:** hard. **Decisive?** yes for the EM claim.

### Gate 6 — Full nonlinear branch-selection closure  ☐  ⭐ THE DEEPEST LAYER / likely WALL
- **Tests:** everything in Gates 2–5 is **linearized about a stationary background**. The genuinely
  **nonlinear** branch-selection closure — does the true completed PDE preserve the three exact quotient
  coordinates `(𝔠_tr, 𝔠_nt, ε_η)` (translation-dictionary §15.3)? — sits *underneath* all of it.
- **Roadmap caveat:** the roadmap deliberately says **don't lead** with the full nonlinear theorem; the
  linearized response is what produces the data. This gate is the maximal version and **the most likely place
  we hit the wall** (user's "push until a wall").
- **Difficulty:** maximal. **Decisive?** yes (the ultimate open item).

---

## Two standing caveats on what "complete" means
1. **Linearized first, nonlinear last.** Gates 2–5 are linear-response about a stationary background; the
   nonlinear closure (Gate 6) is the deepest layer. Don't mistake a clean linear gate for the whole PDE.
2. **Form, not magnitude.** The quadrupole "hit" (Gate 4) is the *form/branch* of `m̂₀²P₀` + isotropy;
   normalization (`G`) is a calibrated knob (`GENUINE_BLOCKED`). "Complete" = the structure exists and is
   self-consistent and matches the PN ladder — not a first-principles derivation of `G`.

## Pointers
- Roadmap: `moving_throat_pde_roadmap.md` (Phase 0→6). Scaffold: `moving_throat_pde_phase1_linearized_scaffold.md`
  (§9 unit test, §11 ordered calcs). Engine + targets: `moving_throat_pde_handoff_full.md`.
- Cross-ℓ gate + the two open-item framings: `moving_throat_pde_open_item_reconciliation.md`.
- PN ladder (don't re-derive): `research/4d_*pn*`, `research/1pn_orbital_dynamics`.
- Memory: `project-moving-throat-pde-push` (resume), `project-calibrated-pde-goal`, `project-pn-gravity-ladder`.
