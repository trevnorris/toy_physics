# Branch-Realization Execution Plan

**Companion to:** [`docs/branch_realization_brief.md`](branch_realization_brief.md) (currently at `notes/branch_realization_brief.md`; user plans to move it to `docs/`).

**Status:** Forward plan, written 2026-05-27, based on a working-conversation review of the brief, prior-attempt lessons, and the current PDE-ledger red-team. Pre-execution: the parent operator is not yet frozen (waiting on the active red-team to land its remaining batches). This doc fixes the *shape* of the work so a future session can pick it up cold.

**Read this if you are:** a future session that needs to know how Stage 1 will actually be run, what tooling to use, why the framing is what it is, and what *not* to repeat from prior failed attempts.

---

## 1. What this plan does and does not cover

This plan is the engineering-side companion to the scientific brief. The brief (see §0 above) specifies *what* must be solved and *what* must be reported. This plan specifies *how*: the discretization shape, the hardware tier, the validation sequence, and the relationship to other parts of the program.

This plan does **not** alter, override, or weaken the brief. In particular, every §5 validation requirement in the brief is mandatory regardless of what's in this plan, and the pre-registration rule (§3.4 of the brief) is non-negotiable.

This plan assumes the user's stated project goal (see [[project-analog-framework-goal]] in auto-memory): the deliverable is a working *mathematical analog* between EM and gravity, not an ontology claim. Stage 1 scoring is "the analog operationally works as a bridge," not "evidence the model is true."

---

## 2. Computational shape of the problem

### 2.1 Why the working domain is 2D

The brief asks for two solves:

1. **WP1 — Stationary isotropic moving-throat branch** (see `notes/5pn/5pn_stage354_355_computational_handoff.md` §3 and `notes/moving_throat_pde_engine_handoff.md` §9–11).
2. **WP3 — Grouped real P2 weak-axisymmetric tangent** linearized around that branch.

The parent system is (4+1)D — fields `ψ(X,t)`, `A_M(X,t)`, `R(Ω, w, t)`. Two symmetry facts collapse the working domain to 2D:

- **Stationary** drops `t`.
- **Isotropic** (spherical symmetry in `(x,y,z)`) drops `Ω`, collapsing the three brane spatial dimensions to a single radial coordinate `r`.

So on the isotropic branch the fields are functions of `(r, w)` only. The P2 tangent has angular structure fixed by `Y_{2m}`, so it too lives on `(r, w)` as a linear PDE around the WP1 solution.

This is not an approximation. The 2D mesh **parameterizes 4D field configurations that happen to have spherical symmetry** — the angular content is identically present, just constant along the angular directions because of the symmetry. The same maneuver is what reduces the hydrogen `1s` to a 1D radial ODE, FRW cosmology to a 1D `a(t)` ODE, axisymmetric CFD jets to 2D `(r, z)` meshes, etc. The 4D measure factor `r²` appears explicitly in the discrete operators; the discretization is 3D physics on a 2D mesh, not 2D physics.

### 2.2 What 2D does and does not capture

| Captures | Does not capture |
|---|---|
| Full radial and axial flux of fluid through and along the throat | Angular flux (identically zero by symmetry on the isotropic branch) |
| All coefficients `D0, D2, D4, N0, N2, N4, P0, P2, P4` defined on the isotropic branch | Whether the isotropic branch is dynamically stable against general angular perturbations |
| The grouped P2 response operator's low-frequency expansion (`R_pole`, `R_norm`) | Two-throat / multi-body interactions |
| The "induced P_{22} response" mechanism from `notes/lepton_work.md` §10 (it is precisely the linearized P2 tangent) | The full spin-1/2 chain (needs §14 autonomous-eigenmode closure separately) |
| Outgoing DtN data via Robin impedance / PML on the linear tangent | Long-time dynamic radiation propagation |

The "does not capture" column is intentionally *not* in the brief's scope. Those are downstream tests gated on Stage 1 passing.

### 2.3 Field count and grid scale

Unknowns on the `(r, w)` mesh (post-gauge-fixing):

- `ψ_real, ψ_imag` (or amplitude/phase — Madelung form)
- `A_0, A_r, A_w` (Maxwell components in the axisymmetric reduction)
- `R_0(w)` (throat shape — 1D auxiliary)
- One gauge-fixing scalar or Lagrange multiplier

Roughly 5–7 coupled fields per cell. Estimated grid scales for Stage 1:

| Tier | Cells `(N_r × N_w)` | Total DOFs | Hardware | Wall clock per Newton solve |
|---|---|---|---|---|
| Smoke | 512 × 256 | ~0.8M | laptop CPU | minutes |
| Production | 2048 × 1024 | ~12M | A100/A6000 class, 64–128 GB RAM | hours |
| Refined | 8192 × 4096 | ~200M | H100 or dual-A100, 256 GB RAM | ~1 day |

The §5.2 convergence study needs 3–5 levels; the WP3 tangent on the same mesh is 2–5× cheaper than WP1 because it is linear. Total compute for a complete Stage-1 attempt with validation: **1–4 weeks of intermittent runpods GPU time**, ballpark **$500–$5000** of cloud GPU spend depending on iteration count.

---

## 3. Hardware and tooling plan

### 3.1 Hardware tier

**Sweet spot:** runpods.io API-triggered runs on single high-end GPU instances (A100 80 GB, A6000, or H100), spun up per resolution level and torn down after artifact extraction. The user already plans to drive runs remotely from this working repo via the runpods API.

**Not required:** multi-GPU, multi-node, supercomputer allocation, GPU farm. The Newton-Krylov + tangent structure is single-GPU-friendly. Multi-GPU only becomes relevant if the program later moves to 3D (isotropy stability test) or 4+1D dynamic (multi-body or radiation tests).

### 3.2 Software stack (to be chosen at solver-build time, not pre-committed)

Candidates in order of suitability:

- **Dedalus** (spectral, native support for nonlinear PDE + Newton-Krylov + continuation; widely used in astrophysical-fluid community). Strong default.
- **Firedrake** or **FEniCS** (FEM, axisymmetric reduction is straightforward; less idiomatic for spectral GPE work).
- **Custom JAX or PyTorch solver** with finite differences + autodiff Jacobian for JFNK. Best fit for GPU; worst fit for "off-the-shelf preconditioner."
- **PETSc/Trilinos** if going CPU-cluster with serious multigrid preconditioning.

The deciding question at solver-build time is preconditioner availability for the stiff `n=5` polytropic Jacobian (see §6 below). Whichever stack gives the cleanest path to multigrid on the linearized operator wins.

### 3.3 Run-orchestration pattern

Spin up GPU instance → run solver with frozen config → checkpoint artifacts back to working repo → tear down. No long-running stateful instances. Each Newton solve produces:

- the converged stationary branch (or a clean non-convergence diagnostic),
- the tangent response packet,
- the §5 validation report (convergence orders, conservation drifts, noise floor),
- a frozen-hash manifest matching `NONLINEAR_PROTOCOL_V2.md` §"Freeze Boundary".

---

## 4. Why this is the right shape of problem (vs. the prior failed attempt)

The prior project at `/var/projects/4d_1pn_sim` was a long-time-evolved, two-defect, 4+1D PDE simulation aiming to extract 1PN orbital precession (`β_eff`) from time-series fits. It hit five structural walls (recorded in this conversation):

1. **Gravity behaves differently in 2D vs 3D** → forced full 3+w dimensionality.
2. **Two bodies must both move freely for perihelion** → couldn't lock one to save compute.
3. **Mercury-vs-Sun astronomical scale ratio blew out float64** → measurements lost to error.
4. **Box couldn't be made big enough to hold orbits without wave reflections** → boundary noise dominated.
5. **Vortex inflow depleted box, no clean refill scheme** → mass conservation drifted over long runs.

The branch-realization brief is structurally a *different shape* of problem:

| Prior attempt | Branch realization |
|---|---|
| Two defects | One throat |
| Full 3+w spatial mesh | (r, w) 2D from symmetry |
| Time-evolve thousands of steps | Newton solve, no time stepping |
| Extract `β_eff` from orbit fit | Extract scalar coefficients from static + linearized solve |
| Cost `O(N⁴ × N_t)` | Cost `O(N² × Newton iters × continuation steps)` |
| No noise-floor gating | §5 validation suite is the deliverable |

Each of the five issues maps as follows:

| Issue | Status in branch realization |
|---|---|
| Gravity 2D vs 3D | 2D mesh discretizes **3D physics with axial symmetry**, full `r²` measure factor. Not "2D gravity." |
| Two-body wobble | Single defect — does not arise. |
| Float64 dynamic range | All targets are dimensionless ratios; pick natural units (`a=c_s=c=G=1`); everything `O(1)–O(10²)`. |
| Box reflections | Stationary main solve has no traveling waves. PML on the linear tangent is single-frequency, textbook. |
| Vortex inflow / mass refill | Stationary solve has `∂_t ρ = 0`; no time-progressive depletion. For dynamic validation runs, grand-canonical / chemical-potential boundaries are standard cold-atom-GPE practice. |

A future session reading this should not interpret "we're going simpler" as "we're cutting corners." The brief asks an easier *and still decisive* question than the prior project asked; the easier shape is precisely what makes it tractable. The prior project was a harder problem, not a more rigorous one.

---

## 5. Relationship to spin-1/2 / lepton_work / autonomous-eigenmode closure

The brief's framing intersects `notes/lepton_work.md` in a non-obvious way that future sessions need to understand.

`lepton_work.md` §10 documents that the spin-1/2 mechanism rides on an **induced area-preserving `P_{22}` deformation** of the throat mouth, and §14 documents that the same-charge half-integer quantizer is conditional on the isolated particle being an **autonomous self-reproducing standing-wave GNLS soliton**.

The relationship to Stage 1:

- The **stationary isotropic branch** (WP1) **is** the autonomous-eigenmode candidate. A Newton solve that converges to a stationary isotropic configuration with stable support modes is implicitly testing whether the parent PDE admits such an autonomous soliton at all. If WP1 has no solution, the §14 closure is unreachable in this model.
- The **grouped P2 tangent** (WP3) **is** the induced-deformation response. The P_{22} mechanism in `lepton_work.md` §10 is a special case of the broader `Y_{2m}` quadrupole linearized response that the brief tests. The coefficients the brief extracts (`D_2, D_4, N_2, N_4`) characterize the response operator whose poles include the resonance the user described ("bumping a particle causes it to resonate").
- **Stage 1 does not directly test spin-1/2.** That's downstream of `lepton_work.md` §14 + §12 (central-sign holonomy) + the recirculation/plumbing closure. Stage 1 only tests whether the foundational autonomous-eigenmode and P2-response structures exist with the right gravity-side coefficients.

The user explicitly flagged: *"never assume the simple solution works"* in this system. The brief honors that — its targets are coefficients of the **response operator's low-frequency expansion**, not naïve properties of the static configuration. The feedback structure (inflow pressure ↔ wall energy ↔ standing-wave support ↔ brane closure) is *what determines those coefficients* on the realized branch. A successful Stage 1 means the feedback closes self-consistently and produces the right numbers; a failure tells you which sector of the feedback didn't close.

---

## 6. The real failure mode to watch for

Not in the five prior-project issues. Specific to this problem:

**Stiff `P = K ρ^5` polytropic conditioning.** The enthalpy `h(ρ) = (5K/4)ρ^4` becomes near-singular as `ρ → 0` inside the throat or near the wall. This makes the Newton Jacobian poorly conditioned (condition number ~ 10⁶ regime, not catastrophic but unfriendly without a good preconditioner). Mitigations:

- **Geometric or algebraic multigrid** on the linearized operator (this is why the software stack choice matters).
- **Change of variables** — e.g., work in `log ρ` or in a regularized density, so the near-singular region maps to a well-behaved one.
- **Continuation in `K`** — start with a milder EOS, continue to the stiff target value.

Poorly conditioned Newton solves fail *loudly* (residual stops dropping, GMRES iteration count blows up), which is the opposite of how the prior project's issues manifested. That's a feature: the diagnostic is unambiguous.

---

## 7. Sequencing

Implementation order, with the answer each step gives:

1. **GPE-soliton or GPE-vortex benchmark on a 2D grid.** Cheap (single-GPU hour or two). Proves discretization + boundary treatment work on a known-analytic problem. *Brief §5.1 requirement.*
2. **Manufactured-solution tests per operator block** (matter, Maxwell, wall). Cheap. Proves discrete operators match continuum operators.
3. **Stationary isotropic branch Newton solve at coarse grid.** Few GPU-hours. End-to-end pipeline smoke test producing a smoke `D0, P0, …` packet.
4. **Convergence study at 3–5 grid levels.** Order-of-day GPU work per level. Establishes noise floor and convergence order for each observable. *Brief §5.2.*
5. **PML / sponge characterization on the tangent.** Demonstrates boundary reflections sit below the target signal. *Brief §5.3.*
6. **Conservation diagnostics.** Mass / charge / energy drift over each run. *Brief §5.4.*
7. **Error budget statement.** Quantitative uncertainty per extracted observable. *Brief §5.5.*
8. **WP3 — P2 tangent on the converged stationary branch.** Day or two. Produces `d ln R_tr, d ln R_target, d ln ε_η`.
9. **Stage-1 verdict.** Compare extracted observables to brief §4 targets within validated error. Honest either way.

If Stage 1 passes, **Stage 2** (divergence extraction, brief §8) becomes worth running. If it fails trustworthy, the result is publishable and diagnoses *which* sector of the parent theory needs revisiting.

---

## 8. Prerequisites

These must be true before Stage 1 work begins:

1. **Parent operator frozen.** `research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md` records that `S_eta`/`S_Sigma` promotion and the full coupled physical residual equations are not yet frozen. The active red-team is part of what will land them.
2. **PDE-ledger red-team complete.** Currently 60/253 stages verified (see [[project-moving-throat-verification]] in auto-memory). Stage 1 work should not start until the red-team's first end-to-end pass is done; the user has also planned a second cross-check pass before downstream computational work (see [[project-full-second-pass]]).
3. **A pre-registration document.** Before any compute is spent, freeze: which branch, which boundary protocol, which gauge convention, which outgoing-port convention, which observable is primary. Per brief §3.4 — once written, do not edit after seeing data.

---

## 9. What success means in the analog-framework frame

Stage 1 success **is not** "evidence the model describes physical reality." (User has not made and does not intend to make that claim — see [[project-analog-framework-goal]].)

Stage 1 success **is** "the analog framework is self-consistent and operationally demonstrable: one set of equations, reduced two ways, returning the right coefficients in both reductions, with a validated noise floor." That delivers the project's stated deliverable — a working mathematical bridge between EM and gravity formalisms — at the strongest level a single solo-researcher-plus-AI program can deliver.

A trustworthy Stage 1 miss is also a valid outcome. It says: the analog as currently written does not close all the way through; here is which coefficient does not match; here is the validated error budget around that statement. That is publishable, citable, and useful — particularly because the prior 4d_1pn_sim attempt could not produce a trustworthy verdict either direction.

Stage 2 (novel-prediction extraction) is the test of whether the analog also *says something new and right*. Stage 1 is the gate to running Stage 2 honestly.

---

## 10. References

**Primary scope (read first):**
- `docs/branch_realization_brief.md` — the scientific brief itself (currently `notes/branch_realization_brief.md`)
- `notes/5pn/5pn_stage354_355_computational_handoff.md` — work-package decomposition (WP0–WP5)
- `notes/moving_throat_pde_engine_handoff.md` — exact parent action, modes, conservative bundle formulas

**Compact reference:**
- `notes/moving_throat_pde_program_compact.md` §7A.6 (solver pipeline card), §7A.7 (V2-23 prior reduced run), §13.5 (solver card), §12, §13 (full coefficient and target definitions)

**Spin / lepton chain (relationship to WP1 / WP3):**
- `notes/lepton_work.md` §10 (area-preserving P_{22}), §14 (autonomous-eigenmode closure)
- `notes/sponge_mode_leakage_exploration.md` — exploratory note on `S_leak` and distributed charge

**Prior-attempt artifacts (do not reuse; reference for context only):**
- `/var/projects/4d_1pn_sim/` — failed 4+1D two-body orbital-precession attempt; `docs/results_journal.md` documents the brick wall
- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — reduced 1D-FEM surrogate, NOT a real GPE+Maxwell solve; useful only as a fixture
- `research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md` — records the parent-operator freezing requirement
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — frozen-input checklist (gauge, boundary, stability gates, source hashes)

**Red-team state:**
- `notes/STAGE_VERIFICATION_COVERAGE.md` (and other tracker docs under `research/pde_ledger/notes/`)
- See [[project-moving-throat-verification]] in auto-memory for current batch state

**Methodology context:**
- `docs/methodology_paper_outline.md` (companion doc)
- `docs/redteam_thoroughness.md`
- `docs/adversarial_audit_directive.md`
