# BUILD DIRECTIVE — Step 5 (CPU): boundary-control characterization of the coupled branch

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 5 ("PML / sponge characterization …
demonstrates boundary reflections sit below the target signal"); brief §5.3 ("Boundary control. Use
absorbing / PML / sponge boundary layers; demonstrate quantitatively that boundary reflections sit
below the target signal"). Pre-reg `docs/stage1_preregistration.md` **§J.3 (non-negotiable):**
"Boundary control — absorbing/PML/sponge layers; quantitatively demonstrate boundary reflections sit
below the target signal" + **§D** (Exit BC = "Open impedance / Robin condition with finite
`R_exit > 0`"; D/N ladder = "AC impedance/reflection, **not** a hard cap"; "Forbidden: hard caps").
**Stack:** `decisions/01_stack_choice.md`. **Builds on:** step 3a (`dcd6a1b`) coupled residual +
open-impedance Robin exits, step 3b (`affc745`) preconditioner, step 4 (`c9b8b2c`) convergence study
+ the target-blind surrogate observables. REUSE all of it.

**Methodology scope (read carefully — this is the load-bearing framing).** The solver currently
realizes the open exits as **Robin impedance BCs** (`matter_exit_impedance_alpha`,
`a0_radial_impedance_alpha`, `a0_exit_impedance_alpha`); there is no sponge/PML layer yet. For the
**stationary isotropic branch there are no propagating waves to reflect**, so the meaningful §J.3
demonstration *at this step* is that the **finite-domain boundary treatment does not contaminate the
interior** where observables live — i.e. the interior solution is insensitive to (a) where the domain
is truncated and (b) the impedance-coefficient values, and a sponge/absorber (if added) only reduces
any residual boundary error. The genuine **outgoing-wave reflection coefficient** (a true PML
reflection number) is a property of the **WP3 linearized tangent** (the wave-like sector), which is
step 8 — the plan explicitly says PML characterization is "on the tangent." So that wave-reflection
characterization is **DEFERRED to step 8**; this step does the stationary-branch boundary-contamination
study and explicitly flags the wave piece as deferred. If you find the stationary study genuinely
cannot be done without the tangent, STOP and flag it rather than inventing a wave problem.

**Contract.** Codex DESIGNS and WRITES the study (the sweeps, the interior/boundary decomposition,
any sponge layer, the metrics, the report). Claude REVIEWS code + output with clean agents + an
independent arbiter re-run and iterates with you until acceptance. This directive states requirements
+ acceptance criteria only; the mechanism is your design call. Run everything LOCALLY ON CPU; no GPU,
no network. Never wrap your own session in a shell `timeout`.

---

## 1. Objective

Quantitatively demonstrate that the finite-domain boundary treatment of the coupled stationary
isotropic branch does **not contaminate the interior** — i.e. boundary-induced error in the interior
(where the target-blind surrogate observables are evaluated) sits **below a stated numerical
threshold** (a self-referential / target-blind reference — see §2-D), across reasonable variation of
the truncation placement and the impedance coefficients. Add a sponge/absorbing layer only if the
Robin impedance exits alone do not achieve this, and characterize its effect. This is the §J.3 gate
for the stationary branch; the outgoing-wave reflection coefficient is deferred to the WP3 tangent
(step 8). Still the **engineering smoke** regime: placeholder params, NO physical packet, export guard
untouched, **target-blind**.

---

## 2. Scope

### IN — mandatory

**A. Truncation-placement sensitivity.** Vary the domain truncation (e.g. the outer radial extent
and/or the axial extent / `R_exit` placement) across **≥3 settings** and measure how much the
**interior** solution and the step-4 target-blind surrogate observables change. The interior must be
defined as a region held away from the boundary zone (so moving the boundary tests contamination of a
fixed interior). Demonstrate the interior change shrinks / stays below threshold as the boundary moves
out — the signature that the truncation is not contaminating the interior.

**B. Impedance-coefficient sensitivity.** Vary the Robin open-impedance coefficients (the matter and
A0 exit/radial alphas) across a defensible range (**≥3 settings**) and show the interior is
insensitive — the signature that the impedance exit is acting transparently, not reflecting energy
back or pinning the interior. (These alphas are placeholder `free_choice`, not targets.)

**C. Sponge/absorber — only if needed (honest either way).** If A/B reveal non-negligible boundary
contamination, add an absorbing/sponge layer near the exit and show it drives the interior boundary
error below threshold; report its parameters. If the Robin impedance exits are already adequate,
**document that explicitly** (a sponge is then not required for the stationary branch) — do NOT add a
sponge for show. State which outcome held and the evidence.

**D. Quantitative boundary-error metric vs a target-blind reference.** Report a measure-consistent
metric of boundary-induced interior error (e.g. interior-restricted L2 change between boundary
settings, and/or an interior-vs-boundary-zone energy/flux decomposition). The "below the signal"
demonstration must use a **target-blind reference scale** — the interior signal magnitude of the
solution itself, and/or the step-4 numerical floor / discretization error — **NOT** any §H target
signal. State the reference and the achieved ratio.

### OUT — do not touch this step

The reduced impedance / DtN extraction and the outgoing-port observable (`chi_Q`, `N_Q`,
`m̂0`, `S_port` — §D outgoing-port convention; that is EXTRACTION/TARGETS, firewalled); the genuine
outgoing-wave reflection coefficient (WP3 tangent, step 8); ANY §H target/observable; the
field→coefficient extraction map; the frozen physical run / `free_choice` value freeze / export-guard
flip (gated); the conservation budget (step 6) and error budget (step 7) as such; GPU / PETSc port.
**Hard caps are FORBIDDEN** (§D) — do not "fix" a boundary by clamping the field to a hard value.

---

## 3. Requirements

**R1 — Reuse, don't fork.** Build on `coupled_branch.py` (residual + Robin exits + R_0(w) geometry),
`newton.py` + the step-3b preconditioner, `convergence.py`'s surrogate observables + interior norms,
`boundaries.py`, `grid.py`, `manifest.py`. Add the sweep / interior-decomposition / sponge as new
functions; do not duplicate the residual, the Newton core, or the preconditioner.

**R2 — Same problem except the boundary under test.** In each sweep, change ONLY the variable being
characterized (truncation placement, OR an impedance coefficient, OR the sponge) — everything else
(placeholder params, gauge, continuation schedule, the other BCs) stays fixed, so the measured
interior change is attributable to the boundary variable. Use a fixed, sufficiently-resolved grid
(informed by step 4) so discretization error does not masquerade as boundary error.

**R3 — Defensible interior metric.** The interior region and the difference norm must be
measure-consistent (conservative r² cell volumes) and identical across settings; the interior must
exclude the boundary/sponge zone. State the interior definition and why it is a fair fixed window as
the boundary moves.

**R4 — No fabricated success (the integrity axis).** Do NOT pick the sweep range, interior window, or
metric that happens to show insensitivity and hide the rest. Report the honest sensitivity even if it
reveals contamination — that is a real finding to diagnose (and may motivate the sponge in C). Do not
loosen Newton/GMRES tolerances. If a setting fails to converge, report it, don't drop it silently.

**R5 — Determinism & manifest.** Deterministic (fixed seed, single-thread, float64). Every solve
emits a manifest (grid, solver controls incl. preconditioner, boundary settings). Emit a
machine-readable boundary-characterization table.

**R6 — TARGET-BLINDNESS (the #1 review axis).** No §H target value anywhere — not in code, config,
params, sweep ranges, the interior metric, thresholds, or the reference scale. The "signal" reference
is target-blind (interior solution magnitude / step-4 floor). Impedance alphas, truncation extents,
and sponge params are conditioning/boundary-driven, never target-driven. Verifiable by grep.

**R7 — NON-FROZEN smoke discipline / guard firewall intact.** Placeholder parameters only, labelled
engineering-smoke / target-blind. No physical packet; no impedance/DtN/chi_Q extraction; do NOT write
under `research/pde_audit/simulation/output/`; do NOT import/modify/satisfy/trip the export guard
scripts; no frozen-packet schema. No hard caps (§D forbidden list).

**R8 — Environment & deps.** Self-contained on the existing stack (torch/numpy/scipy). If you believe
a new dependency is genuinely needed, STOP and FLAG it (reason + tradeoff); do not silently install.

---

## 4. Acceptance criteria

1. A **truncation-placement sweep (≥3 settings)** and an **impedance-coefficient sweep (≥3 settings)**,
   each changing only the variable under test (R2), on a fixed sufficiently-resolved grid.
2. A **quantitative interior boundary-error metric** (R3/D) with a stated **target-blind reference**
   and the achieved ratio — demonstrating interior boundary contamination sits below threshold, OR an
   honest diagnosis that it does not + a sponge (C) that brings it below threshold.
3. An explicit statement of which outcome held (Robin-adequate vs sponge-required) with evidence (R4).
4. Determinism + manifest + machine-readable table (R5).
5. Target-blind (R6, grep-verifiable); guard firewall intact; non-frozen smoke discipline; no hard
   caps (R7); no new dependency (or STOP-flagged, R8).
6. An explicit **deferral note** that the outgoing-wave reflection coefficient is a WP3-tangent
   (step 8) property, not done here, and why.
7. A concise **report** (markdown/YAML, KB-scale) under `reports/` capturing all of the above.

Claude reviews code + report (clean agent — boundary-contamination integrity / no fabricated
insensitivity / target-blindness / firewall / no hard caps) + an independent arbiter re-run
reproducing the boundary-characterization table; we iterate with you until all hold.

---

## 5. Repo hygiene (firewall is live)

- **Code** → `software/stage1_solver/src/`, **tests** → `software/stage1_solver/tests/` (tracked).
- **All run output** → `software/stage1_solver/{runs,data,figures}/` — **gitignored**.
- **The one tracked artifact** → a small report under `software/stage1_solver/reports/`.
- Do **not** commit data; do **not** `git add` anything (orchestrator stages + commits after
  review); **never** `git add -A`; **never** write under `research/pde_audit/simulation/`. If you
  update `README.md`, keep it purely additive (step-5 row + run command) and flag it for review.

---

## 6. Deliverable for review

Report back with: the boundary-characterization design (the two sweeps + ranges, the interior
definition + difference norm, the fixed grid + why it is resolved enough); the truncation-placement
sweep results; the impedance-coefficient sweep results; the quantitative interior boundary-error
metric + target-blind reference + achieved ratio; the Robin-adequate-vs-sponge-required verdict (+
sponge params if added); the deferral note for the wave-reflection coefficient; the report path; and
any STOP-and-flag items. Then Claude reviews before anything is committed.
