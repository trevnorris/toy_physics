# BUILD DIRECTIVE — Step 2: manufactured-solution tests per operator block

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 2 ("Manufactured-solution
tests per operator block (matter, Maxwell, wall) — proves discrete operators match continuum
operators"); pre-registration `docs/stage1_preregistration.md` §J validation gates. **Stack:**
`software/stage1_solver/decisions/01_stack_choice.md` (torch, conservative FD/FV, float64).
**Builds on:** step 1 (`directives/step1_gpe_benchmark.md`, committed `240fa40`) — the solver
skeleton, conservative-FV radial/tensor operators, boundary framework, Newton/JVP core, manifest
and determinism discipline are already in place and reviewed clean. This step certifies the
**discrete operators against their continuum operators** and brings the **production sectors**
online via manufactured solutions.

**Contract.** Codex DESIGNS and WRITES all code (the MMS harness, the discrete operator forms for
each sector, the forcing terms, the convergence study, the tests). Claude REVIEWS code + output
with clean agents + an independent re-run, and iterates with you until the acceptance criteria are
met. **This directive states requirements and acceptance criteria only — it does NOT prescribe
the discretization route or the discrete operator stencils; those are your design call** (mirrors
the audit's "Claude reviews, Codex codes"). Run everything **locally on CPU** (free); no GPU, no
runpods, no network. Iterate until the harness exits 0 and every criterion below holds.

---

## 1. Objective

Use the **method of manufactured solutions (MMS)** to certify, for **each operator block of the
production residual**, that the *discrete* operator converges to the *continuum* operator at the
scheme's formal order. MMS is the validation-first vehicle for the transition from the step-1
cubic-GPE stand-in to the **production physics**: we pick an analytic field, push it through the
continuum operator analytically to get an exact forcing, then confirm the discrete operator
reproduces that forcing under refinement. This certifies operators we will *not* otherwise have a
closed-form answer for (the stiff quintic matter residual, the localized-Maxwell operator, the
effective-closure wall operator), **before** any of them is trusted inside the coupled stationary
branch solve (step 3).

This is still validation-only: **no coupled nonlinear branch solve, no §H target, no observable
from §G of the pre-registration.** We are proving the machinery, one operator at a time.

The operators built and certified here are **production machinery** (stack-decision audit risk #2):
the discrete forms, measures, boundary operators, and forcing/residual layer MUST be the same ones
WP1/WP3 will use later. MMS is exactly the test that they are.

---

## 2. Scope

### IN — mandatory deliverables

**A. Reusable MMS harness.** A backend-neutral harness that, given (i) an analytic manufactured
field, (ii) the continuum operator applied to it analytically (the exact forcing), and (iii) the
discrete operator, measures the weighted-norm error vs grid refinement and reports the **observed
order of convergence**. It must work for scalar radial, 2D `(r,w)` tensor, complex, and
vector/multi-component fields. Reuse the step-1 grid/measure/boundary/manifest layers; do not fork
them.

**B. Matter sector (gauged GNLS, production quintic EOS).**
  Canonical source: `notes/moving_throat_pde_program_compact.md` §2.4.1 (matter Lagrangian +
  covariant derivatives), §2.5.1 (gauged GNLS field equation), with the **frozen stiff EOS**
  `P(rho)=K rho^5`, `U(rho)=(K/4) rho^5`, `h(rho)=dU/drho=(5K/4) rho^4`, `c_s^2=(5K/m) rho^4`.
  - Build the **production quintic matter residual operator** (the `h(rho)` nonlinearity replacing
    the step-1 cubic `g|psi|^2`), in the same conservative-FV form. MMS-certify it at formal order
    on a manufactured density/amplitude profile that actually exercises the `rho^4` stiffness
    (non-trivial `rho(r)`, not a constant). Keep the step-1 cubic-GPE benchmark intact as a
    **regression** (it must still pass).
  - **Follow-up (a) — genuine 2D `(r,w)` MMS convergence.** The tensor `(r,w)` Laplacian is
    currently only null-space tested (constant → 0). Add a manufactured solution with nontrivial
    `r`- and `w`-dependence so BOTH the radial flux (with full `4*pi*r^2` measure) and the `w`
    flux, plus the boundary operators on all faces, are exercised, and report a measured 2D order
    of convergence driving to a noise floor.
  - **Follow-up (b) — complex / phase-winding field exercising the current operator.** The current
    diagnostic `j = (hbar/m) Im(psi^* D_i psi)` (compact §2.5.2) is today checked only on real
    `l=0` fields, where it is identically zero — a tautology. Add a complex manufactured field with
    a known **nonzero analytic current** (e.g. a phase-winding / vortex-like `psi = f(r) e^{i*S}`
    with prescribed `f`, `S`), and MMS-certify the discrete current operator against that nonzero
    analytic current at formal order. (Gauge-covariant `D_i` if the gauge field is in scope for
    this benchmark; otherwise the plain current with `A=0`, stated explicitly.)

**C. Maxwell sector (localized, gauge-fixed H=Z).**
  Canonical source: compact §2.4.2 (localized Maxwell kinetic term `-(Z(w)/4 mu_0) F_{MN}F^{MN}`
  + weighted gauge-fixing `-(H(w)/2 xi mu_0)(d.A)^2`, **structural choice H=Z**, `xi_4=xi`) and
  §2.5.3 (exact localized Maxwell equation `d_M(Z(w) F^{MN}) + (1/xi) d^N(H(w) d.A) = mu_0 J^N`).
  - Build the **discrete localized-Maxwell operator** (the `A_0, A_r, A_w` field equation with the
    `Z(w)` weight and the H=Z gauge-fixing term) in conservative-FV form on the `(r,w)` grid, and
    MMS-certify it at formal order on a manufactured gauge potential with a known analytic RHS.
    The gauge-fixing convention is **frozen to H=Z** (prereg §D, M-1) — do not introduce a
    different gauge.
  - The matter current must **not** be double-counted as an explicit Maxwell source in the same
    term (compact §2.4.2 bookkeeping rule); for this MMS the source is the manufactured forcing.

**D. Wall sector (effective-closure quadratic operator `S_eta^(2)`).**
  Canonical source: compact §2.4 / §13.1 (`parent_action_status = effective_closure`, the wall is a
  fixed effective linear closure for Stage 1) and prereg §E (the quadratic wall action `S_eta^(2)`
  with constitutive packet `mu_eta, T_w, T_Omega, K_eta` held as **fixed posited inputs** — frozen,
  not derived, not refit). Provenance:
  `provenance/_synthesis/batch_01/fit_stage001_wall_action_constitutive_coefficients__mu_eta.yaml`,
  `__k_eta.yaml`.
  - Build the **discrete wall operator** (the linear operator from `S_eta^(2)` acting on the
    wall/support profile `chi_eta(s)`, with the frozen constitutive coefficients as explicit config
    inputs) and MMS-certify it at formal order on a manufactured wall profile.
  - The constitutive coefficients are **frozen inputs**, quoted from the provenance files — they are
    NOT tuned, fit, or chosen here. State their values and cite the source in the report.

**E. config_hash hygiene fix (follow-up c).** Restrict `HarnessConfig.config_hash()` to the
  reproducibility-relevant content — physics parameters, numerics, mesh spec, tolerances, and the
  source revision — and **exclude cosmetic output paths** (`run_root`, `report_path`). A move of the
  output directory must not change the config hash. Keep the full config still recorded verbatim in
  the manifest; only the *hash input* is narrowed.

### Per-sector "underspecified → STOP and flag" rule

The continuum forms for B/C/D are believed pinned in the cited sources at discrete-ready fidelity.
If, while building any operator, you find a continuum form is **not** pinned well enough to
manufacture a faithful solution (an ambiguous sign/factor/weight, an undefined boundary on a sector
field, a missing constitutive value), **STOP and flag it in your response with the specific
ambiguity and the source line** — do **not** guess an operator form. Deliver the sectors that are
unambiguous + the harness + the follow-ups, and flag the rest. A wrong operator silently passing
MMS against its own wrong forcing is the worst outcome; the harness only proves discrete==continuum
*for the continuum form you encoded*, so the form must be transliterated faithfully, not invented.

### OUT — do not touch this step

The coupled stationary branch Newton solve (step 3); any multi-sector *coupled* solve; PML/sponge
on the tangent (step 5); the WP3 P2 tangent (step 8); any §H target constant; any §G observable;
GPU/runpods/network; refitting any frozen constitutive/branch input.

---

## 3. Requirements

**R1 — Continuum fidelity (the #1 substantive axis).** Each discrete operator must transliterate
its **frozen continuum operator** from the cited source (B: §2.4.1/§2.5.1+EOS; C: §2.4.2/§2.5.3,
H=Z; D: prereg §E `S_eta^(2)`). No sign, factor, power, weight, or measure may drift from the
source. The MMS forcing must be derived from the **same continuum operator** the discrete operator
targets — and the analytic forcing derivation must be shown (so review can check the operator is
not being tested against a forcing that secretly encodes the discrete scheme). Where the analytic
forcing is algebraically heavy, deriving it with `sympy` is encouraged; the symbolic derivation is
part of the deliverable.

**R2 — Production-transferable, conservative form.** Operators stay in conservative FD/FV form with
the full `r^2` / `(r,w)` measure and `r=0` regularity via zero inner-face area (audit risk #2 — the
machinery here MUST be the machinery WP1/WP3 uses). Reuse step-1 primitives; extend, do not fork.

**R3 — Convergence per operator.** ≥3 grid levels per manufactured benchmark; report the
**observed order of convergence per operator** in a weighted norm, and show error driving to a
noise floor. Observed order must match the scheme's formal order; a shortfall is a finding to
resolve, not to paper over. For the 2D benchmark, refine in both `r` and `w`.

**R4 — Newton/Jacobian integrity carries forward.** For any sector whose MMS uses a nonlinear
solve (e.g. the quintic matter residual solved to its manufactured field), the JVP/autodiff
Jacobian is cross-checked against an independent finite-difference directional-derivative probe and
the agreement reported (audit risk #6). All solver controls explicit in config and logged (audit
risk #3). Pure-linear operators (Maxwell, wall) may be certified by direct operator-application
error without a Newton solve — state which path each sector uses and why.

**R5 — Determinism & manifest.** Deterministic runs; every run emits a manifest (dtype,
tolerances, mesh, config hash, git revision, library versions) — the `NONLINEAR_PROTOCOL_V2.md`
freeze discipline, **minus any target fields** (there are none). Apply the §2-E config_hash fix.

**R6 — TARGET-BLINDNESS (the #1 review axis).** No §H target constant — `54 G c_s^5/5 a^5 c^5`,
`chi_Q=1`, `R_pole/R_norm/P_2/P_4`, the GR quadrupole, any benchmark target value — may appear in
code, config, manufactured fields, forcing, stopping criteria, or tuning. Manufactured solutions
are arbitrary smooth analytic fields with **no relation** to any physical target; choosing one that
encodes a target is a target-blindness violation. The constitutive wall coefficients (§2-D) are the
only externally-fixed numbers and they are frozen *inputs* quoted from provenance, not targets.

**R7 — Frozen inputs are not refit.** The quintic `K`, the gauge convention H=Z, and the wall
constitutive packet (`mu_eta, T_w, T_Omega, K_eta`) are **frozen**. They are config inputs quoted
from their sources; they are never tuned, fit, or adjusted to make a benchmark pass.

**R8 — Environment & deps.** Use the installed env (numpy/scipy/torch/sympy/h5py/matplotlib). If a
new dependency seems genuinely needed, **STOP and flag it** with the reason — do not silently
install.

**R9 — No fabricated success.** If an operator cannot reach formal order, report the shortfall and
diagnose it; do not loosen tolerances or hand-pick grids to manufacture a pass. A reduced order
that is *understood and documented* (e.g. a known boundary-stencil order drop) is acceptable; an
unexplained one is a finding.

---

## 4. Acceptance criteria (the gate to step 3)

1. **MMS harness** exists, is backend-neutral, reuses the step-1 grid/measure/boundary/manifest
   layers, and is driven by config (grid levels, manufactured-field choice, tolerances explicit).
2. **Matter quintic operator** MMS-certified at formal order, error → floor; step-1 cubic-GPE
   benchmark still passes as a regression.
3. **Follow-up (a):** genuine 2D `(r,w)` manufactured-solution convergence at formal order in both
   directions, error → floor (not the null-space-only check).
4. **Follow-up (b):** complex/phase-winding manufactured field with nonzero analytic current; the
   discrete current operator certified against it at formal order (no longer a `j≡0` tautology).
5. **Maxwell (H=Z)** and **wall (`S_eta^(2)`)** operators each MMS-certified at formal order —
   OR a precise STOP-and-flag per §2 if a continuum form is found underspecified.
6. **Follow-up (c):** `config_hash` excludes cosmetic output paths; verified that moving
   `run_root`/`report_path` does not change the hash, while full config is still in the manifest.
7. Jacobian cross-check (R4) passes for any nonlinear-solve sector; conservation/symmetry drift at
   noise floor where applicable.
8. Runs deterministic; manifest emitted per run.
9. A concise **validation report** (markdown/YAML, KB-scale) states, per operator block: the
   continuum source + cited lines, the manufactured field, the analytic forcing (and how derived),
   grid levels, observed order, final error vs the manufactured answer, the solver/operator config
   used, and (matter) the Jacobian-check residual + conservation drift.
10. Target-blindness (R6) holds — verifiable by `grep`.

Claude reviews code + report (clean agent), checks transliteration fidelity against the cited
sources, runs an independent arbiter re-run reproducing the convergence orders, and we iterate with
you until all hold.

---

## 5. Repo hygiene (the firewall is live — do not blow up the repo)

- **Code** → `software/stage1_solver/src/` and tests → `software/stage1_solver/tests/` (both
  tracked). Extend the existing module layout cleanly (e.g. an `mms/` or `manufactured.py` module,
  per-sector operator modules); do not duplicate step-1 code.
- **All run output** (fields, convergence arrays, plots, checkpoints) → `software/stage1_solver/`
  `runs/` / `figures/` / `data/` — **gitignored**; nothing there gets committed.
- **The one tracked result artifact** is the small validation report →
  `software/stage1_solver/reports/` (tracked; KB-scale, markdown/YAML).
- **Do NOT commit** `*.npz/*.h5/*.pt/*.vtk/*.pkl/...` or any large data — the `.gitignore` firewall
  blocks these; do not `git add -f` around it. **Never `git add -A`.**
- You apply + run + iterate; the orchestrator (Claude) handles staging/commits with explicit paths
  after review.

---

## 6. Deliverable for review

When the acceptance criteria are met, report back with: the module layout (what you added vs
reused), the discrete operator form you chose for each sector and the continuum source line it
transliterates, the manufactured fields + how each analytic forcing was derived, the per-operator
convergence tables (observed order, final error, noise floor), the Jacobian-check result for the
matter sector, the config_hash change, the path to the validation report, and any STOP-and-flag
items. Then Claude reviews before anything is committed.
