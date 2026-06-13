# BUILD DIRECTIVE — Step 6 (CPU): conservation diagnostics of the coupled branch

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 6 ("Conservation diagnostics. Mass /
charge / energy drift over each run. *Brief §5.4.*"); brief `docs/branch_realization_brief.md` §5.4
("Conservation diagnostics. Report mass / charge / energy drift (as appropriate to the sector) over
the run."). Pre-reg `docs/stage1_preregistration.md` **§J.4 (non-negotiable):** "Conservation
diagnostics — report mass/charge/energy drift over each run (plan §7 step 6)." Exact conservation
laws live in `notes/moving_throat_pde_engine_handoff.md` §3.2 (number current + continuity), §1.3 /
§2.1 (charge ontology + matter Lagrangian), §2.2 / §3.3 (localized Maxwell + current conservation).
**Stack:** `decisions/01_stack_choice.md`. **Builds on:** step 3a (`dcd6a1b`) coupled residual, step 3b
(`affc745`) preconditioner, step 4 (`c9b8b2c`) convergence study + target-blind surrogate observables +
the null/floor vocabulary, step 5 (`63cd885`) boundary characterization + the real interior-zero sponge.
REUSE all of it.

---

## Methodology scope (read carefully — this is the load-bearing framing)

**This branch is a STATIONARY Newton fixed point, not a time-integration.** There is no timestep, so
"drift over each run" in the prior-attempt sense (a conserved quantity wandering as a PDE is marched in
time — brief §6 / the §112 failure mode) **does not exist here by construction**. The conservation
*content* of a stationary state is instead:

1. the **conserved-current divergences vanish in the bulk** — exact continuity `∂_i j^i = 0` (number),
   `∂_M J^M_tot = 0` (charge), and the stationary energy-flux balance; and
2. the **global budgets close** — net flux through a closed interior surface balances the explicit
   sponge sink plus the discretization residual.

Realize the diagnostics in that stationary form. Do **not** invent a time-marching loop to manufacture a
drift number.

**Expect null/floor diagnostics — report them honestly, do not manufacture signal.** The step-4 study
already found the isotropic stationary branch carries **essentially zero net spatial current and zero
spatial gauge field by symmetry** (a real-amplitude isotropic ground state has no phase winding and no
transverse vector potential, so `j = (ℏ/m)Im(ψ*Dψ) ≈ 0` to the solver floor). Therefore the local
*current-divergence* continuity residuals (mass, charge spatial-current) and the interior energy-*flux*
divergence will sit at the numerical floor not because the scheme is perfectly conservative but because
**there is no transport to conserve in this branch**. That is the correct, honest result. Label these
exactly that way using the **existing null/floor vocabulary** from `convergence.py` (`"null diagnostic"`,
`"solver-floor diagnostic"`). The genuinely current-*carrying* conservation test belongs to the WP3
linearized tangent (the wave-like sector) and is therefore **exercised at step 8, not here** — say so.

**The genuinely informative, non-trivial conservation diagnostics for THIS branch are §3-B (Gauss-law
closure) and §3-C (global budgets with explicit sponge accounting).** Those carry the step. Gauss-law
closure is an *independent* test of the Maxwell solve's divergence consistency (not a restatement of the
Newton residual); the sponge budget cross-checks step 5's interior-zero claim from the conservation side.

**Contract.** Codex DESIGNS and WRITES the diagnostics module, harness, report, and tests. Claude REVIEWS
code + output with clean agents + an independent arbiter re-run, and the standing transliteration-fidelity
audit runs on the new math. This directive states requirements + acceptance criteria only; the mechanism
is your design call. Run everything LOCALLY ON CPU; no GPU, no network. Never wrap your own session in a
shell `timeout` (the `timeout 600` cap is only for the short scripts you invoke). Still the **engineering
smoke** regime: placeholder params, **NO physical packet**, export guard untouched, **target-blind**.

---

## 1. Objective

For the converged coupled stationary isotropic branch, quantitatively report — to a stated numerical
floor and under grid refinement — that:

- **(mass)** the discrete number-continuity residual `∂_i j^i` is at the floor in the interior, and the
  total number is accounted for (interior + sponge-absorbed), with explicit sponge bookkeeping;
- **(charge)** the localized **Gauss law** closes — the surface flux of `Z F^{i0}` over closed interior
  surfaces equals `μ₀` times the enclosed charge `q∫ρ dV` to a stated relative residual that **shrinks
  under refinement** — and the charge-current divergence residual is reported (expected null/floor);
- **(energy)** the stationary energy budget closes — interior energy-flux divergence at the floor, with
  the sponge-dissipated energy quantified and the interior/sponge split reported.

Each non-null diagnostic must carry an observed **order of convergence** against the step-4 grid ladder,
and an honest **null/floor** label where the branch's symmetry makes the diagnostic trivial.

---

## 2. Conserved quantities and their exact laws (the spec — take exact forms from the source, not here)

From `notes/moving_throat_pde_engine_handoff.md` (cross-check against `notes/moving_throat_pde_program_compact.md`):

- **Number / mass.** Bulk number current `j^i = (ℏ/m) Im(ψ* D_i ψ)` with `D_i = ∂_i − i(q⋆/ℏ)A_i`
  (§3.2). Exact continuity `∂_t ρ + ∂_i j^i = 0`; stationary ⇒ `∂_i j^i = 0`. `ρ = |ψ|²`.
- **Charge.** Electric charge density `J^0 = q⋆ ρ`, spatial charge current `J^i = q⋆ j^i`
  (§1.3 / §2.1 charge ontology; `q⋆ = η_Q e⋆`). Current conservation `∂_M J^M_tot = 0` (§3.3, required by
  gauge invariance). **Localized Gauss law** = the `N=0` component of `∂_M(Z(w) F^{MN}) + (1/ξ)∂^N(∂·A)
  = μ₀ J^N_tot` (§2.2 / §3.3): `∂_i(Z F^{i0}) = μ₀ J^0_tot` (in the stationary gauge with `∂·A` handled
  by the existing gauge-fixing term).
- **Energy.** Energy density from the matter Lagrangian `𝓛_ψ` (§2.1: kinetic `(ℏ²/2m)|D_iψ|²`,
  confinement `V_conf ρ`, internal `U(ρ)=(K/4)ρ⁵`) plus the EM field energy from `𝓛_EM`
  (§2.2: `(Z/4μ₀)F_{MN}F^{MN}`). Stationary ⇒ energy-flux divergence balances sources/sinks.

These are PDE structure derived from the parent action — they are **not** §H targets and **not** §F
extraction coefficients. Diagnosing them is firewall-safe (see §5).

---

## 3. Scope — IN (mandatory)

**A. Local conservation residuals (honest null/floor expected).**
Compute, on the converged step-4-resolved coupled state, the discrete divergence of each conserved
current over the **interior window** (same interior region step 4/5 use for surrogate evaluation,
excluding the sponge/boundary layer):
- number-continuity residual `∂_i j^i`,
- charge-current divergence `∂_i J^i = q⋆ ∂_i j^i`,
- interior energy-flux divergence.
Report each as interior L2 **and** L∞, normalized by a **target-blind** reference scale (an aggregate
norm — e.g. total interior number `N=∫ρ dV`, total charge, total energy, or `‖j‖` — **never** an
extraction coefficient). Apply the existing `convergence.py` null/floor labelling: if the current itself
is ≤ ~1e-14 the residual is a `"null diagnostic"` (no transport to conserve — isotropic-branch symmetry);
otherwise it is a `"solver-floor diagnostic"`. State plainly which sectors are null for this branch and
why, and that the current-carrying test is deferred to the WP3 tangent (step 8).

**B. Gauss-law closure (the genuinely non-trivial charge diagnostic — REQUIRED).**
Over **≥2 nested closed interior surfaces** (e.g. constant-`r` shells inside the interior window),
verify `∮ Z F^{i0} dA_i = μ₀ q⋆ ∫ρ dV` (enclosed charge). Use the **same conservative face-flux /
measure machinery** the operators already use (so the discrete divergence theorem holds cell-by-cell —
see §4). Report the absolute and **relative** closure residual `(LHS − RHS)/|RHS|` for each surface, and
its **observed order of convergence** across ≥3 grid levels from the step-4 ladder. This residual must
**decrease under refinement**; a closure residual that does not shrink is a red flag — report it as a
failure, do not hide it.

**C. Global budgets with explicit sponge accounting (REQUIRED — cross-checks step 5).**
For mass, charge, and energy, report the closure:
> (net outward flux through the interior-window boundary) = (sponge-absorbed amount in the outer layer)
> + (interior local residual).
Quantify the **sponge sink** explicitly: how much number/charge/energy the `σ·χ·field` term removes,
and confirm from the conservation side that it is **interior-zero** (the sponge contributes nothing
inside the interior window — the conservation-side restatement of step 5's interior-zero property). Run
this **both with the sponge ON and OFF** so the budget difference isolates exactly the sponge
contribution. Report each closure absolutely and relative to the corresponding total.

**D. Convergence of every non-null diagnostic.**
Reuse the step-4 grid ladder (≥3 levels). For each non-null diagnostic (Gauss closure; any non-null
local residual; the budget closures), report the value at each level and the observed order via the
existing `observed_order_from_three` helper. Tie the floor to the step-4 numerical floor.

**E. Machine-readable table + tracked report** (see §6).

---

## OUT — explicitly deferred / forbidden

- **No time integration / no drift-over-time loop.** Stationary branch only (see Methodology scope).
- **No current-carrying conservation test** (that is the WP3 tangent — step 8). Flag it as deferred.
- **A virial / Pohozaev energy identity is OPTIONAL, not required.** Include it *only* if you can derive
  it faithfully from the parent action `𝓛_ψ + 𝓛_EM` and verify it term-by-term; if there is any doubt,
  **omit it and say so** rather than ship a plausible-but-wrong identity. (A faithful-but-wrong identity
  is exactly the failure the transliteration-fidelity audit exists to catch — don't create one under
  time pressure.)
- **No physical packet, no export.** Export guard untouched (`physical_export_permitted` stays False).
- **No new dependencies.** torch + numpy + sympy (for any symbolic check) only — already in the stack.
- **Never write under `research/pde_audit/simulation/`.** All run output under the gitignored
  `runs/`/`data/`/`figures/`/`_scratch/`; only code + tests + the small report + this directive are tracked.

---

## 4. Reuse requirements (fidelity — diagnose, do NOT reimplement physics)

- **Number current:** use the existing `_matter_number_current(fields, grid, cfg)` in
  `coupled_branch.py` (it already implements `j = (ℏ/m)Im(ψ*Dψ)` with the covariant `−(q/m)A` term). Do
  **not** introduce a second, parallel current formula.
- **Divergence / flux:** use the **same conservative finite-volume face-flux construction** the
  `tensor_laplacian` / `tensor_weighted_laplacian` operators use (area-weighted face fluxes differenced
  and divided by `grid.cell_volumes`), so the discrete divergence theorem holds exactly cell-by-cell and
  the surface-flux ↔ volume-integral identity in §3-B/§3-C is exact at the discrete level. Do **not** use
  a naive non-conservative centered finite difference for any divergence in a budget/closure check.
- **Measure tensors:** integrate with `grid.cell_volumes` / `grid.face_areas` / `grid.radial_shell_volumes`
  (the explicit `r²` measure — 3D-on-2D). Do not hand-roll a measure.
- **Sponge:** use the existing `boundary_sponge_profile_torch` + the `sponge_*_strength` config fields;
  the sponge sink terms are `σ·χ·field` already present in `coupled_pde_residual`. Read the absorbed
  amount off those exact terms, not a re-derivation.
- **Null/floor vocabulary + order helper:** reuse `convergence.py`'s null/floor labels and
  `observed_order_from_three` (do not invent parallel verdict strings).
- **EOS:** `U(ρ)=(K/4)ρ⁵`, `h(ρ)=(5K/4)ρ⁴` (stiff, `P=Kρ⁵`) — matches `physics.quintic_enthalpy`. Do
  **not** "simplify" to a soft polytrope or a `|ψ|⁴ψ` quintic.

---

## 5. Firewall / target-blindness (the #1 review axis)

- Conservation diagnostics use only conserved-current structure (continuity, Gauss, energy) derived from
  the parent action. They must touch **no** §H target (R_norm, R_pole, P2=P4=0, chi_Q, GR quadrupole,
  Λ₀) and **no** §F extraction coefficient (D0/D2/D4, N0/N2/N4, P0/P2/P4, chi_Q, N_Q, m̂0, S_port, Xi_1,
  varrho). Do not import, read, or compare against the extraction map or `benchmarks.py` target values.
- All normalization scales must be **target-blind aggregate norms** (total N / charge / energy / ‖j‖),
  never an extraction coefficient.
- Export guard stays False; emit no physical packet.

---

## 6. Deliverables

1. A diagnostics module (your naming; e.g. `conservation_diagnostics.py`) + a thin harness
   (e.g. `conservation_harness.py`) mirroring the step-4/5 harness pattern, writing a machine-readable
   table (JSON or markdown table) under `runs/…` (gitignored).
2. `software/stage1_solver/reports/step6_conservation_diagnostics.md` — the tracked result report (§7).
3. Tests in `tests/test_stage1_solver.py` (see §8).
4. Do **NOT** `git add` or commit — the orchestrator stages + commits after re-review.

---

## 7. Report requirements (`reports/step6_conservation_diagnostics.md`)

- State the **stationary-vs-time-drift framing** up front (no time-marching; conservation = bulk
  divergence + budget closure).
- A per-sector table (mass / charge / energy): the local residual (interior L2/L∞ + null/floor label),
  the budget closure (with/without sponge, absolute + relative), and — for non-null diagnostics — the
  observed order of convergence.
- The **Gauss-law closure** table: per nested surface, relative residual at each grid level + observed
  order; an explicit statement that this is an *independent* test of the Maxwell solve, not a restatement
  of the Newton residual.
- The **sponge accounting**: absorbed number/charge/energy, and the conservation-side confirmation that
  the sponge is interior-zero (cross-reference step 5).
- An explicit **honest-limitations** paragraph: which diagnostics are null/floor for this branch and why
  (isotropic symmetry → no spatial transport), and that the current-carrying conservation test is
  deferred to the WP3 tangent (step 8). State whether the optional virial identity was included or
  omitted and why.
- Reproduction command(s); confirm target-blind + export guard untouched.

---

## 8. Acceptance criteria

1. **Physics-faithful & non-reimplemented:** number current via `_matter_number_current`; all
   divergences via the conservative FV face-flux machinery; measure via the grid measure tensors; EOS
   `U=(K/4)ρ⁵`. (The transliteration-fidelity audit will check this term-by-term.)
2. **Gauss-law closure decreases under refinement** with a sane observed order (≈ the step-4 order);
   reported honestly if it does not.
3. **Budgets close** to the floor with the sponge accounted for, ON and OFF; the sponge confirmed
   interior-zero from the conservation side.
4. **Null/floor diagnostics labelled honestly** (reusing `convergence.py` vocabulary); no manufactured
   signal; the WP3/step-8 deferral stated.
5. **Target-blind:** no target / extraction-coefficient access; export guard untouched; no new deps.
6. **`pytest` passes** (expect the prior count + the new test(s)); the conservation harness **exits 0**;
   results reproduce to the digit on a re-run (only wall-clock may differ).
7. No script runs > 10 min (`timeout 600` on each script you invoke); a timeout is a failure to
   reformulate, not a cap to raise.

---

## 9. Tests (§8.3 detail)

Add focused, **non-tautological** tests to `tests/test_stage1_solver.py`:
- **Discrete divergence theorem (the load-bearing test):** on a hand-constructed field, assert the
  conservative FV volume-integral of `∂_i j^i` over a region equals the face-flux `∮ j·dA` over that
  region's boundary to ~machine precision (this is what makes §3-B/§3-C exact and non-tautological).
- **Gauss closure on a known field:** construct an `A0` whose `−∂_i(Z F^{i0})` is a known charge
  distribution and assert the closure recovers it.
- **Sponge interior-zero from the conservation side:** assert the sponge contributes exactly zero to the
  interior budget and a positive, expected amount in the outer layer (sign + magnitude), cross-checking
  the step-5 wiring test.
- A **null-diagnostic guard:** assert the labelling correctly flags a zero-current sector as
  `"null diagnostic"`.

**When done**, report: the new module + harness names; each diagnostic's headline number (local
residuals + labels, Gauss closure residuals + orders, budget closures, sponge-absorbed amounts); the new
test names + what each asserts; confirmation `pytest` passes and the harness exits 0; and confirmation
target-blindness + export guard are intact.
