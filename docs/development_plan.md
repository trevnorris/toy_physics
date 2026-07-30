# Development Plan — completing the PDE to "simulation-ready" (all sectors)

**Why this doc exists.** On **2026-06-26** we made a scope decision (memory
`project-simulation-deferred-complete-pde-strategy`): the full nonlinear 4D self-consistent fluid **simulation is out of reach**
(solver tractability, not hardware — hundreds of hours of multi-AI effort could barely hold a stable single throat in a box, and
the `S_leak` return throws spurious waves that corrupt measurements). So instead of stalling on the simulation, we **pivot wide**:
bring **every sector** of the PDE up to **"simulation-ready,"** and leave the actual simulation as delineated future work.

This decision **added weeks of work** — it converts "finish the gravity ladder" into "complete *all* sectors." This doc is the
master checklist of that full distance so no session forgets it. It is the superset above the gravity-only
`research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`.

**Front door:** `STATUS.md`. **Conceptual vision:** `docs/conceptual_foundation.md`. **This doc = the work breakdown.**

---

## 0. Definition of done — what "simulation-ready" means (per sector)

A sector is **simulation-ready** when:
1. its action / constitutive form is **fixed** (postulated where not derivable, and **labeled** as a postulate);
2. every coupling + calibration constant is **identified with a value or range** (the calibration budget);
3. every consistency gate checkable **algebraically/perturbatively** has been run — **passed, or a no-go recorded**;
4. **boundary / initial conditions** are specified;
5. the genuinely **simulation-dependent** quantities are listed with the **exact equations they must satisfy** and the
   **measurement that would confirm them**.

When gravity, EM/light, brane, and throat **all** hit (1)–(5), the **assembled PDE is ready to hand off** to whoever can run the
simulation later.

### The guardrail (load-bearing — never let this drift)
This plan **completes the specification; it does NOT prove the theory.** The simulation-dependent questions (does the
self-selected throat exist with the right `L/a`; does the cross-ℓ return close consistently; is the whole thing stable) stay
**genuinely OPEN — precisely posed, not answered.** And the EM/brane sectors can still hit a real **no-go** — which is welcome
(`feedback-falsification-is-the-goal`), not a failure of the plan.

---

## 1. Gravity sector — most complete (form delivered, magnitude calibrated)

**Status:** the PN ladder is GR-matched (calibrated; `research/4d_*pn*`, DON'T re-derive). The moving-throat realization is a
~6-gate ladder; Gates 1–4 DONE & CALIBRATED.

| Item | State | Artifact |
|---|---|---|
| Gate 1 — frozen-wall D/N resonator | ✅ `DN_UNITTEST_BC_DEPENDENT` | `pathA_30` |
| Gate 2 — breathing → `(a,L)` closure | ✅ `BREATHING_CALIBRATED` | `pathA_31` |
| Gate 3 — ℓ=2 isotropy | ✅ `ISOTROPY_CALIBRATED` | `pathA_32` |
| Gate 4 — `54/5` quadrupole | ✅ `QUAD_CALIBRATED` | `pathA_33` |
| **Gate 5 — scalar/dipole + cross-ℓ unification** | ✅ **DONE = `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (EARNED)** — linear theory can't pin the ℓ=0/1 return (`Z0_ret,Z1_ret` are a pathA_29 premise); rank audit nullity 2, flips to 0 under a real selector; Gate-6 must supply that selector. v1 rejected (rigged audit) → remediated → re-verified EARNED | `pathA_34` |
| Gate 6 — nonlinear closure | ☐ **algebra-tractable part only** (reduced/collective closure, Derrick existence, write the balance equations); the full 4D field solve is **SIM-DEFERRED** | — |
| PN match-back (threaded) | ☐ decisive falsifier — match extracted low-freq data to `research/4d_*pn*` | — |

**Remaining to sim-ready:** finish Gate 5; do Gate-6 Tier-A (reduced nonlinear closure + existence + posed balance equations);
run the PN match-back. Then the gravity sector's sim-dependent residue is just `{L/a, ε, stability}`, precisely posed.

---

## 2. EM / Light sector — largely OPEN (this is new bulk work)

**Status:** light = brane in-plane MacCullagh shear is **POSTULATED, conditional** (`pathA_23` Stage-2 CRUX =
`FAIL_UNSPECIFIED_SUBSTRUCTURE`: not derivable medium-natively). The **EM half of the ladder has not been run.**

**Remaining to sim-ready:**
- **Maxwell mixed channels** (the ladder's "EM half"): keep `A_w / J^w / F_{μw}` alive (do NOT zero them as ontology); derive
  the brane-effective EM reduction that rides the same level-set lift **without breaking Magnus (magnetism)**.
- **MacCullagh constitutive closure:** specify the postulated rotational-elastic form precisely; run the gates that CAN be
  checked algebraically — **Gate L** (2-transverse + no-longitudinal + bounded-below / inertially anchored + leak-free — THE
  CRUX), **Gate S** (magnetism: in-plane stiffness confined to the layer, bulk shear-free).
- **Cone-lock `λγ` (Gate K):** `c_γ` vs `c_s` — expected a **calibration gap**, recorded as such.
- **Output:** EM sector specified + calibrated + algebraically-gated (or a recorded no-go). Riding on the brane (§3), so it is
  partly **gated on the brane sector**.

---

## 3. Brane sector — OPEN HYPOTHESIS (the find-one-consistent-structure program)

**Status:** the density route is **CLOSED** (`FAIL_NOT_CODIM1` + `RC_DENSITY_SMECTIC_LIGHT_NOGO`); the brane must be a
**light-confining shear surface**, currently the **4D-throat-light hypothesis** (`project-light-is-4d-throat-hypothesis`). Not yet
a derived/consistent structure.

**Remaining to sim-ready:**
- Specify the postulated brane structure precisely (the shear surface + its embedding), labeled as postulate.
- Run the consistency gates that are algebraically checkable: **Gate L** (light on it — crux, shared with §2), **Gate S**
  (magnetism), **Gate B** (brane↔gravity compat — must not disturb long-wavelength `c_s`/flow or the GR-quadrupole bundle),
  **Gate Q** (two `±w` charge signs), **Gate T** (a defect punctures it → finite trapped-wave throat).
- **This is the most likely place for a genuine NO-GO** (two requirements incompatible in any structure). That is a first-class
  result, not a failure.
- **Output:** a specified, consistency-checked brane structure — or a recorded no-go that kills/redirects the program.

---

## 4. Throat / defect sector

**Status:** `pathA_26` Derrick + open-stability = `THROAT_DRAIN_DESTABILIZED` (NOT a kill — a conservative throat-soliton exists
generically). The drain-sector derivation is **PAUSED**.

**Remaining to sim-ready:**
- **Drain-sector derivation** (resume the paused work): pin `g_phys ≪ g_crit`; resolve the `J_w=0` mechanism.
- Write out the **self-selection balance equations** (trapped-wave pressure + drain backpressure + multi-brane return → `R*(E)`,
  the mass–radius relation). The actual **self-selected solve is SIM-DEFERRED** — but the equations + the existence argument
  (Derrick) are in-scope.
- **Output:** throat specified, existence supported, self-selection posed as equations + the measurement that confirms `L/a`.

---

## 5. Integration — assemble the simulation-ready PDE

**Remaining:**
- **Assemble** all four sectors into the single frozen parent action + BCs + a unified **calibration map** (every postulate and
  calibration constant labeled, with values/ranges).
- **Whole-system dimensional consistency** under the dimensional firewall (no LLM prose; dual-engine `dim_of`/`dimOf` on real
  expressions — `feedback-dimensional-consistency-check`).
- **Tracked debt to clear:** (a) ✅ **#98 DONE (2026-06-26)** — Gate-1–3 dimensional checks retrofitted (real-expression +
  able-to-fail + dual-engine; checks out); (b) the **χ_Q reconciliation** (`pathA_22b` numeric `≈0.712` vs `pathA_33` derived `=1`
  — different computations, same name) — still open.
- **The simulation hand-off spec:** the consolidated list of all sim-dependent quantities `{L/a, ε (per ℓ), stability, …}`, each
  with the exact equation it must satisfy and the measurement that confirms it. **This is the literal deliverable.**

---

## 6. Suggested sequencing (dependency- and risk-ordered; adjustable)

1. **Gravity arm — at its sim-ready boundary.** Gates 1–5 DONE (Gate 5 = `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`, EARNED). **Gate-6
   Tier-A is DEFERRED as an optional cap** (decision 2026-06-26): it would most likely just re-confirm that the `Z0_ret,Z1_ret`
   selector / `L/a` self-selection need the sim-deferred full solve (pathA_29 `Z_is_premise`; self-selection "requires dynamics") —
   a tidy-up, not a discovery; its 4D field solve is sim-deferred regardless. PN match-back stays the threaded falsifier.
2. **⭐ Brane sector (§3) — THE CHOSEN NEXT FRONT (2026-06-26).** The crux; it gates light/charge (light rides the brane); a no-go
   here would save all downstream EM work — highest-information step, AND algebra-tractable / NOT sim-gated. First: specify the
   light-confining shear-surface candidate, then **Gate L** (the crux).
3. **EM/Light sector (§2)** — on the brane that survives §3: MacCullagh closure + Maxwell mixed channels + `λγ`.
4. **Throat sector (§4)** — drain-sector + self-selection equations (can run in parallel with §2/§3; existence already supported).
5. **Integration (§5)** — assemble + whole-system dim check + clear debt + write the sim hand-off spec.

Each sector item runs through the **standard process** (directive → **ONE** design-review pass on a fresh reviewer → fold →
dual-engine execute → arbiter re-run → **TWO mutually independent fresh review legs** doing fidelity + per-tooth ablation, because a
sector item's math is **physics-bearing** → the **blocking physics leg**), per `docs/development_pipeline.md`. ⛔ **Retired
2026-07-29/30:** the Codex→GLM→fold→Codex bookend, the fixed three-leg tri-review, and the per-gate user gate — stop for the user at
a decision, a blocking finding or a no-go, not at every item. ⚠ The second leg (restored 2026-07-30) is a second **independent
reader**, not a second model in a sequence; prose and tooling still get one.

---

## 7. Honest scale note

This is **weeks of work**, deliberately. The trade we accepted on 2026-06-26: we cannot run the simulation, so the win is a
**complete, internally-consistent, calibrated, simulation-ready PDE** — every sector specified to the point a solver could take
over — with the genuinely-nonlinear questions precisely posed rather than answered, and a real chance of a no-go (especially in
§3/§2) along the way. That is a finishable, enumerable program; this doc is its ledger.

*Cross-refs:* memories `project-simulation-deferred-complete-pde-strategy`, `project-calibrated-pde-goal`,
`project-moving-throat-pde-push`, `project-brane-existence-defect-structure`, `project-light-is-4d-throat-hypothesis`,
`project-em-brane-shear-picture`, `project-pn-gravity-ladder`. Docs: `docs/conceptual_foundation.md`,
`docs/medium_requirements_and_prior_art.md`, `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`.
