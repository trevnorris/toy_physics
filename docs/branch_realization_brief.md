# Branch-Realization Test — Problem Brief

**For:** a computational physicist, or an AI system capable of validated numerical PDE work.
**From:** Trevor Norris (independent researcher).
**Status of the work being tested:** a speculative research program. This brief asks for a *decisive numerical test*, not an endorsement of the underlying model.

---

## 1. What I'm asking for, in one paragraph

I have a parent field theory and a set of reduced target observables that the theory is *supposed* to reproduce. I need someone to numerically realize the relevant branch of the parent PDE and determine — **target-blind, with a proper error budget** — whether the realized branch returns those targets, or whether the apparent agreement is numerical artifact. A trustworthy *no* is as valuable to me as a trustworthy *yes*. What I cannot use is another run that produces numbers without establishing that the numbers are scientifically valid.

---

## 2. The system (decoupled from interpretation)

You do **not** need to assess the physical interpretation of the model to do this work. Treat it as a well-defined PDE plus a list of observables plus a validation requirement.

- **Parent system:** a (4+1)D gauged nonlinear Schrödinger / Gross–Pitaevskii field coupled to a localized Maxwell sector, with a stiff polytropic equation of state (exponent n = 5), plus a moving-throat geometry sector.
- **Exact action, gauge handling, boundary protocol, and the throat-field promotion** are specified in the source ledger (see §7). Use those definitions verbatim; do not reconstruct them.
- **The relevant object** is the moving-throat branch: solve it, then extract the reduced impedance / Dirichlet-to-Neumann (DtN) data and the grouped-coefficient observables from the solution.

---

## 3. The task

1. **Realize the branch.** Numerically solve the appropriate stationary (and, where specified, dynamic) moving-throat branch of the parent PDE.
2. **Extract the reduced observables** as defined in the source ledger — at minimum the coefficient bundle (`D0, D2, D4, N0, N2, N4, P0, P2, P4`), the outgoing-normalization scalar `χ_Q`, the prefactor slope `Ξ1 = P1/P0`, the rigid-mouth orbit-lock chart `(U, V)`, and the selected support-placement coordinate. Pull exact definitions from §13.2 and §12 of the source ledger.
3. **Keep the extraction target-blind.** This is the single most important constraint. The target values (§4) must **not** be inserted into the solve, used to tune any free parameter, used to choose a branch, or used to set a stopping criterion. Solve first, extract second, compare to targets **last**.
4. **Pre-register the branch.** Before running, fix and write down *which* branch, boundary conditions, and closure are being realized, and record that choice. The program admits many branches and relaxed/open-system companions; if a miss is allowed to send you to the next branch in search of a pass, the test stops being a falsification and becomes an unbounded search that no negative result can ever end. Name the branch first, then let the result stand.

---

## 4. The targets (used only for comparison, after extraction)

From the source ledger's target card (§13.3 / §13.4). Check the realized branch against these; do not steer toward them.

- `R_pole = D0(B4+Z4) − 3(M+B2+Z2)^2 = 0`
- `R_norm = m̂0² · S_port · N0/D0 − 54 G c_s^5 / (5 a^5 c^5) = 0`
- `P2 = P4 = 0`
- `R_tail = Θ_tail (c/c_s)^3 − 1 = 0` (if the tail sector is active)
- `χ_Q = 1`, and the 5PN actual-branch conditions `d ln R_tr = d ln R_target = d ln ε_η = 0`

Take the exact, current form of each from the ledger, not from this summary.

---

## 5. Validation requirements (non-negotiable — this is where the prior attempt failed)

The previous attempt produced a stable, open branch but the data was noise-dominated and could not be certified as scientifically valid. The following are mandatory, and a result without them does not count as an answer.

1. **Validate before trusting.** Demonstrate the solver on at least one limit with a known analytic or well-benchmarked answer (e.g. a free/linear limit, a known GPE vortex or soliton solution, or a standard GPE benchmark) *before* extracting any unknown observable.
2. **Convergence study.** Refine spatial grid and timestep; show each extracted observable converges and report the observed order of convergence. An observable that drifts under refinement is not a measurement.
3. **Boundary control.** Use absorbing / PML / sponge boundary layers; demonstrate quantitatively that boundary reflections sit below the target signal.
4. **Conservation diagnostics.** Report mass / charge / energy drift (as appropriate to the sector) over the run.
5. **Error budget and noise floor.** State the numerical noise floor explicitly. Every extracted observable must carry a quantitative uncertainty. A match or mismatch with a target is only meaningful relative to that floor.

---

## 6. Deliverables and acceptance criteria

**Deliverables**
1. The validated solver and a reproducible configuration.
2. A validation + convergence report establishing the solver is trustworthy *before* any physics claim.
3. Target-blind extracted values for each observable in §3.2, each with a stated uncertainty.
4. A per-target verdict: does the realized branch meet each target within validated error, and if not, by how much relative to the noise floor?

**What counts as an answer**
- **Trustworthy NO** — the branch misses one or more targets beyond validated error, or the apparent signal is shown to be numerical artifact. *This is a valid, valuable result.*
- **Trustworthy YES** — the branch returns at least one target, target-blind, within validated error.
- **Not an answer** — numbers produced without the §5 validation report. (The whole point of the exercise is to avoid exactly this.)

---

## 7. Interpreting the result (what a pass and a fail actually mean)

So that the outcome is read correctly by whoever runs this:

- **A validated, target-blind pass** on one or more targets means the unfitted dynamics reproduce that piece of known physics in a way that could not have been rigged. That is strong evidence the model is a **viable candidate** — it is *consistent with* reality without having been told the answer. It is **not** proof the model is the correct description of reality, because reproducing known physics is non-unique: more than one theory reproduces 1PN gravity, and inverse-square attraction follows from almost any sink model. A pass earns the model a seat among viable candidates; it does not single it out as true.
- **A validated miss, on a pre-registered branch (§3.4)**, falsifies that branch's claim. If the branch was the program's committed prediction, that is a clean negative on the program's core premise.
- **Numbers without the §5 validation report** settle nothing either way.

The step that would move the model from "viable candidate" to "true claim about the world" is **not** another reproduced target — it is a *novel, confirmed prediction* that diverges from GR/QED and turns out right. That is the purpose of Stage 2.

---

## 8. Stage 2 — divergence extraction (conditional; only after a Stage-1 pass)

Run this only if Stage 1 produced at least one validated target-blind pass. A novel prediction from a branch that cannot reproduce known physics is worthless, so this stage is gated on Stage 1.

1. **Use the same branch, nothing re-tuned.** All divergence predictions must come from the identical validated, unfitted branch realized in Stage 1. Any new free parameter introduced here voids the result.
2. **Find where the model diverges from GR/QED.** Identify regimes or observables where the model's prediction differs from the incumbent theory, and compute the size of that difference from the branch.
3. **Classify every divergence by testability, honestly:**
   - *Measurable now* — the divergence exceeds current or near-future experimental sensitivity in a regime where data exists or is obtainable. **These are the only ones that matter.** Prioritize them.
   - *Real but below the noise floor* — the divergence is smaller than the §5 numerical noise floor. Not a prediction yet; report as such.
   - *Unmeasurable in principle* — the divergence lives only where no experiment can reach. This is a promissory note, not a test, and must be labeled that way rather than presented as a prediction.
4. **Validate the difference to the same standard.** A claimed divergence must clear the same convergence and error-budget bar as the Stage-1 targets. A "prediction" buried in numerical noise is not a prediction.

The deliverable of Stage 2 is a short list of *measurable-now* divergences with magnitudes and the specific measurement that would see them — or an honest statement that the model diverges from GR/QED only in regimes nothing can currently probe.

---

## 9. Source documents (exact equations live here, not in this brief)

`moving_throat_pde_program_compact.md`:
- §1.4 — strict parent action and throat-field promotion status
- §13.1 / §13.1A — parent / gauge / boundary and projection-first EM cards
- §13.2 — full-bundle coefficient definitions
- §13.3 — target card
- §13.4 — 5PN actual-branch card
- §13.5 — solver card and prior-run status
- §12.11–12.12 — branch / orbit / support / open-system variable definitions
- V2-19 through V2-23 — the executable branch-realization protocol and the first target-blind open-throat residual extraction

**Prior-run note for whoever picks this up:** the first reduced open branch (protocol V2-23) came out *stable and open but target-failing*, and the runs were noise-dominated. The validation requirements in §5 exist specifically to determine whether that negative result is real or an artifact of insufficient numerical fidelity. Settling that question — either way — is the deliverable.
