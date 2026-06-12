# Decision 01 — numerical stack for the Stage-1 solver

**Status:** DECIDED 2026-06-12 (Claude+Codex methodology consult; the user delegated the
numerical-stack choice to us). Nothing here changes what the model physically claims, so no
user gate was required.

**Provenance:** consult prompt `decisions/01_stack_choice_consult.md`; Codex session
`019ebdb9-4b01-7722-90e5-d519c6f500f6` (gpt-5.5, read-only, xhigh); raw transcript
`software/stage1_solver/_scratch/01_stack_choice_raw.txt` (gitignored — the distilled decision
is *this* file). Codex read the frozen pre-registration + the engineering plan before answering.

---

## D1 (immediate, blocks step 1) — RATIFIED

Build steps 1–2 as a **custom, structured-grid finite-difference/finite-volume solver in
torch**, with:

- **Conservative FD/FV form** for the weighted divergence terms — carries the full `r^2`
  measure through face/cell weights, supports `r=0` regularity, Robin/open-exit conditions,
  and PML/sponge layers without a heavy mesh stack.
- **float64 from the start** — the production Jacobian conditions at ~1e6, which leaves no
  margin for float32 in any Newton/observable verdict. (Mixed precision allowed *only* for
  preconditioning/exploration, never for final extraction.)
- **Clean separation of the physics-residual layer from the backend array layer**, so the
  same residual code can later move to another backend with no rewrite.
- A **backend-neutral validation harness** (convergence, conservation, manufactured-solution,
  Jacobian checks) that does not depend on torch internals.

**Why torch (vs Dedalus / JAX / FEM / PETSc) for steps 1–2:** already installed (zero install,
runs on this CPU today); same code runs CPU-now / GPU-later; autodiff gives JFNK Jacobians;
fully transparent → best surface for the target-blindness audit. Dedalus is the *worst* fit for
this geometry (open boundaries + PML + low-density wall regions undermine spectral's advantages,
plus install friction and a poorer audit surface). JAX just swaps one custom-autodiff stack for
another and doesn't address the real D2 risk.

## D2 (deferred, informed by step-2/3 conditioning data) — direction set

**Defer the production-stack lock-in** (honors plan §3.2/§6). Stay in torch *only if* step-2/3
data show that matrix-free JFNK + a simple geometric-multigrid / block preconditioner gives
stable Newton convergence with roughly **grid-independent Krylov counts**.

**Escape hatch (biased to PETSc, NOT Dedalus):** if GMRES iteration counts grow strongly under
refinement, or the target-`K` Newton residual stalls, port the linearized-operator +
preconditioner path to PETSc-class structured-grid multigrid (hypre) / FEM. The
physics-residual, discrete-operator, boundary, extraction, and manifest layers carry forward
regardless — only the linear-algebra core moves.

**Decision trigger (the thing that settles D2):** *"Does the linearized target-`K` solve show
scalable, target-blind Newton–Krylov behavior under refinement?"* — NOT "does torch run."

---

## Load-bearing refinements from Codex (fold into every build directive)

- **The real wall is scalable linear algebra, not the discretization.** torch gives
  autodiff/JVPs but not mature AMG / Schur-block preconditioners. A 5–7-field/cell
  variable-coefficient Jacobian at production scale becomes the trap if we assume "autodiff
  gives Newton." **Treat preconditioning as the main project at step 3.**
- **Rewrite-trap is avoided iff** the torch code *owns* the discrete operators, validation
  harness, boundary operators, extraction pipeline, and frozen-manifest discipline. Those all
  carry forward; only the linear-algebra core is at risk of a port.
- **Discretization family:** conservative FD/FV is least-fighting for `(r,w)`. FEM is the
  fallback if weak-form/multigrid maturity dominates. Spectral is last.
- **First `n=5` preconditioner move (step 3):** (1) predeclared **continuation in `K`** with
  ordinary block/diagonal or simple geometric preconditioning on the coarse grid; (2)
  regularized density/amplitude variable *only if* near-vacuum positivity is the observed
  failure (do not reflexively switch to `log rho` if `psi_real,psi_imag` is smoother); (3)
  multigrid for production scalability.

## Target-leakage audit checklist (Codex Q6 — bake into reviews)

A solver passes review only if NONE of these occur:
1. Any §H target constant (`54 G c_s^5/5 a^5 c^5`, `chi_Q=1`, `R_pole/R_norm/P_2/P_4`) leaking
   into initialization, stopping criteria, PML tuning, continuation schedule, branch choice, or
   adaptive-mesh choices.
2. A toy GPE benchmark whose operators **diverge from the production operators** → validation
   non-transferable. (The discrete differential/measure/boundary/Newton machinery built at
   step 1 MUST be the machinery WP1 uses.)
3. Backend defaults hiding tolerances, dtype, nondeterminism, line-search, or stopping rules.
4. PML/sponge parameters tuned against `R_norm` instead of reflection benchmarks.
5. Continuation path / mesh family / branch gates changed after seeing residuals.
6. Autodiff Jacobians trusted without manufactured/JVP finite-difference cross-checks.
