# Directive — M0: Mathematica 15 environment / toolchain check

You are the **CODER** (Codex codes / Claude reviews). Build and run a Mathematica 15 environment-check
script for a new second-engine workstream: a Wolfram-Language **moving-throat branch producer** that will
solve a stationary open-throat profile + a linearized wall/BdG/Maxwell modal problem and export a frozen
branch packet for an existing Python residual chain.

This step does **one** thing: **discover and honestly report which Mathematica 15 capabilities the route
actually has.** A plan handed to us assumes several features that I am NOT sure exist (e.g.
`LinearSolve[..., Method -> "MUMPS"]`, `Method -> "HybridCPUGPU"`, and some matrix decompositions). Your
job is to find out the truth, not to make anything pass.

## Deliverables

1. `software/stage1_solver/mathematica/mt15_00_environment_check.wls` — a reproducible WolframScript.
2. `software/stage1_solver/mathematica/environment_report.md` — a **markdown** report (NOT JSON; this file
   is read by an LLM/human — see the project's YAML/markdown-not-JSON rule) summarizing every probe result.

## What to probe — for EACH item report one of: `PRESENT` / `PRESENT (caveat: ...)` / `ABSENT`

Do NOT abort on an absent feature — `Quiet`/`Check`/`FailureQ` around each probe and record the outcome.
Map each probe to the route need in parentheses.

1. **Version & license.** `$VersionNumber` (confirm ≥ 15.0), `$Version`, `$MachineName`,
   `$ProcessorCount`, `$MaxLicenseProcesses` (we have a 2-seat cap — confirm it). (gating prerequisite)
2. **FEM / NDSolve** (stationary open-throat + linearized BVPs):
   - Load `Needs["NDSolve`FEM`"]`; report success.
   - Solve a trivial 1-D BVP (`u''[x]==-1`, `u[0]==u[1]==0`) with `Method -> {"FiniteElement"}` and report
     it returned a function. Then a trivial 2-D Poisson on a unit square via FEM. Report success/failure.
3. **Coordinate / curvilinear support** (keep throat coords natural): does `CoordinateChartData` exist and
   return data for `"Spherical"`? Can `Laplacian`/`Div`/`Grad` take a curvilinear coordinate spec
   (e.g. `Laplacian[f[r,θ,φ], {r,θ,φ}, "Spherical"]`)? Report.
4. **Eigensolvers** (BdG / mixed modal extraction): `Eigensystem` on a small dense matrix; a **generalized**
   eigenproblem `Eigenvalues[{A, B}]`; and a sparse Arnoldi probe
   `Eigenvalues[sparseA, 3, Method -> {"Arnoldi"}]`. Report which work.
5. **Sparse LinearSolve methods** (the real bottleneck — scalable linear algebra): build a small
   `SparseArray`, then for EACH method string probe whether `LinearSolve[A, b, Method -> m]` is accepted and
   returns a correct solution: `"Multifrontal"`, `"Pardiso"`, `"Krylov"`, `"Banded"`, **`"MUMPS"`**. For
   `"MUMPS"` and any unknown method, explicitly report PRESENT vs ABSENT — **do not assume MUMPS exists.**
6. **Matrix decompositions** (stability classification + pole extraction). For EACH symbol, report whether it
   is a *defined* function (`Head` is not `Symbol`-unknown / no `::shdw` / not undefined) AND runs on a tiny
   matrix: `LUDecomposition`, `CholeskyDecomposition`, `SchurDecomposition`, `JordanDecomposition`,
   `QRDecomposition`, `SingularValueDecomposition`, and the **uncertain** ones the plan named —
   `LDLDecomposition`, `BunchKaufmanDecomposition`, `OrderedSchurDecomposition`, `PolarDecomposition`,
   `PopovDecomposition`, `JordanReduce`, `FrobeniusReduce`, `PartialFractions`, `PolynomialReduction`. For
   each uncertain one, state plainly whether it actually exists in this install.
7. **GPU** (nice-to-have, NOT required): is `CUDALink`/`GPUArray`/`TargetDevice -> "GPU"` available? Probe
   `Method -> "HybridCPUGPU"` for LinearSolve and report PRESENT/ABSENT. Do not require GPU.
8. **Export** (the V2-22B packet handoff): `Export` to `.json` works (round-trip a small assoc through
   `Export`/`Import`), and report whether any YAML export is available.
9. **ModelFit** (post-freeze diagnostic only): present? quick `LinearModelFit`/`NonlinearModelFit` smoke.

## Working rules

- Invoke via `wolframscript` if present, else `math -script` / `MathKernel`. Find the working invocation;
  the `.wls` must run clean to **exit 0** (probes that find a feature ABSENT are recorded, not errors).
- Wrap each script RUN in `timeout 600`. Iterate until the script exits 0 and the report is complete.
- ≤ 2 concurrent `math`/`wolframscript` kernels (2-seat license) — this step uses one.
- Create files only under `software/stage1_solver/mathematica/`. Do **not** `git add` or commit.
- No network. No GPU requirement.

## Report back to the reviewer

- The exact working kernel invocation (`wolframscript` vs `math -script`) and `$VersionNumber`.
- A compact PRESENT/ABSENT table for items 2–9.
- **Loudest:** which capabilities the route was hoping for that are ABSENT — especially `"MUMPS"`,
  `"HybridCPUGPU"`, and any of the uncertain decompositions — so we can re-plan around what truly exists.
- Anything surprising (e.g. a better sparse method that IS present, GPU available, etc.).
