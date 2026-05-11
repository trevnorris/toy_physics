# Moving-Throat Mathematica Verification Plan

## Goal

Build `research/pde_ledger/mathematica/` as a second, independent symbolic
verification layer for the moving-throat PDE derivation.

This is **not** a replacement for the existing SymPy audits and **not** a
line-by-line translation project. The Mathematica layer should instead:

- re-derive the same claims from the note equations independently,
- use a different symbolic engine to stress-test assumptions and branch logic,
- provide a second saved audit output for the most important stages, and
- make disagreements actionable by localizing whether the problem is in the
  notes, the SymPy script, the Mathematica script, or engine-specific
  simplification behavior.

---

## Operating Principles

1. **Independent derivation, not transliteration**

   A Mathematica script should not simply copy the SymPy script structure.
   It should recompute the claim from the note equations in a way that would
   still be useful if the SymPy script were wrong.

2. **Checkpoint-first coverage**

   Port the load-bearing stages first: checkpoints, final closure stages, and
   stages that previously had weak or tautological audits.

3. **Assumptions must be explicit**

   Each `.wl` script should state its assumptions clearly and use them only
   where justified. Mathematica can over-simplify under aggressive assumptions
   just as SymPy can.

4. **No tautological checks**

   Avoid tests of the form “define `x = expr`, then verify `x == expr`”.
   Each script must derive at least one key result from a different starting
   point than the statement being checked.

5. **Saved outputs are part of the audit trail**

   Every script we keep should have a corresponding saved output under
   `research/pde_ledger/mathematica/output/`.

---

## Proposed Layout

```text
research/pde_ledger/mathematica/
  VERIFICATION_PLAN.md
  run_all_audits.sh
  common/
    assertions.wl
    formatting.wl
    assumptions.wl
  output/
    _summary.txt
    moving_throat_pde_stageNNN_*.txt
  moving_throat_pde_stageNNN_*.wl
```

### Notes on layout

- Keep the file naming parallel to `scripts/moving_throat/` where practical.
- Use `.wl` rather than notebooks so the scripts are runnable noninteractively.
- Put shared helpers in `common/` only after at least 3-5 scripts justify the
  abstraction. Do not overbuild the framework up front.

---

## Script Contract

Each Mathematica stage script should satisfy the following contract.

### Input contract

- It should be runnable via:
  ```bash
  math -script research/pde_ledger/mathematica/<script>.wl
  ```
- It should not require manual notebook interaction.
- It should not depend on hidden kernel session state.

### Behavior contract

- Clear or localize symbols at the top.
- Print a stage banner and the key identities being checked.
- Evaluate at least one nontrivial symbolic check using Mathematica
  primitives such as:
  - `FullSimplify`
  - `Simplify`
  - `Reduce`
  - `Solve`
  - `Eliminate`
  - `Series`
  - `D`
  - `Integrate`
  - `Limit`
- Exit nonzero on failure.

### Output contract

- The corresponding saved output file should include:
  - script name,
  - date,
  - exit code,
  - `PASS` or `FAIL`,
  - the printed symbolic checks.

---

## Porting Strategy

Port in four phases.

### Phase 0: Scaffold

Create the minimal infrastructure:

- `research/pde_ledger/mathematica/`
- `research/pde_ledger/mathematica/output/`
- `research/pde_ledger/mathematica/run_all_audits.sh`
- one tiny helper layer only if repeated patterns appear

Deliverable:

- one reproducible runner that can execute a single stage or all currently
  ported stages and save outputs.

### Phase 1: Pilot High-Value Slice

Port a small set of high-value stages first:

- `058` coupled support-source operator
- `082` master quadrupole residual
- `083` Family-1 direct operator window
- `099` reduced finish line
- `106` canonical outgoing reduced closure
- `125` positive source theorem
- `134` Family-1 mouth fixedpoint
- `155` frozen traction fixedpoint
- `156` renormalized canonical branch
- `164` microscopic log channels
- `170` linear grouped outlet map
- `185` microscopic monomials
- `187` orbit quotient closure

Deliverable:

- a representative end-to-end Mathematica verification spine covering the
  main structural checkpoints and the final closure chain.

Success criterion:

- these scripts run cleanly,
- outputs are saved,
- no unresolved disagreement remains against the notes/SymPy audits.

### Phase 2: Previously Weak-Audit Slice

Port the stages that historically needed audit hardening:

- `022`, `023`, `024`
- `028`, `029`
- `033`, `034`, `035`, `036`
- `045`, `046`, `048`
- `064`
- `071`
- `135`

Deliverable:

- Mathematica confirms the places where we previously had the highest risk of
  tautological or assumption-fragile checking.

### Phase 3: Batch Expansion

Expand coverage batch by batch, following the existing review order in
`notes/moving_throat/review/REVIEW_PLAN.md`.

Suggested order after the pilot:

1. Batches 8-11 (`058-089`)
2. Batches 13-17 (`100-157`)
3. Batches 18-19 (`158-187`)
4. Early foundation batches (`003-057`)

Rationale:

- this keeps the first Mathematica coverage concentrated on the stages most
  likely to influence current reasoning and future work.

### Phase 4: Optional Full Coverage

Only after the previous phases are stable should we consider porting the
entire SymPy audit corpus.

At that point, decide whether the value is high enough to justify maintaining
near-total dual-engine coverage.

---

## How To Choose What To Port

A stage is a strong Mathematica candidate if it satisfies one or more of:

- it is a checkpoint or final-stage closure,
- it carries a theorem gate used by later stages,
- it contains a nontrivial solve/elimination/integral/limit,
- it previously needed audit hardening,
- it is numerically delicate or branch-sensitive,
- it compresses many earlier results.

A stage is lower priority if:

- it is purely a status/consolidation note with no independent computation,
- it is already trivial once adjacent stages are verified,
- its claims are fully subsumed by a stronger checkpoint script.

For status notes such as `084`, `149`, or `157`, prefer a Mathematica script
only if it independently verifies the consolidated formulas or quoted carried
values. Do not create a script that merely restates prose.

---

## Verification Policy For Each Ported Stage

For each stage we port:

1. Read the note.
2. Read the SymPy script.
3. Identify the minimum set of claims worth independently checking in
   Mathematica.
4. Write the `.wl` script from the note equations, not from SymPy variable
   choreography.
5. Run it with `math -script`.
6. Save the output under `research/pde_ledger/mathematica/output/`.
7. Compare the Mathematica result against:
   - the note,
   - the SymPy audit source,
   - the saved SymPy output.
8. If there is a discrepancy, classify it before fixing anything.

Discrepancy classes:

- `notes issue`
- `SymPy issue`
- `Mathematica issue`
- `assumption mismatch`
- `equivalent but differently normalized expression`

---

## Runner Plan

Create a shell runner analogous to the SymPy one:

```bash
bash research/pde_ledger/mathematica/run_all_audits.sh
bash research/pde_ledger/mathematica/run_all_audits.sh --force
bash research/pde_ledger/mathematica/run_all_audits.sh stage106
```

Responsibilities:

- discover `moving_throat_pde_stage*.wl`,
- run each via `math -script`,
- capture stdout/stderr,
- prepend a metadata header,
- write per-script outputs,
- maintain `output/_summary.txt`.

Important:

- use elevated execution when needed for `math -script`,
- but keep all writes inside the repo tree.

---

## Coding Standards For `.wl` Files

1. Use ASCII filenames and comments.
2. Keep symbols and stage names close to the notes where possible.
3. Prefer `Module`/`With`/local scoping over global symbol leakage.
4. Print both the expression being checked and the reduced residual.
5. Use helper assertions only if they stay transparent.
6. Avoid notebook-generated noise or front-end formatting constructs.

---

## First Concrete Implementation Slice

When implementation starts, do this exact sequence:

1. Create `research/pde_ledger/mathematica/run_all_audits.sh`.
2. Port `106`, `125`, and `170`.
3. Make sure all three run cleanly through the runner.
4. Port `058`, `082`, and `083`.
5. Port `134`, `155`, and `156`.
6. Port `164`, `185`, and `187`.

Why this slice:

- it gives one early reduced-closure chain,
- one positivity/source branch check,
- one grouped outlet map,
- one support/source operator chain,
- one mouth/core chain,
- and one final closure chain.

That is enough coverage to tell us whether the Mathematica effort is yielding
real independent confidence before we scale it.

---

## Milestones

### Milestone 1

Scaffold plus first 3 scripts:

- runner exists,
- outputs save correctly,
- `106`, `125`, `170` pass.

### Milestone 2

Pilot spine complete:

- all Phase 1 scripts exist,
- all Phase 1 outputs are saved,
- disagreements are triaged and resolved.

### Milestone 3

Historically weak-audit stages covered:

- all Phase 2 scripts exist,
- Mathematica independently confirms the previously hardened checks.

### Milestone 4

Decision point:

- either continue to full staged coverage,
- or stop with a high-value dual-engine verification layer.

---

## What Counts As Success

This project is successful if:

- Mathematica reproduces the main theorem gates and closure relations
  independently,
- it catches at least some real assumption/normalization risks or confirms
  their absence,
- the saved outputs become a usable second audit trail,
- and the maintenance burden stays reasonable.

It is **not** necessary for success that every SymPy script be ported.

---

## Recommendation

Proceed with the Mathematica layer.

But do it as a staged verification program:

- scaffold first,
- pilot high-value stages second,
- expand only after the pilot proves useful.

That gives the strongest confidence gain per unit of effort and avoids
building a second large symbolic codebase before we know which parts are
actually worth maintaining.
