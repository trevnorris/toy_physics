# Moving-Throat PDE — Proof Hardening Plan

## Purpose

This document defines the next verification pass after the full stage-by-stage
review and the complete SymPy/Mathematica dual-audit buildout.

The current state is strong:

- derivation review is complete through Stage `187`,
- the SymPy audit stack has been reviewed and cleaned up,
- the Mathematica verification layer now covers every SymPy stage audit,
- the current Mathematica summary is `167/167 PASS`.

That does **not** mean the proof is beyond further improvement.
It means the dominant remaining risk is no longer algebra-engine failure.
The dominant remaining risk is now **shared conceptual error**:

- a wrong premise copied into both script stacks,
- an unstated branch assumption,
- a theorem statement slightly stronger than what the derivation proves,
- a hidden dependency where later stages use more hypotheses than earlier stages recorded,
- or a boundary/degenerate regime that was never stress-tested.

This plan is for hardening against those risks.

---

## Current Verification Baseline

As of `2026-04-04`:

- review coverage: foundations plus Stages `001-187`
- SymPy audit coverage: complete for all existing audit stages
- Mathematica audit coverage: complete for all SymPy audit stages
- ledger status: no known blocker in the present derivation chain

What this baseline already gives us:

- strong protection against sign errors, dropped factors, stale formulas, and weak symbolic checks
- independent confirmation of the main algebraic identities in two CAS systems
- much better confidence in the derivation as an internal theorem/program document

What it does **not** fully guarantee:

- that every theorem has exactly the right hypotheses
- that every later stage uses earlier results under valid assumptions
- that branch/limit statements are complete at the edges of admissible parameter space
- that no shared modeling assumption is wrong in both script stacks

---

## Risk Model

The remaining proof-risk categories are:

1. Hypothesis drift
   A stage proves a statement under assumptions `A`, but a later stage uses it as if it held under weaker assumptions `B`.

2. Silent branch selection
   A computation implicitly chooses a positive root, a stable branch, a nonvanishing denominator, or a one-sided limit without recording it clearly.

3. Statement inflation
   Notes sometimes describe a result as exact, global, unique, or constructive when the derivation only proves a weaker claim.

4. Shared premise error
   SymPy and Mathematica both verify the same wrong reduced formula because that wrong formula was copied into both.

5. Boundary fragility
   Identities may hold generically but become false, ambiguous, or under-justified on threshold surfaces or degenerate limits.

6. Dependency opacity
   The chain is long enough that a correct local stage can still be globally fragile if its inputs are not tracked explicitly.

---

## Workstreams

## 1. Assumption Ledger

### Goal
Produce a stage-by-stage ledger of every assumption actually used.

### For each stage, record

- positivity assumptions
- reality assumptions
- nonvanishing denominators
- discriminant nonnegativity assumptions
- branch/root choices
- limit directions
- perturbative small-parameter assumptions
- asymptotic regime assumptions
- monotonicity-domain assumptions

### Deliverable

- `notes/moving_throat/review/ASSUMPTION_LEDGER.md`

### Acceptance criteria

- every theorem/proposition stage has an explicit hypothesis block
- every script-specific simplification assumption is recorded if it matters mathematically
- every root/branch choice is identified and justified
- every later use of a stage result can be checked against the ledger

### Why this matters

This is the single most likely place the current dual-CAS setup can still miss a real issue.

---

## 2. Dependency Ledger

### Goal
Make the proof chain explicit enough that we can see exactly what each stage consumes and what it exports.

### For each stage, record

- direct inputs: prior stages actually used
- local assumptions
- exact outputs: formulas/theorems produced
- output type:
  - exact identity
  - existence statement
  - asymptotic statement
  - perturbative statement
  - status/consolidation summary

### Deliverable

- `notes/moving_throat/review/DEPENDENCY_LEDGER.md`

### Acceptance criteria

- every non-status stage has a concise dependency entry
- every checkpoint/status stage cites the precise stages it summarizes
- there is no silent leap where a stage uses an object not previously defined or justified

### Why this matters

The derivation is now long enough that global correctness depends on explicit dependency control, not just local symbolic correctness.

---

## 3. Theorem Statement Tightening

### Goal
Audit the notes for places where the prose may overstate what the derivation actually proves.

### For each theorem-style stage, check

- whether the conclusion is exact or only perturbative
- whether uniqueness is actually proved
- whether existence is conditional or unconditional
- whether positivity/nonnegativity is global or only on a branch
- whether a bound is attained, strict, one-sided, or only a supremum/infimum
- whether "if and only if" is really proved in both directions

### Deliverable

- direct edits to the stage notes
- optional summary file:
  - `notes/moving_throat/review/THEOREM_TIGHTENING_NOTES.md`

### Acceptance criteria

- no theorem statement is stronger than its proof
- no asymptotic statement is labeled exact
- no branch-local result is phrased globally without the branch

### Why this matters

At this stage, prose drift is a more realistic failure mode than raw algebra drift.

---

## 4. Adversarial Numerical Validation

### Goal
Stress-test the critical formulas at representative and near-boundary parameter points.

### Required test regions

- zero-load limits
- threshold surfaces
- near-degenerate discriminants
- small/large parameter asymptotics
- denominator pinch surfaces
- branch-transition neighborhoods
- isotropic vs anisotropic specializations

### Required comparisons

- note formula vs direct residual substitution
- SymPy output vs Mathematica output
- exact formula vs numerical solve when applicable
- asymptotic approximation vs exact expression in its intended regime

### Deliverables

- `research/pde_ledger/scripts/numerical/` for Python numerical spot checks
- `research/pde_ledger/mathematica/numerical/` for Mathematica spot checks
- `notes/moving_throat/review/NUMERICAL_STRESS_PLAN.md`

### Acceptance criteria

- every choke-point stage has at least a small adversarial sample set
- every sample records assumptions and whether it is near a boundary
- discrepancies, if any, are triaged into:
  - script bug
  - note bug
  - domain/assumption issue
  - numerical instability only

### Why this matters

This is the best remaining defense against branch and boundary mistakes.

---

## 5. Choke-Point Re-Derivation

### Goal
Do not re-review all `187` stages uniformly again. Re-derive the load-bearing joints from scratch enough to test whether the global chain would fail if one of them were wrong.

### Priority stages

- `003-004`
- `037-048`
- `058`
- `106`
- `125`
- `142-157`
- `170`
- `185-187`

### What "re-derive" means here

- start from the stated inputs, not from the final target formula
- avoid reusing the exact same intermediate parameterization when possible
- verify the theorem through a different route if practical
- confirm the stated hypotheses are sufficient

### Deliverable

- `notes/moving_throat/review/CHOKE_POINT_REDERIVATION_PLAN.md`
- per-stage notes appended to the existing review files or stored in a dedicated addendum file

### Acceptance criteria

- each choke point has an explicit independent argument path
- each choke point has a short statement of what would break downstream if it were wrong
- each choke point has a signed-off status after re-derivation

### Why this matters

The proof does not need another uniform pass. It needs concentrated scrutiny where downstream propagation risk is highest.

---

## 6. Script and Output Regression Harness

### Goal
Turn the current audit stacks into a reproducible regression system.

### Required components

- one command to run all SymPy audits in a controlled way
- one command to run all Mathematica audits sequentially
- one command to rebuild summaries
- output diff checks against saved artifacts

### Deliverables

- a lightweight regression script or Make-style entry point
- documented run order and licensing note for Mathematica
- optional CI config if the environment permits it

### Acceptance criteria

- drift in saved audit outputs is intentional and reviewable
- it is easy to rerun the full verification spine after edits
- the proof no longer depends on manual memory of which scripts need to be rerun

### Why this matters

Once the notes continue evolving, regression control becomes part of proof safety.

---

## 7. Status/Consolidation Stage Audit

### Goal
Re-check the non-script status stages specifically as theorem-summary documents rather than symbolic derivations.

### Focus stages

- `049`
- `084`
- `090`
- `093`
- `096`
- `103`
- `113`
- `117`
- `120`
- `128`
- `132`
- `136`
- `141`
- `145`
- `149`
- `153`
- `157`

### What to check

- whether they summarize prior stages exactly
- whether they introduce new claims not present earlier
- whether any theorem wording is stronger than the underlying stage support
- whether cited stage numbers and filenames are correct

### Deliverable

- direct note cleanup where needed

### Why this matters

Consolidation notes are where local truths often become globally overstated statements.

---

## Suggested Execution Order

The recommended order is:

1. Assumption ledger
2. Dependency ledger
3. Theorem statement tightening
4. Choke-point re-derivation
5. Adversarial numerical validation
6. Regression harness
7. Final status/consolidation audit

Reason:

- steps `1-3` clarify what the proof actually claims
- step `4` tests the most dangerous joints
- step `5` stresses those joints near boundaries
- step `6` locks in the improved verification workflow
- step `7` ensures the summary documents do not overclaim

---

## Practical First Batch

If we want to start immediately with the highest-value work, the first batch should be:

1. Build `ASSUMPTION_LEDGER.md`
   Start with stages:
   - `003-004`
   - `037-048`
   - `058`
   - `106`
   - `125`
   - `142-157`
   - `170`
   - `185-187`

2. Build `DEPENDENCY_LEDGER.md`
   At minimum for the same choke-point block.

3. Tighten theorem wording in those same stages before doing any more symbolic work.

This gives the biggest reduction in residual proof-risk per unit effort.

---

## Definition Of Done

This hardening pass is complete when all of the following are true:

- every theorem stage has an explicit hypothesis record
- every choke-point stage has a dependency entry and an independent re-derivation note
- every theorem statement has been checked for exactness/branch scope
- critical boundary regimes have numerical adversarial coverage
- the SymPy and Mathematica stacks are runnable as a regression harness
- status/consolidation notes have been audited for statement drift
- no open blocker or unresolved `ISSUE` remains in the hardening tracker

At that point the proof will still not be a refereed external publication, but it will be in a materially stronger internal state than "dual-CAS checked." It will be **hypothesis-tracked, dependency-tracked, branch-audited, numerically stressed, and regression-protected**.

---

## Working Rule

From this point onward, every new stage should ideally add:

- a short explicit hypothesis block in the note,
- a clear statement of exact vs asymptotic status,
- an entry in the dependency ledger,
- and an audit script that derives rather than merely restates the claimed formula.

That prevents us from needing another large retrospective cleanup later.
