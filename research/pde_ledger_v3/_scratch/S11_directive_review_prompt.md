# REVIEW — two build directives, before their scripts are trusted

**Read-only review.** ⛔ Do not edit any file. ⛔ Do not commit. Report findings as text.

## What you are reviewing

Two directives that an AI builder follows to produce the two engines for one derivation step:

1. `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11_wl_directive.md` — builds a
   **blind Mathematica** audit (imports nothing, reads no registry).
2. `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11_sympy_directive.md` — builds a
   **SymPy** audit that **does** import the project's shared registry, plus registry insertions.

## ⭐⭐ Why this review exists — the specific failure it must catch

The project's guard against a wrong result is that **two independently written engines must agree**. The
two directives are the **one thing both engines share**. An error in a directive propagates into **both**
engines, they agree, and the dual-engine check reports success on wrong physics. Nothing downstream can
catch that. **You are the only check on it.**

## ⛔ Do not read these

- `research/pde_ledger_v3/V3_STEP_PLAN.md`, and anything under `research/pde_ledger_v3/steps/`
- `research/pde_ledger_v2/scripts/ledger_stage003_*`,
  `research/pde_ledger_v2/mathematica/ledger_stage003_*`

These contain a related calculation and the step's own narrative. Reading them would make you check the
directives against a remembered answer instead of checking them on their merits — which is the whole
point of the exercise. Everything you need is in the two directives.

⚠ You will notice the directives do not state expected results anywhere. That is deliberate. ⛔ Do not
report it as an omission, and ⛔ do not supply expected results in your report.

---

## What to check

**C1 · Is the physical setup correct as stated?** Both directives define a Lagrangian density for an
in-plane displacement field on a `D`-dimensional elastic sheet, with a curl-squared term and a
divergence-squared term. Check the definitions themselves:
- Is `Curl2[u] ≡ ½ Σ_{i,j} (∂_i u_j − ∂_j u_i)²` a well-defined rotational invariant in **every** `D ≥ 2`,
  or does it silently assume `D = 3` (where the cross-product curl exists)?
- Is each term's placement and sign consistent with a Lagrangian `L = T − V` with a stable potential?
- Does the plane-wave/Euler–Lagrange route the directives ask for actually yield the eigenproblem they
  describe, or is there a step being skipped that would change the operator?

**C2 · Is the bulk-matching setup right?** Both directives state that an unbounded medium off the sheet
supports a scalar sound mode `ω² = c_s0²(k_in·k_in + k_w²)`, and that a sheet mode couples to it sharing
**both** `ω` and `k_in`. Check that dispersion relation and that matching condition. Are they the correct
statement of the problem, or is a boundary condition, a mode-conversion channel, or a geometric factor
missing that would change whether the answer comes out real or imaginary?

**C3 · Is the dimensional derivation well-posed?** Both ask the builder to derive each coefficient's
`[L,T,M]` exponent vector as a closed function of `D` from the single premise that the Lagrangian is an
energy density on a `D`-dimensional spatial brane. Is that premise sufficient to fix them all uniquely?
If any is underdetermined by it, say which and what else is needed.

**C4 · ⭐ Is every task ABLE TO FAIL?** For each numbered task in each directive, ask whether there is a
possible state of the physics in which it would report something other than what it reports for the
intended physics. ⛔ Flag any task whose phrasing admits only one answer regardless of the underlying
physics — that task is decoration, and worse, it will read as confirmation.
Pay particular attention to the two **control** tasks, which are meant to distinguish a control that
tests the *shape* of the answer from one that only tests arithmetic. Does each actually do what it is
meant to?

**C5 · ⭐⭐ DO THE TWO DIRECTIVES SPECIFY THE SAME PHYSICS?** Compare them term by term: the Lagrangian,
the invariant definitions, the bulk dispersion, the matching conditions, the controls. **Any divergence
is a blocking finding** — if the two engines are set different problems, their agreement or disagreement
means nothing. Also flag anything one asks for that the other should ask for and does not.

**C6 · Does either directive LEAK an answer?** Check for any expected value, target, comparison, or
phrasing that tells the builder what to find rather than what to compute. Include subtle forms: a task
that names the property it expects to hold, a variable name that encodes an outcome, an ordering that
implies a result.

**C7 · Is anything MISSING** that the step needs in order to be trustworthy, or that one engine needs in
order to be a genuine check on the other?

---

## ⭐⭐ THE FILTER — apply it before you write anything down

Report a finding **only if it catches a way the PHYSICS or the SPECIFICATION could be wrong.**

⛔ Do **not** report: style, formatting, wording, file organisation, "this script would be wrong on a
different input", missing error handling, or process suggestions.

⚠ This project has repeatedly lost days to reviews that generated many findings of which nearly none
were about the physics. A short report with one real finding is worth far more than a long one. If the
directives are sound, **say so plainly and stop** — "no blocking findings" is a complete and welcome
result.

## Output format

```
VERDICT: SOUND | BLOCKING FINDINGS
BLOCKING   (physics or specification is wrong — list, or "none")
DIVERGENCE (C5: the two directives disagree — list, or "none")
LEAKS      (C6 — list, or "none")
NOT-ABLE-TO-FAIL (C4 — list, or "none")
MISSING    (C7 — list, or "none")
WHAT IS SOUND (brief — what you checked and found correct)
```
