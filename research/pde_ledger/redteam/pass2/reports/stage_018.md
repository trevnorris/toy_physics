---
unit_id: 018
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 018 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_018.tex`
- notes: `(none)` — no `notes/stages/moving_throat_pde_stage018_*.md` exists
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 58)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.txt`

## What the paper claims

Stage 018 ("Parent throat action bundle master") collects the parent-wall inputs into the grouped-bundle language and exports a two-equation "parent-wall bundle bridge." The card's `\stagefield{Output}` reads: "Stage~018 exports the parent-wall bundle bridge \eqref{eq:stage018-msigma}--\eqref{eq:stage018-ksigma}." Those two branch integrals are the wall inertia `M_\Sigma = \int dw\,\mu_\eta\beta_2^2` (eq:stage018-msigma, lines 14-17) and the wall stiffness `K_\Sigma = \int dw\,[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2]` (eq:stage018-ksigma, lines 18-22). The card states these "replace abstract wall knobs in the grouped bundle." The appendix row 58 summarizes the stage as "Bundle-level parent-action identities used by the projected electromagnetic response," status `\StatusExactClosure{}`. The card explicitly claims an algebraic-identity status only "inside the declared parent/projection data" (claimstatus, line 5) — i.e. it asserts the *form* of the two branch integrals, not numeric values of the coefficient fields `\mu_\eta, T_w, K_\eta, T_\Omega`.

## What the script claims to verify

Both scripts verify a single concrete Gaussian realization of the two branch-integral structures. With the parent-wall profile `\beta = \exp(-w^2/2)` and the coefficient fields set to unity (`\mu_\eta = T_w = K_\eta + 6T_\Omega = 1`), they evaluate the definite integrals and compare to closed forms: the inertia integral `\int_{-\infty}^{\infty}\beta^2\,dw = \sqrt\pi`, and the stiffness integral `\int_{-\infty}^{\infty}[(\beta')^2 + \beta^2]\,dw = 3\sqrt\pi/2`. The SymPy print (`output .txt:2`) and the WL comment ("M1 Gaussian wall inertia and stiffness integrals," `.wl:4`) both describe this as a concrete-example sanity check of the inertia and stiffness branch integrals — the script does not claim to verify the fully general symbolic identity with all coefficient fields retained.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `M_\Sigma = \int dw\,\mu_\eta\beta_2^2` (inertia branch integral, eq:stage018-msigma) | SymPy line 23/25 + WL line 17/19-22: `\int\beta^2 = \sqrt\pi` with `\beta=e^{-w^2/2}`, `\mu_\eta\to1` | partial (structure of inertia integral exercised on one concrete profile; coefficient field `\mu_\eta` set to unity) |
| `K_\Sigma = \int dw\,[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2]` (stiffness branch integral, eq:stage018-ksigma) | SymPy line 24/26 + WL line 18/24-27: `\int[(\beta')^2+\beta^2] = 3\sqrt\pi/2`, `T_w\to1`, `(K_\eta+6T_\Omega)\to1` | partial (two-term stiffness structure `(\beta')^2 + \beta^2` exercised on one concrete profile; coefficient fields set to unity) |

The script faithfully exercises the *structure* of both branch integrals declared in the card (an inertia term `\propto\int\beta^2` and a two-term stiffness `\propto\int(\beta')^2 + \int\beta^2`) on a concrete Gaussian, with coefficient fields collapsed to 1. The card's `claimstatus` scopes its claim to "displayed algebraic identities inside the declared parent/projection data," and the script's own output labels itself a concrete-example check — so the concrete-profile scope is mutually consistent and disclosed on both sides. No deliverable is unaddressed; nothing in the script lacks a paper counterpart. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 25 | `assert_zero(MSigma_example - sqrt(pi))` where `MSigma_example = integrate(beta**2,(w,-oo,oo))` | M_\Sigma inertia branch integral (eq:stage018-msigma), concrete Gaussian | yes (genuine integral evaluation vs independent closed form) |
| A2 | sympy | 26 | `assert_zero(KSigma_example - 3*sqrt(pi)/2)` where `KSigma_example = integrate(diff(beta,w)**2 + beta**2,(w,-oo,oo))` | K_\Sigma stiffness branch integral (eq:stage018-ksigma), concrete Gaussian | yes |
| A3 | mathematica | 21-22 | `If[FullSimplify[massIntegral - Sqrt[Pi]] =!= 0, Exit[1]]`; `massIntegral = Integrate[beta^2,...]` | M_\Sigma inertia branch integral, concrete Gaussian | yes |
| A4 | mathematica | 26-27 | `If[FullSimplify[stiffnessIntegral - 3*Sqrt[Pi]/2] =!= 0, Exit[1]]`; `stiffnessIntegral = Integrate[D[beta,w]^2+beta^2,...]` | K_\Sigma stiffness branch integral, concrete Gaussian | yes |

All four assertions are non-tautological: each independently evaluates a definite integral with the engine's native integrator and compares against a closed form that is *not* defined as that integral. If the integrand structure or the closed form were wrong, the residual would be nonzero and the script would fail (`raise AssertionError` / `Exit[1]`).

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration of the SymPy algebra in any defect sense. The stage's claim reduces to evaluating two definite Gaussian integrals; there is no multi-step symbolic choreography to copy. Each engine evaluates the integrals through its own, entirely distinct backend:

- SymPy (`.py:23-24`): `sp.integrate(beta**2, (w, -sp.oo, sp.oo))` and `sp.integrate(sp.diff(beta, w)**2 + beta**2, (w, -sp.oo, sp.oo))` — SymPy's Risch/Meijer-G definite-integration path.
- Mathematica (`.wl:17-18`): `Integrate[beta^2, {w, -Infinity, Infinity}]` and `Integrate[D[beta, w]^2 + beta^2, {w, -Infinity, Infinity}]` — Wolfram's independent definite-integration engine.

Both then subtract the *same* hand-stated closed forms (`Sqrt[Pi]`, `3 Sqrt[Pi]/2`). Sharing the target closed form is correct here — it is the claim being verified, not borrowed algebra. The two integral evaluations are genuinely independent. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines assert the identical residuals and both report zero:

- SymPy output (`.txt:3`): `STATUS: PASS` (assertions at lines 25-26 did not raise).
- Mathematica output (`.txt:2-4`): `M1 Gaussian inertia integral residual = 0`, `M1 Gaussian stiffness integral residual = 0`, `STAGE 018 MATHEMATICA AUDIT PASS`.

Independent hand-check confirms both: for `\beta=e^{-w^2/2}`, `\int e^{-w^2}dw = \sqrt\pi` (inertia ✓); and `(\beta')^2 = w^2 e^{-w^2}` gives `\int w^2 e^{-w^2}dw = \sqrt\pi/2`, plus `\int\beta^2 = \sqrt\pi`, summing to `3\sqrt\pi/2` (stiffness ✓). Engines agree; both agree with the analytic result. `engines_agree: true`.

## Verdict justification

`clean`. I read the paper card, confirmed there is no notes file, read the appendix row 58, then attacked both scripts. The two paper deliverables — the inertia branch integral (eq:stage018-msigma) and the two-term stiffness branch integral (eq:stage018-ksigma) — each have a faithful, non-tautological script-side check that genuinely evaluates the corresponding integral structure on a concrete Gaussian and compares to an independently-stated closed form. Attacks tried and failed: (1) tautology — no, the integrals are computed, not assumed; (2) wrong closed form — no, `\sqrt\pi` and `3\sqrt\pi/2` are both correct by hand; (3) symbol-domain error — no, `w` real, Gaussian integrals converge, no positivity assumption needed; (4) transliteration — no, each engine integrates natively; (5) hardcoded result — the closed forms are the verification targets, legitimately stated, and the LHS integrals are derived. The one scope limitation (coefficient fields `\mu_\eta, T_w, K_\eta, T_\Omega` collapsed to unity, single Gaussian profile) is disclosed by the script's own output label and is consistent with the card's `claimstatus` scoping the claim to "displayed algebraic identities inside the declared parent/projection data" — so it is not a `paper_misalignment` or `insufficient_verification` defect, just an honest concrete-example confirmation of the exported bundle bridge.

Minor cosmetic note (not a finding, pre-existing known drift): the SymPy docstring (`.py:2`) and its print banner (`.py:28`, output `.txt:1`) say "step_16"/"STEP 16" while this is stage 018 — the known +N pre-renumber label drift documented in project memory. It is a label-only artifact in prose strings, does not touch any math or assertion, and the WL banner correctly says "STAGE 018." No script-side math finding.

## Self-test notes

I checked: (1) variable independence — both `sp.diff(beta, w)` and `D[beta, w]` differentiate w.r.t. the genuine integration variable `w` on which `beta` depends, so no identically-zero-derivative trap. (2) Parity/symmetry on the unbounded domain — `beta^2 = e^{-w^2}` is even (integral nonzero, `\sqrt\pi`); `(\beta')^2 = w^2 e^{-w^2}` is even (integral nonzero, `\sqrt\pi/2`); the assertions claim nonzero closed forms, consistent with even integrands, no false "this vanishes" claim. (3) Trivial-case substitution — hand-evaluating both Gaussian integrals reproduces `\sqrt\pi` and `3\sqrt\pi/2` exactly, matching both saved outputs. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 2 deliverable-level values checked, 0 misaligned.

The two scripts emit exactly two result values, both of which are the closed forms of the *concrete Gaussian example* used to sanity-check the exported branch-integral structures. These example values are verification scaffolding (the sanity-check targets), not stage deliverables that the card is meant to report — the card's deliverables are the two *symbolic* branch-integral forms, which are present and correct in the `.tex`. I list both example values plus the symbolic deliverables for completeness:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_\Sigma = \int dw\,\mu_\eta\beta_2^2` (symbolic inertia branch integral — the actual deliverable) | structure exercised by py:23 / wl:17; deliverable form | `paper/stages/stage_018.tex:14-17` (eq:stage018-msigma) | MATCH |
| `K_\Sigma = \int dw\,[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2]` (symbolic stiffness branch integral — the actual deliverable) | structure exercised by py:24 / wl:18; deliverable form | `paper/stages/stage_018.tex:18-22` (eq:stage018-ksigma) | MATCH |
| `\sqrt\pi` (Gaussian inertia example numeric value) | py:25 `MSigma_example - sp.sqrt(sp.pi)`; sympy out:3 PASS; wl:21 / math out:2 `residual = 0` | not in `.tex`/`.md` (no notes file) | INTERNAL (concrete-example sanity-check target, not a card deliverable) |
| `3\sqrt\pi/2` (Gaussian stiffness example numeric value) | py:26 `KSigma_example - 3*sp.sqrt(sp.pi)/2`; sympy out:3 PASS; wl:26 / math out:3 `residual = 0` | not in `.tex`/`.md` (no notes file) | INTERNAL (concrete-example sanity-check target, not a card deliverable) |

INTERNAL items (accounted for, no finding): the Gaussian profile `\beta = e^{-w^2/2}`; the intermediate `MSigma_example`/`massIntegral` and `KSigma_example`/`stiffnessIntegral` objects; the pass/fail flags and `residual = 0` strings; the unity collapse of the coefficient fields `\mu_\eta, T_w, K_\eta, T_\Omega`.

Both stage deliverables (the symbolic branch-integral bridge) reconcile to the `.tex` card verbatim. The two emitted numeric values are example sanity-check targets, legitimately absent from the terse card (and there is no notes file), so per the augmentation guards they are MATCH/INTERNAL, not MISSING. Zero MISMATCH, zero MISSING-DELIVERABLE. No `paper_misalignment` raised from reconciliation.
