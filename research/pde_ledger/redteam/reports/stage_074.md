---
unit_id: 074
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md
  paper_appendix: present
---

# Audit unit 074 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_074.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`

## What the paper claims

The paper card (`paper/stages/stage_074.tex`) declares the Family--1 healing-length lock with `Inputs` of `L/ell = 37` and closure `ell = hbar / (2 m c_{s,w})`. The `Output` line states the deliverables verbatim: "Family--1 support values `chi_s = 37/2`, `kappa = 12321/5`, `eta = 37`." The body of the card derives `chi_s = m c_{s,w} L / hbar = 37/2` via the healing-length substitution and gives the closed forms `kappa = (9/5) Lambda_ell^2 = 12321/5` and `alpha = sqrt(kappa) = 128/sqrt(5)` in a boxed equation. The notes (`notes/stages/moving_throat_pde_stage074_family1_healing_lock.md`) lay out the derivation: identify the active wall width with the GNLS healing length, deduce `chi_s = Lambda_ell/2`, then evaluate `kappa = (1 + 4/5) Lambda_ell^2 = 9/5 * 37^2 = 12321/5` and `alpha = sqrt(12321/5) ~ 49.6407091`. The appendix row (table line for stage 074) summarises: "`ell = hbar/(2 m c_{s,w})`, `chi_s = 37/2`, `kappa = 12321/5`."

## What the script claims to verify

The (now-updated) SymPy script (`scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`) introduces positive symbols `Lambda_ell, hbar, m_psi, c_s, ell, L`, builds the physical definition `chi_def = m_psi * c_s * L / hbar`, substitutes `c_s -> hbar/(2*m_psi*ell)` to obtain `chi_in_ell = L/(2*ell)`, and then substitutes `L -> Lambda_ell * ell` to obtain `chi_lock = Lambda_ell/2`. It then forms `kappa_lock = 4 * chi_lock**2 + (4/5) Lambda_ell^2` (with an inline comment anchoring the `4` and `4/5` to the Family-1 Euler-Lagrange branch carry-forward) and asserts `chi_lock - Lambda_ell/2 == 0`, `kappa_lock - (9/5) Lambda_ell^2 == 0`. Reference-branch assertions specialize `Lambda_ell = 37` and check `chi_ref - 37/2 == 0`, `kappa_ref - 12321/5 == 0`. `alpha = sqrt(kappa_ref)` is computed and printed (output: `111*sqrt(5)/5 ~ 49.6407091...`) but not asserted to a literal. The Mathematica script mirrors the substitution chain.

## Paper -- script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `ell = hbar/(2 m c_{s,w})` (healing-length closure) | sympy line 37, mathematica line 33: `c_s -> hbar/(2 m_psi ell)` substitution | match |
| `chi_s = 37/2` (Output line) | sympy line 53/68 `chi_lock - Lambda_ell/2`, `chi_ref - 37/2`; mathematica line 42/57 | match |
| `kappa = 12321/5` (Output line) | sympy line 54/69 `kappa_lock - (9/5) Lambda_ell^2`, `kappa_ref - 12321/5`; mathematica line 43/58 | match |
| `eta = 37` (Output line) | Carry-forward from Stage 56/073 (`Lambda_ell = 37` referenced and used as the numeric specialization). Not separately asserted in this unit but is the same scalar as `Lambda_ell`. | partial (carry-forward via `Lambda_ref = 37`; no separate `eta` assertion) |
| `alpha = sqrt(kappa) = 128/sqrt(5)` (boxed body equation, paper line 31) | sympy line 59-66 computes/prints `alpha_ref = 111*sqrt(5)/5 = 111/sqrt(5) ~ 49.6407091`; mathematica line 48-55 prints `alpha = 111/Sqrt[5]` | mismatch (paper literal `128/sqrt(5)` is arithmetically wrong; `sqrt(12321/5) = 111/sqrt(5)` since `111^2 = 12321`. The scripts produce the arithmetically correct value but do not assert it; the paper card states a different literal that disagrees with its own `kappa`.) |
| `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2` (carry-forward from Stage 54) | sympy line 48 / mathematica line 38: literal `4` and `4/5` with inline comment provenance | partial (provenance is a comment, not a script-internal derivation; acceptable as a carry-forward) |

Dominant pattern: most deliverables match; one paper-side literal (`alpha = 128/sqrt(5)`) is internally inconsistent with the paper's own `kappa = 12321/5`. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero("chi_s - Lambda_ell/2", chi_lock - Lambda_ell/2)` where `chi_lock = (m c_s L / hbar) /. c_s -> hbar/(2 m ell) /. L -> Lambda_ell * ell` | `chi_s = 37/2` (general form) | yes |
| A2 | sympy | 54 | `expect_zero("kappa - (9/5) Lambda_ell^2", kappa_lock - (9/5) Lambda_ell^2)` | `kappa = 12321/5` (general form) | yes |
| A3 | sympy | 68 | `expect_zero("chi_ref - 37/2", chi_ref - 37/2)` | `chi_s = 37/2` (numeric) | yes (numeric specialization of A1) |
| A4 | sympy | 69 | `expect_zero("kappa_ref - 12321/5", kappa_ref - 12321/5)` | `kappa = 12321/5` (numeric) | yes |
| A5 | mathematica | 42 | `expectZero[chi_s - Lambda_ell/2, ...]` with `chiLock` from `m c_s L / hbar` -> healing-length substitution -> `L/ell -> Lambda_ell` | same as A1 | yes |
| A6 | mathematica | 43 | `expectZero[kappa - (9/5) Lambda_ell^2, ...]` | same as A2 | yes |
| A7 | mathematica | 57 | `expectZero[chi_ref - 37/2]` | same as A3 | yes |
| A8 | mathematica | 58 | `expectZero[kappa_ref - 12321/5]` | same as A4 | yes |
| - | sympy/mathematica | 59-66 / 48-55 | `alpha_ref` is computed and printed but never asserted to a literal value | `alpha = 128/sqrt(5)` (paper body equation) | no (no assertion; print only) |

The `alpha` quantity from the boxed body equation has no script-side assertion. The scripts print the arithmetically correct value `111/sqrt(5)` while the paper boxes `128/sqrt(5)`; the lack of an assertion means the inconsistency went undetected by the existing pipeline.

## Findings

### F1 -- paper_misalignment

**Subtype:** value_mismatch
**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_074.tex:25-32` (paper-side)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage074_family1_healing_lock.md:113-117` (notes-side, related typo)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:59-66` (script-side)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.wl:48-55` (script-side)

**What's wrong:**

The paper card boxes the following equation (paper line 25-31):

    kappa = (9/5) Lambda_ell^2 = 12321/5,    alpha = sqrt(kappa) = 128/sqrt(5)

But `sqrt(12321/5) = sqrt(12321)/sqrt(5) = 111/sqrt(5)`, because `111^2 = 12321`. The literal `128` is incorrect: `128^2 = 16384`, so `(128/sqrt(5))^2 = 16384/5 != 12321/5`. The paper card therefore contains an internal inconsistency: it asserts `kappa = 12321/5` and `alpha = sqrt(kappa)` and `alpha = 128/sqrt(5)`, and these three statements cannot all be true simultaneously.

The SymPy script computes (line 59) `alpha_ref = sp.simplify(sp.sqrt(kappa_ref))` and prints `111*sqrt(5)/5` (= `111/sqrt(5)`). The Mathematica script computes (line 48) `alphaRef = FullSimplify[Sqrt[kappaRef]]` and prints `111/Sqrt[5]`. Both engines produce the arithmetically correct value (`111/sqrt(5) ~ 49.6407091`), in agreement with each other and with the numerical value `~ 49.6407091` that the notes also state.

The notes file (`notes/stages/moving_throat_pde_stage074_family1_healing_lock.md`, section 5) state:

    alpha = sqrt(12321/5) = 179/sqrt(5) ~ 49.6407091

The `179/sqrt(5)` literal is also wrong (`179^2 = 32041`); however, the numerical value `49.6407091` in the notes does match `111/sqrt(5)`. Both prose docs have inconsistent literal values for alpha; the scripts have the arithmetically correct value.

Quotes:

paper/stages/stage_074.tex lines 26-31:
> `\kappa=\frac95\Lambda_\ell^2=\frac{12321}{5}, \qquad \alpha=\sqrt\kappa=\frac{128}{\sqrt5}`

notes section 5 (lines 113-117):
> `alpha = sqrt(12321/5) = 179/sqrt(5) ~ 49.6407091.`

script output (sympy and mathematica):
> `alpha = 111*sqrt(5)/5` / `alpha = 111/Sqrt[5]`, `alpha (numeric) = 49.640709100495331260`

**Why this matters:**

The boxed `alpha = sqrt(kappa)` value is a load-bearing carry-forward into Stage 41/42-style kernel formulas (the notes explicitly call out "alpha is the support-decay scale entering the exact Stage-41/42 kernel formulas"). A reader who quotes the paper's `128/sqrt(5)` and inserts it into a downstream kernel formula will produce a numerically wrong result (`~57.24` instead of `~49.64`). The fact that both engines agree on `111/sqrt(5)` and the numerical value in the notes agrees with `111/sqrt(5)` makes it overwhelmingly likely the paper's `128` and the notes' `179` are typos for `111`. Resolution direction is the user's call (per audit prompt rules for paper_misalignment).

**Required change:**

This is a paper_misalignment finding. Per audit-framework policy, Codex does not auto-edit paper.tex or notes/. The orchestrator must surface the question to the user before any change is applied. See the `## Resolve before fix_loop` block in the directive.

**Verification:**

Once the user picks a direction:
- If the user confirms the paper should match the scripts: paper line 31 `\alpha=\sqrt\kappa=\frac{128}{\sqrt5}` becomes `\alpha=\sqrt\kappa=\frac{111}{\sqrt5}`, and notes section 5 line `alpha = sqrt(12321/5) = 179/sqrt(5)` becomes `alpha = sqrt(12321/5) = 111/sqrt(5)`. No script change. The verifier then re-reads the paper card and confirms the literal matches `111*sqrt(5)/5` in the script output (these forms are equal: `111/sqrt(5) = 111*sqrt(5)/5`).
- If the user chooses to add an explicit script-side assertion `alpha_ref == 111*sqrt(5)/5` to harden the check (optional follow-up), the verifier runs the script and confirms the assertion passes.

## Independent-derivation check (Mathematica)

The Mathematica `.wl` is now structurally parallel to the (updated) SymPy script, but the parallel emerged because both scripts independently chose the same physically meaningful derivation chain: build `chi_s` from `m c_s L / hbar`, substitute `c_s -> hbar/(2 m ell)`, substitute `L/ell -> Lambda_ell`, form `kappa = 4 chi^2 + (4/5) Lambda^2`. The Mathematica script (lines 33-38) and the SymPy script (lines 33-48) share the substitution choreography because the underlying physics has exactly one natural derivation path here, not because one was line-by-line ported from the other. The Mathematica script also has its own assumption-block (`$Assumptions` declarations at lines 30-31) that the SymPy script does not mirror, and uses `FullSimplify` + `Together[Expand[...]]` chaining that is genuinely Mathematica-idiomatic rather than a SymPy port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on all assertions and printed values:

| Quantity | SymPy | Mathematica |
|---|---|---|
| chi (after healing-length substitution) | `L/(2*ell)` | (same, intermediate, implicit) |
| chi_s (locked) | `Lambda_ell/2` | `lambdaEll/2` |
| kappa (locked) | `9*Lambda_ell**2/5` | `(9*lambdaEll^2)/5` |
| chi_ref | `37/2` | `37/2` |
| kappa_ref | `12321/5` | `12321/5` |
| alpha_ref (symbolic) | `111*sqrt(5)/5` | `111/Sqrt[5]` |
| alpha_ref (numeric) | `49.640709100495331260` | `49.64070910049533126028365544...` |

`111/sqrt(5) = 111*sqrt(5)/5`; the symbolic forms are equal, the numerics match to 20 digits. No `engine_disagreement`.

## Verdict justification

The math itself holds up: with `ell = hbar/(2 m c_{s,w})`, the dimensionless support scale `chi_s = m c_{s,w} L / hbar` reduces to `Lambda_ell/2`, and with `Lambda_ell = 37` and the Family-1 EL-branch carry-forward `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`, one gets `chi_s = 37/2`, `kappa = 12321/5`. Both engines now derive these non-tautologically (the v1 SymPy tautology has been fixed; v1 F1 and F2 are resolved). The remaining v2 finding is a paper-side internal inconsistency: the paper card's `alpha = sqrt(kappa)` is stated as `128/sqrt(5)`, but the only value consistent with the paper's own `kappa = 12321/5` is `111/sqrt(5)`, which the scripts compute correctly. The notes have a similar typo (`179/sqrt(5)`) but with the numerically correct decimal. No `stop_cold`: the discrepancy is local to the paper/notes prose and does not change any constant the scripts use or any downstream carry-forward (the scripts compute alpha correctly and stage 075 will inherit the script value, not the paper literal). The user must choose whether to correct the paper or the notes; neither auto-fix is appropriate.

## Self-test notes

I checked: (1) Variable independence not applicable -- this unit has no `sp.diff` or `D[...]` operations. (2) Symmetry/parity not applicable -- no integrals on symmetric domains. (3) Trivial-case substitution: manually verified `sqrt(12321/5) = sqrt(12321)/sqrt(5) = 111/sqrt(5)` since `111^2 = 12321` (`100^2 + 2*100*11 + 11^2 = 10000 + 2200 + 121 = 12321`), confirming `111/sqrt(5)` is correct and `128/sqrt(5)` and `179/sqrt(5)` are both wrong as literals (though the notes' numerical `49.6407091` is right). (4) Path specifications n/a -- no missing-script findings. (5) Paper round-trip: the finding is paper-side; no script edits are prescribed, so there is no risk of introducing a fresh `paper_misalignment` from script changes. Outputs are fresh (script-py mtime 1779513065 < output txt mtime 1779513120; script-wl mtime 1778522212 < output txt mtime 1779513124).
