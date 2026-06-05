---
unit_id: 078
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md]
  paper_appendix: present
---

# Audit unit 078 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_078.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 134, 132, 130, 321 referencing this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.txt`

## What the paper claims

Stage 078 is a comparison/verdict stage: it inserts the explicit Family-1 wall-depth datum `Theta_w` (extracted in Stage 077) into the Stage-075 threshold window and computes the resulting Péclet-demand success/failure bands. `\stagefield{Output}` reads verbatim: *"The first explicit Family--1 support/source verdict: wall-depth is not the leading bottleneck for moderate demand."* The concrete deliverables are the two boxed equation blocks: for the natural shell-weighted (chi) datum, `Pe_suff^(chi) ≈ 96.5285247264386 λ_μ²` and `Pe_fail^(chi) ≈ 11220.5441626259 λ_μ²`; for the conservative Jensen floor (J), `Pe_suff^(J) ≈ 22.0062226330754 λ_μ²` and `Pe_fail^(J) ≈ 2558.01892349205 λ_μ²`. The notes add the input provenance (Stage-077 `Theta_w^(chi) ≈ 4.06863235008162 λ_μ²`, `Theta_w^(J) ≈ 0.927552032539308 λ_μ²`; Stage-075 coefficients `Theta_fail/Pe_req ≈ 3.62605617972939e-4`, `Theta_suff/Pe_req ≈ 4.21495341569977e-2`) and the qualitative verdict that the branch succeeds for modest demand and fails only for very large demand. The appendix row 134 is a status-only one-liner (`\StatusExactClosure`, "Wall-depth side is not the bottleneck except for anomalously large demand").

## What the script claims to verify

The SymPy script imports the two Stage-077 `Theta_w` numbers and the two Stage-075 threshold coefficients as cited literals, forms the four ratios `Pe = Theta_w / coeff`, prints them, and asserts (i) the success band lies below the failure band in both routes (`Pe_suff < Pe_fail`), and (ii) the Jensen-floor windows nest below the chi windows with overlap (`Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, `Pe_suff_chi < Pe_fail_J`). The Mathematica script independently re-derives the two threshold coefficients from their Stage-075 symbolic cosh/sinh closed forms (`thetaFailSym`, `thetaSuffSym`, α = 111√5/5), adopts the chi/J `Theta_w` coefficients as extended-precision Stage-077 numerics (explicitly stating their independent re-derivation is Stage-077's job, not this one's), prints the same four Pe values, and asserts the identical ordering/overlap inequalities plus numeric self-consistency checks. Both scripts' bottom line is the four Pe windows and their ordering.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Pe_suff^(chi) ≈ 96.5285247264386` | sympy line 46 print = `96.528524726438575954`; mma line 68 print = `96.5285247264384…` | match |
| `Pe_fail^(chi) ≈ 11220.5441626259` | sympy line 47 print = `11220.544162625905301`; mma line 69 print = `11220.5441626259…` | match |
| `Pe_suff^(J) ≈ 22.0062226330754` | sympy line 48 print = `22.006222633075413597`; mma line 70 print = `22.0062226330754…` | match |
| `Pe_fail^(J) ≈ 2558.01892349205` | sympy line 49 print = `2558.0189234920526360`; mma line 71 print = `2558.0189234920537…` | match |
| Verdict: success band below failure band (success for modest demand) | sympy lines 51–54 `Pe_suff < Pe_fail` (both routes); mma lines 83–84 | match |
| Notes: chi/J windows nest with overlap | sympy lines 60–74; mma lines 87–89 | match (notes-level, not boxed; substantive extra check) |
| Inputs: Stage-077 `Theta_w^(chi)=4.06863235008162`, `Theta_w^(J)=0.927552032539308` | sympy lines 30–31 literals (cited); mma lines 51–54 extended-precision | match |
| Inputs: Stage-075 coeffs `3.62605617972939e-4`, `4.21495341569977e-2` | sympy lines 38–39 literals (cited); mma lines 39–47 symbolic re-derivation | match |

`paper_alignment: aligned`. Every boxed deliverable in the card and every value in the notes is reproduced by both engines to the precision the paper quotes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 51–52 | `Pe_suff_chi < Pe_fail_chi` (else raise) | chi success-below-failure verdict | yes |
| A2 | sympy | 53–54 | `Pe_suff_J < Pe_fail_J` (else raise) | J success-below-failure verdict | yes |
| A3 | sympy | 60–64 | `Pe_suff_J < Pe_suff_chi` (else raise) | notes window-nesting | yes |
| A4 | sympy | 65–69 | `Pe_fail_J < Pe_fail_chi` (else raise) | notes window-nesting | yes |
| A5 | sympy | 70–74 | `Pe_suff_chi < Pe_fail_J` (else raise) | notes window-overlap | yes |
| A6 | mma | 79–82 | `expectApprox[peSuff/Fail vs target, tol]` | four Pe values (numeric self-consistency) | partial (self-referential; see Verdict) |
| A7 | mma | 83–84 | `expectTrue[Pe_suff < Pe_fail]` (both routes) | success-below-failure verdict | yes |
| A8 | mma | 87–89 | `expectTrue` nesting + overlap | notes window structure | yes |

The genuinely independent cross-engine verification is the agreement of the four PRINTED Pe values: SymPy produces them from hardcoded Stage-077/075 floats, Mathematica produces them from a *symbolic re-derivation* of the Stage-075 threshold coefficients, and the two agree to ~13 significant figures (and to the paper box). The ordering/overlap assertions (A1–A5, A7–A8) can genuinely fail (a sign error, a swapped numerator/denominator, or `Theta_J ≥ Theta_chi` would flip them) and are the load-bearing verdict checks.

## Findings

None. Zero math findings and zero `paper_misalignment` findings.

(One non-finding numbering observation is recorded below under "Self-test notes" / the deferred SCRIPT/OUTPUT-band numbering pass — it is a stale docstring cross-reference, not a value/target/math defect, and is explicitly out of scope for this per-stage math audit per the documented numbering-drift remediation plan.)

## Independent-derivation check (Mathematica)

Not a transliteration. The SymPy script obtains the Stage-075 threshold coefficients as plain decimal literals (`Theta_fail = 3.62605617972939e-4 * Pe_req`, line 38; `Theta_suff = 4.21495341569977e-2 * Pe_req`, line 39). The Mathematica script instead carries the *symbolic closed forms*:

- mma 39–43: `thetaFailSym = (37 Cosh[111√5/5] + (111√5/5) Sinh[111√5/5]) / (136900 (-1 + (√5/3) Sinh[111√5/5] + Cosh[111√5/5]))`
- mma 47: `thetaSuffSym = -(45 Cosh[111√5/5] + 27√5 Sinh[111√5/5]) / (2500 - 2500 Cosh[111√5/5])`

and only then numericizes (lines 55–56). This is a structurally different route on the threshold side (symbolic cosh/sinh evaluation vs. a pre-baked float), so the two engines reach the four Pe numbers by independent paths. Hand-check of the large-α limit confirms both symbolic forms reproduce the SymPy literals: `thetaSuffSym → (45 + 27√5)/2500 ≈ 0.042148` (matches `4.21495341569977e-2` after the `cosh−1` and finite-α corrections); `thetaFailSym → (37 + 111√5/5)/(136900 (1+√5/3)) ≈ 3.627e-4` (matches `3.62605617972939e-4`). The chi/J `Theta_w` coefficients are *adopted* numerically in both engines (correctly — their re-derivation is Stage-077's audit, as the mma comment lines 48–50 state).

## Engine cross-check

The four printed Pe values agree:

| value | SymPy (output) | Mathematica (output) |
|---|---|---|
| Pe_suff^(chi) | `96.528524726438575954` | `96.52852472643840784595…` |
| Pe_fail^(chi) | `11220.544162625905301` | `11220.54416262589764…` |
| Pe_suff^(J) | `22.006222633075413597` | `22.00622263307539771…` |
| Pe_fail^(J) | `2558.0189234920526360` | `2558.0189234920537145…` |

Agreement to ~13 significant figures (the divergence in the last digits is the expected float-vs-30-digit-symbolic rounding difference, well within the `expectApprox` tolerances `10^-12`…`10^-9`). All ordering/overlap booleans agree (both engines: `Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, `Pe_suff_chi < Pe_fail_J`, all `True`). No `engine_disagreement`.

## Verdict justification

Clean. I attacked the following and each held: (1) **Hardcoded-numbers tautology** — the SymPy script does divide literals, but the inputs are cited to Stage-077/075 outputs and, crucially, the Mathematica engine re-derives the threshold coefficients *symbolically* and reaches the same four Pe values, so the deliverables are independently anchored, not confirmed against themselves. (2) **Ordering checks failability** — A1–A5/A7–A8 are genuine inequalities that would fail under a sign flip, a swapped ratio, or `Theta_J ≥ Theta_chi`; they are not algebraically guaranteed. (3) **Mathematica `expectApprox` self-reference** — these checks (A6) ARE weak: `peSuffChi = N[thetaChiCoeff/N[thetaSuffSym,30],30]` is compared against `peSuffChiTarget = N[thetaChiCoeff/thetaSuffSym,30]`, i.e. both sides descend from the same `thetaSuffSym` and differ only in whether the symbolic coefficient is rounded to 30 digits before or after the division — so they can essentially only fail on a >`10^-12` rounding artifact. However, A6 is redundant scaffolding, not the load-bearing verification: the real cross-engine check is that the independently-(symbolically-)derived printed values match SymPy and the paper box, which they do. A6 being near-self-consistent does not let the audit pass falsely, so it stays below the finding threshold (no `insufficient_verification` filed). (4) **Transliteration** — ruled out; the threshold side is derived by genuinely different routes. (5) **Paper alignment** — all four boxed Pe values, the input `Theta_w` data, and the threshold coefficients reconcile exactly across `.tex`, `.md`, and both script outputs (see Value Reconciliation). (6) **Stale output** — both `.txt` files post-date their scripts and their content matches the current scripts. I read the card, the notes, and the appendix rows; the script's verified claim matches the paper's claim.

## Self-test notes

- **Variable independence / derivatives**: no `sp.diff`/`D[]` in either script — not applicable.
- **Symmetry/parity**: no integrals over unbounded domains — not applicable.
- **Trivial-case / sign sanity**: confirmed `Theta_J (0.9276) < Theta_chi (4.069)` and both threshold denominators positive, so the nesting/overlap inequalities are the correct sign and genuinely failable; large-α hand-checks of `thetaSuffSym`/`thetaFailSym` reproduce the SymPy literals (no sign or factor error in the cosh/sinh closed forms).
- **Numbering observation (non-finding, deferred)**: the SymPy docstring line 6 ("Insert the explicit Stage-60 Theta values into the Stage-58 threshold window") carries STALE pre-renumber cross-references — Stage-60→77 and Stage-58→75 (the documented +17 EM-extension renumber drift). The script's own inline comments two lines below (lines 26, 29, 32, 37) already use the CORRECT Stage-77/Stage-75 and the correct source filenames (`stage077…`, `stage075…`), and every self-label ("STAGE 078"/"Stage 078") is correct in both scripts and both outputs. This is an unambiguous label-only cross-reference, not a value/target/math defect and not a `paper_misalignment` (the actual numbers and the paper agree). Per the documented remediation policy, mechanical script/output-band label fixes are owned by the dedicated DEFERRED `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` pass, not this per-stage math audit; recorded here so the deferred pass can pick it up. No directive written.

## Value Reconciliation (pass-2 augmentation)

**Deliverable-level table** (every RESULT value the scripts emit, located in the docs):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Pe_suff^(chi) = 96.5285247264386…` | py line 46 / out line 5; wl line 68 / out line 5 | `.tex:17` boxed `96.5285247264386`; `.md:38` `96.5285247264386` | MATCH |
| `Pe_fail^(chi) = 11220.5441626259…` | py line 47 / out line 6; wl line 69 / out line 6 | `.tex:19` boxed `11220.5441626259`; `.md:43` `11220.5441626259` | MATCH |
| `Pe_suff^(J) = 22.0062226330754…` | py line 48 / out line 7; wl line 70 / out line 7 | `.tex:25` boxed `22.0062226330754`; `.md:54` `22.0062226330754` | MATCH |
| `Pe_fail^(J) = 2558.01892349205…` | py line 49 / out line 8; wl line 71 / out line 8 | `.tex:27` boxed `2558.01892349205`; `.md:59` `2558.01892349205` | MATCH |
| input `Theta_w^(chi) = 4.06863235008162` | py line 30; wl line 51 | `.md:14` `4.06863235008162` (Stage-077 carry; `.tex` row terse) | MATCH |
| input `Theta_w^(J) = 0.927552032539308` | py line 31; wl line 52 | `.md:17` `0.927552032539308` | MATCH |
| input `Theta_fail/Pe_req = 3.62605617972939e-4` | py line 38; wl `thetaFailSym` lines 39–43 | `.md:8` `3.62605617972939e-4` | MATCH |
| input `Theta_suff/Pe_req = 4.21495341569977e-2` | py line 39; wl `thetaSuffSym` line 47 | `.md:10`/`.md:37` `4.21495341569977e-2` | MATCH |

**Internal / scaffolding values (no prose expected, no finding):** the three SymPy ordering booleans printed at out lines 9–11 (`Pe_suff_J < Pe_suff_chi` etc.); the Mathematica `expectApprox` diff residuals (`0\`\`27.7…` etc., out lines 9–16) and the `expectTrue` booleans (out lines 17–26); the symbolic intermediates `thetaFailSym`/`thetaSuffSym`/`α = 111√5/5`; the per-check tolerances (`10^-12`…`10^-9`); the PASS flags and "audit passed" banner.

`reconciliation: complete; 8 deliverable/input values checked, 0 misaligned`
