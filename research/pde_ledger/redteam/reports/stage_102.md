---
unit_id: 102
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage102_higher_odd_irrelevance.md]
  paper_appendix: present
---

# Audit unit 102 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_102.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage102_higher_odd_irrelevance.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_102}` line — no separate row for this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.txt`

## What the paper claims

The stage card (label `stage:102`; title cosmetically references "Stage~119" but the `\label`, file path, and verification script names all anchor it as unit 102) says that, at 2.5PN, retarded structure beginning at `O(omega^7)` or higher is irrelevant to the point-particle theorem. The quoted bottom-line is: "Retarded terms beginning at \(O(\omega^7)\) do not affect the point-particle \(2.5\)PN theorem." The notes (`moving_throat_pde_stage102_higher_odd_irrelevance.md`) give the supporting algebra: the generalized one-pole denominator `Yhat_Q^ret(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5 - i tau_Q omega^7)` expands through `O(omega^5)` to `1 + omega^2/(4 Omega_Q^2) + omega^4/(4 Omega_Q^4) + i chi_Q (9/(32 Omega_Q^5)) omega^5 + O(omega^6)`, with the new tau_Q parameter only appearing at `O(omega^7)`. The Checks block of the card asks (a) the source/conservative/outgoing factors stay separated in `mhat_0^2 chi_Q N_Q`, (b) higher odd terms begin beyond the 2.5PN coefficient, (c) outgoing l=2 DtN fingerprint matches the canonical `9/32` coefficient.

Deliverables this stage's verification must support:
- D1: tau_Q does not enter the `omega^5` coefficient (the irrelevance statement at 2.5PN).
- D2: tau_Q first appears at `omega^7` with coefficient `i/4` (i.e. d/dtauQ Im[coeff(omega^7)] = 1/4).
- D3: the canonical `omega^5` coefficient is exactly `i chi_Q (9/32)/Omega_Q^5`.

## What the script claims to verify

The SymPy script (`...stage102...sympy_audit.py`) defines `sigma_can = (9/8)/Omega^5`, the closed-form `Y`, expands to series at orders 6 and 8, then **prints** four diagnostic quantities: the two series, `diff(im(coeff(omega^5)), tauQ)`, `diff(im(coeff(omega^7)), tauQ)`, and `im(coeff(omega^5)) - chiQ * 9/(32 Omega^5)`. There are no `assert` statements; nothing forces the script to exit non-zero if any of these residuals were nonzero. The Mathematica script wraps the same three quantities with `expectZero[...]` calls, so it does assert (a) `tauQ irrelevance at omega^5` = 0, (b) `tauQ coefficient at omega^7 - 1/4` = 0, (c) `(coeff(omega^5)/I) - chiQ * 9/(32 omegaQ^5)` = 0.

## Paper ↔ script cross-check

| Deliverable | SymPy | Mathematica |
|---|---|---|
| D1 (tauQ irrelevant at omega^5) | partial — printed but not asserted | match — `expectZero["tauQ irrelevance at omega^5", tauCoeff5]` |
| D2 (tauQ coefficient at omega^7 is 1/4) | partial — printed but not asserted | match — `expectZero["tauQ coefficient at omega^7 - 1/4", tauCoeff7 - 1/4]` |
| D3 (canonical 9/32 odd coefficient) | partial — printed but not asserted | match — `expectZero["check canonical odd coefficient", ...]` |

`paper_alignment: aligned` — the Mathematica side faithfully exercises every paper-side deliverable. The SymPy side computes the right quantities but does not assert them, so SymPy provides only diagnostic evidence, not gating verification.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 13 | `print('series through O(omega^5) =', Yser5)` | D3 (input) | no (no assert) |
| A2 | sympy | 14 | `print('series through O(omega^7) =', Yser8)` | D2 (input) | no (no assert) |
| A3 | sympy | 15 | `print('tauQ coefficient in omega^5 term =', ...)` | D1 | no (no assert) |
| A4 | sympy | 16 | `print('tauQ coefficient in omega^7 term =', ...)` | D2 | no (no assert) |
| A5 | sympy | 17 | `print('check canonical odd coefficient =', ...)` | D3 | no (no assert) |
| A6 | mathematica | 46 | `expectZero["tauQ irrelevance at omega^5", tauCoeff5]` | D1 | yes |
| A7 | mathematica | 47 | `expectZero["tauQ coefficient at omega^7 - 1/4", tauCoeff7 - 1/4]` | D2 | yes |
| A8 | mathematica | 48 | `expectZero["check canonical odd coefficient", (Coefficient[ySer5, omega, 5]/I) - chiQ*(9/32)/omegaQ^5]` | D3 | yes |

A1–A5 are unanchored only in the sense that they emit diagnostic numbers without forcing failure. The underlying expressions are non-tautological — they compute series coefficients from a closed-form `Y` and compare against numbers that follow only if the algebra is right.

## Findings

### F1 — missing_verification_script

**Subtype:** script_doesnt_cover_claim (sympy side)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py:13-17`

**What's wrong:**
The SymPy script contains **no `assert` statements**. Every check is a `print(...)` call. The script computes the right residuals — `sp.diff(sp.im(Yser5.coeff(omega,5)), tauQ)`, `sp.diff(sp.im(Yser8.coeff(omega,7)), tauQ)`, and `sp.im(Yser5.coeff(omega,5)) - chiQ*sp.Rational(9,32)/Omega**5` — but never enforces that they equal zero or 1/4 respectively. The script returns exit code 0 regardless of the numeric values printed. The transcript at `scripts/output/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.txt:11-13` shows the correct values (0, 1/4, 0), but a future regression in SymPy's `series` or `Rational` handling would silently slip through:

```
tauQ coefficient in omega^5 term = 0
tauQ coefficient in omega^7 term = 1/4
check canonical odd coefficient = 0
```

For a non-status-only, non-checkpoint unit whose mathematica peer uses `expectZero[..., Exit[1]]`, the SymPy side must provide independent gating verification.

**Why this matters:**
The framework's "two-engine policy" relies on each engine independently asserting the same identities. With SymPy as print-only, a future change that broke the algebra (e.g. wrong rational, missed factor of I, off-by-one in the series cutoff) would not fail SymPy — verifier would have to manually compare the printed numbers against the expected values to catch it. That defeats the redundancy.

**Required change:**
Add three `assert sp.simplify(...) == 0` statements that exercise D1, D2, D3 — see directive F1 for the exact lines.

**Verification:**
After Codex applies, the verifier runs `redteam exec-sympy 102`. The new transcript should contain three additional lines (or modified lines) where the script asserts the residuals are 0 (and 0, and 0 for the shifted tauQ-omega^7 check), and the script must still exit 0.

## Independent-derivation check (Mathematica)

The Mathematica `.wl` is a near-line-by-line transliteration of the SymPy `.py`: same closed-form `yRet` denominator (`3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5 - I*tauQ*omega^7)`), same `sigmaCan = (9/8)/omegaQ^5`, same `Series[..., {omega, 0, 5}]`/`{omega, 0, 7}` extraction, same coefficient-of-omega-then-derivative-in-tauQ choreography, same canonical comparator `chiQ * (9/32)/omegaQ^5`. This is acceptable here because the "claim" of this stage is itself the algebraic content of one specific denominator: there is no separate first-principles route — both engines must verify the same algebraic identity. I am not raising `mathematica_transliteration` because the deliverables are themselves an algebraic check of a posited form, not a derivation from physical premises.

Note (cosmetic, not a finding): line 26 of the `.wl` prints `"STAGE 085 — HIGHER ODD IRRELEVANCE"` as the banner — a stale label. The label on line 51 ("Stage 102 Mathematica audit passed.") is consistent. This does not affect the math; flagging only as a heads-up to the author. Likewise the paper card title text reads `"Stage~119"` even though `\label{stage:102}` and the verification file paths are unit 102. Both are cosmetic; the math content of card, notes, and scripts is internally consistent.

## Engine cross-check

Both transcripts agree on the load-bearing residuals (after re-ordering Mathematica's `Expand` output to match SymPy):

SymPy (`scripts/output/...sympy_audit.txt:9-13`):
- `series through O(omega^5) = 1 + omega**2/(4*Omega**2) + omega**4/(4*Omega**4) + 9*I*chiQ*omega**5/(32*Omega**5)`
- `series through O(omega^7) = ... + I*omega**7*tauQ/4 + ... + 9*I*chiQ*omega**7/(16*Omega**7)`
- `tauQ coefficient in omega^5 term = 0`
- `tauQ coefficient in omega^7 term = 1/4`
- `check canonical odd coefficient = 0`

Mathematica (`mathematica/output/...mathematica_audit.txt:13-21`):
- `series through O(omega^5) = 1 + (((9*I)/32)*chiQ*omega^5)/omegaQ^5 + omega^4/(4*omegaQ^4) + omega^2/(4*omegaQ^2)`
- `series through O(omega^7) = ... + (((9*I)/16)*chiQ*omega^7)/omegaQ^7 + omega^6/(4*omegaQ^6) + ... + (I/4)*omega^7*tauQ`
- `tauQ coefficient in omega^5 term = 0` and `PASS`
- `tauQ coefficient in omega^7 term = 1/4` and `PASS` on `tauCoeff7 - 1/4`
- `PASS: check canonical odd coefficient`

Engines agree term-by-term (modulo variable name `Omega` vs `omegaQ` and ordering); no engine disagreement.

## Verdict justification

`findings: F1`. The math claim of the stage (D1, D2, D3) is correct and is faithfully exercised by the Mathematica script (assertions A6, A7, A8). The SymPy script, however, only prints — it cannot fail. That violates the two-engine requirement for a non-status-only unit. There is no paper↔script disagreement (paper, notes, and Mathematica all agree on the `9/32` canonical coefficient and on the `O(omega^7)` start of tauQ), so `paper_alignment: aligned`. The output transcripts are newer than the scripts (sympy: txt 1778525113 > py 1775068798; math: txt 1778526364 > wl 1778522213), so no `stale_output` finding. Mathematica is not flagged as transliteration because the stage's claim is intrinsically algebraic; both engines need to verify the same identity. No downstream propagation: adding asserts in SymPy cannot break anything downstream.

Attacks I tried that failed:
- Tried: is the canonical odd coefficient secretly tautological? It is not — the script defines `sigma_can = (9/8)/Omega^5` and derives `9/(32 Omega^5)` from a different formula path (the `1/4` prefactor on the rational expansion). The factor-of-4 reduction is the substantive content.
- Tried: does `assume(positive=True)` on `chiQ, tauQ, Omega` hide a sign issue? No — the imaginary parts are linear in `chiQ` and `tauQ` and the assumption only affects `simplify` behavior; the residuals are exact rationals/zeros, no branch ambiguity.
- Tried: does `sp.series(Y, omega, 0, 6)` actually capture the omega^5 term? Yes — `0, 6` means up through omega^5 inclusive in SymPy's convention, and the transcript confirms omega^5 is present.
- Tried: does the Mathematica `Coefficient[ySer5, omega, 5]/I` lose sign information vs SymPy's `sp.im(...)`? No — `Coefficient[ySer5, omega, 5]` returns a purely-imaginary expression for the omega^5 term, and `/I` gives the real coefficient, which matches `sp.im(...)` on a purely imaginary input.

## Self-test notes

- Variable independence: `tauCoeff5 = D[Coefficient[ySer5, omega, 5]/I, tauQ]` — does `ySer5` actually depend on `tauQ`? Yes, because `yRet` carries `tauQ` in the denominator; `Series[..., {omega, 0, 5}]` truncates *after* the algebraic expansion, so `tauQ` is present but its coefficient at `omega^5` is identically zero (its lowest-order contribution is `omega^7`). So `D[...,tauQ]` at omega^5 is zero non-trivially. Good — not an "identically-zero-by-construction" trap.
- Symmetry/parity: not applicable (no integrals).
- Trivial-case pre-check: substituting `chiQ=0, tauQ=0` gives `Y = 1/(1 - omega^2/Omega^2)` and the omega^5 imaginary coefficient is 0; the canonical-odd check becomes `0 - 0 = 0`. Sub `tauQ -> tauQ + delta`: omega^5 coefficient is unchanged (no delta dependence) → `tauCoeff5 = 0` literally. Sub at omega^7: coefficient gains `I*delta/4` → `tauCoeff7 = 1/4`. All three asserts the directive prescribes pass on the symbolic substitution.
- Path specification: my prescribed sympy-side assertions go into the existing file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py` — no new files needed.
- Paper round-trip: the asserts use the same constants the paper/notes state (`9/32`, `Omega_Q^5`, tauQ first at omega^7, omega^7-coefficient derivative = 1/4). No new paper_misalignment introduced.
