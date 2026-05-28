---
unit_id: 173
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage173_axisymmetric_loading_mismatch.md]
  paper_appendix: present
---

# Audit unit 173 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_173.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage173_axisymmetric_loading_mismatch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 77 table row; detail block lines 489-554)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.txt`

## What the paper claims

Stage 173 ("Weak-axisymmetric loading mismatch") computes the physical grouped slopes directly from the weak-axisymmetric expansion of the actual grouped moving-throat response moments. From the expansions `D_{A,0}=D0+eps*lam_A*D01`, `D_{A,2}=D2+eps*lam_A*D21`, `D_{A,4}=D4+eps*lam_A*D41`, `N_{A,0}=N0+eps*lam_A*N01` (signature `lam_20=1, lam_21=1/2, lam_22=-1`), and the definitions `u2^(A)=-D_{A,2}/D_{A,0}`, `u4^(A)=(D_{A,2}^2-D_{A,0}D_{A,4})/D_{A,0}^2`, `P_0^(A)=N_{A,0}/D_{A,0}`, the stage proves: (1) `u2^(1)=-(D21+u2 D01)/D0`; (2) `u4^(1)=-(5 D01+18 D21+81 D41)/(81 D0)` on the canonical branch `(u2,u4)=(1/9,4/81)`; (3) `P1/P0=N01/N0-D01/D0`; (4) the hidden-even relation `u4^(1)=(8/9)u2^(1)` is equivalent to `D41=(2/3)D21+(1/27)D01`; (5) imposing `u2^(1)=0` gives `D21=-D01/9` and hence `D41=-D01/27`. The verbatim `\stagefield{Output}`: "On the even-preserving branch, the remaining grouped defect is the scalar \(\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0\)." The lane form is `Delta_Q^(20)=eps*Xi_load`, `Delta_Q^(21)=(eps/2)*Xi_load`, `Delta_Q^(22)=-eps*Xi_load` (appendix eq. `app-part05-Xiload-lanes`).

## What the script claims to verify

The docstring lists four checks: the weak-axisymmetric expansions for u2/u4/P0; the hidden-even operator identity; the even-preserving collapse `D21=-D01/9, D41=-D01/27`; and the residual defect `Xi_load=N01/N0-D01/D0`. The assertions: build `u2A,u4A,P0A` by first-order series in `eps`, extract first-order coefficients `u21,u41,P1` via `diff/.eps->0 /lam`, then `expect_zero` four identities: (i) `u21` equals the paper `u2^(1)` formula; (ii) canonical `u41` equals the paper `u4^(1)` formula; (iii) `P1/P0` equals `N01/N0-D01/D0`; (iv) hidden-even residual is zero. It then solves `u21_can=0` for D21 and the hidden-even equation for D41, and `expect_zero` checks `D21+D01/9=0` and `D41_even+D01/27=0`. Xi_load and the three lane forms are computed and printed (not asserted).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `u2^(1)=-(D21+u2 D01)/D0` | `expect_zero("u2 slope identity", ...)` sympy:54 / math:49 | match |
| `u4^(1)=-(5D01+18D21+81D41)/(81D0)` canonical | `expect_zero("u4 canonical formula", ...)` sympy:68 / math:60 | match |
| `P1/P0=N01/N0-D01/D0` | `expect_zero("P1/P0 formula", ...)` sympy:72 / math:61 | match |
| hidden-even `D41=(2/3)D21+(1/27)D01` | `expect_zero("hidden-even residual", ...)` sympy:82 / math:65 + solve sympy:84 | match |
| even-preserving `D21=-D01/9` | `expect_zero("D21 + D01/9", ...)` sympy:94 / math:77 | match |
| even-preserving `D41=-D01/27` | `expect_zero("D41 + D01/27", ...)` sympy:95 / math:78 | match |
| `Xi_load:=N01/N0-D01/D0` (Output) | defined+printed sympy:98 / math:81 — backed by P1/P0 assertion (Xi_load is a *definition*, identical to the asserted P1/P0 form) | match |
| lane form `Delta_Q^(20,21,22)=(1,1/2,-1)*eps*Xi_load` | computed+printed sympy:101-111 / math:84-89 (trivial scalar multiples by checked signature constants) | match |

All paper deliverables are faithfully exercised. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54-57 | `expect_zero(u21 - paper u2^(1))` | claim 1 (u2 slope) | yes |
| A2 | sympy | 68-71 | `expect_zero(u41_can + (5D01+18D21+81D41)/(81D0))` | claim 2 (u4 canonical) | yes |
| A3 | sympy | 72-75 | `expect_zero(P1_ratio - (N01/N0 - D01/D0))` | claim 3 (P1/P0) | yes |
| A4 | sympy | 77-82 | `expect_zero(hidden-even residual)` | claim 4 (hidden-even D41 law) | yes |
| A5 | sympy | 94 | `expect_zero(u21_zero_D21 + D01/9)` | claim 5 (D21=-D01/9) | yes |
| A6 | sympy | 95 | `expect_zero(D41_even + D01/27)` | claim 5 (D41=-D01/27) | yes |
| A7 | math | 49 | `expectZero[u21 - paper u2^(1)]` | claim 1 | yes |
| A8 | math | 60 | `expectZero[u41Can + (5d01+18d21+81d41)/(81d0)]` | claim 2 | yes |
| A9 | math | 61 | `expectZero[p1Ratio - (n01/n0 - d01/d0)]` | claim 3 | yes |
| A10 | math | 65 | `expectZero[hidden-even residual]` | claim 4 | yes |
| A11 | math | 77 | `expectZero[u21ZeroD21 + d01/9]` | claim 5 | yes |
| A12 | math | 78 | `expectZero[d41Even + d01/27]` | claim 5 | yes |

Every assertion traces to a specific paper-side deliverable and is non-tautological: each LHS quantity (`u21`, `u41_can`, `P1_ratio`, hidden-even residual, solved `D21`/`D41`) is derived independently by series expansion / Solve, then compared to the paper's closed form. None assign-then-assert the same expression.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:32-78`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections, in identical order with identical intermediate-variable choreography:

- Premise encoding — sympy:37-40 `D0A = D0 + eps*lam*D01` (and d2A, d4A, n0A) vs math:32-35 `d0A = d0 + eps*lam*d01` (and d2A, d4A, n0A): identical.
- Slope extraction — sympy:42-48 `u2A = ...series(-D2A/D0A...); u21 = ...diff(u2A,eps).subs(eps,0)/lam` vs math:37-43 `u2A = ...Series[-d2A/d0A...]; u21 = (D[u2A,eps]/.eps->0)/lam`: same expansion, same derivative-at-zero / lam coefficient extraction, same variable names (u2A,u4A,p0A,u21,u41).
- Even-preserving solve — sympy:84,88,91 `solve(Eq(u41_can,(8/9)u21_can),D41); solve(Eq(u21_can,0),D21); D41_hidden.subs(...)` vs math:67,71,74 `Solve[u41Can==8 u21Can/9,d41]; Solve[u21Can==0,d21]; d41Hidden /. d21->u21ZeroD21`: same solve targets, same order, same substitution.

The two engines do not derive the result by independent routes; the `.wl` echoes the `.py`'s algebra. Per the second-engine policy, both engines should reach each result from the physical premises along distinct paths, so a hidden bug in one is not silently mirrored in the other.

**Why this matters:**
A transliterated second engine provides no genuine cross-check: any algebraic mistake or symbol-domain slip in the SymPy choreography would be reproduced verbatim in Mathematica and both would still "PASS." The independent-engine guarantee is the entire point of carrying a `.wl`.

**Required change:**
Restructure at least the load-bearing slope-extraction steps in the `.wl` so the result is reached by a route distinct from the SymPy series-then-`D[...]/.eps->0` choreography — e.g., obtain the first-order coefficients via `Coefficient[Series[...],eps,1]` directly (or `Limit[(expr-base)/eps, eps->0]`) and verify the canonical-branch identities by substituting the explicit canonical operator values rather than mirroring the SymPy `Solve` sequence. The final `expectZero` targets (the paper closed forms) must stay identical; only the derivation path changes. See directive F1.

**Verification:**
Verifier confirms the `.wl` no longer uses the `D[u2A,eps]/.eps->0)/lam` + sequential `Solve` choreography that mirrors sympy:42-48,84-91, that the six `expectZero` checks still target the paper forms, and that the script still exits 0 with all checks PASS.

### F2 — stale_output (informational — cosmetic banner mislabel)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:26`

**What's wrong:**
Both scripts print a section banner labelled for the wrong stage: `banner("STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT")`. This unit is Stage 173. The mislabel propagates into both saved transcripts (sympy output line 11, mathematica output line 11: "STAGE 156 — ..."). The `.wl`'s trailing line correctly says "Stage 173 Mathematica audit passed." (math:101) and the `.py` docstring correctly names stage173, so the banner is an isolated copy-paste error, not a sign the wrong physics was run. No assertion or computed quantity is affected.

**Why this matters:**
Purely a transcript-labeling defect: a reader scanning the output sees "STAGE 156" for a Stage 173 audit, which undermines traceability of the saved transcript. No mathematical content is wrong.

**Required change:**
Change the banner string in both scripts from `"STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"` to `"STAGE 173 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT"`. See directive F2.

**Verification:**
Verifier confirms both scripts' banner line reads "STAGE 173 ..." and the regenerated transcripts' line 11 reads "STAGE 173 ...".

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py` — see F1. The premise encoding (math:32-35), the series + derivative-at-zero coefficient extraction (math:37-43), the canonical substitutions (math:52-53), the hidden-even residual (math:64), and the two `Solve` calls (math:67,71) are one-to-one with sympy:37-40, 42-48, 60-61, 78-81, 84,88. It is not an independent re-derivation.

## Engine cross-check

Both engines pass identically. Final symbolic forms agree (modulo CAS ordering):
- `u2^(1) general`: sympy `(-D0*D21 + D01*D2)/D0**2` (out:13) vs math `(d01*d2 - d0*d21)/d0^2` (out:13) — same.
- `u4^(1) general`: sympy `(-D0**2*D41 + D0*(D01*D4 + 2*D2*D21) - 2*D01*D2**2)/D0**3` vs math `(-2*d01*d2^2 + 2*d0*d2*d21 + d0*d01*d4 - d0^2*d41)/d0^3` — same.
- `P1/P0`, `Xi_load`, lane forms, `D21=-D01/9`, `D41=-D01/27`, hidden-even `D41=(d01+18 d21)/27`: identical in both transcripts.
All six `expect_zero`/`expectZero` residuals are 0 in both engines; both exit 0. `engines_agree: true`. (The agreement is expected and carries no cross-check weight precisely because the `.wl` is a transliteration — see F1.)

## Verdict justification

The mathematics is correct and the paper alignment is exact. I independently re-derived all five paper deliverables: the first-order Taylor coefficient `u2^(1)=-D21/D0+D2 D01/D0^2` matches the paper formula under `u2=-D2/D0`; the canonical `u4^(1)` reduces to `-(5D01+18D21+81D41)/(81D0)` after substituting `D2=-D0/9, D4=-D0/27` (and I verified the canonical values `u2=1/9, u4=4/81` are consistent with those operator values); `P1/P0=N01/N0-D01/D0`; the hidden-even residual vanishes (the `2/3, 1/27` coefficients are correct: `54/81=2/3`, `3/81=1/27`); and the even-preserving `D21=-D01/9 → D41=-D01/27` follows. The Output scalar `Xi_load` is a definition equal to the asserted `P1/P0` form, and the lane signature `(1,1/2,-1)` matches the paper. Attacks that failed: I checked for a tautology in each `expect_zero` (each LHS is independently series-derived or Solve-derived, none is assign-then-assert); I checked the symbol domains (only `D0`,`N0` divisors, both `nonzero`; no `positive`/branch-cut hazards in a fully rational expression); I checked the series order (`eps,0,2` removeO and `{eps,0,1}` both retain through eps^1); and I checked the `Solve[...][0]`/`First[Solve]` picks (linear equations, unique solution). The verdict is `findings` (not `clean`) for two low-severity script-side issues: (F1) the `.wl` is a line-by-line transliteration of the `.py`, violating the independent-second-engine policy; and (F2) both scripts print a wrong-stage banner ("STAGE 156" for a Stage 173 audit). Neither affects the verified physics, so no stop-cold.

## Self-test notes

Checked the four traps. Variable-independence: every `diff(EXPR,eps)`/`D[expr,eps]` is on an expression that genuinely depends on `eps` (each of d0A,d2A,d4A,n0A contains `eps*lam*...`), so no identically-zero derivative; the `/lam` extracts the first-order coefficient correctly. No unbounded integrals, so no parity trap. Trivial-case: substituting the canonical operator values `D2=-D0/9, D4=-D0/27` into the general `u4^(1)` reproduces `4/81`-consistent canonical slopes and the asserted residuals reduce to 0 by hand. Paper round-trip: the F1 restructure only changes the `.wl` derivation path and keeps the six `expectZero` targets identical to the paper closed forms (no new constant introduced); the F2 banner edit is a pure string change touching no math — neither fix introduces a new paper_misalignment.
