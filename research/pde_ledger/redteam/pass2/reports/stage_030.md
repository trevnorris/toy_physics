---
unit_id: 030
batch: II.1
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md]
  paper_appendix: present
---

# Audit unit 030 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_030.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage030_selected_mode_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 50)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage030_selected_mode_normalization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage030_selected_mode_normalization_mathematica_audit.txt`

## What the paper claims

Stage 030 translates the Stage-029 selected lower wall mode into the normalized-response convention of Stages 022/025. Its `\stagefield{Output}` states: "Stage~030 outputs the lower eigenvalue \eqref{eq:app-stage030-lambda-minus}, the selected overlap \eqref{eq:app-stage030-s-minus}, the Hellmann--Feynman identity \eqref{eq:app-stage030-HF}, the selected prefactor \eqref{eq:app-stage030-P0minus}, and the target \eqref{eq:app-stage030-selected-target}." Concretely, the five deliverables are: (1) the exact lower eigenvalue `λ_- = (A+B − α₀σ − R)/2` with `R = √((ΔK_ax+α₀δ_κ)² + 4α₀²Π_κ)` and `A+B = 2A+ΔK_ax`; (2) the closed selected overlap `s_- = ½[σ + ((ΔK_ax+α₀δ_κ)δ_κ + 4α₀Π_κ)/R]`; (3) the Hellmann–Feynman identity `s_- = −dλ_-/dα₀`; (4) the exact selected static prefactor `P_{0,-} = β₀s_-/λ_- = −β₀ d(ln λ_-)/dα₀`; (5) the selected-branch target `m̂_-² P_{0,-} = N_Q^target`. The card body also gives the normalized-response coefficients `u_{2,-} = −D_{-2}/D_{-0}`, `u_{4,-} = (D_{-2}² − D_{-0}D_{-4})/D_{-0}²`, `Γ_{5,-} = C_{5,-}/D_{-0}` (eq:app-stage030-selected-even), and the selected Γ form `Γ_{5,-} = (a⁵/27c_s⁵)·β₀s_-/λ_-` (eq:app-stage030-Gamma-selected). The notes add a determinant identity `λ_-λ_+ = AB − α₀(Bκ₀²+Aκ₁²)`, the constants `Γ₂^port = a⁵/(27c_s⁵)` and `N_Q^target = 54Gc_s⁵/(5a⁵c⁵)`, and two limit checks (`s_-→κ₀²` as α→0, `s_-→σ` as α→∞).

## What the script claims to verify

Both engines verify the same five-part program (Parts I–IV). Part I: a generic normalized-response series expansion `D_-/D_{-0}` confirms `u_{2,-}=−D_2/D_0`, `u_{4,-}=(D_2²−D_0D_4)/D_0²`, and `Γ_{5,-}=C_5/D_0` (imaginary part of the w⁵ coefficient). Part II: builds the explicit loaded 2×2 wall matrix, computes `λ_-,λ_+` from it, and checks (a) the HF identity `−dλ_-/dα` equals the typed closed overlap form, (b) the eigenvector-derived `(v·e_-)²` equals the same closed form, and (c) the α→0 weak-loading limit gives `κ₀²` (=`x0`). Part III: confirms `P_{0,-} = −β₀ d(ln λ_-)/dα` and the exact determinant identity `λ_-λ_+ = A(A+ΔK) − α·T0`. Part IV: confirms the two target forms are equivalent (`cond1 − G5_phys·cond2 = 0`, which tests the numeric constant stack 27/54/5/2). The SymPy and Mathematica scripts are independent at the matrix level (both route eigenvalues through the actual matrix rather than typing the closed form).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) λ_- closed form (eq lambda-minus) | py:73 `lam_minus`, validated by eigenvalues of `M_eff` (HF/det checks py:108–109,147) | match |
| (2) s_- closed form (eq s-minus) | py:103 `s_minus_closed`, cross-checked vs HF + eigenvector (py:108–109) | match |
| (3) HF identity s_-=−dλ_-/dα (eq HF) | py:102/108 `s_minus_hf - s_minus_closed`; py:140–143 `P0_sel + β₀ d(log λ_-)/dα` | match |
| (4) P_{0,-}=β₀s_-/λ_-=−β₀ d(ln λ_-)/dα (eq P0minus) | py:126 def + py:140–143 assertion | match |
| (5) target m̂_-²P_{0,-}=N_Q^target (eq selected-target) | py:161–167 `cond1 − G5_phys·cond2` (equivalence of both target forms) | match |
| body: u_{2,-},u_{4,-},Γ_{5,-} (eq selected-even) | py:53–55 series-coefficient checks | match |
| body: Γ_{5,-}=Γ₂^port·P_{0,-} (eq Gamma-selected) | constants Γ₂^port, N_Q^target used in Part IV (py:157–167); the Γ=C5/D0 piece via Part I | match (constants verified in equivalence) |
| notes: det identity λ_-λ_+=AB−α₀(Bκ₀²+Aκ₁²) | py:146–147 `det identity` | match |
| notes: α→0 limit s_-→κ₀² | py:113 `weak-loading overlap limit` | match |
| notes: α→∞ limit s_-→σ | (not checked) | partial (notes-only secondary check, not an Output deliverable) |

`paper_alignment: aligned`. The only uncovered item is the notes' secondary α→∞ limit (`s_-→σ`), which is an illustrative "expected check" in the notes prose, not a stage Output deliverable and not a load-bearing identity. It does not rise to a `paper_misalignment` (the deliverable values are all present and verified). Noted as `partial` for completeness only.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `coeff(w,2) − u2 == 0` | u_{2,-} (eq selected-even) | yes |
| A2 | sympy | 54 | `coeff(w,4) − u4 == 0` | u_{4,-} (eq selected-even) | yes |
| A3 | sympy | 55 | `im(Ysel).coeff(w,5) − Gamma5 == 0` | Γ_{5,-}=C5/D0 (eq selected-even) | yes |
| A4 | sympy | 108 | `s_minus_hf − s_minus_closed == 0` | HF identity + s_- (claims 2,3) | yes |
| A5 | sympy | 109 | `s_minus_eig − s_minus_closed == 0` | s_-=(v·e_-)² (claim 2) | yes |
| A6 | sympy | 113 | `s_minus_closed(α=0) − x0 == 0` | notes α→0 limit | yes |
| A7 | sympy | 140–143 | `P0_sel + β₀ d(log λ_-)/dα == 0` | P_{0,-} / HF (claims 3,4) | yes |
| A8 | sympy | 146–147 | `λ_-λ_+ − (A(A+DK) − α·T0) == 0` | notes det identity; validates λ_-,λ_+ | yes |
| A9 | sympy | 161–167 | `cond1 − G5_phys·cond2 == 0` | target equivalence + constants (claim 5) | yes |
| B1–B9 | mathematica | 39,40,41,97,98,100,122,125,139 | mirror of A1–A9 via `Eigenvalues`/`Eigensystem`/`Series` | same claims | yes |

All load-bearing rows are "yes". The `gamma5Sel − g5*p0Sel` and `lam_- − lambda_req` relations are intentionally NOT asserted (the scripts document at py:135–139 / py:170–177 and wl:117–121 / wl:142–145 that these are definitional-by-construction and that the physical content is exercised elsewhere — a correct and honest call, not a tautological assertion).

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent derivation at the load-bearing step, not a transliteration. The matrix `mMat` (wl:55–56) equals the SymPy `M_eff` (py:82–85), but both engines then compute the spectrum FROM the matrix rather than typing the closed form: the `.wl` uses `Eigenvalues[mMat]` (wl:63) and `Eigensystem[mMat]` (wl:80) and `Select`s the branch matching the typed reference (wl:64–73), while the `.py` uses `M_eff.eigenvects()` (py:88) for the eigenvector and types `lam_minus` directly (py:73) but cross-validates it through the determinant identity (py:147) and HF derivative (py:102). The Part I series checks both use the natural `Series`/`series` route (wl:31 vs py:44). The typed closed forms (u2/u4/gamma5, s_minus_closed) are identical between engines BECAUSE they are the paper's claimed forms being tested — that is correct, not an echo. The comment at wl:48–54 explicitly documents the de-transliteration intent ("routing lamMinus/lamPlus through Eigenvalues[mMat] rather than through a typed closed form"). No `mathematica_transliteration` finding.

## Engine cross-check

The two outputs agree at every check:
- Part I: SymPy `Y_-(omega)` = `1 − D2ω²/D0 + D2²ω⁴/D0² − D4ω⁴/D0 + iC5ω⁵/D0` (sympy txt 5–9) matches Mathematica `1 − (d2*w^2)/d0 + (d2^2*w^4)/d0^2 − (d4*w^4)/d0 + (I*c5*w^5)/d0` (math txt 5). All three coefficient checks = 0 in both.
- Part II: `λ_-` (sympy txt 18–23) matches `(2*a + dK − alpha*(x0+x1) − Sqrt[dK^2 + 2*alpha*dK*(x0−x1) + alpha^2*(x0+x1)^2])/2` (math txt 16); note `Sqrt[dK²+2α·dK·δκ+α²σ²]` is the algebraic expansion of `√((dK+α·δκ)²+4α²·x0·x1)`, identical residual. `s_-` forms agree (sympy txt 32–47 ↔ math txt 22). HF, eigenvector, and weak-limit checks all = 0 in both.
- Parts III–IV: `det identity`, `P0_sel + β₀ d(log λ_-)/dα`, and `normalization equivalence` all = 0 in both transcripts (sympy txt 101–102,107; math txt 32–35,40–41).

`engines_agree: true`.

## Verdict justification

Clean. I attacked every assertion for tautology, hidden definitional collapse, branch coverage, and constant mismatch, and all held. The λ_- and s_- closed forms are not merely re-stated: they are cross-validated through a genuinely independent route — the explicit loaded 2×2 matrix `M_eff`/`mMat`, whose trace `2A+DK−α·σ` and determinant `A(A+DK)−α·T0` I verified by hand reproduce exactly the gap `R=√((DK+α·δκ)²+4α²·x0·x1)` and the det identity, so the eigenvalue/eigenvector checks (A5,A8) genuinely constrain the typed forms. The HF check (A7) differentiates the eigenvalue independently of the typed overlap. The Part IV equivalence (A9) is a real constant-stack check (27·5 vs 54·… for 2G/5c⁵) that would fail if any of N_Q^target's literals were wrong. The two skipped relations are correctly documented as definitional, not silently dropped. The constants Γ₂^port=a⁵/27c_s⁵ and N_Q^target=54Gc_s⁵/5a⁵c⁵ match the notes verbatim. I read the paper card, notes, and appendix row; the script's verified claim matches the paper's stated Output exactly. The only blemishes are cosmetic, non-blocking numbering-drift residuals (see below), which belong to the numbering-reconciliation track, not script-math.

## Self-test notes

(1) Variable independence: the only derivatives are `d/dα` of `λ_-`/`ln λ_-`, and `λ_-` genuinely depends on α (via `α·σ` and `R`), so these derivatives are non-trivial and the HF/P0 checks are not trivially-zero. (2) Parity/integrals: no unbounded-domain integrals in this unit — Part I uses finite series truncation only. (3) Trivial-case: the det identity and HF checks reduce correctly under hand-substitution (trace/det of `M_eff` verified by hand to match `2A+DK−α·σ` and `A(A+DK)−α·T0`; eigenvalue gap verified `= R`); the weak-loading limit `s_-(α=0)=x0` checks out from the closed form. No directive is written (zero script-side findings).

## Value Reconciliation (pass-2 augmentation)

I enumerated every RESULT/deliverable value the scripts emit (from script source + committed `.txt` outputs; nothing executed) and located each in the `.tex` card and/or `.md` notes. The committed `.txt` outputs are stale ONLY in one cosmetic banner line ("STAGE 13 AUDIT COMPLETE", sympy txt 126 / math txt 45) — the June-3 numbering reconciliation (e2a4780) renamed the banner `STAGE 13`→`STAGE 30` in the scripts but the outputs were captured May 26 (before that edit). The substantive math output is unchanged by that label-only edit, so I base the reconciliation on the committed outputs + the script source. All deliverable values reconcile.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `λ_- = (2A+DK − α·σ − R)/2`, `R=√((DK+α·δκ)²+4α²·x0·x1)` | py:73 / sympy txt 18–23; math txt 16 | tex:33–37 (eq lambda-minus); md:60–66 | MATCH |
| `λ_+ = (2A+DK − α·σ + R)/2` | py:74 / sympy txt 24–29; math txt 17 | md:68–70 (λ_+) | MATCH |
| `s_- = ½[σ + ((DK+α·δκ)δκ + 4α·x0x1)/R]` | py:103 / sympy txt 32–47; math txt 22 | tex:43–50 (eq s-minus); md:92–95 | MATCH |
| HF identity `s_- = −dλ_-/dα` (verified residual 0) | py:108 / sympy txt 30; math txt 18–19 | tex:51–56 (eq HF); md:84 | MATCH |
| weak-loading limit `s_-(α→0)=κ₀²` | py:113 / sympy txt 48; math txt 23 | md:98–99 | MATCH |
| `u_{2,-} = −D_{-2}/D_{-0}` | py:49,53 / sympy txt 11 | tex:73 (eq selected-even); md:125 | MATCH |
| `u_{4,-} = (D_{-2}²−D_{-0}D_{-4})/D_{-0}²` | py:50,54 / sympy txt 12 | tex:74–76 (eq selected-even); md:127 | MATCH |
| `Γ_{5,-} = C_{5,-}/D_{-0}` | py:51,55 / sympy txt 13 | tex:77 (eq selected-even); md:131 | MATCH |
| `C_{5,-} = β₅ s_- = G5·β₀·s_-` | py:124 / sympy txt 53–68; math txt 29 | md:116 (`C_{5,-}=β_5 s_-`) | MATCH |
| `Γ_{5,-}(sel) = β₅ s_-/λ_-` | py:125 / sympy txt 69–84; math txt 30 | tex:80–83 (eq Gamma-selected); md:132 | MATCH |
| `P_{0,-} = β₀ s_-/λ_-` | py:126 / sympy txt 85–100; math txt 31 | tex:86–90 (eq P0minus); md:152 | MATCH |
| identity `P_{0,-} = −β₀ d(ln λ_-)/dα` (residual 0) | py:140–143 / sympy txt 101; math txt 32–33 | tex:88–89 (eq P0minus); md:158 | MATCH |
| det identity `λ_-λ_+ = A(A+DK) − α·((A+DK)x0+A x1)` (residual 0) | py:146–147 / sympy txt 102; math txt 34–35 | md:211 | MATCH |
| `Γ₂^port = a⁵/(27 c_s⁵)` | py:157 / used Part IV | md:177 | MATCH |
| `N_Q^target = 54 G c_s⁵/(5 a⁵ c⁵)` | py:158 / used Part IV | tex:95 (eq selected-target, N_Q^target); md:181,191 | MATCH |
| target equivalence `m̂²Γ₂^port P_{0,-} = 2G/(5c⁵)` ⇔ `m̂²P_{0,-}=N_Q^target` (residual 0) | py:161–167 / sympy txt 107; math txt 40–41 | tex:91–96 (eq selected-target); md:169–183 | MATCH |
| `λ_req = m̂²β₀s_-/N_Q^target` (printed, not asserted) | py:169 / sympy txt 108–123; math txt 42 | md:187,199 (λ_- = β₀ s_- · 5a⁵c⁵/(54Gc_s⁵)) | MATCH |

INTERNAL scaffolding (no prose expected, no finding): `D_-(ω)` template symbols `D0,D2,D4,C5`; the per-check residual=0 print values; the generic `Y_-(ω)` printed series form; `β₅=G5·β₀` intermediate; `T0`/`A,B,σ,δκ,Π_κ,R` notation symbols.

reconciliation: complete; 18 deliverable values checked, 0 misaligned.

### Note — non-blocking numbering-drift residuals (not script-math findings)

These are cosmetic stale stage-NUMBER tokens belonging to the numbering-reconciliation track (doc/label-only; never equation/value), recorded here for the orchestrator, not raised as math findings:
- `scripts/...stage030...sympy_audit.py:3` docstring still reads "Stage 13 SymPy audit." (the June-3 reconciliation fixed the banner at py:182 but missed this header line).
- Both committed outputs show "STAGE 13 AUDIT COMPLETE" (sympy txt 126 / math txt 45) — stale banner only; the current scripts emit "STAGE 30". A fresh capture (which the orchestrator's re-run performs) clears this.
These do not affect any verified identity or value.
