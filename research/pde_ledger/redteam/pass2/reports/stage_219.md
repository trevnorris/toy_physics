---
unit_id: 219
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 219 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_219.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (Stage-219 narrative lines 233–364; index row line 50)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.txt`

## What the paper claims

Stage 219 inserts the isotropic one-port wall/BdG/Maxwell/mixed bundle into the completed same-charge audit and asks whether the static mixed sector creates a new long-range attractive family. `\stagefield{Output}` states verbatim: "Exact square-law suppression theorem: the one-port static mixed bundle produces only $x^{-6}$, $e^{-2\kappa x}x^{-4}$, and $e^{-4\kappa x}x^{-2}$ product families, and no new $-1/x$ or Yukawa-strength long-range same-charge attraction." The derivation ledger and notes enumerate the distinct deliverables: (1) the reduced static $3\times3$ kernel $\mathcal K_{\rm red}$ with $\det\mathcal K_{\rm red}=\Delta D_0$; (2) the Schur-complement form of $D_0$ on the admissible $(U,W)$ block; (3) the six exact inverse-susceptibility entries $\chi_{qq},\chi_{qU},\chi_{qW},\chi_{UU},\chi_{UW},\chi_{WW}$, all over denominator $\Delta D_0$; (4) the on-shell shift $\delta V_{\rm mix}=-\tfrac12 J^T\mathcal K_{\rm red}^{-1}J$ and the collinear-source factorization $\delta V_{\rm mix}=-\tfrac12\chi_s\mathcal S^2$ with the explicit $\mathcal N_s$; (5) the outgoing-prefactor bridge $\Lambda=P/\Delta$, $N_0=\Lambda^2$, $P_0=\Lambda^2/D_0$, $\chi_{qW}=\Lambda/D_0$, and $\chi_{qW}^2=P_0/D_0$; (6) the primitive-source product-kernel theorem giving only the three families $x^{-6}$, $e^{-2\kappa x}/x^4$, $e^{-4\kappa x}/x^2$ (Theorem "Static one-port square-law suppression", appendix line 359). The admissibility side is $\Omega_U^2>0,\Delta>0,D_0>0$. The card's `\stagefield{Verification}` field (line 11) reads "Mathematica audit: none yet."

## What the script claims to verify

Both scripts build $\mathcal K_{\rm red}$ from the paper's matrix and verify, in sequence: (M1/det) $\det\mathcal K_{\rm red}=\Delta D_0$; (M2/Schur) the $(U,W)$-block Schur complement equals $D_0$; (M3) all six $\mathcal K_{\rm red}^{-1}$ entries equal their paper closed forms; (M4) the collinear-source factorization $-\tfrac12 J^T\mathcal K_{\rm red}^{-1}J=-\tfrac12\chi_s S^2$ with $\chi_s=\mathcal N_s/(\Delta D_0)$; (M5) the bridge identities $\chi_{qW}=\Lambda/D_0$ and $\chi_{qW}^2=P_0/D_0$; (M6) the primitive-source product kernel; and (M7) a constructive admissible numeric slice ($K_*=11,\Omega_U=3,\Omega_W=4,R=2,G_U=1,G_W=2$) confirming positive-definiteness. The `.wl`'s M6 is the load-bearing square-law check and is implemented as a *structural family extraction* (not a coefficient compare): it scales $-2\,\delta V_{\rm p}\,x^6 e^{4\kappa x}$, substitutes $e^{2\kappa x}\to z$, then uses `Collect`/`Coefficient`/`CoefficientList`/`Exponent` to prove the result is a polynomial in $(x,z)$ of $z$-degree exactly 2, $x$-degree $\le 4$, with the three $z$-coefficients sitting at $x$-degrees exactly $\{4,2,0\}$, no $x^5$ cross term, and all three families nonzero.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| $\det\mathcal K_{\rm red}=\Delta D_0$ | py L46 / wl M1 | match |
| $D_0$ = Schur complement of $(U,W)$ block | py L59 / wl M2 | match |
| six $\chi$ inverse entries over $\Delta D_0$ | py L83-85 / wl M3 | match |
| $\delta V_{\rm mix}=-\tfrac12 J^T K^{-1}J$ + collinear $\chi_s,\mathcal N_s$ | py L110 / wl M4 | match |
| $\Lambda,N_0,P_0,\chi_{qW}=\Lambda/D_0,\chi_{qW}^2=P_0/D_0$ | py L130-131 / wl M5 | match |
| product kernel: only $x^{-6}, e^{-2\kappa x}/x^4, e^{-4\kappa x}/x^2$ (square-law suppression, Output) | py L164 / wl M6 (structural) | match |
| admissibility $\Omega_U^2>0,\Delta>0,D_0>0$ | py L199-202 / wl M7 | match |
| card `Verification`: "Mathematica audit: none yet" | a `.wl` now exists and passes | mismatch (card understates engine coverage) |

`paper_alignment: partial` — all seven mathematical deliverables match exactly; the only discrepancy is the card's stale `Verification` field, which now mis-states that no Mathematica audit exists.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `simplify(det_identity - Delta*D0) == 0` | det identity | yes |
| A2 | sympy | 59 | `simplify(schur - D0) == 0` | Schur $D_0$ | yes |
| A3 | sympy | 85 | `simplify(together(actual - expected)) == 0` ×6 | six $\chi$ entries | yes |
| A4 | sympy | 110 | `simplify(together(delta_V_col - delta_V_col_expected)) == 0` | collinear $\chi_s$ | yes |
| A5 | sympy | 130-131 | `chi_qW - Lambda/D0 == 0`, `chi_qW^2 - P0/D0 == 0` | bridge identities | yes |
| A6 | sympy | 164 | `simplify(together(delta_V_primitive - delta_V_primitive_expected)) == 0` | product kernel | yes |
| A7 | sympy | 199-202 | `Delta_sample>0`, `D0_sample>0`, minors>0, eigvals real & >0 | admissibility | yes |
| M1 | mathematica | 22 | `zeroQ[Det[Kred]-Delta*D0]` | det identity | yes |
| M2 | mathematica | 32 | `zeroQ[schur-D0]` | Schur $D_0$ | yes |
| M3 | mathematica | 44 | `And@@(zeroQ/@m3Residuals)` (6) | six $\chi$ entries | yes |
| M4 | mathematica | 59 | `zeroQ[dV - (-1/2 chiS S^2)]` | collinear $\chi_s$ | yes |
| M5 | mathematica | 70 | `chiqW-Lambda/D0`, `chiqW^2-P0/D0` | bridge identities | yes |
| M6 | mathematica | 108-112 | `zeroQ[m6Residual] && m6StructuralOK` (degree/x5/nonzero family checks) | square-law suppression (Output) | yes |
| M7 | mathematica | 143 | residuals + `PositiveDefiniteMatrixQ` + eig real/positive | admissibility | yes |

Every script-side check traces to a specific paper deliverable. No orphaned scaffolding; M7's numeric values are an internal positive-definiteness sanity slice, not a paper deliverable.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (engine-coverage statement, not a value)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_219.tex:11` quote: "SymPy audit: \StageFile{...stage219..._sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage219_..._sympy_audit.md:462-486` (section 9 "Script-backed status" lists only the SymPy file)
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage219_..._mathematica_audit.wl` exists and its committed output passes all of M1–M7.

**What's wrong:**
The card's `\stagefield{Verification}` field asserts "Mathematica audit: none yet." A Mathematica audit `.wl` now exists (added in the VII.1 dual-engine forward pass, commit `1dfc3fe`), passes all seven checks, and is part of the saved-output record. The card and the notes' section 9 both still describe the stage as SymPy-only. This is a documentation staleness in the *direction* of understating coverage (the card claims less verification than exists), so it does not create a false verification claim — but the card's stated `Verification` is now factually inaccurate.

**Why this matters:**
A reader citing the card would believe Stage 219 is single-engine when it is in fact dual-engine verified. Direction of resolution (update the card/notes to reference the `.wl`) is a paper-side edit and therefore the user's call, not Codex's. No script change is implied.

**Required change:**
`## Resolve before fix_loop` (see directive). The user decides whether to update `stage_219.tex:11` and notes §9 to reference the existing `.wl`. Codex applies nothing to the scripts for this finding.

**Verification:**
After user resolution, the card's `Verification` field references the `.wl` path; no script edit; both engines continue to pass.

## Independent-derivation check (Mathematica)

The `.wl` is a **genuinely independent verification route**, not a transliteration of the `.py`.

The two scripts necessarily share the *physical premises* — `Kred`, `Delta`, `Q`, `P`, `PU`, `D0` — because those are the paper's definitions (notes §1–2, appendix eqs. `app-part07-static-delta-qp`/`app-part07-static-d0`/`app-part07-static-kred-matrix`), not derivation choices. Sharing premises is required, not a port signature. The verification *route* diverges:

- M1–M5: each engine computes its own `Inverse[Kred]` / `K_red.inv()` and `Det[Kred]` and compares to the paper's claimed closed forms. The expected right-hand sides come from the paper card, not from copying the other engine's intermediate algebra. Acceptable independent second-engine practice.
- M6 (the load-bearing square-law theorem) is *more* independent than the `.py`. The `.py` (L155–164) compares `delta_V_primitive` against a hand-assembled `C6/C4/C2` closed form. The `.wl` (L86–102) does NOT compare to a pre-baked answer; it scales the residual `-2 dVp x^6 Exp[4 kappa x]`, substitutes `Exp[2 kappa x] -> z`, then *structurally proves* the family content via `Collect`/`Coefficient`/`CoefficientList`/`Exponent`:
  - `Exponent[familyZ, z] == 2`, `Exponent[familyZ, x] <= 4`,
  - `zDegrees === {4, 2, 0}` (the $z^0/z^1/z^2$ coefficients sit at $x$-degree $4/2/0$),
  - `zeroQ[x5Coefficients]` (no spurious $x^5$ cross term),
  - all three `Coefficient[familyZ, z, {0,1,2}]` nonzero.

This is a structurally different attack on the same claim. Mapping back: $C_6/x^6 \to C_6 z^2$ ($x$-deg 0), $2C_4 e^{-2\kappa x}/x^4 \to 2C_4 x^2 z$ ($x$-deg 2), $C_2 e^{-4\kappa x}/x^2 \to C_2 x^4$ ($x$-deg 4) — exactly the `{4,2,0}` by-z-power $x$-degree pattern the script asserts. The combination (z-deg = 2, x-deg ≤ 4, by-power x-degrees exactly {4,2,0}, x^5 coeff zero) is precisely what excludes any longer-range term: a $-1/x$ contribution would scale to $x^5 z^2$ (forbidden by both `x<=4` and `x5Coefficients==0`), and an $e^{-2\kappa x}/x$ term to $x^5 z$ (same exclusion). So the `.wl` independently establishes the Output's *suppression* clause, not merely the presence of the three families.

**Independence call: independent.** No shared-operation port signature beyond the unavoidable shared premises; the load-bearing check uses a distinct extraction method.

## Engine cross-check

Both saved outputs agree. SymPy emits the closed-form symbolic results (det, six χ, δV_mix, χ_s, Λ/N0/P0/χ_qW, C6/C4/C2, product kernel) and the numeric slice; Mathematica reports `M1…M7 = PASS` with all residuals `0`, the M6 scaled family in `(x,z)` matching the SymPy product kernel term-by-term, `zDegrees = {4,2,0}`, `CoefficientList` lengths `{5,3,1}`, and `PositiveDefiniteMatrixQ = True`. The M7 numeric slice agrees exactly across engines: `Delta=140`, `D0=74/7`, `det=1480`, minors `{11,98,1480}`, eigenvalues `{17.0188, 11.2538, 7.7274}`. I independently re-derived the slice by hand (Q=60, D0=11−60/140=74/7, det=140·74/7=1480, minor_22=99−1=98) — all consistent. No engine disagreement.

## Verdict justification

Attacks attempted and failed: (a) tautology hunt — every assertion compares an engine-computed object (`Inverse`, `Det`, the matrix product `J^T K^{-1} J`) against an independently-written paper closed form, so none is true-by-construction; the `.py` M6 builds `delta_V_primitive` (matrix product) and `delta_V_primitive_expected` (coefficient assembly) by different paths. (b) port check — the `.wl` shares only premises and uses a structurally different M6 extraction; it is independent. (c) suppression-clause check — the M6 degree/x5/nonzero checks genuinely exclude the forbidden $-1/x$ and $e^{-2\kappa x}/x$ long-range families, so the Output is actually exercised, not merely its presence half. (d) value reconciliation — all seven mathematical deliverables match the card/notes/appendix verbatim. The single finding is a low-severity paper-side staleness: the card and notes §9 still say "Mathematica audit: none yet" while a passing `.wl` exists. Per the prompt, paper↔script disagreements route to the user, never to Codex. Verdict: `findings` (one paper_misalignment, user resolution; no script-side defect). The residual stale stage-number labels flagged in the first pass were not re-encountered as new blockers here and remain the known deferred numbering-cleanup class.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table (symbolic results are the deliverables here; the numeric slice is internal sanity scaffolding):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| $\det\mathcal K_{\rm red}=\Delta D_0$ | py out L11 / wl out L2-3 (M1) | tex L13; md L126-130; appx L288-291 | MATCH |
| $D_0=K_*-Q/\Delta$ (Schur form) | py out L8,L15-16 / wl M2 | md L64-68; appx L248-251 | MATCH |
| $\chi_{qq}=1/D_0$ | py out L19 / wl M3 | md L177; (notes) | MATCH |
| $\chi_{qU}=P_U/(\Delta D_0)$ | py out L20 / wl M3 | md L180-183 | MATCH |
| $\chi_{qW}=P/(\Delta D_0)$ | py out L21 / wl M3 | md L185; appx L311-315 | MATCH |
| $\chi_{UU}=(K_*\Omega_W^2-G_W^2)/(\Delta D_0)$ | py out L22 / wl M3 | md L188-191 | MATCH |
| $\chi_{UW}=(K_*R+G_UG_W)/(\Delta D_0)$ | py out L23 / wl M3 | md L193-196 | MATCH |
| $\chi_{WW}=(K_*\Omega_U^2-G_U^2)/(\Delta D_0)$ | py out L24 / wl M3 | md L198-201 | MATCH |
| $\delta V_{\rm mix}=-\tfrac12 J^T K^{-1}J$ | py out L26-27 / wl M4 | md L155-159; appx L302-307 | MATCH |
| $\chi_s=\mathcal N_s/(\Delta D_0)$ (collinear) | py out L30 / wl M4 | md L216-234 | MATCH |
| $\Lambda=P/\Delta$ | py out L34 / wl M5 | md L243-245; appx L317 | MATCH |
| $N_0=\Lambda^2$ | py out L35 / wl M5 | md L250; appx L320 | MATCH |
| $P_0=\Lambda^2/D_0$ | py out L36 / wl M5 | md L251; appx L322 | MATCH |
| $\chi_{qW}=\Lambda/D_0$ | py out L37-38 / wl M5 | md L256-258; appx L324 | MATCH |
| $\chi_{qW}^2=P_0/D_0$ | py out L38 / wl M5 | md L263-265; appx L326 | MATCH |
| $\mathcal C_6=\chi_{qq}\beta_Q^2$ | py out L41 / wl M6 | md L320-321 | MATCH |
| $\mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W$ | py out L42 / wl M6 | md L323-327 | MATCH |
| $\mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2$ | py out L43 / wl M6 | md L329-335 | MATCH |
| product kernel: $x^{-6},e^{-2\kappa x}/x^4,e^{-4\kappa x}/x^2$ (square-law) | py out L44-45 / wl out L13-16 (M6) | tex L15 (Output); md L304-316; appx L348-362 | MATCH |

INTERNAL items (raise no finding): numeric admissible slice `Delta(sample)=140`, `D0(sample)=74/7`, `det(sample)=1480`, leading principal minors `{11,98,1480}`, sample eigenvalues `{17.0188, 11.2538, 7.7274}`, `PositiveDefiniteMatrixQ=True`, and all per-check pass/fail flags & zero residuals. These are positive-definiteness sanity scaffolding for a chosen parameter point, not stated card/notes deliverables.

reconciliation: complete; 19 deliverable values checked, 0 value-misaligned. (The one finding F1 is an engine-coverage staleness in the card's `Verification` field, not an emitted-value mismatch.)

## Self-test notes

I checked: (1) variable independence — no zero-derivative trap here; the only `diff`-style operations are SymPy/Mathematica `Inverse`/`Det` of `Kred`, which depend on all six bundle symbols, so no identically-zero residual masquerading as a pass. (2) Parity/symmetry — n/a (no unbounded-domain integrals; the kernel work is algebraic). (3) Trivial-case pre-check — I hand-verified the M7 numeric slice (Q=60, D0=74/7, det=1480, minors 11/98/1480) and the M6 family-degree mapping ({4,2,0} ↔ the three product families) reduce correctly. (4)/(5) Path/round-trip — the sole finding is paper-side (no script edit prescribed), so no new paper_misalignment is introduced; the directive routes it to the user.
