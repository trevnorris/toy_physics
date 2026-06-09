---
unit_id: 183
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage183_triangular_normal_form.md]
  paper_appendix: present
---

# Audit unit 183 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_183.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage183_triangular_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 97 status row, lines 819-916 narrative, line 1543 \input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_183.tex:15): *"Compresses the coherent weak-axisymmetric problem to \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\) with triangular map to \((\Theta_1,\Xi_1,\mathcal R_1)\)."* The notes flesh this out into a small set of boxed deliverables: (D1) the definition of the nontracking transfer-shape slippage `Sigma_nt` (notes 36-46); (D2) the exact triangular observable form `Theta_1 = -C_tr*Sigma_tr`, `Xi_1 = A_tr*Sigma_tr + Sigma_nt`, `R_1 + Xi_1 = -(eps_eta/(1-eps_eta))*Sigma_eta` with `C_tr`, `A_tr` as defined (notes 49-70, 118, 151-153); (D3) the three exact inverse reconstructions for `Sigma_tr`, `Sigma_nt`, `Sigma_eta`, including the simplified ratio `A_tr/C_tr = 2(1+chi0+deltaU)/deltaU` (notes 207-257); and (D4) the triple-rigidity theorem `Theta_1=Xi_1=R_1=0 <=> Sigma_tr=Sigma_nt=Sigma_eta=0` on the constructive branch `chi0>0, deltaU>0, 0<eps_eta<1` (notes 334-343). The appendix (part05:824-879) carries the same `Sigma_nt` definition, triangular form, `C_tr`/`A_tr` defs, and zero-defect statement. This is an algebraic normal-form repackaging of upstream (Stage 182) results, not a new physical derivation; claim status is ExactClosure inside the coherent local D/N branch.

## What the script claims to verify

The docstring (py:8-15) lists four checks: (1) the Stage-182 formulas reduce to `(Sigma_tr, Sigma_nt, Sigma_eta)`; (2) the observable ledger is triangular; (3) the inverse reconstructions are correct; (4) triple-rigidity reduces to the three scalars vanishing. Concretely, the script hardcodes the Stage-182 expressions for `Theta1`, `Xi1`, `R1` (with `Sigma_tr` carried as an opaque symbol), defines `Sigma_nt` by deleting the tracking term from `Xi1`, then asserts: `Xi1 - (A_tr*Sigma_tr + Sigma_nt) == 0`; `R1 + Xi1 + eps_eta/(1-eps_eta)*Sigma_eta == 0`; the three inverse round-trips `== 0`; the ratio identity `A_tr/C_tr - 2(1+chi0+deltaU)/deltaU == 0`; and three `expect_nonzero` prefactor checks. The `.wl` performs the identical sequence.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 `Sigma_nt` definition (notes 36-46) | `SigmaNT_def` py:74-81 / wl:53-58 (verbatim same coefficients) | match |
| D2 `Xi_1 = A_tr*Sigma_tr + Sigma_nt` (notes 165) | `expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)")` py:88 / wl:63 | match |
| D2 `R_1+Xi_1 = -(eps_eta/(1-eps_eta))Sigma_eta` (notes 179-182) | py:89 / wl:64 | match |
| D2 `Theta_1 = -C_tr*Sigma_tr`, `C_tr`/`A_tr` defs (notes 151-158, 118) | `C_tr`,`A_tr` defs py:85-86 / wl:61-62; `Theta1` py:57-59 | match (used in inverse/ratio checks) |
| D3 `Sigma_tr` inverse (notes 209-213) | `expect_zero("Sigma_tr inverse")` py:101 / wl:73 | match |
| D3 `A_tr/C_tr = 2(1+chi0+deltaU)/deltaU` (notes 229-232) | py:105 / wl:77 | match |
| D3 `Sigma_nt` inverse (notes 235-239) | py:108 / wl:80 | match |
| D3 `Sigma_eta` inverse (notes 252-256) | py:111 / wl:83 | match |
| D4 triple-rigidity on branch (notes 334-343) | three `expect_nonzero` prefactor-on-branch checks py:120-122 / wl:90-92 | partial — only verifies diagonal prefactors are nonzero; relies on the inverse round-trips above for invertibility. The iff "<=>" is not directly exercised as a logical biconditional, but the inverse formulas + nonzero prefactors are a faithful constructive proof of it. |

`paper_alignment: aligned` — every boxed deliverable maps to a script check; no value mismatch, no missing deliverable, no extra. The only defects are method-quality (transliteration + decomposition tautology) and stale captured output, not paper misalignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 88 | `expect_zero(Xi1 - (A_tr*SigmaTr + SigmaNT_def))` | D1/D2 | partial — true by construction (SigmaNT_def = Xi1 minus the A_tr*Sigma_tr term, term-for-term); the residual cancels identically. Confirms the split is consistent but cannot fail. |
| A2 | sympy | 89 | `expect_zero(R1 + Xi1 + eps_eta/(1-eps_eta)*SigmaEta)` | D2 | partial — `R1` is *defined* (py:71) as `-eps_eta/(1-eps_eta)*SigmaEta - Xi1`; so `R1+Xi1+eps_eta/(1-eps_eta)*SigmaEta` is identically 0 by construction. Tautological. |
| A3 | sympy | 101 | `expect_zero(SigmaTr_inv - SigmaTr)` | D3 | yes — exercises that the `C_tr` prefactor inverts correctly; fails if `C_tr` is wrong. |
| A4 | sympy | 105 | `expect_zero(A_tr/C_tr - 2(1+chi0+deltaU)/deltaU)` | D3 | yes — genuine non-trivial cancellation; the single strongest real check. |
| A5 | sympy | 108 | `expect_zero(SigmaNT_inv - SigmaNT_def)` | D3 | yes — depends on `ratio = A_tr/C_tr` AND `Theta1 = -C_tr*SigmaTr`; non-trivial. |
| A6 | sympy | 111 | `expect_zero(SigmaEta_inv - SigmaEta)` | D3 | yes — prefactor-inverse round-trip; mild but can fail. |
| A7-A9 | sympy | 120-122 | `expect_nonzero(C_tr / A_tr / dressing_pref)` | D4 | yes — confirms diagonal prefactors nonzero on the positive branch; supports invertibility/rigidity. |
| B1-B9 | mathematica | 63,64,73,77,80,83,90,91,92 | same nine checks, same names/order | same | same anchoring profile as A1-A9 |

A1 and A2 are tautological-by-construction (the residual cancels identically regardless of any physical input). A4/A5 and the inverse round-trips A3/A6 are the substantive checks that could actually fail.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.txt`

**What's wrong:**
Both saved transcripts are older than their scripts:
- sympy `.py` mtime 2026-06-03 15:59:11; sympy `.txt` mtime 2026-05-30 01:40:47.
- mathematica `.wl` mtime 2026-06-03 15:59:11; mathematica `.txt` mtime 2026-05-30 01:40:57.

The content confirms staleness independently: both `.txt` files print the banner `STAGE 166 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT` (sympy .txt:3, mathematica .txt:3), whereas the current scripts print `STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT` (py:42, wl:32). The body checks otherwise match the current scripts (all PASS / `= 0`), so only the banner string is out of date in the captured output; the math result is unchanged.

**Why this matters:**
The committed transcripts do not reflect the current script banner, so a reader trusting the `.txt` would see the wrong stage number. Low severity because the asserted residuals/PASS lines all match what the current script would emit.

**Required change:**
Re-run both scripts and recommit the transcripts so the banner reads `STAGE 183`. No source-code edit is required for this finding (the script banners are already correct).

**Verification:**
After re-run, both `.txt` line 3 should read `STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT` and mtimes should be newer than the scripts.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl:34-92`
- (corresponding) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py:45-122`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:

1. eps definition — py:46 `eps = epsW*(1 - sp.Rational(2,11)*deltaU/(1+deltaU))` ↔ wl:38 `eps = FullSimplify[epsW*(1 - (2/11)*deltaU/(1 + deltaU)), ...]`. Same literal.
2. The Stage-182 observable expressions are *hardcoded identically* in both — py:57-71 builds `Theta1`, `Xi1`, `R1` from the same literal coefficient strings as wl:39-50; neither engine re-derives them from Stage-182 physical premises (`Sigma_tr` is an opaque symbol in both). `Sigma_nt` is then obtained in both by deleting the same `2*chi0/((1+chi0)(1+deltaU))*Sigma_tr` term (py:74-81 ↔ wl:53-58).
3. The assertion sequence is identical in count, order, content, AND label strings: `"Xi_1 - (A_tr Sigma_tr + Sigma_nt)"`, `"R_1 + Xi_1 + eps_eta/(1-eps_eta) Sigma_eta"`, `"Sigma_tr inverse"`, `"A_tr/C_tr - 2(1+chi0+deltaU)/deltaU"`, `"Sigma_nt inverse"`, `"Sigma_eta inverse"`, and the three `expect_nonzero` prefactor checks (py:88-122 ↔ wl:63-92).

Because both engines start from the SAME hardcoded Stage-182 forms and run the SAME algebraic decomposition, the second engine is not an independent route — it cannot catch an error in the shared hardcoded inputs (e.g., a wrong coefficient transcribed identically into both). The only genuine independence is that the two CAS use different simplifiers on the same algebra.

**Why this matters:**
The dual-engine policy requires each engine to verify the result independently of the other's algebra. Here a transcription error in the Stage-182 inputs (the load-bearing content of this stage) would pass both engines identically.

**Required change:**
This is a method-quality finding. The pragmatic correction that adds real independence without rewriting the stage: in the `.wl`, replace the opaque `sigmaTr` with its Stage-182 microscopic definition `Sigma_tr := (1+chi0)*Sigma_del + (1+deltaU)*Sigma_chi` (stated in both scripts' carry-forward print block, py:125 / wl:96) and re-derive `Theta1`/`Xi1` from the underlying microscopic slippages, OR add a concrete numerical spot-check (substitute distinct positive numbers for `chi0, epsW, eps_eta, deltaU` and all five `Sigma_*` and confirm the triangular map + inverses numerically). Either gives the Mathematica engine a route not identical to the SymPy line ordering. Note this is a quality improvement, not a math error — the current checks are correct as far as they go.

**Verification:**
The `.wl` should no longer mirror the `.py` assertion-for-assertion on the identical hardcoded symbolic objects; a numeric branch substitution or microscopic re-expansion appears and exits 0.

## Independent-derivation check (Mathematica)

The `.wl` does NOT derive the claim from first principles independently. Quoted corresponding triples (see F2): the eps definition (py:46/wl:38), the hardcoded `Xi1` (py:61-69/wl:43-49), and the identical `expect_zero` assertion list with identical name strings (py:88-122/wl:63-92). Both engines treat `Sigma_tr` as an opaque carried symbol and operate on the SAME hardcoded Stage-182 expressions. This is a transliteration → `mathematica_transliteration` finding F2.

## Engine cross-check

Both engines agree fully. Every shared check yields `= 0` / `PASS` in both transcripts, and the non-trivial identity matches: sympy `.txt:44` `A_tr/C_tr = 2*(chi0 + deltaU + 1)/deltaU` ↔ mathematica `.txt:27` `A_tr/C_tr = (2*(1 + chi0 + deltaU))/deltaU`; the three nonzero prefactors print the same symbolic forms (sympy .txt:52-54 ↔ mathematica .txt:38-42). The simplified `Sigma_nt` forms differ only in surface presentation (the SymPy form folds `(1-eps)` to its expanded `(-11(1+deltaU)+(11+9deltaU)epsW)/(11(1+deltaU))` denominator; the Mathematica form keeps an equivalent factorization) — equivalent, not a disagreement. `engines_agree: true`.

## Verdict justification

The scripts faithfully and correctly verify every boxed deliverable the paper card and notes state — `paper_alignment: aligned`, no value mismatch, no missing deliverable. Attacks tried: (a) checked the `eps` definition against Stage-182 notes (matches, line 14); (b) checked that `Xi1` retains the `+Sigma_Z` term (it does, py:62); (c) verified `eps_eta` is a distinct symbol from `eps` (it is); (d) verified the `A_tr/C_tr` cancellation and the inverse round-trips are genuinely non-trivial (they are). What does NOT hold up to the dual-engine bar: the `.wl` is a transliteration of the `.py` (F2), and two of the nine assertions (A1, A2) are tautological-by-construction (the split-consistency check and the `R1`-definition round-trip cannot fail) — though the substantive checks A3-A6 carry the real verification weight, so this is `insufficient`-flavored rather than fatal. Verdict is `findings` driven by F2 (transliteration, medium) and F1 (stale output, low). No `paper_misalignment`, no stop-cold.

## Self-test notes

Checked: (1) Variable independence — no `diff`/`D` in either script, so the zero-derivative trap does not apply. (2) No integrals, so parity trap N/A. (3) Trivial-case pre-check on A1/A2 — confirmed both residuals cancel identically by construction (flagged as tautological, not as a fix to "make pass"); on A4 confirmed the ratio cancellation is a real non-trivial identity. (4) Path specs — F2's suggested route stays in `mathematica/`; no new-script path needed. (5) Paper round-trip — the F2 quality improvement (numeric branch substitution or microscopic re-expansion of `Sigma_tr`) introduces no new constant and uses only values already in the notes/scripts, so it cannot create a new paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

The scripts emit no figure-of-merit numerics; every emitted deliverable is a symbolic closed form. Output `.txt` is stale (F1) but the printed forms match the current script source, so reconciliation is based on script source cross-checked against the (consistent) transcripts.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Sigma_nt` def (= Sigma_Z + [2 epsW/(1-eps)][(11+9 deltaU)/(11(1+deltaU))] Sigma_eps - [2 chi0/(1+deltaU)+4 epsW deltaU/(11(1-eps)(1+deltaU)^2)] Sigma_del) | py:74-81 / wl:53-58; sympy .txt:9-20, mathematica .txt:9 | notes 36-46 & 99-109; appendix part05:824-838 (eq:app-part05-Sigma-nt-def) | MATCH |
| `Theta_1 = -C_tr*Sigma_tr`, `C_tr = chi0 deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU))` | py:57-59,86 / wl:39-42,62; sympy .txt:28-30 | notes 146,151-159; appendix part05:847,855-858 | MATCH |
| `Xi_1 = A_tr*Sigma_tr + Sigma_nt`, `A_tr = 2 chi0/((1+chi0)(1+deltaU))` | py:88,85 / wl:63,61; sympy .txt:31-34 | notes 115-118,165; appendix part05:849,860-862 | MATCH |
| `R_1 + Xi_1 = -(eps_eta/(1-eps_eta)) Sigma_eta` | py:89 / wl:64; sympy .txt:35-38 | notes 179-182; appendix part05:851 | MATCH |
| `Sigma_tr` inverse `= -((1+chi0)(1+deltaU)(1+chi0+deltaU)/(chi0 deltaU)) Theta_1` | py:100-101,132 / wl:72-73,103 | notes 209-213, 364-366 | MATCH |
| `A_tr/C_tr = 2(1+chi0+deltaU)/deltaU` | py:105 / wl:77; sympy .txt:44 / mathematica .txt:27 | notes 229-232 | MATCH |
| `Sigma_nt` inverse `= Xi_1 + 2(1+chi0+deltaU)/deltaU * Theta_1` | py:107-108,133 / wl:79-80,104 | notes 235-239, 369-371 | MATCH |
| `Sigma_eta` inverse `= -((1-eps_eta)/eps_eta)(R_1+Xi_1)` | py:110-111,134 / wl:82-83,105 | notes 252-256, 373-377 | MATCH |
| triple-rigidity: diagonal prefactors `C_tr`, `A_tr`, `eps_eta/(1-eps_eta)` nonzero on branch | py:120-122 / wl:90-92; sympy .txt:52-54 | notes 334-343 (iff statement); appendix part05:874-879 | MATCH |
| carry-forward `Sigma_tr = (1+chi0)Sigma_del + (1+deltaU)Sigma_chi` | py:125 / wl:96 (print only) | notes (Stage-182 carry); appendix part05 context | MATCH (display, not asserted) |

INTERNAL (scaffolding, no finding): pass/fail flags, the `expect_zero`/`expect_nonzero` residual values (all 0 / nonzero prefactor strings), the `eps` substitution helper, the banner string, and the opaque carried symbols `Sigma_Z/Sigma_chi/Sigma_eps/Sigma_del/Sigma_eta/Sigma_tr/Sigma_nt`.

reconciliation: complete; 10 values checked, 0 misaligned.
