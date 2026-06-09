---
unit_id: 187
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage187_orbit_quotient_closure.md]
  paper_appendix: present
---

# Audit unit 187 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_187.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage187_orbit_quotient_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 105 table row; line 265 part summary; lines 1087-1136 "Orbit--quotient closure" subsection; line 1551 input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.txt`

## What the paper claims

The `.tex` `\stagefield{Output}` states: *"Proves finite level sets of \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) are exactly \(\mathcal G_*\)-orbits."* The notes upgrade Stage 186's tangent statement to a finite one. Concretely the stage proves: (D1) the finite log-ratio equality of the three branch monomials \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) between two positive states is exactly the linear system \(M_*\Delta\mathbf x=0\), with the *same* rank-3 matrix \(M_*\) found in Stage 186 (notes §3, boxed matrix); (D2) the selected 3×3 minor over the dependent triple \((\Delta_\eta,\Delta_T,\Delta_\mu)\) has determinant \(1+\chi_{0,*}>0\), so the solve is unique (notes §4); (D3) solving the fibre equations reproduces the three Stage 186 finite orbit laws \(\Delta_\eta=2\Delta_c-\Delta_U\), \(\Delta_T=\Delta_U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U)\), and the longer \(\Delta_\mu\) law (notes §4, three boxed formulas); (D4) therefore the invariant-map fibres coincide exactly with the \(\mathcal G_*\)-orbits, giving \(\mathcal M_+/\mathcal G_*\cong(\mathbb R_{>0})^3\) with quotient coordinates \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\) (notes §§5-6, appendix eq:app-part05-orbit-fibre-theorem / eq:app-part05-quotient-space). The notes §7 observable map (\(\Theta_1,\Xi_1,\mathcal R_1\)) is explicitly carried forward from Stages 183-185, not a Stage-187-original deliverable.

## What the script claims to verify

The docstring enumerates four checks: (1) the finite log-ratio equalities of the three monomial invariants equal \(M_*\Delta\mathbf x=0\); (2) solving the fibre equations for \((\Delta_\eta,\Delta_T,\Delta_\mu)\) reproduces the Stage 186 orbit laws; (3) substituting those solved laws annihilates the invariant log-ratio equations; (4) hence the three monomials are complete finite orbit invariants. The assertions match: the script (a) hand-writes the three rows, then independently builds each from the *actual* Stage-187 monomial ratios via `expand_log`/`PowerExpand` and checks they coincide (py:75-86, wl:57-59); (b) builds \(M_*\) and checks each matrix row equals the corresponding hand-written row (py:97-99, wl:69-71); (c) prints the minor determinant — and the `.wl` *additionally* asserts it equals \(1+\chi_{0,*}\) (wl:79); (d) `solve`s the fibre system and checks each solved component equals the hand-written Stage-186 expected law (py:122-124, wl:98-100); (e) round-trip: substitutes the solution back and checks the rows vanish (py:127-129, wl:102-104).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 — three monomial fibre eqs = \(M_*\Delta\mathbf x=0\), same matrix as S186 | `log C_* ratio - row_*`==0 (rows derived from real monomials) + `matrix row k - exact row_k`==0; matrix `M` literal matches notes §3 boxed matrix exactly | match |
| D2 — selected minor det = \(1+\chi_{0,*}>0\) | `.py` prints det only (py:106); `.wl` asserts `Det[minor]-(1+chiStar)==0` (wl:79) | partial (sympy prints, does not assert; see F1) |
| D3 — three finite orbit laws reproduced | `Delta_eta/Delta_T/Delta_mu finite law`==0 vs hand-written expected forms (py:114-124, wl:89-100) | match |
| D4 — fibres = orbits ⇒ quotient ≅ (ℝ_>0)³ | round-trip `row_* after solve`==0 + prose conclusion; the equivalence is a logical consequence of D1+D3, exercised by uniqueness (D2) + solve (D3) | match (proof-by-print of the logical step, but every algebraic premise is asserted) |

`paper_alignment: aligned`. Every load-bearing identity in the notes/appendix is asserted on the script side; the matrix, minor, and three orbit laws match the boxed forms verbatim. The only gap is that SymPy prints the minor determinant without asserting it (F1, `insufficient_verification`), which the Mathematica engine does assert.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75-78 | `expand_log(log(Ctr_ratio)) - row_tr == 0` | D1 | yes (row built from real \(\mathfrak C_{{\rm tr},*}\) ratio) |
| A2 | sympy | 79-82 | `expand_log(log(Cnt_ratio)) - row_nt == 0` | D1 | yes |
| A3 | sympy | 83-86 | `expand_log(log(eps_ratio)) - row_eta == 0` | D1 | yes |
| A4-A6 | sympy | 97-99 | `M*Dx[k] - row_k == 0` | D1 (matrix encodes rows) | partial (two hand encodings cross-checked) |
| — | sympy | 106 | `print(minor.det())` — NO assert | D2 | no (print only → F1) |
| A7-A9 | sympy | 122-124 | `sol[k] - k_expected == 0` | D3 | yes (solve independent of expected forms) |
| A10-A12 | sympy | 127-129 | `row_k.subs(sol) == 0` | D4 (round-trip) | partial (consistency of solve) |
| B1-B3 | mathematica | 57-59 | `PowerExpand[Log[*Ratio]] - row_* == 0` | D1 | yes |
| B4-B6 | mathematica | 69-71 | `m.dx[[k]] - row_k == 0` | D1 | partial |
| B7 | mathematica | 79 | `Det[minor] - (1+chiStar) == 0` | D2 | yes (this is the assert .py lacks) |
| B8-B10 | mathematica | 98-100 | `(k /. sol) - k_expected == 0` | D3 | yes |
| B11-B13 | mathematica | 102-104 | `row_k /. sol == 0` | D4 | partial |

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py:101-106`

**What's wrong:**
The minor-determinant claim (D2 — that the selected \((\Delta_\eta,\Delta_\mu,\Delta_T)\) minor has determinant \(1+\chi_{0,*}>0\), which is what guarantees the unique exact solve in notes §4: *"Because the selected minor has determinant \(\det M_*^{(\Delta_T,\Delta_\eta,\Delta_\mu)}=1+\chi_{0,*}>0\), this solve is exact and unique"*) is only **printed** in the SymPy script, not asserted:

```
print("det selected minor (Delta_eta, Delta_mu, Delta_T) =", sp.simplify(minor.det()))
```

There is no `expect_zero("selected minor determinant", minor.det() - (1 + chi))`. By contrast the Mathematica script *does* assert it (wl:79: `expectZero["selected minor determinant", Det[minor] - (1 + chiStar)]`). So if `minor` were ever mis-transcribed in the `.py`, the SymPy run would still exit 0; the non-singularity premise underpinning the uniqueness of the fibre solve is unguarded in one engine.

**Why this matters:**
D2 is a load-bearing premise: without \(\det\neq 0\) the "exact and unique" solve in §4 (and hence the fibre = orbit equivalence in §5) does not follow. A printed-but-unasserted value is not a verification; the two engines are asymmetric on this one claim.

**Required change:**
Add, after `scripts/...stage187...py:106`, an assertion mirroring the Mathematica engine:
```python
expect_zero("selected minor determinant", minor.det() - (1 + chi))
```
This brings the SymPy engine to parity with `wl:79` and converts the printed D2 value into a guarded check.

**Verification:**
After the fix, the SymPy transcript should contain a new line `selected minor determinant = 0` (and the script still exits 0). The new assertion appears at the line immediately following the existing `print("det selected minor ...")`.

## Independent-derivation check (Mathematica)

The `.wl` is a **line-by-line transliteration** of the `.py`, not an independent re-derivation. Evidence (corresponding sections):

1. Identical hand-coded rows — `.py:46-48` defines `row_tr/row_nt/row_eta`; `.wl:34-36` defines `rowTr/rowNt/rowEta` with byte-for-byte equivalent algebra. Both start from the *same hand-written* linear forms rather than each deriving the row from the monomial independently first.
2. Identical monomial-ratio definitions — `.py:70-74` `Ctr_ratio/Cnt_ratio/eps_ratio`; `.wl:48-50` `ctrRatio/cntRatio/epsEtaRatio`, same factor structure, same exponents, same `log_subs`/`logSubs` substitution dictionaries (`.py:60-69` ↔ `.wl:44-47`).
3. Identical matrix, minor, solve, and expected laws — `.py:89-93` matrix `M` ↔ `.wl:61-65` `m`; `.py:101-105` `minor` ↔ `.wl:73-77` `minor`; `.py:109` `sp.solve([row_tr,row_nt,row_eta],[DEta,DT,DM])` ↔ `.wl:81` `Solve[{rowTr==0,rowNt==0,rowEta==0},{deta,dt,dm}]`; `.py:114-120` `DEta_expected/DT_expected/DM_expected` ↔ `.wl:89-96` `detaExpected/dtExpected/dmExpected`. The check ordering, names, and printed banners are 1:1.

This is genuinely a port. However: for this stage the verification is pure symbolic linear algebra over rational functions of the four parameters — there is essentially one route (build the monomials, take logs, assemble the matrix, solve the 3×3 minor). A `Series`/`Coefficient` black-box is not involved; the "independence" that matters here is whether the second engine's CAS (Solve + FullSimplify) independently confirms the same rational-function identities. It does. I record `mathematica_transliteration` as a structural observation in the cross-check, but do **not** raise it as a blocking finding because (a) the per-engine algebra kernels (SymPy `solve`/`simplify` vs Mathematica `Solve`/`FullSimplify`) are genuinely different solvers cross-checking the same identities, and (b) the `.wl` actually adds the D2 assertion the `.py` omits, so it is not a strict subset echo. The orchestrator may upgrade this to a finding per the dual-engine policy if desired; flagging it here as a NOTE.

## Engine cross-check

Both transcripts agree. SymPy (`...sympy_audit.txt`): `Delta_eta = -Delta_U + 2*Delta_c`, `Delta_T = (Delta_U*chi0_star + ... )/(chi0_star+1)`, `det selected minor = chi0_star + 1`, and all `... law = 0` / `... after solve = 0`. Mathematica (`...mathematica_audit.txt`): `Delta_eta = 2*dc - du`, `Delta_T = (-((1+deltaStar)*(dc+dg)) + (2+chiStar+deltaStar)*du)/(1+chiStar)` (algebraically identical to the SymPy form), `det selected minor = 1 + chiStar`, plus `PASS:` on every `expectZero` including `selected minor determinant`. The solved \(\Delta_T\) forms differ only in unexpanded vs expanded presentation; both reduce to the boxed notes §4 law (confirmed by the `Delta_T finite law = 0` assertion in each). No residual, sign, or factor disagreement. `engines_agree: true`.

## Verdict justification

The scripts faithfully and non-tautologically exercise every load-bearing identity the notes/appendix require: the three monomial fibre equations are derived from the *actual* Stage-187 monomials (not asserted by fiat), the matrix `M_*` literal matches the notes §3 boxed matrix exactly, the 3×3 solve is compared to independently hand-written Stage-186 orbit laws, and the round-trip confirms consistency. Both engines agree. Attacks tried and failed: (i) the first three checks are NOT tautological — `row_tr` is hand-written but checked against `expand_log(log(Ctr_ratio))` built from the real monomial, so a wrong coefficient would fail; (ii) the solve-vs-expected checks are not round-trips — `sol` comes from `solve`, `*_expected` is hand-written, so they could disagree; (iii) symbol domains are correct (`positive=True` on the four \(\chi_{0,*},\delta_{U,*},E_*,F_*\) parameters justifies `1+chi != 0` for the minor; the `r_*` ratios are `positive` justifying the `expand_log(force=True)`/`PowerExpand` log expansion). The single defect is asymmetric coverage of D2: SymPy prints the minor determinant but does not assert it, whereas Mathematica does (F1, low). Verdict `findings` for that one engine-parity gap. I read the paper card, the notes, and the appendix; the script's verified claim matches the paper's `\stagefield{Output}`.

## Value Reconciliation (pass-2 augmentation)

The stage's deliverables are symbolic (no numeric constants are pinned beyond the structural `1+\chi_{0,*}`). Reconciling each emitted closed-form result against the notes/`.tex`/appendix:

| value | source (py/wl + output line) | .tex/.md/appendix location | status |
|---|---|---|---|
| `row_tr = (1+δ)(Δ_γ+Δ_c−Δ_U)+(1+χ)(Δ_T−Δ_U)` | py:46 / wl:34 / out lines 6/6 | notes:160-164 (C_tr fibre eq) | MATCH |
| `row_nt = 2(1+E)Δ_λ+2E Δ_γ+(F−E)Δ_U−Δ_η−(2+E)Δ_W+Δ_μ−F Δ_T` | py:47 / wl:35 / out 7/7 | notes:171-179 | MATCH |
| `row_eta = 2Δ_c−Δ_U−Δ_η` | py:48 / wl:36 / out 8/8 | notes:187-188 | MATCH |
| Matrix `M_*` (3×8) | py:89-93 / wl:61-65 | notes:199-206 boxed `M_*`; appendix eq:app-part05-finite-fibre-equations | MATCH |
| Selected minor det `= 1+χ_{0,*}` | py:106 print / wl:78-79 / out 15/21-23 | notes:226-229; appendix:1165 `\det P=1+\chi_{0,*}>0` | MATCH (asserted only in .wl; see F1) |
| `Δ_η = 2Δ_c−Δ_U` | py:114/122 / wl:89/98 / out 18/26 | notes:235 boxed; appendix eq:app-part05-Keta-orbit-law | MATCH |
| `Δ_T = Δ_U−((1+δ)/(1+χ))(Δ_γ+Δ_c−Δ_U)` | py:115/123 / wl:90/99 / out 19/27 | notes:240-245 boxed | MATCH |
| `Δ_μ = 2Δ_c−Δ_U+2Δ_W−2Δ_λ−E(2Δ_γ+2Δ_λ−Δ_U−Δ_W)−F((1+δ)/(1+χ))(Δ_γ+Δ_c−Δ_U)` | py:116-120/124 / wl:91-96/100 / out 20/28 | notes:248-256 boxed; appendix:1079-1082 | MATCH |
| Quotient result `M_+/G_* ≅ (ℝ_>0)³` (printed conclusion) | py:143-146 / wl:119-122 | notes:343-353 boxed; .tex Output line 15; appendix eq:app-part05-quotient-space | MATCH (prose/print) |

INTERNAL scaffolding (no finding): the `log_subs`/`logSubs` substitution dictionaries; the `r_*` primitive-ratio symbols; `Mx`/`mx` intermediate matrix product; `sol` solution dict; all `expect_zero`/`expectZero` residuals (all `= 0`); `PASS:`/banner lines; `det selected minor` is a deliverable (listed above), not internal.

reconciliation: complete; 9 deliverable values checked, 0 misaligned.

## Self-test notes

Variable-independence: no `diff`/`D` derivatives in this stage — all checks are linear-algebra residuals, so the "derivative-of-the-wrong-variable" trap does not apply. Trivial-case: each `expect_zero` residual is a rational function of the four parameters; the solve-vs-expected pairs would be nonzero if either side's coefficients were wrong (confirmed they are independently sourced), so the zeros are substantive, not vacuous. Symbol domains: `positive=True` on (χ₀,δ_U,E,F) legitimately makes `1+chi ≠ 0` (minor non-singular) and the `r_*>0` assumptions justify the log-expansion — both physically warranted by the positive-coupling state space. Round-trip checks (`row.subs(sol)`) are mildly tautological but are not load-bearing (D1/D3 carry the proof). F1's proposed `expect_zero("selected minor determinant", minor.det()-(1+chi))` was self-tested: `minor.det()` over the literal `[[0,0,1+chi],[-1,1,-F],[-1,0,0]]` expands by the first column to `(-1)·det([[0,1+chi],[1,-F]])·(cofactor sign)` → `1+chi`, so the residual is identically 0 — the assertion passes for the right reason and does not introduce a new paper_misalignment (it matches notes:228 and appendix:1165).
