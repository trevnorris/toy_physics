---
unit_id: 023
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md]
  paper_appendix: present
---

# Audit unit 023 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_023.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage023_full_grouped_bundle.md`
- part appendix: `/var/projects/toy_physics/../paper/appendices/stage_appendix_part01.tex` → `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 68 + prose lines 9, 131)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.txt`

## What the paper claims

Stage 023 is the Part-I checkpoint that writes the entire grouped real 20/21/22 bundle in one place. Its `\stagefield{Output}` (line 290) enumerates seven distinct deliverables: (1) the exact weighted projectors `P_x̄, P_a, P_b` on the grouped metric `G_grp = diag(1,2,2)` with the orthogonal basis `e_x̄=(1,1,1), e_a=(4,-1,-1), e_b=(0,1,-1)` (norms 5, 20, 4) and the identities `P_x̄+P_a+P_b=I, P_iP_j=δ_ij P_i`; (2) the BdG support moments `B_{An}`; (3) the conservative Maxwell/mixed moments `Z_{An}^{(r)}`; (4) the outgoing-transfer moments `N_{An}^{(r)}`; (5) the assembled coefficients `D_{A0}=K_A−B_{A0}−Z_{A0}, D_{A2}=−(M_A+B_{A2}+Z_{A2}), D_{A4}=−(B_{A4}+Z_{A4})` with their exact grouped `(bar,a,b)` decomposition; (6) the isotropic normalization ratio `m̂_0² N_0/(K−B_0−Z_0) = 54 G c_s⁵/(5 a⁵ c⁵)` together with the isotropic `u_2,u_4,P_0,P_2,P_4` formulas and the constant-prefactor conditions `N_2=2D_2N_0/D_0, N_4=[2D_0(D_2N_2+D_4N_0)−3D_2²N_0]/D_0²`; and (7) the first-order anisotropy transport laws `a_{u,2}=−(a_{D,2}+u_2 a_{D,0})/D_0` and `a_{P,0}=(a_{N,0}−P_0 a_{D,0})/D_0` (and their `b` partners). The notes add the underlying one-port `Δ,S,Q,G,P` definitions and the monotonicity reading of `P_0=N_0/D_0`. The card's status is `\StatusExactClosure{}` for the algebra / `\StatusOpen{}` for microscopic branch values. The appendix row (line 68) summarizes: "Weighted grouped projectors, full D,N,P coefficient ledger, isotropic branch test, and first-order anisotropy transport laws."

## What the script claims to verify

Both engines verify, symbolically and (in the `.wl`) also numerically: (I) `G_grp`-orthogonality and norms 5/20/4, the projector outer-product forms, idempotency `P_i²=P_i`, mutual orthogonality `P_iP_j=0`, completeness `P_x̄+P_a+P_b=I`, and the exact decomposition of an arbitrary grouped vector; (II) a Schur-complement derivation of the one-port denominator `Δ−Sω²+ω⁴` and a series expansion confirming the closed `Z_n, N_n` forms, plus the grouped `(bar,a,b)` decompositions of `D_{An}` and an additivity (linearity) check on the `N` bundle; (III) the isotropic `u_2,u_4,P_0,P_2,P_4` formulas, the constant-prefactor `N_2/N_4` conditions (cross-checked against independent closed forms), and an independent derivation of the normalization RHS `54 G c_s⁵/(5 a⁵ c⁵)` from the Stage-5 DtN anchor `Γ_5=a⁵/(27 c_s⁵)` and `γ_GR=2G/(5c⁵)`; (IV) the first-order transport laws via an explicit `ε`-series; and (V) monotonicity derivatives of `P_0=N_0/(K−B_0−Z_0)` in `N_0,B_0,Z_0,K`. The `.wl` adds a numerical-substitution cross-check of `Z_n/N_n` and a direct Bessel small-z expansion for `Γ_5`, independent of the SymPy choreography.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| Projectors + basis + norms (5,20,4) | I.1, I.2 (py 75–100; wl 48–65) | match |
| `P_iP_j=δ_ij P_i`, completeness, decomposition | I.2, I.3 (py 94–112; wl 59–71) | match |
| `B_{An}` moments | carried as Stage-003 symbols; assembled into `D` (py 199–223; wl 131–144) | match (carry-forward) |
| `Z_{An}^{(r)}` closed forms | II.0 Schur + series (py 137–167; wl 88–103, 124–129) | match |
| `N_{An}^{(r)}` closed forms | II.0 series (py 169–193; wl 96, 101–103, 127–129) | match |
| `D_{A0/A2/A4}` assembly + grouped (bar,a,b) | II.1 (py 213–261; wl 136–167) | match |
| `u_2,u_4,P_0,P_2,P_4` isotropic | III.1 (py 295–313; wl 182–195) | match |
| Normalization `m̂²N_0/(K−B_0−Z_0)=54Gc_s⁵/(5a⁵c⁵)` | III.3 (py 331–367; wl 205–244) | match |
| Constant-prefactor `N_2,N_4` conditions | III.2 (py 315–329; wl 197–203) | match |
| Transport laws `a_{u,2}, a_{P,0}` (+b) | IV (py 383–419; wl 246–261) | match |
| Monotonicity reading of `P_0` (notes §8) | V (py 426–451; wl 263–270) | match (extra-but-notes-grounded) |

`paper_alignment: aligned` — every Output deliverable maps to a faithful, non-tautological script check; the only "extra" (monotonicity) is explicitly anchored in notes §8.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 75–80 | `expect_zero(orthogonality, norms−{5,20,4})` | projector basis | yes |
| A2 | sympy | 94–100 | `P_i²−P_i, P_iP_j, ΣP−I` | projector identities | yes |
| A3 | sympy | 109–112 | exact decomposition residuals | decomposition | yes |
| A4 | sympy | 137–148 | Schur denom + `Z` rational form | Z derivation | yes |
| A5 | sympy | 158–193 | `Z_n,N_n` series-coeff − closed form | Z/N moments | yes |
| A6 | sympy | 251–270 | grouped `(bar,a,b)` of `D`; `N` additivity | D decomposition | yes |
| A7 | sympy | 306–313 | `u_2,u_4,P_0,P_2,P_4` residuals | isotropic formulas | yes |
| A8 | sympy | 328–329 | `N2/N4_target − closed form` | const-prefactor | yes |
| A9 | sympy | 341–367 | normalization abstract==explicit; `Γ_5` anchor; ratio=54... | normalization | yes |
| A10 | sympy | 400–419 | `du2,dP0` + grouped transport residuals | transport laws | yes |
| A11 | sympy | 448–451 | monotonicity derivatives | notes §8 | yes |
| A12 | math | 48–71 | mirror of A1–A3 | projectors | yes |
| A13 | math | 88–129 | A4/A5 + numerical cross-check (indep) | Z/N moments | yes |
| A14 | math | 159–174 | A6 mirror | D decomposition | yes |
| A15 | math | 191–203 | A7/A8 mirror | isotropic/const-pref | yes |
| A16 | math | 211–244 | A9 + direct Bessel small-z `Γ_5` (indep) | normalization | yes |
| A17 | math | 256–270 | A10/A11 mirror | transport/monotonicity | yes |

No tautological or unanchored rows.

## Findings

None. (`findings_count: 0`.)

## Independent-derivation check (Mathematica)

The `.wl` is not a line-by-line transliteration. It adds two genuinely independent verification paths absent from the SymPy script: (1) lines 105–129, a fixed-parameter numerical substitution `{Ω_U→2, Ω_W→3, R→1, g_U→1, g_W→2}` whose `Z_n,N_n` closed-form values are compared against `SeriesCoefficient[...]` of the rational function — this breaks the symbolic-series correspondence; and (2) lines 231–241, a direct small-z Taylor expansion of `j2 + I·y2` and its derivative (`h2SmallZ`, `h2DerivSmallZ`) feeding `Γ_5`, deliberately bypassing the `ω·D[h2,z]/h2` ratio path used at lines 220–230. For the shared symbolic identities the two engines use the same algebraic targets (as they must — both must verify the same paper identity), but the Mathematica engine reaches the load-bearing constants (`Z_n,N_n`, `Γ_5`) by an independent route. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts show all-zero residuals / PASS for every shared check. SymPy output (lines 9–249) reports `= 0` for all `expect_zero` lines and prints the expected closed forms (`2D₂N₀/D₀`, `N₀(2D₀D₄+D₂²)/D₀²`, `Γ_5 = a⁵/(27c_s⁵)`, required `m̂²P_0 = 54Gc_s⁵/(5a⁵c⁵)`). Mathematica output (lines 9–155) reports `PASS` for the same set plus the six numerical cross-checks (lines 61–72) and the Bessel-path `Γ_5` (line 119–120). Engines agree.

## Verdict justification

`clean`. I attacked: (a) the Section V derivatives for the variable-independence trap — `P_0=N_0/(K−B_0−Z_0)` genuinely depends on all four differentiation variables, so none of the four monotonicity derivatives is identically zero; the assertions are real. (b) The constant-prefactor closed forms — I verified by hand that substituting `N_2=2D_2N_0/D_0` into the paper's general `N_4` condition (tex line 251) yields exactly the script's `N_4_target_closed = N_0(2D_0D_4+D_2²)/D_0²` (py line 327); consistent, and the script independently re-derives `N2_target/N4_target` via `solve` then checks against the closed form, so the substitution is not tautological. (c) The load-bearing normalization constant 54 — the script derives it independently from `γ_GR/(m̂²·Γ_5) = [2G/(5c⁵)]/[a⁵/(27c_s⁵)] = 54Gc_s⁵/(5a⁵c⁵)`, which I confirmed by hand and which matches the paper box (tex line 231) exactly. (d) Sign/factor conventions on `D_{A2},D_{A4}` (both negative) — match tex lines 165–167. (e) The `.wl` numerical-substitution re-check `nNum2 .../. zRule` at lines 119–120 applies `/. zRule` to an already-numeric expression (a harmless no-op, since `pNum/deltaNum/sNum` carry no remaining symbols) — not a defect, the residual cross-check still compares the correct number against `SeriesCoefficient`, confirmed PASS in the transcript (lines 69–71). I read the paper card, the notes, and the appendix row; the script's claim matches the paper's stated Output on all seven deliverables. The only flagged signal is `stale_output` (informational only — see Self-test notes), which does not change the verdict.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 13 deliverable values checked, 0 misaligned.

The stage's substantive deliverables are symbolic closed forms plus the small set of numeric/closed constants the scripts pin. Every one reconciles against the `.tex` card and/or the `.md` notes.

| value | source (py / wl + output) | .tex / .md location | status |
|---|---|---|---|
| `G_grp = diag(1,2,2)` | py 69; wl 43; out — | tex 19; md 32 | MATCH |
| basis `e_x̄=(1,1,1), e_a=(4,−1,−1), e_b=(0,1,−1)` | py 71–72; wl 44–46 | tex 24–28; md 46–50 | MATCH |
| norms `5, 20, 4` | py 78–80; wl 51–53 | tex 33–37; md 62–66 | MATCH |
| `P_x̄,P_a,P_b` matrices | py 82–84 / out 19–36; wl 55–57 | tex 47–61; md 78–82 | MATCH |
| `Z_{A0/A2/A4}^{(r)}` closed forms | py 158–167; wl 98–100 | tex 120–122; md 166–170 | MATCH |
| `N_{A0/A2/A4}^{(r)}` closed forms | py 177–193; wl 101–103 | tex 141–143; md 182–190 | MATCH |
| `D_{A0}=K−B_0−Z_0, D_{A2}=−(M+B_2+Z_2), D_{A4}=−(B_4+Z_4)` | py 213–223 / out 259–263; wl 136–144 | tex 165–167; md 210–214 | MATCH |
| `u_2=−D_2/D_0, u_4=(D_2²−D_0D_4)/D_0²` | py 306–307; wl 191–192 | tex 211–214; md 292–294 | MATCH |
| `P_0=N_0/D_0, P_2, P_4` | py 308–313; wl 193–195 | tex 218–224; md 296–300 | MATCH |
| const-prefactor `N_2=2D_2N_0/D_0`; `N_4=N_0(2D_0D_4+D_2²)/D_0²` (P2=0 branch) ↔ general `N_4=[2D_0(D_2N_2+D_4N_0)−3D_2²N_0]/D_0²` | py 326–327 / out 160–168 + 484; wl 200–201 | tex 250–251; md 332–334 | MATCH |
| normalization RHS `54 G c_s⁵/(5 a⁵ c⁵)` | py 339,366 / out 184–190; wl 209,244 | tex 231; md 304,308 | MATCH |
| transport `a_{u,2}=−(a_{D,2}+u_2a_{D,0})/D_0`, `a_{P,0}=(a_{N,0}−P_0a_{D,0})/D_0` (+b) | py 416–419; wl 258–261 | tex 270–286; md 368–388 | MATCH |
| monotonicity `∂P_0/∂N_0=1/D_0, ∂P_0/∂{B_0,Z_0}=N_0/D_0², ∂P_0/∂K=−N_0/D_0²` | py 448–451 / out 226–249; wl 267–270 | md §8 (lines 405–417) | MATCH (notes-carried) |

INTERNAL (genuine scaffolding, no prose expected, no finding): `Δ,S,Q,G(H),P` one-port intermediate definitions; the Schur block matrix `M(ω)` and `det_M`; the `ε`-series intermediates `u2_eps/P0_eps`; the Stage-5 DtN intermediates `j2,y2,h2,Λ2,Y2,Y2_static,Y2_hat`; the upstream anchors `Γ_5=a⁵/(27c_s⁵)` and `γ_GR=2G/(5c⁵)` (Stage-5/022 carry-forwards used only to *produce* the 54 constant — legitimately not stage-023 deliverables, so their absence from the 023 card is correct); the `.wl` numerical-substitution sample values `{Ω_U→2,Ω_W→3,R→1,g_U→1,g_W→2}`; all residual `= 0` / PASS flags.

Every emitted deliverable value reconciles. No `value_mismatch` or `script_missing_paper_claim` raised.

## Self-test notes

Checked all five traps. (1) Variable-independence: every `sp.diff`/`D[...]` in Section V differentiates `P_0` w.r.t. a symbol it actually contains, so no derivative is spuriously zero. (2) Symmetry/parity: no unbounded-domain integrals in this unit — the only integral-flavored object is the small-z series of the spherical-Hankel ratio, evaluated by Taylor coefficient extraction, not an even/odd integral. (3) Trivial-case: the constant-prefactor and normalization residuals reduce to 0 only via the genuine algebraic identities (hand-verified `N_4` specialization and the `54` constant). (4) Path specs: n/a (no missing-script finding). (5) Paper round-trip: no fix prescribed. One informational note: `stale_output` — script mtimes (2026-06-03 15:59) post-date output mtimes (2026-05-25 22:24), but `git show e2a4780` confirms the only post-capture change was a single docstring label edit ("Stage 6" → "Stage 23"; py line 4, and the analogous `.wl` banner), with zero effect on any computed value; the captured transcripts therefore still faithfully reflect current script behavior. Flagging it as the standard freshness signal only — not a content finding, not blocking.
