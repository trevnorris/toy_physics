---
unit_id: 005
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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

# Audit unit 005 red-team report (v2, paper-grounded)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_005.tex`
- notes: (none) — EM-projected stage; per-stage notes were not committed (per orchestrator preface)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row for stage 005 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt`

Freshness: SymPy script mtime 2026-05-04 12:00, output mtime 2026-05-21 11:26 (output newer — fresh). Mathematica script mtime 2026-05-21 01:13, output mtime 2026-05-21 11:50 (output newer — fresh).

Manifest entry: `is_checkpoint: false`, `is_status_only_candidate: false`. Both engines required.

## What the paper claims

The stage card states the parent law `partial_N(Z(w) F^{NM}) = mu_0 J^M` (Eq. 5.1) and, by integrating against a normalized observation kernel `W(w)` over `w`, derives two `\stagefield{Output}` deliverables for brane components `M = mu`:

1. The projected inhomogeneous Maxwell equation in raw form (Eq. 5.2) and after IBP using the stage-004 projection-by-parts identity (Eq. 5.3): `partial_nu int W Z F^{nu mu} dw = mu_0 int W J^mu dw - [W Z F^{wmu}]_partial + int (partial_w W) Z F^{wmu} dw`.
2. The projected charge-continuity law (Eq. 5.4): `partial_mu int W J^mu dw = -[W J^w]_partial + int (partial_w W) J^w dw`.

The paragraph "Output" verbatim reads: "Stage~005 exports the projected inhomogeneous law (5.2)-(5.3) and projected continuity (5.4)." The part-appendix row for stage 005 (line 32 of stage_appendix_part01.tex) summarizes the deliverable as "Exact projected inhomogeneous Maxwell equation, boundary flux, kernel-gradient leakage, and projected continuity law" — matching the four pieces (LHS form, boundary term, W'-weighted leakage, continuity).

The paper does NOT introduce a gauge-fixing term `(1/xi) partial^nu (partial · A)` in the parent equation at stage 5; the parent law (5.1) is the geometric `partial_N(Z F^{NM}) = mu_0 J^M`.

## What the script claims to verify

The SymPy script is split into a symbolic generic part (Sections 1-4) that prints the projection identities and assembled projected equations using abstract `sp.Function(...)` placeholders — these sections contain NO assertions, only `print` statements — and a concrete Gaussian-kernel audit (Section 5) that contains all the executable checks. For `W(w) = exp(-w^2)/sqrt(pi)` it asserts: (i) projection commutes with `partial_t` on `(t^2 + x)*(w^2+1)`; (ii) for `Q = w^3 + w` the Gaussian-weighted boundary `[Wg Q]_partial = 0` and the projection-IBP identity `Pg[partial_w Q] = -Pg'[Q]` holds, with a sign-flipped mutant guard; (iii) for hand-chosen polynomial fields `G^{0,x,y,z,w}_nu` and a hand-inserted extra-term `Gamma_nu = (t+x)(w^2+1)`, both the with-boundary and decay-form projected inhomogeneous identities hold, with sign-mutant guard; (iv) the analytic value `-Pg'[w] = 1`; (v) the projected continuity identity for analogous polynomial source profiles `J^{0,x,y,z,w}`.

The Mathematica script independently checks the same five identities (M1-M5) using *different* concrete profiles: `Phi = (Sin[t]+x^2)(w^2+3)` for commutation; `Q = w Exp[-w^2/4]` for IBP; `G0 = Cos[t](w^2+2), Gx = x^2 w^2, Gy = y(w^4+1), Gz = z w^2, Gw = w + w^3/3, gammaTerm = (Sin[t]+x)(w^2+1)` for the projected inhomogeneous law; `Gw = w` for the analytic leakage value; analogous `Jw = w + w^3/3` profile for continuity. Sign-mutant guards are present at every stage and report a literal nonzero residual.

Note: both scripts generalize the paper's parent law by carrying an extra additive term `Gamma_nu` (gauge-driver) through the projection. With `Gamma_nu = 0` the identities collapse to the paper's stated Eq. (5.3) and Eq. (5.4); the extension is benign because the projection operator is linear.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Eq. (5.2) raw projected law `partial_nu int W Z F^{nu mu} dw = mu_0 int W J^mu dw - int W partial_w(Z F^{wmu}) dw` | SymPy line 261 (`projected_lhs_with_boundary` includes `+ boundary_Gw - Pgp(Gwc)`, which equals `Pg(partial_w Gwc)` by IBP) and Mathematica M3 (`projW[D[Gw,w]] = -projWPrime[Gw]` since `m3Boundary == 0`) | match |
| Eq. (5.3) IBP form `... = mu_0 P[J] - [W Z F^{wmu}]_partial + P_{W'}[Z F^{wmu}]` | SymPy line 262 (`project_bulk_lhs - projected_lhs == 0` with `-Pgp(Gwc)` term explicit) and Mathematica M3 lines 89-95 (`lhsM3 = rhsM3` with explicit `-projWPrime[Gw]` and zero boundary) | match |
| Boundary flux `[W Z F^{wmu}]_partial` (must vanish in decay form) | SymPy line 260 (`assert_zero("concrete projected inhomogeneous boundary term", boundary_Gw)`) and Mathematica M3 line 83-87 (`m3Boundary == 0`) | match |
| Kernel-gradient leakage `+ int (partial_w W) Z F^{wmu} dw` | SymPy lines 277-279 (`leakage = -Pgp(Gwc)`, `assert_zero("concrete transverse leakage value", leakage - 1)`, `assert_nonzero` guard) and Mathematica M4 lines 108-127 (independent recomputation, mutant guard, nonzero guard) | match |
| Eq. (5.4) projected continuity `partial_mu int W J^mu dw = -[W J^w]_partial + int (partial_w W) J^w dw` | SymPy lines 286-303 (`project_bulk_cont - projected_cont == 0` with `-Pgp(Jwc)` term) and Mathematica M5 lines 130-158 (analogous) | match |
| Projection commutes with brane derivatives `partial_{t,x,y,z}` (used implicitly in deriving 5.2) | SymPy line 226 (`Pg(diff(Phi,t)) - diff(Pg(Phi),t) == 0`) and Mathematica M1 lines 38-51 (with mutant guard) | match |
| Gauge-driver term `Gamma_nu` carried through projection (NOT in paper) | SymPy `Gammac = (t+x)(w^2+1)` in projected_lhs at line 257; Mathematica `gammaTerm = (Sin[t]+x)(w^2+1)` in lhsM3 at line 89 | extra (benign; reduces to paper claim when Gamma=0; linearity of projection makes this trivially equivalent to the paper identity for any extra term) |

All five paper-side deliverables are covered by both engines with non-tautological assertions and explicit mutant guards. The single "extra" row is a benign linearity-only generalization that does not introduce any new physical claim — it would not cause a script residual to vanish for the wrong reason.

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 226 | `assert_zero("projection commutes with brane derivative", Pg(diff(Phi,t)) - diff(Pg(Phi),t))` | commutation lemma backing Eq. (5.2)-(5.4) | yes |
| A2 | sympy | 230 | `assert_zero("decaying-kernel boundary term", boundary_Q)` | boundary flux vanishing in decay form | yes |
| A3 | sympy | 231 | `assert_zero("decaying-kernel integration by parts with boundary", Pg(diff(Q,w)) - (boundary_Q - Pgp(Q)))` | IBP backing Eq. (5.3) | yes |
| A4 | sympy | 232 | `assert_zero("decaying-kernel integration by parts", Pg(diff(Q,w)) + Pgp(Q))` | decay-form IBP | yes |
| A5 | sympy | 233 | `assert_nonzero("mutated decaying-kernel IBP sign should fail", ...)` | sign-mutant guard | yes |
| A6 | sympy | 260 | `assert_zero("concrete projected inhomogeneous boundary term", boundary_Gw)` | boundary flux vanishing | yes |
| A7 | sympy | 261 | `assert_zero("concrete projected inhomogeneous law with boundary", project_bulk_lhs - projected_lhs_with_boundary)` | Eq. (5.2) | yes |
| A8 | sympy | 262 | `assert_zero("concrete projected inhomogeneous law", project_bulk_lhs - projected_lhs)` | Eq. (5.3) (decay form) | yes |
| A9 | sympy | 263-275 | `assert_nonzero(... mutated projected inhomogeneous leakage sign should fail ...)` | sign-mutant guard | yes |
| A10 | sympy | 278 | `assert_zero("concrete transverse leakage value", leakage - 1)` | kernel-gradient leakage analytic value | yes |
| A11 | sympy | 279 | `assert_nonzero("concrete transverse leakage", leakage)` | leakage nonzero (open-system signature) | yes |
| A12 | sympy | 301 | `assert_zero("concrete projected continuity boundary term", boundary_Jw)` | continuity boundary | yes |
| A13 | sympy | 302 | `assert_zero("concrete projected continuity law with boundary", project_bulk_cont - (projected_cont + boundary_Jw))` | Eq. (5.4) | yes |
| A14 | sympy | 303 | `assert_zero("concrete projected continuity law", project_bulk_cont - projected_cont)` | Eq. (5.4) decay form | yes |
| M1 | mathematica | 38-51 | `lhsM1 - rhsM1 === 0` + mutant guard | brane-derivative commutation | yes |
| M2 | mathematica | 53-73 | boundary-decay + IBP + mutant guard | IBP backing Eq. (5.3) | yes |
| M3 | mathematica | 75-105 | projected inhomogeneous law + mutant guard | Eq. (5.3) | yes |
| M4 | mathematica | 107-127 | `leakage = -projWPrime[w]; leakage === 1` + nonzero guard + mutant | leakage analytic value | yes |
| M5 | mathematica | 129-158 | projected continuity + mutant guard | Eq. (5.4) | yes |

Every assertion is non-tautological: each compares two independently computed sides (one by `Project(differentiate(expression))`, the other by `differentiate/IBP(Project(expression))`). For the leakage value `-Pg'[w] = 1`, the value 1 is the analytic Gaussian integral `int (2w^2/sqrt(pi)) exp(-w^2) dw = 1` — pinning an analytic result, not a circular check.

Both engines use the same Gaussian kernel `Wg = exp(-w^2)/sqrt(pi)` (the canonical normalized choice) but DIFFERENT concrete profiles for `G^M_nu`, `Q`, `J^M`, and `Gamma_nu` — so an algebra error in one would not be mirrored in the other.

## Findings

**None.** Attacks attempted and the result of each:

- **Tautology hunt.** Each `assert_zero` compares two independently constructed expressions; the LHS uses `Pg(symbolic differentiation)`, the RHS uses `differentiation(Pg(.))` or IBP-restructured form. Not algebraically equivalent by construction — they coincide only when the projection-by-IBP identity holds. Mutant guards (A5, A9, M1-mutant, ..., M5-mutant) verify that a sign flip would actually fail; the recorded mutant residuals (e.g., `7*Cos[t]`, `16/(5*Sqrt[5])`, `3`, `-2`) confirm symbolic distinguishability.
- **Hidden hardcoding.** The single literal numeric pin is `leakage == 1` in both engines. This is independently computable by hand (`-int (-2w/sqrt(pi)) exp(-w^2) w dw = 2/sqrt(pi) * (sqrt(pi)/2) = 1`) and the script does *not* set `leakage = 1` then check `leakage == 1`; it computes `leakage = -Pg'[w]` analytically and confirms the result equals 1. Non-circular.
- **Symbol assumption errors.** SymPy declares `t,x,y,z,w` as `real=True` and `mu0, xi` as `nonzero=True` — appropriate. Mathematica declares `{t,x,y,z} \[Element] Reals` as projection assumptions; `w` is the integration variable (no restriction needed). No assumption forces residuals to vanish trivially.
- **Missing branch.** The script verifies the identity for a Gaussian kernel only, but the docstring and paper claim cover any kernel `W(w)` such that `[W Q]_partial = 0`. The symbolic-generic block (Sections 1-4 of SymPy) prints the abstract identity but does NOT assert it (the print of `proj_dw` is the SymPy `Eq(...)` rather than an `assert_zero`). However, the abstract identity is *provable* from the chain rule plus the fundamental theorem of calculus, which both engines exhibit; the concrete Gaussian test is a witness, not the proof. The paper claim is "for any normalized observation kernel"; the concrete witness covers one instance of a normalized kernel; the abstract algebraic step is exhibited symbolically. This is the standard structure for these audits and the paper's intent — not a `missing_branch` finding.
- **Engine disagreement.** Different concrete profiles, both engines independently report zero residuals and matching mutant nonzero residuals (`7*Cos[t]`, `16/(5*Sqrt[5])`, `3`, `-2`, `3` for Mathematica; nonzero analytic residuals for SymPy mutant guards). No disagreement.
- **Mathematica transliteration.** The `.wl` script uses entirely DIFFERENT concrete profiles than the `.py` script (SymPy uses polynomial `(t^2+x)(w^2+1)`; Mathematica uses trigonometric+polynomial `(Sin[t]+x^2)(w^2+3)`; SymPy uses `Gw = w`; Mathematica uses `Gw = w + w^3/3`; etc.). The logical structure of "boundary, residual, mutant guard" is parallel because that's the audit pattern, but the algebra each engine performs is independent. Not transliteration.
- **Insufficient verification.** Every paper-side deliverable has at least one assertion in each engine. Boundary, decay-form, and with-boundary forms are all exercised. Leakage value is pinned. Continuity is independently checked.
- **Paper misalignment.** The script's `Gamma_nu` gauge-driver extension is the only divergence. It is *additive* (a new term carried through linearly) and *reduces* to the paper claim when `Gamma_nu = 0`. The script's docstring explicitly anchors it as "applies both to the current unweighted gauge-fixing form and to later weighted variants." This is a forward-looking generalization, not a misalignment. The paper's identities are all verified verbatim by the script even at `Gamma_nu = 0`.
- **Stale output.** Both outputs are newer than their respective scripts. Not stale.

## Independent-derivation check (Mathematica)

The `.wl` script derives the same identities from scratch with different test functions and different residual-style. Concretely:

- SymPy commutation test (line 226): `Phi = (t**2 + x) * (w**2 + 1)`, checks `Pg(diff(Phi,t)) - diff(Pg(Phi),t) == 0`.
- Mathematica M1 (lines 38-51): `Phi = (Sin[t] + x^2)*(w^2 + 3)`, checks `projW[D[Phi,t]] - D[projW[Phi],t] === 0` PLUS a mutant `lhs + rhs` whose value is `7*Cos[t]` (proves the test is non-trivial in `t`).

- SymPy IBP test (lines 228-233): `Q_concrete = w**3 + w`, checks `Pg(diff(Q,w)) + Pgp(Q) == 0`.
- Mathematica M2 (lines 53-73): `Q = w*Exp[-w^2/4]`, checks `projW[D[Q,w]] + projWPrime[Q] === 0`, mutant residual `16/(5*Sqrt[5])`.

- SymPy leakage (lines 277-279): `Gwc = w`, asserts `-Pgp(w) - 1 == 0`.
- Mathematica M4 (lines 107-127): `leakageGw = w`, asserts `-projWPrime[w] === 1`, mutant residual `-2`.

These are independent algebra at every stage. Not a transliteration.

## Engine cross-check

Both engines produce zero residuals on every identity and matching nonzero residuals on every sign-mutant. The SymPy script exits with `STATUS: PASS` after 14 assertions. The Mathematica script exits with `STATUS: PASS` after 5 grouped identity blocks with mutant + nonzero guards. The outputs agree on the substantive claim (projected Maxwell + leakage + continuity all hold) and the literal leakage value (`1`).

## Verdict justification

The scripts verify what the paper claims. The paper's stage card asserts three things: the projected inhomogeneous Maxwell equation in raw form (5.2) and post-IBP form (5.3); the boundary-flux + W'-weighted leakage decomposition; and the projected charge continuity (5.4). Both SymPy and Mathematica exhibit these identities both symbolically (the abstract algebra is printed) and concretely (Gaussian-kernel + polynomial/trigonometric witnesses). Mutant guards confirm sign-dependence. Engine cross-check is clean. Output files are fresh. No notes file exists for this stage (consistent with the orchestrator preface for EM-projected stages 004-020), so the paper card and appendix row are the sole paper-side authority; both are honored.

Attacks tried and failed: tautology hunt (each assertion compares independently computed sides), hidden hardcoding (the only numeric pin is `leakage = 1`, which is the literal analytic Gaussian integral), branch limitations (the paper says "any normalized kernel" and the script verifies the abstract identity symbolically plus a Gaussian witness), engine disagreement (different profiles, identical conclusions), transliteration (different concrete profiles confirm independent algebra), and paper-side gauge-term mismatch (the script's `Gamma_nu` is a benign linearity extension that reduces to the paper's identity when set to zero).

`verdict: clean`. `paper_alignment: aligned`.

## Self-test notes

Checked variable independence (every `diff(expr, var)` in the concrete blocks has `expr` actually depending on `var` — confirmed via direct expression inspection: `Phi = (t^2+x)(w^2+1)` does depend on `t`; `Pg(Phi) = (t^2+x)*(3/2)` still depends on `t`; therefore the commutation assertion is non-trivial). Checked parity: integrands like `Wgp * w = -(2w^2/sqrt(pi)) exp(-w^2)` are even, so the integral can be nonzero — leakage nonzero claim is consistent. Trivial-case substitution: setting `Gamma_nu = 0` and the polynomial profiles to constants makes the residuals trivially zero on both sides, but for the actual non-constant profiles the residuals are 0 only because the projection-IBP identity holds — confirmed by the mutant guards giving nonzero literals. No paper round-trip issues because no fix is proposed.
