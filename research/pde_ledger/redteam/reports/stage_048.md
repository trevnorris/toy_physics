---
unit_id: 048
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage048_support_compensation_theorem.md"]
  paper_appendix: present
---

# Audit unit 048 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_048.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage048_support_compensation_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 74; `\input` at line 214)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.txt`

## What the paper claims

The stage proves there is no reduced-level support no-go on the coherent tracking branch. The card's `\stagefield{Output}` states verbatim: "The exact threshold \eqref{eq:app-stage048-zeta-req} and success condition \eqref{eq:app-stage048-success-iff}." The two boxed equations are
`zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]` (inverse of `S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)` evaluated at `S = S_req = M_req / M_mix`) and the success condition `zeta_phys >= zeta_req`. The notes elaborate the scaffolding: (i) the tracking law `G_tr(xi,delta;R) = 9 xi (xi+delta) / [9 delta + (9+2 R^2) xi]` is strictly increasing on `0<xi<1` and has critical load `M_crit(delta,R) = 9(1+delta)/(9 delta + 9 + 2 R^2)`; (ii) the normalization function `F_tr` satisfies `F_tr(0,delta;R) = 1` and `lim_{xi->1^-} F_tr = +infty`, so by IVT every finite `R_target>1` is reached by some `xi_req in (0,1)`; (iii) `S(zeta;eps)` is strictly increasing with `S(0;eps)=1`, `lim_{zeta->1/eps^-} S = +infty`, and pole margin `1/eps - zeta_req = (1-eps)/(eps[1+eps(S_req-2)]) > 0`; (iv) branch margin `zeta_crit - zeta_req = (S_crit - S_req)(1-eps) / ((1+eps(S_crit-2))(1+eps(S_req-2))) > 0`; (v) implicit derivative `dxi_phys/dzeta = M_mix (dS/dzeta) / (dG_tr/dxi) > 0`. Per `\stagefield{Downstream use}`, stages 049-057 actually derive `zeta_phys`; stage 048 only sets up the threshold, not the realized value. The part appendix row (line 74) summarises the deliverable as "Unique `zeta_req` threshold and no reduced-level support no-go." Inputs are listed as `S(zeta;epsilon)`, `M_mix`, and `M_req` — all carried into the script as symbolic.

## What the script claims to verify

The SymPy docstring lists four checks: (1) exact tracking-branch critical load and monotonicity in `xi`; (2) exact support-enhancement inverse map `zeta(S)`; (3) exact support-feasibility formulas `zeta_req` and `zeta_crit`; (4) exact stability margins `zeta_req < 1/eps` and `zeta_req < zeta_crit` when `S_req < S_crit`. The Mathematica script mirrors these and adds two strengthening endpoint checks: a `(1-xi) F_tr` softening-coefficient identity at `xi -> 1^-`, and a `(1/eps - zeta) S` pole-coefficient identity at `zeta -> 1/eps^-`; it also derives `zeta_req` via `Solve[S == sReq, zeta]` rather than positing the formula. Both engines additionally verify the implicit-derivative formula for `dxi_phys/dzeta` (notes section 5). All assertions take the form `expect_zero(residual)` / `expectZero[residual]`, i.e. they confirm SymPy/Mathematica-computed quantities match canonical hand-written closed forms. Note: both scripts internally label themselves "STAGE 31" / "STAGE 031" (sympy:3, sympy:31, math:26), tracking the notes' historical numbering; the final-line print and the filenames correctly identify the unit as 048.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Inverse-map formula `zeta_req = (S_req-1)/(1+eps(S_req-2))` (eq. app-stage048-zeta-req) | sympy:87 `S(zeta_req)-S_req == 0`; math:86-88 `Solve[S==sReq,zeta]` + math:93 verifies | match |
| Success iff `zeta_phys >= zeta_req` (eq. app-stage048-success-iff) | Indirect: established via `dS/dzeta > 0` formula (sympy:74, math:77) and inverse identity. `zeta_phys` itself is properly downstream per `\stagefield{Downstream use}`. | match (correct scope) |
| `G_tr` strictly increasing on stable branch (notes §1) | sympy:52-56 + math:49-52 verify `dG_tr/dxi` matches canonical sum-of-squares form (positivity manifest from form) | match |
| `M_crit = 9(1+delta)/(9 delta + 9 + 2 R^2)` (notes §1) | sympy:47 via `sp.limit`; math:40-43 via `Limit[..,Direction->"FromBelow"]`; `M_crit - G_tr` identity asserted (sympy:61-66, math:68-72) | match |
| `F_tr(0,delta;R) = 1` (notes §3) | sympy:57, math:53 assert | match |
| `lim_{xi->1^-} F_tr = +infty` (notes §3, load-bearing for IVT existence of `xi_req`) | sympy:58 printed only, NOT asserted; math:55-64 asserts via softening-coefficient identity at `xi=1` | partial (sympy gap → F1) |
| `S(0;eps) = 1` (notes §2) | sympy:73, math:76 assert | match |
| `dS/dzeta = (1-eps)/(1-zeta eps)^2` (notes §2) | sympy:74, math:77 assert formula | match |
| `lim_{zeta->1/eps^-} S = +infty` (notes §2) | sympy:79-82 asserts under physical `eps = 1/(1+nu)`, `nu>0`; math:79-84 asserts via pole-coefficient identity | match |
| Pole margin `1/eps - zeta_req = (1-eps)/(eps(1+eps(S_req-2)))` (notes §2) | sympy:99-102, math:100 assert | match |
| Branch margin `zeta_crit - zeta_req = (S_crit-S_req)(1-eps)/(...)` (notes §4) | sympy:103-106, math:101-104 assert | match |
| Implicit `dxi_phys/dzeta` formula (notes §5) | sympy:110-117, math:106-112 assert formula | match (positivity manifest from form given `0<eps<1`) |

Dominant pattern: every paper-side equation, both Output items, and every notes-elaborated identity is exercised by at least one engine, with both engines agreeing on the algebraic forms. One partial: the SymPy script omits an explicit assertion of `lim_{xi->1^-} F_tr = +infty` even though it computes and prints the limit. Mathematica covers it independently via softening coefficient. Front-matter `paper_alignment: aligned`.

## Assertion inventory

| #   | Script      | Line   | Form | Exercises which paper claim? | Anchored to claim? |
|-----|-------------|--------|------|------------------------------|--------------------|
| A1  | sympy       | 52-56  | `expect_zero(dG_dxi - canonical)` | notes §1 `G_tr` monotonicity | yes |
| A2  | sympy       | 57     | `expect_zero(F_tr(xi=0)-1)` | notes §3 left endpoint | yes |
| A3  | sympy       | 58     | print only of `limit xi->1- F_tr` | notes §3 right endpoint | no (not asserted) |
| A4  | sympy       | 61-66  | `expect_zero(M_crit-G_tr - canonical)` | notes §1 critical-load gap | yes |
| A5  | sympy       | 73     | `expect_zero(S(zeta=0)-1)` | notes §2 S left endpoint | yes |
| A6  | sympy       | 74     | `expect_zero(dS/dzeta - canonical)` | notes §2 S monotonicity | yes |
| A7  | sympy       | 75     | print only of `limit zeta->1/eps- S` | informational (covered by A8) | n/a |
| A8  | sympy       | 79-82  | `if limit_phys != sp.oo: raise AssertionError` under physical eps | notes §2 S pole divergence | yes |
| A9  | sympy       | 87     | `expect_zero(S(zeta_req)-S_req)` | **paper Output 1 inverse-map** | yes |
| A10 | sympy       | 88     | `expect_zero(S(zeta_crit)-S_crit)` | notes §4 zeta_crit inverse | yes |
| A11 | sympy       | 99-102 | `expect_zero(pole margin - canonical)` | notes §2 pole stability margin | yes |
| A12 | sympy       | 103-106 | `expect_zero(branch margin - canonical)` | notes §4 branch stability margin | yes |
| A13 | sympy       | 112-117 | `expect_zero(dxi/dzeta - canonical)` | notes §5 implicit derivative | yes |
| B1  | mathematica | 49-52  | `expectZero[dG_dxi - canonical]` | notes §1 G_tr monotonicity | yes |
| B2  | mathematica | 53     | `expectZero[F_tr(xi=0)-1]` | notes §3 left endpoint | yes |
| B3  | mathematica | 55-64  | `expectZero[(1-xi) F_tr softening coeff - canonical]` | notes §3 right endpoint (strengthened) | yes |
| B4  | mathematica | 68-72  | `expectZero[M_crit-G_tr - canonical]` | notes §1 critical-load gap | yes |
| B5  | mathematica | 76     | `expectZero[S(zeta=0)-1]` | notes §2 S left endpoint | yes |
| B6  | mathematica | 77     | `expectZero[dS/dzeta - canonical]` | notes §2 S monotonicity | yes |
| B7  | mathematica | 79-84  | `expectZero[(1/eps-zeta) S pole coeff - (1-eps)/eps^2]` | notes §2 S pole divergence | yes |
| B8  | mathematica | 86-88  | `Solve[S==sReq,zeta]` uniqueness then extract zetaReq | **paper Output 1 inverse-map (independent derivation)** | yes |
| B9  | mathematica | 93     | `expectZero[S(zeta_req)-S_req]` | **paper Output 1 inverse-map** | yes |
| B10 | mathematica | 94     | `expectZero[S(zeta_crit)-S_crit]` | notes §4 zeta_crit inverse | yes |
| B11 | mathematica | 100    | `expectZero[pole margin - canonical]` | notes §2 pole stability margin | yes |
| B12 | mathematica | 101-104 | `expectZero[branch margin - canonical]` | notes §4 branch stability margin | yes |
| B13 | mathematica | 108-112 | `expectZero[dxi/dzeta - canonical]` | notes §5 implicit derivative | yes |

No tautologies: each `expect_zero` LHS is reached by a non-trivial chain (`sp.diff` + `sp.factor`, `sp.simplify(sp.limit(...))`, `Solve`, etc.) and the RHS is a distinct hand-written canonical form. Replacing the RHS with anything but the true algebraic value would fail the residual check. Every script-side assertion traces to a specific paper-side deliverable (no orphan checks); no paper-side deliverable lacks at least one engine's coverage (the F_tr right-endpoint gap exists only in SymPy, not in Mathematica).

## Findings

### F1 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py:58`

**What's wrong:**
The notes (section 3, lines 124-126) declare two endpoint facts as exact and load-bearing for the IVT existence argument of `xi_req`:

> `F_tr(0,delta;R) = 1,`
> `lim_(xi -> 1^-) F_tr(xi,delta;R) = +infinity.`

The first is asserted in SymPy at line 57:

```python
expect_zero("F_tr(xi=0)-1", sp.simplify(F_tr.subs(xi, 0) - 1))
```

The second is computed and printed only, at line 58:

```python
print("limit xi->1^- of F_tr =", sp.limit(F_tr, xi, 1, dir="-"))
```

There is no `expect_zero`, no `raise AssertionError`, no guard at all on the value of this limit. The output transcript line 22 records `limit xi->1^- of F_tr = oo`, but a typo or future SymPy version change could silently make this finite or `-oo` without failing the script (`Exit code: 0` would still print).

The Mathematica counterpart (lines 55-64) asserts the strengthened softening-coefficient identity `(1-xi) F_tr | xi->1^- == (9 delta + 9 + 2 r^2)^2 (9 delta + 9 + 2 r)^2 / [81 (9 delta^2 + 18 delta + 9 + 2 r^2)^2]`, which simultaneously certifies the divergence of `F_tr` and pins down the leading coefficient. This is the engine asymmetry the v2 cross-check surfaces.

**Why this matters:**
The IVT existence argument for `xi_req in (0,1)` (notes §3) explicitly relies on both endpoints: `F_tr` starts at 1 and goes to `+infinity`. Asserting only the left endpoint leaves the right endpoint as a print statement that a future code change could silently break. The paper card itself (line 24) states "the support factor is strictly increasing up to its softening pole" — this echoes the divergence claim, which the SymPy script does not actually guard.

**Required change:**
Add a guarded check for the divergence in the SymPy script, immediately after line 58. Mirror the strengthening that Mathematica uses (line 55-64), so the SymPy script also verifies the softening coefficient. Concretely, after the existing `print(...)` line, add a `softCoeff` computation and an `expect_zero` against the closed form, e.g.:

```python
soft_coeff = sp.simplify(sp.limit((1 - xi) * F_tr, xi, 1, dir="-"))
print("softening coefficient for F_tr =", soft_coeff)
expect_zero(
    "(1-xi) F_tr softening coefficient",
    soft_coeff
    - (9 * delta + 9 + 2 * R ** 2) ** 2 * (9 * delta + 9 + 2 * R) ** 2
    / (81 * (9 * delta ** 2 + 18 * delta + 9 + 2 * R ** 2) ** 2),
)
```

This matches the Mathematica `softCoeffExpected` formula (math:59-62) exactly. No change required to paper, notes, or Mathematica.

**Verification:**
After Codex applies, the verifier runs `redteam exec-sympy 048`; the new check appears at sympy:~59-67 (immediately after the existing `print` at line 58), the output transcript shows a new line `softening coefficient for F_tr = ...` and `(1-xi) F_tr softening coefficient = 0`, and the script still exits 0. Cross-engine: the SymPy `soft_coeff` value should textually match the Mathematica output line 20.

## Independent-derivation check (Mathematica)

The Mathematica script is not a transliteration of the SymPy script. Three concrete pieces of independent algebra:

1. **Inverse-map derivation.** SymPy posits `zeta_req = (Sreq - 1) / (1 + eps * (Sreq - 2))` directly at line 84 and then asserts `S(zeta_req) - Sreq == 0`. Mathematica instead invokes `zetaSolutions = Solve[sEnhance == sReq, zeta, Reals]` (line 86), checks uniqueness via `Length[zetaSolutions] != 1` (line 87), and only then extracts `zetaReq = FullSimplify[zeta /. First[zetaSolutions], ...]` (line 88). That is an independent derivation of the inverse, not a re-statement.
2. **Endpoint divergence checks.** SymPy prints `sp.limit(F_tr, xi, 1, dir="-")` and `sp.limit(S, zeta, 1/eps, dir="-")` for visual inspection; Mathematica asserts the strengthening `(1-xi) F_tr | xi->1` softening coefficient (lines 55-64) and `(1/eps - zeta) S | zeta->1/eps` pole coefficient `(1-eps)/eps^2` (lines 79-84) as `expectZero` identities. These checks have no SymPy analogue and use independent `Limit[..,Direction->"FromBelow"]` calls.
3. **Assumption choreography.** Mathematica declares the full physical window `$Assumptions = ... && 0 < xi < 1 && delta > 0 && r > 0 && 0 < eps < 1 && zeta > 0 && mMix > 0 && sCrit > sReq && sReq > 1` (lines 29-32); SymPy declares only `positive=True, real=True` for each symbol (lines 33-36). The two engines therefore explore the assumption space differently; their agreement on residuals is a real cross-check, not a copy.

No `mathematica_transliteration` finding warranted.

## Engine cross-check

Both transcripts report `EXIT_CODE: 0`, `Status: PASS`. Where both engines test the same identity, both report residual `0`:

| Identity | SymPy residual | Mathematica residual |
|---|---|---|
| `dG_tr/dxi - canonical` | 0 | 0 (PASS) |
| `F_tr(xi=0) - 1` | 0 | 0 (PASS) |
| `M_crit - G_tr - canonical` | 0 | 0 (PASS) |
| `S(zeta=0) - 1` | 0 | 0 (PASS) |
| `dS/dzeta - (1-eps)/(1-zeta eps)^2` | 0 | 0 (PASS) |
| `S(zeta_req) - S_req` | 0 | 0 (PASS) |
| `S(zeta_crit) - S_crit` | 0 | 0 (PASS) |
| `pole margin - canonical` | 0 | 0 (PASS) |
| `branch margin - canonical` | 0 | 0 (PASS) |
| `dxi/dzeta - canonical` | 0 | 0 (PASS) |

Side-by-side closed forms agree up to trivial sign/order rearrangement:
- `G_tr`: sympy `9*xi*(delta + xi)/(9*delta + xi*(2*R**2 + 9))` (sympy output line 17) ≡ math `(9*xi*(delta + xi))/(9*delta + (9 + 2*r^2)*xi)` (math output line 13).
- `M_crit`: sympy `9*(delta + 1)/(2*R**2 + 9*delta + 9)` (line 18) ≡ math `(9*(1 + delta))/(9 + 9*delta + 2*r^2)` (line 14).
- `S(zeta;eps)`: sympy `(2*eps*zeta - zeta - 1)/(eps*zeta - 1)` (line 29) ≡ math `1 + ((-1 + eps)*zeta)/(-1 + eps*zeta)` (line 26).
- `zeta_req`, `zeta_crit`, pole margin, branch margin, `dxi/dzeta`: identical canonical forms in both transcripts (lines 40-50 vs lines 34-46).

`engines_agree: true`.

Output freshness: sympy script mtime `May 11 12:48`, sympy output mtime `May 11 12:48` (output equal/newer than script); mathematica script mtime `May 11 11:56`, mathematica output mtime `May 11 12:51` (output strictly newer). Both `outputs_fresh: true`.

## Verdict justification

The script faithfully exercises the paper's stated Output: the inverse-map formula for `zeta_req` (paper eq. app-stage048-zeta-req) is non-tautologically verified by both engines (SymPy by substitution; Mathematica also by independent `Solve`), and the surrounding monotonicity / stability-margin structure that justifies the `zeta_phys >= zeta_req` success condition (paper eq. app-stage048-success-iff) is checked across both engines via explicit canonical-form `expect_zero` identities. The success-condition equation `zeta_phys >= zeta_req` is itself a corollary of `S` strictly increasing and the inverse identity, both of which are verified; `zeta_phys` is correctly downstream (notes §6, paper `\stagefield{Downstream use}` stages 049-057). Banner labelling drift ("STAGE 31"/"STAGE 031" while the unit is 048) is internal to the notes' numbering convention, not a math issue and not in any audit category. Adversarial attacks attempted: (i) probed `S(zeta_req) - S_req` for tautology — fails, because `zeta_req` is a hypothesized closed form and the substitution can land off-target if the hypothesis is wrong; (ii) probed the `dxi/dzeta` sign — SymPy outputs `-M_mix*(eps - 1)*(...)^2 / ...` which is `+M_mix*(1-eps)*(...)^2/...` and matches Mathematica's `-(-1 + eps)*mMix*(...)^2/...`; (iii) probed symbol-domain holes — SymPy lacks `eps < 1` and `S_req > 1` declarations, but every SymPy assertion is a pure algebraic identity of residual = 0, so the missing assumption changes nothing (it would only matter if the script were asserting sign, which it does not — sign is left manifest in the canonical-form structure); (iv) probed the inverse-map for missing branch — `Solve` returns a single solution per math:87 `Length == 1` check, and SymPy's posited closed form matches; no other branch exists in `Reals`. The only real gap surfaced is F1: the SymPy script's failure to assert `lim_{xi->1^-} F_tr = +infty`, which Mathematica does assert via softening coefficient. Verdict: `findings` with one low-severity `insufficient_verification`. No stop-cold.

## Self-test notes

I mentally walked through the proposed F1 patch: (1) Variable independence — `(1-xi)*F_tr` depends on `xi`, so `sp.limit(..., xi, 1, dir="-")` is non-trivial and equals the closed form I prescribed (numerator `(9 delta + 9 + 2 R^2)^2 (9 delta + 9 + 2 R)^2` over `81 (9 delta^2 + 18 delta + 9 + 2 R^2)^2`); cross-checked against Mathematica's `softCoeffExpected` (math:59-62) and the captured Mathematica output line 20. (2) Parity/symmetry — no integrals; the limit `xi -> 1^-` is the relevant one-sided approach. (3) Trivial-case sub — at `delta=1, R=1` the prescribed RHS evaluates to `(9+9+2)^2 (9+9+2)^2 / (81 (9+18+9+2)^2) = 20^2 * 20^2 / (81 * 38^2)`, a finite positive rational, confirming `F_tr` itself diverges as `1/(1-xi)`. (4) Path — `.py` lives in `scripts/`, no path ambiguity. (5) Paper round-trip — the patch only adds a check that the SymPy `print` at line 58 already exposes and that the Mathematica side already asserts; introduces no new paper claims, alters no existing paper claim. Self-test passes.
