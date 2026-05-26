---
unit_id: 041
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage041_rank2_support_completion.md
  paper_appendix: present
---

# Audit unit 041 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_041.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage041_rank2_support_completion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 60 references this stage; full body input at line 200)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage041_rank2_support_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage041_rank2_support_mathematica_audit.txt`

## What the paper claims

The paper card's `\stagefield{Output}` reads: "The support-loading theorem \eqref{eq:app-stage041-nreq}, its monotonicity \eqref{eq:app-stage041-nreq-derivative}, and the tracking/source-tied alternatives." Concretely the stage is asserting four deliverables:

1. The loaded-stiffness determinant decomposition `D_sel = xi(delta+xi) - m[delta+(1+q^2)xi] - n[delta+(1+r^2)xi] + m n (q-r)^2` (eq:app-stage041-Dsel). The linearity in `n` is the structural claim.
2. The exact support-loading theorem `n_req = [xi(delta+xi) - m(delta+(1+q^2)xi)] / [delta+(1+r^2)xi - m(q-r)^2]` (eq:app-stage041-nreq).
3. The monotonicity `dn_req/dm = -[delta+(1+qr)xi]^2 / [delta+(1+r^2)xi - m(q-r)^2]^2 < 0` (eq:app-stage041-nreq-derivative).
4. Two specializations: tracking (`r=q ⇒ n_req = G_q(xi,delta) - m`, eq:app-stage041-tracking) and source-tied (`q = t R_U`, `r = t`, `t^2 = lambda_0`, yielding eq:app-stage041-source-tied).

The notes additionally call out two source-tied feasibility windows — regularity-denominator positivity (`m < (delta+(1+lambda_0)xi)/(lambda_0(R_U-1)^2)`) and a positive-support ceiling (`m <= xi(delta+xi)/(delta+(1+lambda_0 R_U^2)xi)`) — neither of which is in `\stagefield{Output}` (they are derived consequences). The Part-III appendix row classifies the stage as "Exact closure: Linear determinant in support load and exact `n_req` theorem."

## What the script claims to verify

Both engines independently:

1. Build the dimensionless 2x2 loaded matrix `M = diag(1, 1+delta) - m(1,q)^T(1,q) - n(1,r)^T(1,r)`, compute `Det(M - (1-xi)I)`, and assert it equals the paper's `D_sel` decomposition.
2. Solve `D_sel = 0` for `n` and assert the solver result matches the paper's closed `n_req`.
3. Differentiate the closed-form `n_req` with respect to `m` and assert it equals `-(delta+(1+qr)xi)^2 / [delta+(1+r^2)xi - m(q-r)^2]^2`.
4. Substitute `r -> q` and assert `n_req(r=q) = G_q - m`.
5. Substitute `q -> t R_U`, `r -> t` into the general `n_expected`, then `t^2 -> lambda_0 = 2/9`, and assert the result matches the source-tied closed form.
6. Print (informational, not asserted) the `regThreshold` and `numZeroThreshold` feasibility-window expressions from the notes.
7. Differentiate the source-tied closed form and assert it equals `-(delta+(1+lambda_0 R_U)xi)^2 / [delta+(1+lambda_0)xi - m lambda_0 (R_U-1)^2]^2`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) Determinant decomposition eq:app-stage041-Dsel | `expect_zero("determinant decomposition", Det - Det_expected)` (sympy line 72; wl line 56) | match |
| (2) `n_req` closed form eq:app-stage041-nreq | `expect_zero("n_req - expected", ...)` (sympy line 82; wl line 65) | match |
| (3) Monotonicity `dn_req/dm` closed form eq:app-stage041-nreq-derivative | `expect_zero("dn/dm - expected", ...)` (sympy line 94; wl line 76) | match (closed form; sign `<= 0` is manifest from `-X^2/Y^2`, strict `<0` follows on the regular branch the paper qualifies) |
| (4a) Tracking eq:app-stage041-tracking | `expect_zero("tracking collapse", n_track - (G_q - m))` (sympy line 102; wl line 84) | match |
| (4b) Source-tied eq:app-stage041-source-tied | `expect_zero("source-tied formula", n_src - n_src_expected)` (sympy line 121; wl line 113) | match |
| Notes regularity window | `reg_threshold` printed (no assert) | partial (informational only — paper card does not include it) |
| Notes positive-support ceiling | `num_zero_threshold` printed (no assert) | partial (informational only — paper card does not include it) |
| (extra) source-tied `dn/dm` form | `expect_zero("source-tied dn/dm", ...)` (sympy line 140; wl line 117) | extra (derived consequence of (3) at the source-tied substitution; consistent with paper's general dn/dm) |

All four paper deliverables map to non-trivial assertions in both engines. The two notes-only feasibility-window inequalities are printed but not asserted; they are strict inequalities (not equalities), so an `expect_zero` is not the natural test, and they do not appear in `\stagefield{Output}`. The source-tied `dn/dm` is extra coverage but is a manifest consequence of paper deliverable (3) under the source-tied substitution, so it does not introduce a new paper-side commitment. Setting `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `expect_zero("determinant decomposition", Det - Det_expected)` | (1) det decomposition | yes |
| A2 | sympy | 82 | `expect_zero("n_req - expected", n_req - n_expected)` | (2) `n_req` closed form | yes |
| A3 | sympy | 94 | `expect_zero("dn/dm - expected", dn_dm - monotone_expected)` | (3) monotonicity form | yes |
| A4 | sympy | 102 | `expect_zero("tracking collapse", n_track - (G_q - m))` | (4a) tracking specialization | yes |
| A5 | sympy | 121 | `expect_zero("source-tied formula", n_src - n_src_expected)` | (4b) source-tied specialization | yes |
| A6 | sympy | 140 | `expect_zero("source-tied dn/dm", dn_dm_src - dn_dm_src_expected)` | extra consequence of (3) | yes |
| B1 | wl | 56 | `expectZero["determinant decomposition", detExpr - detExpected]` | (1) | yes |
| B2 | wl | 65 | `expectZero["n_req - expected", nReq - nExpected]` | (2) | yes |
| B3 | wl | 76 | `expectZero["dn/dm - expected", dndm - monotoneExpected]` | (3) | yes |
| B4 | wl | 84 | `expectZero["tracking collapse", nTrack - (gQ - m)]` | (4a) | yes |
| B5 | wl | 113 | `expectZero["source-tied formula", nSrc - nSrcExpected]` | (4b) | yes |
| B6 | wl | 117 | `expectZero["source-tied dn/dm", dndmSrc - dndmSrcExpected]` | extra | yes |

Non-tautology spot checks:

- A1/B1: `Det` is computed from the explicit matrix `M`; `Det_expected` is the paper's RHS — independent expressions, so the residual is a genuine identity. Verified by hand expansion: `(xi-m-n)(delta+xi-mq^2-nr^2) - (mq+nr)^2 = xi(delta+xi) - m(delta+(1+q^2)xi) - n(delta+(1+r^2)xi) + mn(q-r)^2`.
- A2/B2: `n_req` from `solve(Det_expected==0, n)` vs the paper's closed form — would diverge if the solver returned the other branch or factored differently.
- A3/B3: `dn_dm = diff(n_expected, m)` vs `-(delta+(1+qr)xi)^2 / [...]^2`. I verified the underlying identity by expansion at general `m`: `-(delta+(1+q^2)xi)(delta+(1+r^2)xi) + xi(delta+xi)(q-r)^2 = -delta^2 - 2 delta xi (1+qr) - xi^2 (1+qr)^2 = -(delta + (1+qr)xi)^2`. Non-trivial.
- A5/B5: This was the tautology that a prior auditor (2026-05-22) flagged on this stage. The current scripts have applied the previously-recommended fix — `n_src` is now derived by substituting `q -> t*R_U, r -> t` into `n_expected` and then replacing `t^2 -> lambda_0`, NOT by rewriting the substituted formula in place (sympy lines 109-113; wl lines 88-97, with the comment on sympy lines 106-108 explicitly stating the intent: "Derive n_src by substituting into the general n_expected from section 24.1, not by re-stating the substituted formula"). The substitution `t^2 -> lambda_0` is safe because after `q -> t R_U, r -> t`, all `t` factors enter only as `t^2` (since `q^2 = t^2 R_U^2`, `r^2 = t^2`, and `(q-r)^2 = t^2 (R_U-1)^2` — there are no `t^1` survivors). The check `n_src - n_src_expected == 0` is therefore a genuine cross-check.

## Findings

None.

## Independent-derivation check (Mathematica)

The two scripts share substantial structural parallelism: both build the dimensionless 2x2 matrix `M`, compute `Det(M - (1-xi)I)`, expand and check against the paper's decomposition, solve for `n`, differentiate, and perform the same `r -> q` (tracking) and `{q -> t R_U, r -> t}` then `t^2 -> lambda_0` (source-tied) substitutions. Examples of corresponding sections:

- sympy 56-59 (matrix) vs. wl 44-47 (matrix) — same matrix layout.
- sympy 74 `sp.solve(sp.Eq(Det_expected, 0), n)[0]` vs. wl 58 `n /. First[Solve[detExpected == 0, n]]` — same solver step.
- sympy 109-113 (source-tied substitution) vs. wl 90-97 (source-tied substitution) — same substitution order.

However, each engine performs its own determinant computation, its own `Solve`/`solve`, and its own differentiation — the algebra is not echoed at the string/expression level. The "physical premise" is the loaded 2x2 matrix; the natural and essentially unique verification path for a 2x2 determinant identity is to compute the determinant, solve, and differentiate. I do not file a `mathematica_transliteration` finding: each engine independently re-derives the algebra and compares against the paper's closed form. A bug in either engine's simplification path would cause an independent assertion failure.

## Engine cross-check

Both engines produce equivalent final symbolic results (modulo expansion ordering and sign-of-numerator-vs-denominator):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| det decomposition residual | `= 0` (txt line 15) | `= 0`, PASS (txt line 10-11) |
| `n_req - expected` residual | `= 0` (txt line 22) | `= 0`, PASS (txt line 13-14) |
| `dn/dm - expected` residual | `= 0` (txt line 34) | `= 0`, PASS (txt line 20-21) |
| tracking collapse residual | `= 0` (txt line 45) | `= 0`, PASS (txt line 27-28) |
| source-tied formula residual | `= 0` (txt line 56) | `= 0`, PASS (txt line 34-35) |
| source-tied dn/dm residual | `= 0` (txt line 74) | `= 0`, PASS (txt line 39-40) |

Form-level agreement:

- `n_req`: SymPy `(delta*m - delta*xi + m*q^2*xi + m*xi - xi^2) / (-delta + m*q^2 - 2 m q r + m r^2 - r^2 xi - xi)`. Mathematica `(delta*(-m+xi) + xi*(-(m*(1+q^2)) + xi))/(delta - m*(q-r)^2 + xi + r^2 xi)`. Multiplying SymPy numerator and denominator by `-1` gives Mathematica's form. Equal.
- `dn/dm`: both yield `-(delta + (1+qr)xi)^2 / [delta - m(q-r)^2 + (1+r^2)xi]^2`.
- source-tied: both yield `(9 delta(-m+xi) + xi(-(m(9+2 rU^2)) + 9 xi))/(9 delta - 2 m(-1+rU)^2 + 11 xi)` up to overall-sign rearrangement of numerator and denominator.
- source-tied dn/dm: both yield `-((9 delta + (9 + 2 rU) xi)^2 / (9 delta - 2 m (rU-1)^2 + 11 xi)^2)` after expanding the denominator.

Output freshness: sympy `.txt` mtime `1779474889` > script mtime `1779474831`; Mathematica `.txt` mtime `1779474898` > script mtime `1779474832`. Both outputs are fresh (`outputs_fresh: true`).

## Verdict justification

`clean`. The paper card states four deliverables — determinant decomposition, closed-form `n_req`, monotonicity `dn_req/dm`, and the two specializations (tracking and source-tied) — and both engines exercise each with non-tautological assertions against expressions extracted from the paper. The previously-flagged tautology in the source-tied check (prior audit 2026-05-22) has been remedied: `n_src` is now derived by substitution from the general `n_expected` (sympy lines 109-113; wl lines 90-97), so A5/B5 are a genuine cross-check between the algebraic specialization and the paper's standalone closed form.

Attacks attempted and failed:

- Hand-expanded `Det(M - (1-xi)I)` to verify the decomposition — matches.
- Hand-expanded the monotonicity-quotient numerator to verify `-(delta+(1+q^2)xi)D + N(q-r)^2 = -(delta+(1+qr)xi)^2` — matches.
- Checked whether the `t^2 -> lambda_0` substitution misses residual `t^1` survivors — confirmed `t` enters `n_expected` after `{q -> t R_U, r -> t}` only through `q^2`, `r^2`, and `(q-r)^2`, all of which collect into pure `t^2` factors.
- Checked symbol-assumption integrity: `delta, xi, R_U > 0` and `m, n, q, r, t` real. None of the assumptions hide a branch ambiguity; `simplify`/`FullSimplify` operates on rational expressions with real-positive denominators only on the regular branch (which is exactly the paper's stated scope).
- Checked monotonicity sign anchoring: the script verifies the closed-form `-X^2/Y^2`, not `< 0` directly. Strict `< 0` follows from the closed form precisely on the "regular branch" the paper qualifies (where `delta + (1+qr)xi != 0` and the denominator is nonzero). The form-anchor is therefore sufficient.
- Checked whether the source-tied `dn/dm` (an extra assertion) introduces a new paper claim — no, it is a manifest consequence of paper deliverable (3) at `q = t R_U, r = t, t^2 = lambda_0`.
- Verified by hand at `R_U = 1`: source-tied formula collapses to `xi(delta+xi)/(delta + 11 xi/9) - m`, consistent with the tracking limit at `q = r = t`, `t^2 = lambda_0`. The script's `n_src = n_expected.subs(...)` would reproduce this.

I confirm I read the paper card, the notes, and the appendix row before the scripts. The script's load-bearing claims match the paper's `\stagefield{Output}`.

## Self-test notes

I checked: (1) variable independence — `n_expected` depends explicitly on `m, q, r, xi, delta` and `n_src_expected` depends on `m, xi, delta, R_U`; the differentiations `sp.diff(n_expected, m)` and `D[nSrcExpected, m]` are non-degenerate. (2) Parity/symmetry — n/a; no unbounded integrals in this stage. (3) Trivial-case substitution — at `R_U = 1`, source-tied collapses correctly to the tracking limit at `q = r = t, t^2 = lambda_0` (verified above); at `m = n = q = r = 0`, `n_req` formula gives `xi`, consistent. (4) Path specs — n/a (no missing-script findings). (5) Paper round-trip — the script's closed forms all appear verbatim in eq:app-stage041-Dsel, eq:app-stage041-nreq, eq:app-stage041-nreq-derivative, eq:app-stage041-tracking, and eq:app-stage041-source-tied; the source-tied dn/dm is an extra consequence of the paper's general dn/dm at the source-tied substitution and adds no new paper-side commitment.
