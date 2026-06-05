---
unit_id: 081
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage081_family1_pi_thresholds.md]
  paper_appendix: present
---

# Audit unit 081 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_081.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 140 references this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 081 re-expresses the explicit Family-1 demand thresholds in the branch-product variable `Pi_tr`. The card's `\stagefield{Output}` is verbatim: *"Product-variable conversion \eqref{eq:app-stage081-Pi-zeta}--\eqref{eq:app-stage081-unblocked}."* The two boxed equations are the support-demand law `zeta_req = (Pi_tr - C_mix)/(C_mix - eps_blk(2 C_mix - Pi_tr))` (carried from Stage 052) and its unblocked limit `Pi_tr/C_mix = 1 + zeta_req` at `eps_blk = 0`. The notes elaborate the load-bearing deliverable: inverting the support law exactly to `Pi_tr = C_mix * Q(zeta;eps_blk)` with the closed form `Q = [1 + (1-2 eps_blk) zeta]/[1 - eps_blk zeta]`, its anchor values `Q(0)=1`, `Q(1)=2`, positive derivative, the five numeric `Pi/C` thresholds at `lambda_mu = 1, eps_blk = 0` (`3.46622291347846`, `3.46752913273870`, `3.44257571477179`, `3.46752736855058`, ceiling `3.46752922945601`), and the blocking ceiling `eps_blk < 1/zeta_max^(F1) ≈ 0.405263689711371`. The five threshold values themselves are carry-forwards from Stages 079/080 (the card's `Inputs` field cites "the support-ratio threshold formula"), not derived in this stage.

## What the script claims to verify

Both scripts construct `zeta_expr` = the support-demand law, `solve` it for `Pi`, form `Q = Pi/C_mix`, and check the inversion. The SymPy script asserts (via `expect_zero`, which raises on nonzero) only the two anchor identities `Q(0)-1 == 0` and `Q(1)-2 == 0`; everything else (the `Q` closed form, `dQ/dzeta`, the five `Pi/C` threshold values, the eps-ceiling `0.40526...`) it merely *prints*. The Mathematica script asserts the full closed-form identity `Q - (1 + zeta - 2 eps_blk zeta)/(1 - eps_blk zeta) == 0` (line 54), the two anchors, five `expectApprox` checks that each `Pi/C|eps=0` equals `1 + zeta`, and that `(1/zeta_max)*zeta_max - 1 == 0`. The scripts hardcode the five upstream zeta thresholds as `sp.Float`/`ToExpression` literals — a legitimate carry-forward, consistent with the card's `Inputs` field.

## Paper ↔ script cross-check

| Paper / notes deliverable | Script-side check | Status |
|---|---|---|
| Support-demand law `zeta_req = (Pi_tr-C_mix)/(C_mix-eps_blk(2C_mix-Pi_tr))` (eq Pi-zeta) | `zeta_expr` definition (py L33 / wl L45) used as the inversion source | match |
| Exact inversion `Q = [1+(1-2eps)zeta]/[1-eps zeta]` (notes main result) | wl L54 asserts full identity; **py only pins Q at zeta=0,1** (L40-41) | partial (SymPy) / match (Mathematica) |
| Anchors `Q(0)=1`, `Q(1)=2` | py L40-41, wl L59-60 (asserted) | match |
| Unblocked limit `Pi_tr/C_mix = 1+zeta` (eq unblocked) | wl L80-84 `expectApprox` (asserted); py prints only | match (Mathematica) |
| Five `Pi/C` threshold values at `lambda_mu=1,eps=0` | printed py L65-69 / wl L74-78 (carry-forward, reconcile to notes) | match |
| Blocking ceiling `eps_blk < 0.405263689711371` | py L72-73 print; wl L86-88 asserted (`eps*zeta_max-1==0`) | match |
| `dQ/dzeta = (1-eps)/(1-eps zeta)^2 > 0` | py L42 print only; wl not checked | extra/info (sign not asserted, but card doesn't require it) |

Set `paper_alignment: aligned` — every paper-side deliverable maps to a faithful script-side check; the SymPy gap is a coverage weakness, not a misalignment. No `value_mismatch`, no `target_mismatch`, no orphaned script claims.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 54 | `expectZero[Q - closedform]` | inversion (main deliverable) | yes |
| A2 | sympy | 40 | `expect_zero(Q(0)-1)` | anchor `Q(0)=1` | yes (point only) |
| A3 | sympy | 41 | `expect_zero(Q(1)-2)` | anchor `Q(1)=2` | yes (point only) |
| A4 | mathematica | 59 | `expectZero[(Q/.zeta->0)-1]` | anchor `Q(0)=1` | yes |
| A5 | mathematica | 60 | `expectZero[(Q/.zeta->1)-2]` | anchor `Q(1)=2` | yes |
| A6 | mathematica | 80-84 | `expectApprox[Pi/C|eps=0 - (1+zeta)]` ×5 | unblocked limit | partial (structural, see F-note) |
| A7 | mathematica | 88 | `expectApprox[eps_ceil*zeta_max - 1]` | blocking ceiling reciprocal | yes |
| — | sympy | 51-73 | `print` only (Pi/C values, ceiling) | thresholds / ceiling | no (not asserted) |
| — | sympy | 42 | `print(dQ/dzeta)` | monotonicity | no (not asserted) |

The "Anchored?" no/partial rows on the SymPy side feed F1. A6 is structurally weak: once A1 fixes `Q` to the closed form, `Q|eps=0 = 1+zeta` holds for *any* zeta by construction, so A6 confirms the closed form's eps=0 limit rather than independently validating each threshold value — redundant but not wrong, and it does exercise the paper's unblocked equation, so it is not raised as a separate finding.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:34-42`

**What's wrong:**
The paper's load-bearing deliverable for this stage is the *exact inversion* of the support law to the closed form `Q = [1 + (1-2 eps_blk) zeta]/[1 - eps_blk zeta]` (notes §1 lines 49-50; card eq `Pi-zeta`/`unblocked`). The Mathematica engine asserts this as a full symbolic identity:

```
wl:54  expectZero["Q matches closed form",
         qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)];
```

The SymPy engine, however, only *prints* the derived `Q` (line 39) and then asserts it at **two points**:

```
py:40  expect_zero("Q(0)-1", sp.simplify(Q.subs(zeta, 0) - 1))
py:41  expect_zero("Q(1)-2", sp.simplify(Q.subs(zeta, 1) - 2))
```

`Q` is a one-parameter family of rational functions in `zeta` (numerator and denominator each linear in `zeta` with an `eps_blk` coefficient). Two interpolation points (`zeta=0,1`) do not pin such a rational form. The only thing forcing the *full* SymPy `Q` to be correct is the `sp.solve` step itself — there is no assertion that the solved form equals the paper's closed form across `zeta` and `eps_blk`. The SymPy script therefore does not load-bearingly verify the stage's primary deliverable; it relies on `solve` being right and on a human reading the printed form.

**Why this matters:**
This is a checkpoint-grade stage (checkpoint: False, but the higher-bar dual-engine policy applies): both engines must substantively verify the central claim. Right now the central inversion identity is asserted in only one engine. If `sp.solve` ever returned a different-but-anchor-consistent branch (or a future SymPy version changed normalization), the SymPy checks A2/A3 would still pass while silently failing to validate the conversion formula the paper boxes.

**Required change:**
Add one SymPy assertion mirroring `wl:54`, immediately after `Q` is formed (after line 39):
```python
expect_zero("Q-closedform",
    Q - (1 + zeta - 2*eps_blk*zeta)/(1 - eps_blk*zeta))
```
This is the exact notes/paper closed form (notes lines 49-50). Self-tested: `Q` from output line 6 is `(2 eps z - z - 1)/(eps z - 1)`, which equals `(1 + z - 2 eps z)/(1 - eps z)` (both numerator and denominator negated), so the residual simplifies to `0` identically — `expect_zero` passes, and it is non-tautological because the LHS comes from `solve` while the RHS is a literal closed form.

**Verification:**
After the fix, the SymPy output must contain a new line `Q-closedform = 0` and the script must still exit 0. The verifier re-runs `redteam exec-sympy 081` and confirms the new assertion appears and passes.

## Independent-derivation check (Mathematica)

The two scripts share the same algorithmic skeleton: define `zeta_expr`, `Solve`/`solve` for `Pi`, form `Q = Pi/C_mix`, check anchors, substitute the same five hardcoded zeta values, compute `1/zeta_max`. On its face this looks transliteration-adjacent. However, the operation being verified is a *unique linear-fractional inversion* — there is essentially one correct route (solve a Möbius-type equation), so "independent derivation" has little room to diverge, and shared choreography here is not evidence of echoing. The Mathematica script is not a line-by-line port: it works under explicit `Reals && zeta>=0 && piTr>0 && cMix>0 && 0<=epsBlk<1` assumptions with `ConditionalExpression` stripping (wl L41-53), adds the full closed-form identity assertion (L54) that SymPy lacks, and uses `expectApprox` numeric machinery absent from the SymPy side. I do **not** raise `mathematica_transliteration`: the Mathematica engine asserts strictly *more* than SymPy and via different mechanics, so it is functioning as an independent second engine, not an echo.

## Engine cross-check

Both engines agree at the level they claim:

- `Q` form: SymPy `(2*eps_blk*zeta - zeta - 1)/(eps_blk*zeta - 1)` (out L6) is algebraically identical to Mathematica `(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)` (out L8) — numerator and denominator both negated.
- Anchors `Q(0)-1` and `Q(1)-2`: both `= 0` in both engines (sympy out L7-8; math out L9-12 PASS).
- `Pi/C|eps=0` values: SymPy prints `3.46622291…, 3.46752913…, 3.44257571…, 3.46752736…, 3.46752922…` (out L16-20); Mathematica's symbolic `2 + r/(d - epsBlk)` forms evaluate to the same at `eps=0` and all five `expectApprox` checks PASS (out L18-27).
- Blocking ceiling: SymPy `0.40526368971137149977`, Mathematica `0.40526368971137148173…` — agree to ~16 significant digits.

No `engine_disagreement`.

## Verdict justification

The math is correct and the paper alignment is exact: the support-demand law, its exact inversion to `Q`, the anchor values, the unblocked limit, the five carried-forward `Pi/C` thresholds, and the blocking ceiling all reconcile between the scripts, the saved outputs, the card, and the notes. Attacks that failed: (1) tautology hunt — the inversion comes from `solve`, not a posited form, and `Q(0)/Q(1)` checks are genuine; (2) symbol-domain attack — Mathematica's `0<=epsBlk<1`, `zeta>=0`, `cMix>0` assumptions are justified by the physical setup and don't enable an invalid simplification; (3) value-reconciliation — all five thresholds and the ceiling match the notes to full printed precision; (4) sign/factor attack on `Q` — the two engines' forms are provably identical. The single finding (F1) is a real but contained coverage gap: SymPy verifies the stage's *primary* deliverable (the closed-form inversion) only at two interpolation points rather than as the full identity that Mathematica already asserts. Hence `verdict: findings`, not `clean`; not `stop_cold` (the fix is a one-line additive assertion that strengthens, not changes, the result). Separately, two stale `+17` pre-renumber cross-reference labels in SymPy comments are noted below (comment-only, content-resolvable, in-loop label hygiene); they affect no assertion and are folded into the directive as low-risk comment edits, not as a categorized finding.

## Value Reconciliation (pass-2 augmentation)

`reconciliation: complete; 11 values checked, 0 misaligned`

The card (`.tex`) is terse and intentionally carries only the two boxed conversion equations; the numeric thresholds live in the notes (`.md`), which is the correct carrier — so card-side omissions of intermediate numbers are MATCH-via-notes, not MISSING.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Q = [1+(1-2eps)zeta]/[1-eps zeta]` | py out L6 / wl out L8 | notes L49-50; card eq `Pi-zeta`/`unblocked` (L16,21) | MATCH |
| `Q(0)=1` | py out L7 / wl out L10 | notes L56 | MATCH |
| `Q(1)=2` | py out L8 / wl out L12 | notes L58 | MATCH |
| `dQ/dzeta = (1-eps)/(1-eps zeta)^2` | py out L9 | notes L62 | MATCH |
| `Pi_suff^(chi)/C_mix = 3.46622291347846` | py out L16 | notes L29, L133 | MATCH |
| `Pi_fail^(chi)/C_mix = 3.46752913273870` | py out L17 | notes L136 | MATCH |
| `Pi_suff^(J)/C_mix = 3.44257571477179` | py out L18 | notes L145 | MATCH |
| `Pi_fail^(J)/C_mix = 3.46752736855058` | py out L19 | notes L147 | MATCH |
| `Pi_max^(F1)/C_mix = 3.46752922945601` | py out L20 | notes L29, L141; (zeta_max L104) | MATCH |
| blocking ceiling `eps_blk < 0.405263689711371` | py out L21 / wl out L28 | notes L118 | MATCH |
| `zeta_max^(F1) = 2.46752922945601` (input) | py L49 / wl L66 | notes L104 | MATCH |

Internal scaffolding (no finding): the four other hardcoded `zeta_*` input literals (suff/fail chi and J at `lambda_mu=1`) are upstream carry-forwards that reconcile via their `Pi/C = 1+zeta` images above; pass/fail flags; `expectApprox` diffs (`0`); tolerances `10^-14`; the `2 + r/(d-eps)` partial-fraction display forms in the Mathematica output (algebraically equal to the SymPy rational forms).

Every emitted deliverable value reconciles to the notes (and, for the two boxed equations, the card). No `value_mismatch` or `script_missing_paper_claim` arises from the reconciliation; F1 stands as the sole finding.

## Self-test notes

Variable-independence trap: F1 introduces no `diff`/`D` and no integral — it is a pure rational subtraction `Q - closedform`, so the zero-derivative failure mode does not apply. Parity/symmetry trap: no unbounded-domain integrals present. Trivial-case pre-check: substituted the derived `Q = (2eps z - z -1)/(eps z -1)` into `Q - (1+z-2eps z)/(1-eps z)` and confirmed it reduces to `0` identically (both fraction parts are sign-negated copies), so the proposed `expect_zero` passes non-tautologically. Paper round-trip: the proposed RHS is the verbatim notes closed form (L49-50) and introduces no new constant, so the fix cannot create a fresh `paper_misalignment`.
