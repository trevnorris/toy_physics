---
unit_id: 051
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage051_lowest_twin_criterion.md"]
  paper_appendix: present
---

# Audit unit 051 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_051.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage051_lowest_twin_criterion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at L80; `\input` at L220; statement-list line at L348)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.txt`

## What the paper claims

Stage 051 rewrites the lowest symmetric-twin sufficiency question in continuum branch-product variables. The paper card's `\stagefield{Output}` is the boxed inequality:

> `Pi_tr <= 2 C_mix = 16 Lambda (1-eps) / pi^2`, with `Pi_tr := F_tr G_tr` on the tracking branch and `C_mix = 8 Lambda (1-eps)/pi^2`.

The notes file decomposes the criterion into six deliverables: (a) the exact closed form for `Pi_tr`; (b) the endpoint facts `Pi_tr(xi=0)=0`, `Pi_tr(xi->1-)=+oo`; (c) elimination of `zeta_req` via `S_req = Pi_tr / C_mix` so that `zeta_req<=1 <=> Pi_tr<=2 C_mix`; (d) the equivalent radiative threshold `Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]`; (e) the mixed-baseline threshold `M_mix^(twin,req) = G_tr/2` and the equivalent wall-overlap threshold `Z_W^(twin,req) = pi^2 (1-eps_eta)(1-eps) G_tr / [16 (1+chi0)^2]` (the latter routed through the Stage 030/047 coherent forward map); (f) the closed quadratic root `xi_(2x) = [2 M_mix(9+2R^2) - 9 delta + sqrt((2 M_mix(9+2R^2)-9 delta)^2 + 648 M_mix delta)]/18` for `G_tr = 2 M_mix`.

`\stagefield{Inputs}` names `Pi_tr = F_tr G_tr` and `C_mix = 8 Lambda(1-eps)/pi^2`; the appendix row (L80) and the part-summary line (L348) both restate the boxed criterion.

## What the script claims to verify

The scripts (i) independently construct `G_tr` and `F_tr` from their Stage 045 closed forms and assert that their product collapses to the notes' closed form for `Pi_tr` (cancellation against one copy of `(9 delta + (9+2R^2) xi)` between F_tr's squared numerator factor and G_tr's denominator); (ii) check both endpoint limits; (iii) derive `zeta_req(Pi, C_mix, eps)` and verify the two anchor identities `zeta_req|_(Pi=C_mix) = 0` and `zeta_req|_(Pi=2 C_mix) = 1` (the lowest-twin threshold); (iv) write down `Lambda_twin_req`, `M_mix^(twin,req)`, and `Z_W^(twin,req)`, and cross-check the latter against the Stage 030/047 forward map applied to `M_mix = G_tr/2`; (v) substitute the closed-form `xi_(2x)` into `G_tr - 2 M_mix` and assert vanishing (Mathematica additionally re-derives `xi_(2x)` via `Solve[{gTr==2 mMix, xi>0}, xi, Reals]` and compares to the closed-form claim).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Closed form of `Pi_tr = F_tr G_tr` | `expect_zero("Pi_tr - expected closed form")` in both engines | match |
| `Pi_tr(xi->0+) = 0` | `Pi0` limit in both engines | match |
| `Pi_tr(xi->1-) = +oo` | direct `Limit` (SymPy); `1/pi1 == 0` (Mathematica) | match |
| Boxed criterion `Pi_tr <= 2 C_mix` (Output) | `zeta_req|_(Pi=C_mix)=0`, `zeta_req|_(Pi=2 C_mix)=1` in both engines | match |
| `Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]` | formula constructed and printed; identity is the algebraic inverse of the boxed criterion already verified through `zeta_req(2 C_mix)=1` | match |
| `M_mix^(twin,req) = G_tr/2` | formula constructed; consistency exercised through the Z_W cross-check (forward map of `M_mix=G_tr/2` reproduces the paper's `Z_W^(twin,req)`) | match |
| `Z_W^(twin,req)` closed form | direct cross-check `ZW_twin_req - forward-map(M_mix=G_tr/2) == 0` in both engines | match |
| `xi_(2x)` closed quadratic root | `G_tr(xi_(2x)) - 2 M_mix == 0` in SymPy; Mathematica additionally derives via `Solve[..., Reals]` and compares to the closed-form claim | match |

`paper_alignment: aligned`. Every notes deliverable has a script-side check; no script-side check verifies an identity not also asserted in the paper card or notes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75 | `expect_zero(Pi_tr - Pi_expected)` | closed form of Pi_tr | yes |
| A2 | sympy | 80-81 | `Pi0 != 0 -> raise` | endpoint xi->0 | yes |
| A3 | sympy | 84-85 | `Pi1 is not sp.oo -> raise` | endpoint xi->1- | yes |
| A4 | sympy | 99 | `expect_zero(zeta_req @ Pi=C_mix)` | zeta_req baseline (S_req=1) | yes |
| A5 | sympy | 100 | `expect_zero(zeta_req @ Pi=2 C_mix - 1)` | boxed criterion threshold (zeta_req=1) | yes |
| A6 | sympy | 104-107 | `expect_zero(factor((zeta_req-1) @ Pi=2 C_mix))` | redundant with A5, alternative simplification path | partial (redundant) |
| A7 | sympy | 130-133 | `expect_zero(ZW_twin_req - forward-map(M_mix=G_tr/2))` | Z_W^(twin,req) vs Stage 030/047 forward map | yes |
| A8 | sympy | 147 | `expect_zero(G_tr(xi=xi_2x) - 2 M_mix)` | xi_(2x) closed root | yes |
| M1 | mathematica | 58 | `expectZero(piTrFactored - piExpectedFactored)` | closed form of Pi_tr | yes |
| M2 | mathematica | 64 | `pi0 =!= 0 -> fail` | endpoint xi->0 | yes |
| M3 | mathematica | 68-69 | `Simplify[1/pi1 == 0]` | endpoint xi->1- (pole-robust idiom) | yes |
| M4 | mathematica | 80 | `expectZero(zeta_req @ Pi=C_mix)` | zeta_req baseline | yes |
| M5 | mathematica | 81 | `expectZero(zeta_req @ Pi=2 C_mix - 1)` | boxed criterion threshold | yes |
| M6 | mathematica | 82 | `expectZero((zeta_req - 1) @ Pi=2 C_mix)` | redundant with M5 | partial (redundant) |
| M7 | mathematica | 94 | `expectZero(zWTwinReq - zWThresholdViaMap)` | Z_W^(twin,req) vs forward map | yes |
| M8 | mathematica | 110 | `expectZero(xi2xDerived - xi2xClaim)` | xi_(2x) closed form derived independently via Solve matches claim | yes |
| M9 | mathematica | 111 | `expectZero(G_tr(xi=xi2xDerived) - 2 M_mix)` | xi_(2x) is the positive root | yes |

A6/M6 are algebraically the same identity as A5/M5 but exercise an alternative simplification path (`sp.factor` vs raw subtraction; `(zeta_req-1)@subs` vs `zeta_req@subs - 1`). They are not tautological — they still substitute `Pi = 2 C_mix` into the rational `zeta_req` formula and assert vanishing — but they are redundant. Not promoted to a finding.

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script shares the same physical symbols and the same defining expressions for `G_tr`/`F_tr` as SymPy (these are inputs to the stage, not outputs to be derived independently). Beyond that, three independent moves keep it clear of `mathematica_transliteration`:

1. Canonicalization route differs. SymPy uses `sp.simplify(sp.factor(Ftr*Gtr))` and compares against `sp.simplify(...)`. Mathematica uses `Factor[Together[fTr gTr]]` and compares against `Factor[Together[piExpected]]` (wl L48-51, L58). The two reach the same canonical form by different normalizers, so a coincidental defect would have to survive both `simplify` and `Factor∘Together`.
2. The xi_(2x) derivation is genuinely independent in Mathematica. SymPy only substitutes the notes' closed-form root into `G_tr - 2 M_mix`. Mathematica invokes `Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]`, strips the `ConditionalExpression` wrapper per the project idiom, `FullSimplify`s, and asserts `xi2xDerived - xi2xClaim == 0` (wl L100-110) before separately verifying `G_tr(xi2xDerived) = 2 mMix`. This is a real independent derivation of the most nontrivial closed-form result in the stage.
3. The xi->1- pole check uses a different probe: SymPy directly tests `Pi1 is sp.oo`; Mathematica uses `1/pi1 == 0` to evade `Limit`'s non-deterministic pole behavior (the project-wide Mathematica idiom).

Not a transliteration finding.

## Engine cross-check

Both engines exit 0. Canonical forms differ only by trivial sign-pairing of `(xi-1)` vs `(1-xi)` factors (SymPy denominator `9*(xi-1)*(...)^2` with leading `-`; Mathematica `-1/9 * (xi*(delta+xi)*...)/((-1+xi)*(...)^2)`). Both reduce to the same rational function. Every cross-engine equivalent assertion produces a zero residual. `xi_(2x)`: SymPy prints `Mmix*(2 R^2+9)/9 - delta/2 + sqrt(...)/18`; Mathematica prints `(-9 delta + 18 mMix + 4 mMix r^2 + Sqrt[...])/18`. Distributing the SymPy form yields the Mathematica form term-by-term: `Mmix(2R^2+9)/9 = (18 mMix + 4 mMix R^2)/18`; `-delta/2 = -9 delta/18`; sum matches. Engines agree.

## Verdict justification

I attacked the script on several fronts:

- **Is `Z_W^(twin,req) - forward-map(M_mix=G_tr/2)` tautological?** No. `ZW_twin_req` is the paper's claimed closed form (with prefactor `pi^2 (1-eps_eta)(1-eps)/[16(1+chi0)^2]`), while `ZW_from_Mmix` is the upstream Stage 030/047 forward map (prefactor `pi^2 (1-eps_eta)(1-eps)/[8(1+chi0)^2]`). The factor-of-2 difference exists because `Z_W^(twin,req)` uses `M_mix = G_tr/2`, and the assertion verifies these two independently-stated formulas reconcile correctly. A perturbation of `(1+chi0)^2` to `(1+chi0)` in either formula would not cancel.
- **Is `Pi_tr - expected closed form` tautological?** No. `F_tr` and `G_tr` are independent Stage 045 closed forms; their product is a degree-7 numerator over a degree-5 denominator before cancellation, and the assertion checks that one factor of `(9 delta + (9+2R^2) xi)` cancels cleanly between F_tr's numerator-square and G_tr's denominator to yield the notes' specific closed form. Substantive.
- **Domain over-strengthening?** SymPy declares all dimensionful parameters `positive=True`; Mathematica adds `0<xi<1`, `0<eps<1`, `chi0 > -1`, `epsEta in (0,1)`. These are consistent with the physical setup in the notes and the paper card; none of them are strong enough to mask sign errors in the rational identities checked (which hold identically on any open subset).
- **xi_(2x) sign / missing branch?** Discriminant `(2 M_mix(9+2R^2) - 9 delta)^2 + 648 M_mix delta` is strictly larger than `(2 M_mix(9+2R^2) - 9 delta)^2` for `M_mix, delta > 0`, so only the `+` root is positive. The script's `+` choice is physically forced. Mathematica's `Solve` with explicit `xi > 0` constraint picks the same root.
- **xi->1- direction sign?** SymPy returns `+oo`; Mathematica's `1/pi1 == 0` test passes. Hand-check: numerator at xi=1 is `-(delta+1)(9 delta+9+2R)^2(9 delta+9+2R^2) < 0`; denominator `9*(xi-1)*(...)^2` is negative-small as xi->1-; quotient → positive infinity.
- **Constants match?** Every paper coefficient (8, 16, 18, 81, 648, the powers of `(9+2R)`, `(9+2R^2)`, `(1+chi0)`) appears identically in both scripts and matches the notes and paper card. No `paper_misalignment` (no `value_mismatch`, no `target_mismatch`).
- **Output freshness?** SymPy script mtime 2026-05-22T17:11, output 2026-05-22T17:18 (fresh). Mathematica script mtime 2026-05-22T17:30, output 2026-05-22T17:31 (fresh).
- **Checkpoint bar?** Stage 051 is a checkpoint. Both engines are present. Every paper-side deliverable has a matching script-side substantive assertion. No status-only carve-outs. Paper alignment is exact on constants, formula structure, and the boxed criterion direction.

Soft observations (not findings):

- The SymPy docstring (L4-6) and both banners still say "Stage 34" / "STAGE 34" — a stale-rename artifact from the global Part-III renumbering. Filename, path, and content are correct; only the prose labels lag. Cosmetic, no impact on math or paper alignment, and does not match any of the 10 finding categories.
- A6/M6 are algebraically redundant with A5/M5. Not tautological — they still substitute and simplify — just duplicative.

Nothing in the audit categories applies. The scripts substantively verify the paper's load-bearing criterion and every notes-side deliverable; both engines agree at the level claimed; outputs are fresh; the paper card, notes, and appendix row are mutually consistent. Verdict: clean.

## Self-test notes

Checked: (1) symbol domains do not over-strengthen the simplifier (all assertions are rational-function equalities holding identically); (2) A6/M6 exercise the same nontrivial substitution as A5/M5 with a different factor route, not a self-comparison; (3) `Pi_tr - expected` collapse is non-trivial because the inputs `F_tr` and `G_tr` are independent rational expressions whose product requires factor cancellation; (4) Mathematica's `Solve[..., Reals]` for `xi_(2x)` provides genuine independent derivation; (5) every numeric coefficient (8, 16, 18, 81, 648, `9+2R`, `9+2R^2`, `(1+chi0)^2`) matches across paper card, notes, and both scripts. No variable-independence, parity, or path-spec traps apply (no derivatives over indep-vars and no integrals over unbounded domains in either script).
