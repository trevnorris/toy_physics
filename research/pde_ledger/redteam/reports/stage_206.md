---
unit_id: 206
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md]
  paper_appendix: present
---

# Audit unit 206 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 43, 210, 737-781, 824-831, 918; equation blocks eq:app-part06-oriented-H through eq:app-part06-turning-bracket)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card's `\stagefield{Output}` reads verbatim: "Certified monotone bracket theorem, certified turning-ray bracket theorem, pairwise ray-ordering theorem, and local search-sieve admissibility test." Stage 206 converts the carried Stage-239 local quadratic predictors into *certified* local root brackets by bounding the oriented log residual's curvature with an exact envelope `\(\underline K_1 \le H_{\mathbf s}'' \le \overline K_1\)` on `\([0,T]\)`. The notes enumerate eight deliverables, of which the symbolically-verifiable core is: (3) the monotone-branch forward root map `\(\mathcal T(H_0,K_0;c) = -2H_0/(K_0 + \operatorname{sgn}(K_0)\sqrt{K_0^2-2cH_0})\)`, with zero-curvature limit `\(-H_0/K_0\)` (notes §3.1) and strict positive monotonicity `\(\partial\mathcal T/\partial c = \mathcal T^2/(2\sqrt{K_0^2-2cH_0})>0\)` (§3.2); (4) the certified monotone bracket with the descent-sign fact `\(K_0+\overline K_1\tau_{\rm hi}=-\sqrt{\Delta_+}<0\)` (§4.1); (5) the turning-ray bracket `\(\tau^{\rm(tp)}=\sqrt{-2H_0/c}\)` (§5); (6) the bracket-width law and its small-envelope leading term `\(W \approx \mathcal T(\bar K_1)^2/(2\sqrt{K_0^2-2\bar K_1 H_0})\cdot\Delta K_1\)` (§6.1) and its zero-curvature simplification `\(W\approx \tau_0^2/(2|K_0|)(\overline K_1-\underline K_1)\)` (§6.2); (7) the pairwise ray-ordering theorem `\(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\Rightarrow\tau_*^{(a)}<\tau_*^{(b)}\)` (§7); and (8) the local search-sieve admissibility test (§8). The card is `is_checkpoint: False` but `\claimstatus` is "Mixed: ExactClosure, Numerical," so the predictor *formulas* must be exact.

## What the script claims to verify

The SymPy script (docstring-less; banner says "STAGE 189 — CERTIFIED RAY RANKING AND LOCAL BRACKETING") verifies, on the oriented negative-slope branch written as `\(K_0=-k,\ k>0\)`: (I) the root-map formula solves its defining quadratic, has the correct zero-curvature limit `\(H_0/k\)`, satisfies the implicit `\(\partial\mathcal T/\partial c\)` identity, and the descent-sign identity `\((-k+c\mathcal T)+\sqrt{k^2-2cH_0}=0\)`; (II) the lower/upper bracket endpoints `Tau_lo`, `Tau_hi` each solve their respective quadratics, plus a degenerate-envelope collapse; (III) the small-envelope width series matches the closed-form leading term, including the zero-mean-curvature simplification `\(H_0^2\eta/(2k^3)\)`; (IV) the turning-ray roots `\(\sqrt{2H_0/a}\)`, `\(\sqrt{2H_0/b}\)` solve `\(H_0-\tfrac12 a\tau^2=0\)` and a turning derivative identity; and (V) a so-called "Stage 205/188 collapse" comparing two expressions `TauStage189` and `TauStage188`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Verdict |
|---|---|---|
| (3) monotone root map + zero-curvature limit | Section I, lines 36-41 | match |
| (3) monotonicity `\(\partial\mathcal T/\partial c>0\)` | Section I, lines 42-45 | match (verifies the formula; positivity is manifest from the formula) |
| (4) certified bracket endpoints + descent sign | Sections I (46-49) and II (56-69) | match |
| (5) turning-ray bracket roots | Section IV, lines 111-128 | match |
| (6) small-envelope width law + zero-curvature law | Section III, lines 84-104 | match |
| (3.3) degenerate-envelope collapse to single predictor | Section II, lines 70-77 | partial (collapse is true by construction — `TauL` is `Tau` with `cL`→`c`; weak but defensible) |
| "Stage 205/188 collapse" (§3.3 / Stage-239 recovery flavor) | Section V, lines 134-148 | mismatch — tautological; both sides are the identical expression (see F1) |
| (7) pairwise ray-ordering theorem | (none) | missing (see F3) |
| (8) local search-sieve admissibility test | (none) | missing (folded into F3 narrative) |

Dominant pattern: most algebraic predictor deliverables are faithfully exercised, but two of the four enumerated `\stagefield{Output}` theorems (pairwise ordering, sieve admissibility) have no script-side check, and Section V's "collapse" check is vacuous. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero(H0 - k*Tau + ½c*Tau²)` | claim 3 (root map) | yes |
| A2 | sympy | 41 | `expect_zero(limit(Tau,c,0) - H0/k)` | claim 3 (zero-curv limit) | yes |
| A3 | sympy | 42-45 | `expect_zero(diff(Tau,c) - Tau²/(2√…))` | claim 3 (monotonicity) | yes |
| A4 | sympy | 46-49 | `expect_zero((-k+c*Tau)+√(k²-2cH0))` | claim 4 (descent sign) | yes |
| A5 | sympy | 62-65 | `expect_zero(H0 - k*TauL + ½cL*TauL²)` | claim 4 (lower endpoint) | yes |
| A6 | sympy | 66-69 | `expect_zero(H0 - k*TauU + ½cU*TauU²)` | claim 4 (upper endpoint) | yes |
| A7 | sympy | 70-73 | `expect_zero(TauL.subs(cL,c) - Tau)` | claim 3.3 (collapse) | partial (true by symbol rename) |
| A8 | sympy | 74-77 | `expect_zero(TauU.subs(cU,c) - Tau)` | claim 3.3 (collapse) | partial (true by symbol rename) |
| A9 | sympy | 96 | `expect_zero(series_width - leading_width)` | claim 6 (width law) | yes |
| A10 | sympy | 104 | `expect_zero(series_width0 - leading_width0)` | claim 6.2 (zero-curv width) | yes |
| A11 | sympy | 117-120 | `expect_zero(H0 - ½a*TauTurnA²)` | claim 5 (turning root a) | yes |
| A12 | sympy | 121-124 | `expect_zero(H0 - ½b*TauTurnB²)` | claim 5 (turning root b) | yes |
| A13 | sympy | 125-128 | `expect_zero(diff(TauTurnA,a)+TauTurnA/(2a))` | claim 5 (turning deriv) | yes |
| A14 | sympy | 148 | `expect_zero(TauStage189 - TauStage188)` | claim 3.3 (collapse) | **no — tautological (F1)** |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py:134-148`

**What's wrong:**
Section V purports to verify a "collapse to the Stage 205 quadratic predictor," but the two operands are the *same expression*:
```python
H0_stage205 = sp.log(Phi0)
TauStage189 = sp.simplify(-2 * H0_stage205 / (Lneg + sp.sign(Lneg) * sp.sqrt(Lneg**2 - 2 * L1 * H0_stage205)))
TauStage188 = sp.simplify(-2 * sp.log(Phi0) / (Lneg + sp.sign(Lneg) * sp.sqrt(Lneg**2 - 2 * L1 * sp.log(Phi0))))
expect_zero("Stage 206/188 collapse", sp.simplify(TauStage189 - TauStage188))
```
`H0_stage205` is bound to `sp.log(Phi0)` on line 137, so `TauStage189` and `TauStage188` are character-for-character identical after substitution. Their difference is identically `0` before `simplify` is ever called. The assertion cannot fail for any value of `Phi0`, `Labs`, or `L1`, regardless of whether any "collapse" relation actually holds. The comment on line 131 claims this checks "Collapse to the Stage 205 quadratic predictor under exact curvature," but no independent Stage-205/Stage-238 predictor expression is constructed to compare against — both sides are the *same* curvature-bracket root map.

The paper's actual collapse claim (notes §3.3, lines 196-210) is non-trivial: when the envelope degenerates `\(\underline K_1=\overline K_1=L_1\)`, the certified bracket map `\(\mathcal T(H_0,K_0;L_1)\)` must equal the Stage-239 quadratic log predictor `\(\tau_{\log,2}=-2\ln\Phi_0/(L_0+\operatorname{sgn}(L_0)\sqrt{L_0^2-2L_1\ln\Phi_0})\)` (appendix eq:app-part06-quadratic-log-predictor, line 713). That is a genuine identity between two *separately written* expressions; the script instead compares a formula to itself.

**Why this matters:**
This is the only check intended to anchor Stage 206's predictor to the upstream Stage-239 log predictor. As written it provides zero verification of that anchoring, so a sign or factor error in either the bracket root map or the carried log predictor would pass silently. It also masquerades as a passing test in the output ("Stage 206/188 collapse = 0"), giving false assurance.

**Required change:**
Replace the self-comparison with a genuine two-sided identity: build the certified bracket root map at the degenerate envelope `\(c=L_1\)` from the Section-I oriented form (with `\(K_0=L_0=-\,\)|L_0|`), and independently build the Stage-239 quadratic log predictor `\(\tau_{\log,2}=-2\ln\Phi_0/(L_0+\operatorname{sgn}(L_0)\sqrt{L_0^2-2L_1\ln\Phi_0})\)` from its own appendix formula, then assert their difference simplifies to zero. The two expressions must be written from independent algebra, not one substituted into the other.

**Verification:**
After the fix, line ~148 should assert equality of two expressions that are *not* syntactically identical (one from the Section-I bracket map, one a fresh transcription of eq:app-part06-quadratic-log-predictor); deleting either expression's distinctive subterm must make the assertion fail. Output line "Stage 206/188 collapse = 0" remains, but now backed by a substantive identity.

### F2 — missing_verification_script

**Severity:** high
**Subtype:** missing_mathematica
**Files:**
- (target) `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl`

**What's wrong:**
Stage 206 is non-status-only (`is_status_only_candidate: False`) and `is_checkpoint: False`, but it nonetheless computes a stack of fully symbolic, independently-verifiable results: the closed-form forward root map and the quadratic it solves, the zero-curvature limit, the implicit curvature derivative, the descent-sign identity at the bracket endpoint, the small-envelope width expansion, the zero-mean width law, and the turning-ray roots. Mathematica can derive every one of these from first principles using native primitives (`Solve`/`Reduce` for the quadratic root, `D[]` for the curvature derivative, `Series` + `Coefficient` for the width expansion, `Simplify`/`Refine` under positivity assumptions). No `.wl` exists for this unit, so there is no second-engine cross-check. Per the project dual-engine contract and the auditor prompt (line 118), both engines are required wherever Mathematica *can* verify the stage; it manifestly can here.

**Why this matters:**
A single-engine stage has no guard against a SymPy-specific simplification quirk or an author transcription error that happens to be self-consistent within SymPy. The whole point of the certified-bracket theorem is exactness of the predictor formulas; an unverified second engine leaves the load-bearing root map (carried forward into Stages 207-209 and the realization compiler, appendix lines 824-831, 918) unchecked by an independent symbolic kernel.

**Required change:**
Add a Mathematica audit script at the target path that independently re-derives and asserts the claim manifest M1-M7 below, using a *different decomposition* than the SymPy script. See the directive for the manifest and the anti-transliteration guard.

**Verification:**
`redteam exec-mathematica 206` produces output with each manifest item printing a zero residual (or `True`), and the script exits 0. The verifier additionally confirms the `.wl` does not mirror the `.py`'s variable choreography (independent-derivation check).

### F3 — paper_misalignment

**Severity:** medium
**Subtype:** script_missing_paper_claim
**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex:15`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md:424-469` (§7) and `:473-494` (§8)
- script: (no corresponding assertion anywhere in the SymPy file)

**What's wrong:**
The card's `\stagefield{Output}` (line 15) enumerates four theorems: "Certified monotone bracket theorem, certified turning-ray bracket theorem, **pairwise ray-ordering theorem**, and **local search-sieve admissibility test**." The SymPy script exercises the first two (Sections I-IV) but has *no* check corresponding to the pairwise ray-ordering theorem (notes §7: `\(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\Rightarrow\tau_*^{(a)}<\tau_*^{(b)}\)`) or the search-sieve admissibility test (notes §8). Two of the four enumerated deliverables are unverified.

The pairwise-ordering theorem is the more clearly symbolic of the two: given two certified brackets with `\(\tau_*^{(a)}\in[\tau_{\rm lo}^{(a)},\tau_{\rm hi}^{(a)}]\)`, `\(\tau_*^{(b)}\in[\tau_{\rm lo}^{(b)},\tau_{\rm hi}^{(b)}]\)`, and `\(\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)}\)`, the conclusion `\(\tau_*^{(a)}<\tau_*^{(b)}\)` follows. The sieve admissibility test (§8) is a conjunction of inequality conditions (`\(H_0>0,\ K_0<0\)`, `\(\Delta_\pm\ge0\)`, `\(\tau_{\rm hi}\le T\)`) and is partly a logical convention (§7.1 explicitly calls the lexicographic tie-break "an audit convention rather than a theorem"), so its verifiability is weaker.

**Why this matters:**
The card advertises a "pairwise ray-ordering theorem" as a stage Output; with no script-side check, a reader trusting the verification banner would believe all four Output theorems are machine-checked when two are not. This is the v2 paper↔script gate: each enumerated `\stagefield{Output}` item should map to a check.

**Required change:**
This is a paper_misalignment — direction routed to the user. See `## Resolve before fix_loop` in the directive. The likely resolution is to add a SymPy (and Mathematica) check for the pairwise-ordering implication and, optionally, an admissibility-predicate check; but whether the card's Output is meant to be exhaustively script-verified or partly carries narrative theorems is the user's call.

**Verification:**
After user resolution: if direction (a), a new assertion encoding the interval-ordering implication appears in both engines and the scripts exit 0; if direction (b), the card text is annotated to mark which Output items are narrative-only.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is moot for now — but F2's directive includes an explicit anti-transliteration guard so the *new* script is not a line-by-line port of the SymPy algebra.

## Engine cross-check

Only SymPy is present; no engine cross-check is possible. This is itself F2.

## Verdict justification

Verdict is `findings`, not `clean` and not `stop_cold`. The SymPy script's algebraic predictor checks (Sections I-IV, A1-A6 and A9-A13) hold up under attack: the root map genuinely solves its quadratic, the zero-curvature limit and curvature derivative are non-tautological and match the notes, the width expansion is a real Series/Coefficient identity, and the turning-ray roots are correctly verified — these survive the adversarial pass. Three things do not: (F1) Section V's "collapse" check compares an expression to itself and can never fail, so it verifies nothing; (F2) no Mathematica second engine exists though the stage is fully symbolically verifiable, violating the dual-engine contract; and (F3) two of the four `\stagefield{Output}` theorems (pairwise ordering, sieve admissibility) have no script-side check. None of these mathematically propagates a wrong constant downstream (the root map itself is correctly verified), so no `CRITICAL_DOWNSTREAM`; the premises are consistent, so no `UNFIXABLE`. F3 is a paper_misalignment requiring user resolution, so the directive carries a `## Resolve before fix_loop` block and the orchestrator will halt before invoking Codex.

## Self-test notes

I checked the variable-independence trap on every proposed/audited derivative: `diff(Tau,c)` (A3) genuinely depends on `c` through both the numerator-implicit and the `√(k²-2cH0)` term, so it is not identically zero; `diff(TauTurnA,a)` (A13) depends on `a`; and the F2 manifest's `D[T,c]` likewise depends on `c` — no zero-derivative silent pass. Parity is not in play here (no unbounded-domain integrals; the width "expansion" is a finite Taylor series whose leading-odd term in `\(\Delta K_1\)` is the claimed result, and the script correctly takes the order-3 truncation so the `O(ΔK1³)` remainder is dropped). Trivial-case pre-check: substituting `c→0` into the F2 root-map manifest gives `T=2H0/(2k)=H0/k=-H0/K0`, matching M2; substituting the degenerate envelope `cL=cU=L1` into the bracket endpoints collapses both to the single predictor (the genuine version of F1). Path spec: the F2 target is `mathematica/…_mathematica_audit.wl` (correct directory). Paper round-trip: the F1 fix re-uses the appendix's own log-predictor formula eq:app-part06-quadratic-log-predictor, introducing no new constant, so it creates no new paper_misalignment.
