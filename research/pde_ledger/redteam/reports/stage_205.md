---
unit_id: 205
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md"]
  paper_appendix: present
---

# Audit unit 205 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_205.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 41, 676–729 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card's `\stagefield{Output}` reads verbatim: "Affine and logarithmic quadratic predictors, discriminant tests, curvature-corrected local expansions, and turning-point/tangency theorem." The derivation ledger enumerates: compute `(Φ₀,Φ₁,Φ₂)`, translate between affine and logarithmic quadratic data, derive discriminants, obtain continuous-branch quadratic predictors, and isolate the turning-point criterion `(1−Φ₀)Φ₂>0`. The notes add the full deliverable list: (D1) directional Hessian operators and the bridge identity `Φ₂=Φ₀(L₁+L₀²)` with `L₁=Φ₂/Φ₀−Φ₁²/Φ₀²`; (D2) the quadratic affine predictor `τ_quad=2(1−Φ₀)/(Φ₁+sgn(Φ₁)√Δ_aff)` with affine discriminant `Δ_aff=Φ₁²−2Φ₂(Φ₀−1)`, reducing to `τ_aff=(1−Φ₀)/Φ₁` as `Φ₂→0`; (D3) the quadratic log predictor `τ_log,2=−2lnΦ₀/(L₀+sgn(L₀)√Δ_log)` with `Δ_log=L₀²−2L₁lnΦ₀`, reducing to `τ_log=−lnΦ₀/L₀` as `L₁→0`; (D4) the turning-point theorem — if `Φ₁=0,Φ₀≠1` then real nearby closure points exist **iff** `(1−Φ₀)Φ₂>0`, with predictors `τ_±=±√(2(1−Φ₀)/Φ₂)`, and **no** closure when `(1−Φ₀)Φ₂<0`; the tangency statement `Δ_s(τ)=½Φ₂τ²+O(τ³)` at `Φ₀=1,Φ₁=0`; (D5) the curvature-corrected expansions `τ_quad=τ_aff−(Φ₂/2Φ₁)τ_aff²+O(τ_aff³)` and `τ_log,2=τ_log−(L₁/2L₀)τ_log²+O(τ_log³)`; and (D6) the agreement theorem `τ_log,2−τ_quad=(L₀²+3L₁)/(6L₀³)ε³+O(ε⁴)` with `Φ₀=1+ε`. The appendix (eqs. app-part06-* at lines 688–726) reproduces D1–D4 verbatim, including the `sgn(Φ₁)`/`sgn(L₀)` factors and the strict-inequality turning criterion. All claims are flagged `\StatusExactClosure{}`.

## What the script claims to verify

The SymPy script (docstring-less; banner self-labels "STAGE 188", a renumbering artifact) checks, in five sections: (I) the bridge identities `L₁=Φ₂/Φ₀−Φ₁²/Φ₀²` and `Φ₂=Φ₀(L₁+L₀²)` over fully symbolic `s`, `g`, and a symmetric 5×5 Hessian; (II) the quadratic affine predictor residual is zero and recovers `τ_aff` as `Φ₂→0`, with `Phi0a,Phi1a,Phi2a` declared **positive**; (III) the log predictor residual is zero and recovers `τ_log` as `L₁→0`, with `Phi0l,L0l,L1l` **positive**; (IV) the two turning-point roots `±√(2(1−Φ₀)/Φ₂)` satisfy the reduced quadratic residual, plus a print-only tangency model; (V) the two curvature-correction series vanish through `O(ε²)` and the predictor-gap leading coefficient equals `(L₀²+3L₁)/(6L₀³)`. The bottom-line `expect_zero` assertions all pass (saved output exit code 0).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 bridge identities `Φ₂=Φ₀(L₁+L₀²)`, `L₁=…` | I, lines 77–78 | match |
| D2 affine predictor residual + `Φ₂→0` limit | II, lines 97–98 | partial (only `Φ₁>0` branch; `sgn(Φ₁)` dropped) |
| D2 affine discriminant `Δ_aff` form | line 86 (defined, used in `tau_quad`) | match |
| D3 log predictor residual + `L₁→0` limit | III, lines 116–117 | partial (only `L₀>0` branch; `sgn(L₀)` dropped) |
| D3 log discriminant `Δ_log` form | line 105 | match |
| D4 turning predictors `τ_±` are roots | IV, lines 129–130 | tautological (root-by-construction; reality criterion untested) |
| D4 turning **criterion** `(1−Φ₀)Φ₂>0` iff real | — | missing (headline theorem unverified) |
| D4 tangency `Δ_s=½Φ₂τ²+O(τ³)` | IV, lines 133–135 | partial (print only, no assertion) |
| D5 ordinary/log curvature corrections | V, lines 156, 163 | match |
| D6 agreement theorem `(L₀²+3L₁)/(6L₀³)ε³` | V, lines 168–172 | match |
| Second engine (Mathematica) | — | missing |

Dominant pattern: substantive coverage of D1/D5/D6 and the affine/log discriminant forms, but the headline turning-point criterion (D4) is not exercised, both predictor branches are restricted to positive slope, and there is no second engine. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 77 | `expect_zero(L1 - (Φ₂/Φ₀ - Φ₁²/Φ₀²))` | D1 | yes |
| A2 | sympy | 78 | `expect_zero(Φ₂ - Φ₀(L₁+L₀²))` | D1 | yes |
| A3 | sympy | 97 | `expect_zero(quadratic affine residual)` | D2 (positive branch) | partial |
| A4 | sympy | 98 | `expect_zero(limit Φ₂→0 = τ_aff)` | D2 limit | yes |
| A5 | sympy | 116 | `expect_zero(quadratic log residual)` | D3 (positive branch) | partial |
| A6 | sympy | 117 | `expect_zero(limit L₁→0 = τ_log)` | D3 limit | yes |
| A7 | sympy | 129 | `expect_zero(turning-point root (+))` | D4 predictor | no (tautological) |
| A8 | sympy | 130 | `expect_zero(turning-point root (-))` | D4 predictor | no (tautological) |
| A9 | sympy | 156 | `expect_zero(ordinary quadratic correction)` | D5 | yes |
| A10 | sympy | 163 | `expect_zero(log quadratic correction)` | D5 | yes |
| A11 | sympy | 168–172 | `expect_zero(gap leading coeff)` | D6 | yes |
| — | sympy | 133–135 | print only (tangency model) | D4 tangency | no assertion |

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `(missing)` — target `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl`

**What's wrong:**
The unit's manifest entry has `is_status_only_candidate: False` and the stage is `is_checkpoint: False`, yet only a SymPy script exists; there is no Mathematica `.wl` anywhere under `mathematica/` (203 other stages present, none for 205). The stage computes purely symbolic, machine-verifiable results — quadratic-form contractions of a 5×5 symmetric Hessian, exact predictor residuals from the quadratic formula, two zero-curvature limits, two `O(ε²)` series corrections, and a cubic-order predictor-gap coefficient. Every one of these is independently reproducible in Mathematica via native primitives (matrix algebra, `Solve`/`Reduce`, `Series`+`Coefficient`, `Limit`, `Refine`). Per the project dual-engine contract and the prompt's line ~118 ("both scripts are required, and missing scripts are findings"), the absence is a finding. There is no legitimate "cannot independently verify" carve-out here: the math is exactly the kind a CAS verifies.

**Why this matters:**
A single-engine stage has no cross-check; a transliteration-free second engine is the project's guard against an entire family of SymPy-specific `simplify`/branch artifacts (e.g. the `sgn`-dropping in F2 and the tautology in F3 below would both be caught by an independent re-derivation that does not start from the same positive-slope assumptions).

**Required change:**
Author a new independent Mathematica audit per the claim manifest in the directive (F1). It must independently re-derive each manifest item via a different decomposition than the `.py`, not transliterate the SymPy algebra.

**Verification:**
Verifier runs the new `.wl` via `math -script`; it must exit 0 with every manifest check passing and the file must not be a line-by-line port (independent-derivation check).

### F2 — missing_branch

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage205..._sympy_audit.py:84,88,97` (affine) and `:104,107,116` (log)

**What's wrong:**
The paper's predictors carry an explicit sign factor (appendix line 710: `\tau_{\rm quad}=\frac{2(1-\Phi_0)}{\Phi_1+\operatorname{sgn}(\Phi_1)\sqrt{\Delta_{\rm aff}}}`; line 714 likewise `L_0+\operatorname{sgn}(L_0)\sqrt{\Delta_{\log}}`; notes §5/§6 state "choose the square-root branch continuously … assume `Φ₁≠0`/`L₀≠0`"). The script declares `Phi1a` (line 84) and `L0l` (line 104) as `positive=True`, so `sgn(·)=+1` and the script never exercises the negative-slope branch `sgn(Φ₁)=−1` / `sgn(L₀)=−1`. The subbanner comments themselves label these "positive-slope continuation branch" (lines 81, 100) — i.e. only the easy branch is tested while the paper formula is general in the sign of the slope.

**Why this matters:**
The continuous-branch claim is exactly that the predictor reduces to the Stage-238 first-order predictor for **either** sign of the slope. With `Φ₁<0`, `τ_quad` must use `Φ₁−√Δ_aff` in the denominator to stay finite and to reduce to `τ_aff=(1−Φ₀)/Φ₁`; the positive-only test cannot detect a sign error in the `sgn` placement.

**Required change:**
Add residual + zero-curvature-limit checks for the negative-slope branch (Codex designs the route; state-of-requirement in directive F2): re-declare slope symbols as nonzero-real (not positive) or add a second block with `Φ₁<0`/`L₀<0`, using the `sgn`-correct denominator `Φ₁+sgn(Φ₁)√Δ_aff` (resp. `L₀+sgn(L₀)√Δ_log`), and assert the residual of `(Φ₀−1)+Φ₁τ+½Φ₂τ²` (resp. `lnΦ₀+L₀τ+½L₁τ²`) is zero and that the `Φ₂→0` / `L₁→0` limit recovers `τ_aff`/`τ_log` on that branch too.

**Verification:**
New negative-branch `expect_zero` lines appear and pass; SymPy exits 0.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage205..._sympy_audit.py:123–130`

**What's wrong:**
`tau_tp = sqrt(2*(1-Phi0t)/Phi2t)` (line 124), then `res_tp_plus = (Phi0t-1) + (1/2)*Phi2t*tau_tp**2` (line 125). Because `tau_tp**2 = 2*(1-Phi0t)/Phi2t` by construction, the residual collapses to `(Φ₀t−1) + ½·Φ₂t·(2(1−Φ₀t)/Φ₂t) = (Φ₀t−1)+(1−Φ₀t) = 0` identically — the same is true for `res_tp_minus` since `(−tau_tp)²=tau_tp²`. The assertion is algebraically guaranteed: "a constructed root of the reduced quadratic is a root." It cannot fail regardless of whether the root is real, and `Phi0t,Phi2t` are declared plain real (line 123) so `tau_tp` may be imaginary with the check still passing.

**Why this matters:**
The actual headline theorem (paper §7.1 / appendix eq:app-part06-turning-criterion, line 717–721) is the **reality criterion**: real nearby closure points exist **iff** `(1−Φ₀)Φ₂>0`, and **no** closure when `(1−Φ₀)Φ₂<0`. That iff — the substantive content of D4 — is never tested. The current check is decorative and gives false confidence that the turning-point theorem is verified.

**Required change:**
Replace the tautological root residuals with a test of the reality criterion (Codex designs the route): with `Phi0t,Phi2t` real, assert that `tau_tp` is real ⟺ `(1−Φ₀t)Φ₂t≥0`, e.g. verify the radicand `2(1−Φ₀t)/Φ₂t≥0` is equivalent to `(1−Φ₀t)Φ₂t≥0` (via `Reduce`/sign analysis), and confirm that under `(1−Φ₀t)Φ₂t<0` the predictor has no real root (radicand `<0`). Keep one non-tautological confirmation that the symmetric pair `±√(2(1−Φ₀)/Φ₂)` solves `(Φ₀−1)+½Φ₂τ²=0` only where the criterion holds.

**Verification:**
A new assertion exercising `(1−Φ₀)Φ₂>0 ⟺ real root` appears and passes; the old line-129/130 by-construction residuals are removed or downgraded to non-load-bearing.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage205..._sympy_audit.py:132–135`

**What's wrong:**
The tangency case (paper §7.2: `Δ_s(τ)=½Φ₂τ²+O(τ³)` at `Φ₀=1,Φ₁=0`) is only *printed* (`Delta_tangent = ½*Phi2g*tau**2`, lines 133–135) with no `expect_zero`/assertion. It is asserted nowhere that this is the actual second-order model obtained by setting `Φ₀=1,Φ₁=0` in the quadratic residual `(Φ₀−1)+Φ₁τ+½Φ₂τ²`; the line just hard-writes the answer.

**Why this matters:**
A print without an assertion verifies nothing; it is a hardcoded restatement of the claimed tangency form rather than a derivation of it.

**Required change:**
Add `expect_zero` confirming that substituting `Φ₀→1, Φ₁→0` into the ordinary quadratic model `(Φ₀−1)+Φ₁τ+½Φ₂τ²` yields exactly `½Φ₂τ²` (i.e. `model|_{Φ₀=1,Φ₁=0} − ½Φ₂τ² == 0`).

**Verification:**
New `expect_zero("tangency model ...", ...)` line appears and passes.

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration assessment applies yet. The directive's F1 includes an explicit anti-transliteration guard so the new script is an independent re-derivation, not a port.

## Engine cross-check

Only one engine present; no cross-check possible. This is itself F1.

## Verdict justification

The SymPy script faithfully and non-tautologically verifies the bridge identities (D1), the two curvature-correction expansions (D5), and the cubic-order predictor-agreement coefficient (D6) — I tried to break A1/A2 by treating the Hessian as a general symmetric form (it holds) and A11 by checking the series order (the leading `ε³` coefficient is exact). However three substantive gaps remain: (F3) the headline turning-point *reality criterion* `(1−Φ₀)Φ₂>0` is never tested — the root-residual checks are zero by construction; (F2) both predictor branches are silently restricted to positive slope despite the paper's general `sgn(·)` formulas; (F4) the tangency form is printed, not asserted. And (F1) there is no second engine for a non-status, non-checkpoint stage whose math a CAS can plainly verify. Hence `verdict: findings`, no stop-cold (none of the fixes flip a quoted forward or a downstream-consumed constant — the predictor formulas are tools, and the values fed to them are PDE/candidate data per appendix line 729). `paper_alignment: partial`. No `paper_misalignment`: every script check that exists matches the paper's stated form; the gaps are under-coverage, not contradiction.

## Self-test notes

Trap 1 (variable independence): the new `.wl` items are matrix-algebra and `Series` in `eps`; `tau_quad_e`, `tau_log2_e` genuinely depend on `eps` via `Φ₀=1+eps`, so no identically-zero-derivative trap. Trap 2 (parity): no unbounded integrals; N/A. Trap 3 (trivial case): verified the F3 reality criterion by hand — `tau_tp²=2(1−Φ₀)/Φ₂`, real iff `(1−Φ₀)Φ₂≥0`, matching the boxed theorem; and confirmed the current residual is identically zero by construction (the tautology). Trap 4 (paths): `.wl` target placed in `mathematica/` with the full canonical filename. Trap 5 (paper round-trip): the F2 negative-branch fix uses the paper's exact `Φ₁+sgn(Φ₁)√Δ_aff` / `L₀+sgn(L₀)√Δ_log` forms (appendix lines 710, 714), introducing no new constant or misalignment.
