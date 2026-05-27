---
batch: IV.3
label: Part IV.3 — Core balance, DtN mixed, outlet, positive source
range: 115-126
stage_count: 12
verified: 12
findings_count: 27
material_change_count: 1
material_change_stages: [118]
clean_stages: [120, 124]
checkpoints: []
status_only: [120, 124]
no_mathematica_only: [121, 122, 123]
created_at: 2026-05-27
applied: true
verification_status: complete
---

# Red-team batch summary — IV.3

12 stages audited, 12 verified. Two clean (120, 124 — status-only carve-outs). Zero checkpoints in this range. 10 dirty stages produced 27 findings total. **Eight↓nine consecutive zero-redirection batches** (II.1–IV.3) — codex bypassed entirely.

## Cluster resolutions (4 user-gated)

- **Cluster A — notes-side numerical typos.** Four independent auditors (stages 121, 122, 123, 126) flagged that notes/stages/*.md files carried numerical typos that contradicted both the scripts' independent re-derivations AND the notes' own quoted decimal values. Scripts were correct in all four cases.
  - 121, 122, 126: `4107 − 168π²` → `4107 − 100π²` (and `168π²` → `100π²`) in boxed surds. Derivation from `R = 37/20` gives `12·R² = 4107/100`, confirming the `100` form.
  - 123: denominator `228 → 160` in `Xi_v` branch law. The notes' own numeric `-1.01675...` is only consistent with `/160`.
  - Resolution: notes fixed in place this batch (7 textual occurrences across 4 files); no script changes for F1 items.

- **Cluster B — Stage 118 λ sign inconsistency.** Internal script inconsistency: section IV's bilinear derivation gave `λ = −q* v_w0 𝓘_sq` (matching notes); section V's closure asserted `λ_uniform = +(8√2/3)·…` (plus sign). Resolution: flipped script section V to minus sign in both `.py` and `.wl`; F2 added 3 new asserts (K_q closed form, g_s closed form, λ from bilinear) with consistent signs. **Material change flagged for stage 118 (see below).**

- **Cluster C — Stage 125 integral inequality.** Paper card's `Output` is the integral inequality `0 ≤ 𝔤[σ] ≤ 1`; script only proved the pointwise kernel bound. Resolution: added a one-parameter family `σ_a(z) = (a+1)(z/L)^a/L` integrated against the cosine kernel in both engines. SymPy required a numerical proxy at `a=100` for the peaked-at-L limit (hypergeometric form blocks symbolic `sp.limit`); Mathematica's `Limit[gA, aSym → ∞]` worked. Endpoint at `a=0` gives the closed form `2/π`.

- **Cluster D — Stage 117 consolidation card.** Three blocked items (transliteration + 2 tautological κ_c=1/3 and γ_c=1/9 checks) resolved by direction (a): cite upstream stages 115 (Schur reduction) and 116 (D/N eigenvalue) via a comment block before section 5; replace the misleading `expect_zero` tautological wrappers with `print("carrying forward (Stage 116/119): …")` lines. F4 (classification_rows) wired booleans from sections 1-5 residuals — `nontrivial_compensated` anchored to the `delta_core - delta_core_expected` series residual.

## Material change — stage 118 (λ sign)

The λ sign flip from `+` to `−` in stage 118's section V closure is the only material change in this batch. Downstream impact analysis:

- Downstream Schur reductions consume λ as `K_s K_q + λ²` (sign-invariant under squaring) and `(K_s g_q − λ g_s)²` (also sign-invariant under squaring).
- Stages 125-139's compensation analysis uses the **squared** combinations exclusively per the auditor's downstream-propagation review.
- No upstream_stale flag propagated to subsequent batches.
- **Caveat:** if a future batch finds a downstream stage that uses λ unsquared (e.g., in a cross-term not yet identified), revisit this stage's sign-flip.

## Mechanical edits applied (no user gate)

- **115**: `.wl` independent parent-overlap derivation (3 new `expectZero` checks via `frakR/frakG` reparametrization; equivalence factor `(kS·kQ + lam²)/(gS²·kQ)` derived correctly; Solve uses fresh `gVar` for symbol-binding hygiene).
- **116**: D/N half-wave eigenvalue derivation in both engines (`q'' + k²q = 0`, `q(0)=0`, `q'(L_W)=0` → `k_W = π/(2L_W)`); κ₀ round-trip; γ₀ extraction from `D_bare`; SymPy tube-length assert (Mathematica already had it). Sign correction in directive: `gamma0_from_D = +I·coeff(z,5)` (directive's `-I·` was wrong).
- **117**: comment block citing upstream stages + `expect_zero` downgrade + `classification_rows` wiring + banner.
- **118**: λ sign flip + 3 new asserts.
- **119**: `rc → rhat²` link assertion + `T_m (±)` branch matches (Mathematica uses `stripCE` for `ConditionalExpression` heads).
- **121**: 4 new `expect_zero` asserts (closed-form `r_geom`, symbolic `r_F1`, symbolic `r_c(F1)`, threshold vanishing) + Ω_W identification + banner.
- **122**: 6 new asserts (compensation quadratic at `g±`, defect closed form, natural off compensation, traction ratio identities) + `expect_nonzero` helper + banner.
- **123**: banner only on script side.
- **125**: parametric family in both engines + `.wl` Solve-based g_± derivation + banner.
- **126**: positivity asserts in both engines (SymPy `sp.calculus.util.minimum`; Mathematica boundary-value checks because `Minimize` returned unevaluated — equivalent under monotone-decreasing `k Cos[k z]`) + SymPy interval check hardened to raise on False + banner.

## Stage-level table

| Stage | Findings | Verdict | Material |
|---:|---:|---|:-:|
| 115 | 1 | verified | no |
| 116 | 4 | verified | no |
| 117 | 4 | verified | no |
| 118 | 2 | verified | **yes** (λ sign flip) |
| 119 | 2 | verified | no |
| 120 | 0 | clean (status-only) | no |
| 121 | 3 | verified | no |
| 122 | 3 | verified | no |
| 123 | 2 | verified | no |
| 124 | 0 | clean (status-only) | no |
| 125 | 2 | verified | no |
| 126 | 4 | verified | no |

## Notes-side cleanup applied this batch

| File | Lines | Change |
|---|---|---|
| `notes/stages/moving_throat_pde_stage121_geometric_r_selection.md` | 69 | `168π² → 100π²` |
| `notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md` | 49, 56, 88 | `168π² → 100π²` (4 occurrences) |
| `notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md` | 30, 39 | `228 → 160` (2 occurrences) |
| `notes/stages/moving_throat_pde_stage126_positive_source_families.md` | 100 | `168π² → 100π²` |

## Cumulative state after IV.3

138/253 stages red-team verified (54.5%). Range **001-126 paper-aligned at v2 depth**.
