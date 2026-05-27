---
batch: IV.3
created_at: 2026-05-27
total_paper_misalignment_findings: 7
total_blocked_findings: 3
clusters: 4
---

# Batch IV.3 paper-alignment consolidation

12 audit directives produced 27 findings total. Two stages clean (120, 124 — status-only carve-outs). This document consolidates the 7 `paper_misalignment` items + 3 `blocked` methodological items into 4 clusters for user decision before fix_loop.

The remaining 17 findings (transliterations, tautological checks, insufficient verification) on stages 115, 116, 117 F4, 119, 121 F2/F3, 122 F2, 123 F2, 125 F2, 126 F2/F3/F4 are mechanical and will be applied per directive after the user gate clears.

---

## Cluster A — Notes-side numerical typos (stages 121, 122, 123, 126)

**Recurring pattern:** Four independent auditors flagged that prose `notes/stages/*.md` files quote symbolic surds that contradict both the scripts' independent re-derivations AND the notes' own quoted numerical values. The scripts are correct; the notes have transcription typos.

### A.1 — `168π² → 100π²` (stages 121, 122, 126)

Three stages share the same `4107 − 168π²` typo for what should be `4107 − 100π²`. Independent derivation from `R = 37/20` gives `12·R² = 4107/100`, confirming the `100` form. Affected boxed equations:

| Stage | Notes file | Lines | Wrong | Correct |
|---|---|---|---|---|
| 121 | `notes/stages/moving_throat_pde_stage121_geometric_r_selection.md` | 67 | `sqrt(4107 - 168 pi^2)/(10 pi)` | `sqrt(4107 - 100 pi^2)/(10 pi)` |
| 122 | `notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md` | 50, 56, 88 | `4107 − 168π²`, `168π²` | `4107 − 100π²`, `100π²` |
| 126 | `notes/stages/moving_throat_pde_stage126_positive_source_families.md` | 100 | `2*sqrt(4107 - 168*pi^2)` | `2*sqrt(4107 - 100*pi^2)` |

All affected stages' quoted **numerical** values (`r_F1 ≈ 1.778`, `g_-^F1 ≈ 0.758`, `xi_* ≈ 0.184`, `C_nat ≈ 1.740`) match the `100` form. No script change needed.

### A.2 — Denominator `228 → 160` (stage 123)

Same pattern, different number. Stage 123 notes (lines 25-31, 37-41) box `Xi_v = -(3√30·π^(3/2)/228)·r_{F1}` but the script and the notes' own numerical box `-1.01675...` are only consistent with denominator `160`. Manually verified: `-3√30·π^(3/2)·1.778/160 ≈ -1.0168` ✓; with `/228` ≈ -0.713 ✗.

### Resolution direction options

- **(a) Notes-side typo confirmation (Recommended).** Fix the 4 notes files in place. No script changes for these F1s. Defer to a notes-cleanup tracker entry analogous to IV.2's P3-13.
- **(b) Revisit upstream sources.** Re-derive `r_F1`, `g_-^F1`, `Xi_v` from upstream stages 220-223 to independently confirm which form is correct. (Auditors already did informal independent derivations; all yielded the `100` and `160` forms.)
- **(c) Defer the notes fix to PAPER_CLEANUP_TRACKER.** Mark as known notes-side drift; no edits this batch.

**Orchestrator note:** Per `feedback_redteam_priorities` ("script math correctness first, alignment with written docs second; notes/ folder is per-topic anchor when present"), the notes typos do not invalidate any script result. The scripts pass all their assertions; the typos only affect future readers of the notes.

---

## Cluster B — Stage 118 λ sign inconsistency (script-side, high severity)

**Auditor verdict:** Notes (`notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:174-182, 194-199`) box `λ = − q_* v_{w0} 𝓘_{sq}` (minus sign). Script's section IV (lines 71-77) **agrees** with the notes minus sign. But script section V closure asserts `λ_uniform = + (8√2/3) q_* v_{w0} a² ℓ √L_W` (plus sign). This is an **internal script inconsistency**, not a notes-vs-script drift.

Downstream propagation: stage 118 feeds stages 125-139. Downstream `λ²` and `K_s K_q + λ²` Schur terms are sign-invariant, but the cross-term `−2 K_s g_q · λ g_s` in `(K_s g_q − λ g_s)²` is sign-sensitive wherever it appears unsquared.

### Resolution direction options

- **(a) Notes are correct (script section V wrong) (Recommended).** Flip the sign on:
  - SymPy line 88: `lam_uniform = sp.simplify(-qstar * v0 * I_sq_uniform)` (add minus)
  - SymPy line 100-101 target: `lam_uniform + 8*sp.sqrt(2)*qstar*v0*a**2*ell*sp.sqrt(L_W)/3` (flip sign in subtraction)
  - Mirror in Mathematica lines 96 and 108.
  - Then apply F2 (3 new asserts: K_q, g_s, λ general form) with consistent signs.
- **(b) Script section V correct, notes + script section IV wrong.** Unlikely per auditor (Madelung kinetic expansion is explicit in notes and matches script section IV). Would require notes edit + section IV re-derivation.
- **(c) Hidden sign convention.** Document the convention and leave signs as-is. Unlikely per auditor.

---

## Cluster C — Stage 125 integral inequality verification gap (substantive)

**Auditor verdict:** Paper card's only quoted `Output` claim is the **integral** inequality `0 ≤ 𝔤[σ] ≤ 1` for any positive normalized σ, with moment representation `𝔤[σ] = ∫₀^L σ(z) cos(π z/(2L)) dz`. Script only verifies the **pointwise** kernel bound (`0 ≤ cos(π z/(2L)) ≤ 1` on `[0,L]`). The convex-combination step from pointwise to integral — the actual content of the theorem — is absent.

### Resolution direction options

- **(a-i) Parametric family extension (Recommended).** Add a one-parameter positive normalized family `σ_a(z) = (a+1)(z/L)^a / L` (a ≥ 0; nonneg + ∫=1) in both engines. Compute `g_a = ∫ σ_a · cos(πz/(2L)) dz` symbolically. Assert bounds at endpoints (`a=0` → `g = 2/π`; `a → ∞` → `g → 0`).
- **(a-ii) Direct pointwise-to-integral lemma.** Assert integrand-level positivity (`1 − cos(π z/(2L)) ≥ 0` and `cos ≥ 0` on `[0,L]`) as `expect_zero` checks and narrate the convex-combination step.
- **(b) Trim the paper card claim** to be the pointwise bound only. Per auditor, unlikely correct — the notes and appendix call out the integral form as the central deliverable.
- **(c) Document the implicit convex-combination** in narrative prints; no new assertion. Weakest — a future kernel edit could break silently.

Auditor recommends **(a-i)**: exercises the actual integral, gives concrete numerical anchor values, and would catch a half-period typo.

---

## Cluster D — Stage 117 blocked methodological items (consolidation card)

Stage 117 (concrete outlet-core status) consolidates upstream results rather than re-deriving them. All three blocked items (F1 transliteration, F2/F3 tautological-by-construction κ_c=1/3 and γ_c=1/9) share the same root cause: this stage's scripts echo upstream κ₀/γ₀/Schur-reduction algebra without forward derivations.

### Resolution direction options

- **(a) Cite upstream forward expressions (Recommended).**
  - F1: Cite stage 116 (which is being patched in this batch with an independent D/N eigenvalue derivation for κ₀) for the Mathematica side; the parent-overlap reduction in stage 115's directive provides the Schur form. Add a comment block in 117's `.wl` pointing readers to the upstream derivations rather than re-running them.
  - F2/F3: Replace `expect_zero("D/N tube fixes kappa_c = 1/3", …)` with definitional substitutions wrapped in non-assertion `print` lines (e.g., `print("Stage 116 fixes kappa0_bare = (1+r_c)/3 via D/N tube; carrying forward.")`). The substantive κ_c = κ₀/(1+r_c) = 1/3 reduction is then an arithmetic consequence, properly documented as such.
- **(b) Accept transliteration limitation.** Update redteam policy doc to allow consolidation cards to echo upstream derivations. Add `Stage 117 is a status-consolidation card; its scripts echo upstream stage 115/116 algebra by design` to the stage's notes header.
- **(c) Leave as-is with explanatory comments.** Weakest — preserves the misleading `expect_zero` wrappers.

Auditor recommends user direction; (a) is the most defensible and matches IV.2's stage-103/113 status-only carve-out pattern.

---

## Mechanical findings (no user gate; applied per directive)

| Stage | Findings | Type | Action |
|---|---:|---|---|
| 115 | 1 | mathematica_transliteration | Insert independent parent-overlap derivation block in `.wl` |
| 116 | 4 | hardcoded + 2 tautological + insufficient | D/N eigenvalue derivation + κ₀ round-trip + γ₀ extraction + tube-length assert |
| 117 F4 | 1 | tautological (capstone) | Wire `classification_rows` booleans from sections 1-5 residuals |
| 119 | 2 | tautological + insufficient | `rc → rhat²` substitution check + `T_m` branch matches |
| 121 F2, F3 | 2 | insufficient + script_missing_paper_claim | 4 new asserts + Ω_W check |
| 122 F2 | 1 | insufficient | 6 new asserts (compensation quadratic + defect + traction) |
| 123 F2 | 1 | banner mislabel | `STAGE 106` → `STAGE 123` |
| 125 F2 | 1 | mathematica_transliteration | Derive g_± via `Solve` in `.wl` |
| 126 F3, F4 | 2 | hardening + banner | Replace silent interval check with raise + banner sweep |

Banner sweep also implicit in 121, 122, 123, 126 (consistent with IV.2 Cluster C pattern).
