# Claude↔Codex math-direction consult — red-team remediation BATCH 3 (stages 117, 118)

- **Date:** 2026-05-29
- **codex_session_id:** `019e74f7-7d0d-7122-a304-136dde942533` (read-only, ephemeral)
- **Verdict:** Codex **CONCUR on all four** (Q1–Q4). Both caveats discharged by the orchestrator (below).
- Raw transcript was deleted after extracting this clean record (per the transcript-bloat policy; the prompt forbade repo-wide grep, raw stayed ~8.6 KB).

## Q1 — Stage 118 F1: λ sign (paper_misalignment / target_mismatch) → **direction (a)** (notes/section-IV correct; §V is MINUS)

**CONCUR.** Evidence: script **section IV independently derives** the bilinear coefficient with a MINUS (`sympy_audit.py:77` `sq_coeff + qstar*varrho_s*v0*A_q`; `mathematica_audit.wl:81`), not asserted; notes carry MINUS in three boxes (`stage118_…extraction.md` 169-171, 176-180, 196-197); downstream **stage 123** consumes UN-squared λ already with MINUS (`stage123_…sympy_audit.py:28,32`), so (a) makes 118 *consistent* with 123 (old PLUS would have contradicted it).

**Codex caveat:** "flag as CONCEPTUAL if the paper's §V states the opposite physical coupling sign." **DISCHARGED** — `paper/stages/stage_118.tex` states no λ sign (only line 16 "…K_s, K_q, λ, g_s, g_q from a … throat-core ansatz"). No opposite-sign paper claim ⇒ not conceptual, no PAPER_CLEANUP item.

**Status:** F1 + all three F2 asserts are **tainted-applied** (present in both engines with the MINUS form, `applied:false`/unrecorded). Fix step = Codex reconciles + records `## Applied: F1/F2` + runs to exit 0; no math change.

## Q2 — Stage 117 F2: `kappa_c = 1/3` tautological_check → **de-tautologize by importing the upstream forward expression** (script change)

**CONCUR.** κ₀ IS genuinely forward-derived upstream at **stage 116** (D/N tube BVP solved from scratch: `stage116_…sympy_audit.py:27-39` → `k_W=π/(2L_W)`, 43-47 `κ₀==4 L_W²/(π²a²)`, 49-62 tube-length law forcing `κ₀=(1+r_c)/3`). Stage 117 already computes `Lw_required` (`stage117_…sympy_audit.py:164`) then discards it for a hardcoded `(1+r_c)/3` (≈py:174). Fix: build κ₀ from `4*L_W_required**2/(pi**2*a**2)`, `κ_c = κ₀/(1+r_c)`, mirror in `.wl`. Codex: "script coverage/provenance, not conceptual."

## Q3 — Stage 117 F3: `gamma_c = 1/9` tautological_check → **correct the false provenance comment**; check stays as assumption-consistency

**CONCUR.** γ₀=(1+r_c)/9 is **NOT derived anywhere** — a pure-scale **ansatz** (modeling choice), postulated in `stage116_…dn_mixed_tube_realization.md:57-75`, hardcoded `stage115_…sympy_audit.py:47` & `stage116_…sympy_audit.py:73` (latter explicitly: "upstream-carried input, not derived in-stage; an expect_zero here would be tautological"). The stage-117 comment falsely said "derived (Stage 119)" (`stage117_…sympy_audit.py:168`, also py:142 "derived upstream"). Fix: correct the comment to "modeling assumption / pure-scale ansatz" citing the stage-116 note, drop the false "Stage 119"/"derived" wording; the `γ_c=1/9` check is acknowledged consistency-of-assumption (cannot be de-tautologized since γ₀ is assumed).

**Codex caveat:** "Not conceptual unless the project intends γ₀ to be *derived* rather than *postulated*." **DISCHARGED** — the notes + stage-116 explicitly postulate it. → recorded as the first entry in the planned ansatz catalog (memory `ansatz-value-catalog`).

## Q4 — Stage 117 F1: mathematica_transliteration → **accepted policy-mirror** (no independent-route rewrite)

**CONCUR.** The `.wl` is a near-line-by-line port; the load-bearing payload is pure rational series-coefficient classification (which outlet-deformation classes preserve `Ŷ₂` to O(z⁵)). A Padé/residue reformulation adds little adversarial value; cross-engine agreement still has worth. Per the mirror policy, accept as a sanctioned mirror.
