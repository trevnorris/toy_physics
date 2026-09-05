# S11c-c2 SHARED_PHYSICS — decision-gate record (rule 7 + rule 2)

The orchestrator-written c2 spec (the self-energy fold) got **two decision legs before any builder** (rule 7).
Both legs did the ⚙ computation demanded and **converged** on the same headline defects; neither softened. The v1
draft was **not sound to build against**. Folded **once** → v2 (rule 7: one pass, then go — the folded spec is
re-checked at the build-directive gate, not re-legged as a standalone spec). ⛔ No commit before both legs reported
(rule 9); both are in.

## Commands (rule 2)

```
# Leg prompt (identical to both): directives/_legs/S11c_c2_spec_review_prompt.md  (8839 bytes)
# Codex sol xhigh (DOC REVIEW model policy):
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
  "$(<directives/_legs/S11c_c2_spec_review_prompt.md)" > <log> 2>&1        # exit 0, 1.88 MB
# Grok:
grok --prompt-file directives/_legs/S11c_c2_spec_review_prompt.md \
  --cwd .../research/pde_ledger_v3 --model grok-4.6 --effort high \
  --permission-mode bypassPermissions --output-format plain > <log> 2>&1   # exit 0
# Reports committed: directives/_measurements/S11c_c2_spec_review_{grok,codex_sol}.txt
# Leg ⚙ scripts (in the reports): /tmp/s11cc2_q{1,2,5,5e,7}*.py + .stdout ; /tmp/s11cc2_review_sympy.{py,stdout}
```

Orchestrator verification (rule 13) — each headline finding confirmed against the **authoritative source**, not the
legs' toy scripts:

```
# F4 leak: the substituted-slot sign conventions are cross-engine-UNVALIDATED (S11c-b's own step record)
steps/S11c_b_variable_coefficient_operator.md:112-115  ("face generalized-force PY +diff vs WL −linearVirtualVariation"; "#90 closure-fold sign")
# Both legs' SymPy: cross-engine increment with opposite carrier = Matrix([[2*beta*c, 2*alpha*f]])  (carriers c,f survive)
# F1 keys: bare face_response is in c1 IMPORT_KEYS (imported FROM S11c-b); c1 write-key is s11c_c1_face_response
scripts/S11c_c1_exports.py:57-58 (IMPORT_KEYS 'face_response'/'_coeffs')  vs :101 ('s11c_c1_face_response')
# F2 Λ_X: closure_shape_deriv (real fold row) carries Λ_A/Λ_V + delta_p_± symbolic, NO Λ_X; Λ_X is only in traction
c1 §1d S11c_c1_SHARED_PHYSICS.md:157-160
# F5 coupling_kernel: the GATED design bars importing the open kernel as a construction operand
directives/export_ledger_bind_closure_design.md:148-153  ("c2 re-extracts … does not bind the open coupling_kernel")
# F6 API: def load_model(base_path, *delta_paths)
scripts/ledger_fold.py:102-103
# θ-row real structure + density map (python3 load_model fold inspection, this session):
#   closure_shape_deriv value: Λ_A(−delta_p_plus/rho_m + mu_theta_L/rho_br)/(ωτ_A+I) + Λ_V W_0 e_W_t/(ωτ_V+I) + d_w_delta_p_plus terms  → J_s ALREADY eliminated; δp_s symbolic
#   background_density_map: Eq(rho_br_bg_rho4_constant, W_bg*rho_br/W_0)  [live] ; Eq(rho_br_bg_rhobr_constant, rho_br) [const]
```

## Findings (union of both legs; all CONFIRMED against source) and how v2 folds them

| # | sev | Finding (spec v1 → source it contradicts) | Legs | Fold in v2 |
|---|---|---|---|---|
| 1 | MUST | v1 §3c/§7 claim the increment cancels **all** deferred S11c-b terms incl. the two sign conventions. FALSE: the face-generalized-force + #90 closure-fold conventions are the **coefficients of the δp_s slots** — they multiply the increment and leak (`2·carrier·incr`). Only the **bulk/kinetic** base cancels. | Grok F4 + Codex 1 (both ⚙ `[2βc,2αf]`) | §3c rewritten: increment drops the **bulk/kinetic** base only; the two conventions are SURFACED by the comparator (§7) and **adjudicated by the new §3d.4 mechanical-power pairing**. §1c/§3c self-contradiction resolved toward §1c. |
| 2 | MUST | v1 says "substitute closed `J_s`", but the real θ-row has `J_s` **already eliminated** (#90); it carries `Λ_A𝒜+Λ_V V` with `δp_s`/`d_w_delta_p_s` symbolic. And v1 §1c wrongly puts `Λ_X` in `closure_shape_deriv`. | Codex 2/3 + Grok F2 | §1c/§1d/§3a rewritten from the real row: the fold **substitutes closed `δp_s`+w-jets into the symbolic `delta_p_±`/`d_w_delta_p_±` slots** (⛔ no closed-`J_s` add — double-count); `Λ_X` is traction-only, routed via `traction` to the mechanical rows. |
| 3 | MUST | v1 emits the "full operator per `(α,ρ,s)`", but `slab_operator` is a **two-face** object (`J₊+J₋`, both tractions summed). | Codex 3 + Grok F7 | §3a: emit the **assembled two-face operator per `(α,ρ)`** + per-face provenance + parity blocks. §0/§4 updated. |
| 4 | MUST | The "exactly three" re-adjudications omit load-bearing obligations: the **traction-vs-slab pairing** c1 assigned to c2; the **flat-resolvent leg-labeling** (enters via MATERIAL off-diag response); the **`μ_R,bg` form control** c1 reserved for c2. | Codex 4 + Grok F4a/F8 | §3d expanded to **six** re-adjudications (density, t_s, DtN kernel/whole-form, traction-slab power pairing, flat-symbol usage, μ_R,bg form); §5 gains the pairing + μ_R,bg ablation. |
| 5 | MUST | Three inequivalent operators all named "self-energy" (§2 commutator remainder / §3b `extract(close)` / §3c increment). | Grok F3 + Codex 5 | §2 names three distinct objects: `CLOSED_COUPLING_KERNEL`, `ORDERING_COMMUTATOR`, `SELF_ENERGY_INCREMENT`. |
| 6 | MUST | Export wiring: `load_model(base=…, deltas=[…])` TypeErrors; §7 binds bare `face_response`/`_coeffs` (S11b open rows) as if c1's; `w1_profile_hat_transfer` not a write-key; §3c subtracts the open `coupling_kernel` (barred). | Grok F1/F5/F6 + Codex 6 | §7 rewritten: **positional** `load_model(base, c1-delta)`; consume-set named by **F9 write-keys** (`s11c_c1_*`); open `coupling_kernel` regression-only (§5a), re-extract for §3c; guard-passes-on-existence caveat stated. |
| 7 | MUST | §5e falsely equates `Z→0 ≡ Λ_A⁰=0`; impermeability is `Λ_A⁰=Λ_V⁰=0`. | Grok F7 + Codex 7 (both ⚙ the 3 distinct solves) | §5e splits into **three distinct limits** (zero-DtN, zero-affinity, impermeable) + the uniform limit; the "bulk-nonlocal part vanishes at Z→0" is the only Z→0 claim. |
| 8 | MUST | §5c one-sided corruption has **no defined pair of independent routes** — N6 requires two constructed routes first. | Codex 8 | §5c: the **two anchorings L/M are the N6 representation-invariance pair** (construct both, map Eulerian↔material, residual must vanish), **then** one-sided corruption. |
| 9 | SHOULD | The live-density lift lacks executable provenance (the real PY `rho_br_bg_rho4_constant` is a bare Symbol; the live relation is in `background_density_map`), and v1's `∇ρ` jet is the wrong object (no jet; O(η) multiplicative; ⛔ not `ρ_m`). | Codex 9 + Grok F9 | §3d.1/§5d: bind `rho_br_bg_rho4_constant` to `background_density_map` before the fold; the **two density representatives ARE the field-vs-constant pair**; emit `DENSITY_LIVE_MINUS_FROZEN`, ⛔ not a `∇ρ` jet, ⛔ not `ρ_m`. |

## What was already sound (both legs agreed — kept)
Close-then-extract is the correct ordering (both ⚙-verified extract/close non-commutation, and that extract-first
commutes only without cross-sector threading — which is not this problem); positional base+c1-delta fold with the
guards working on the two-parent fold (empty exact-key intersection, no overwrite); "no expected coefficient of a c2
output" leak-clean apart from §5e; the increment framed as an export representation, not a check; adjointness residual
withheld unless independent.

## Author note (rule 15 / decision-list length)
A high finding count on a first-draft decision list is a signal about **my** authoring — v1 carried a genuinely-wrong
central design claim (the isolation) and a fold operation undefined against the real θ-row (I read the spec's
*description* of the θ-row, not the exported ROW showing `J_s` eliminated). This is the **first** review, so the fold
is one pass by the same author (rule 7). ⚠ **Rule-15 trigger armed:** if v2 draws another heavy round at the
build-directive gate, hand the re-author to Codex (S11c-b spec round 3 precedent). [[feedback_decision_list_length_is_the_defect_rate]]
