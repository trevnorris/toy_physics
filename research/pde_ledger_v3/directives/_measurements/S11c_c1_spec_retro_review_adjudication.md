# S11c-c1 SHARED_PHYSICS — retrospective spec review, adjudication (rule 2 + rule 13)

The c1 spec was folded once from its decision review and **never re-reviewed to clear** (the gap the user's
"full correct path on c1" closes). This retrospective gate re-reviews the spec's physics decisions for an error
both engines would have executed (invisible to the fidelity legs). Reports:
`_measurements/S11c_c1_spec_retro_review_{grok,codex_sol}.txt`; prompt `_legs/S11c_c1_spec_retro_review_prompt.md`.

## Commands

```
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
  "$(<directives/_legs/S11c_c1_spec_retro_review_prompt.md)"      # → codex_sol report
grok --prompt-file directives/_legs/S11c_c1_spec_retro_review_prompt.md --cwd .../research/pde_ledger_v3 \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain   # → grok report
```

## The two legs DISAGREED on the verdict (rule 6 — the disagreement is the finding)

- **Grok — CLEAR (no MUST-FIX).** Independently re-derived the two-momentum DtN (two routes agree, corruptions
  bite), the operator inverse, Λ_X-in-traction, the three distinct dissipation objects, N11a/grazing, and the
  inheritance — all sound. Leaving **seal-5 UNDECIDED is the right spec-gate call**: the spec required live
  density; **PY froze it** (an engine defect, not a spec defect). 3 SHOULD-FIX / 3 NIT, all claim/control-level.
- **Codex sol — BLOCK (3 MUST-FIX).** (1) **reopen c1** — the RHO4_CONSTANT response dropped an O(εη) density
  channel (`d(μ_s)/dη|₀ = −μ_θ·w₁/ρ_br ≠ 0`), so "δp_s/J_s agree" is invalid; regenerate the response. (2)
  grazing claimed outside the inverse's analytic domain. (3) rule-5/rule-3 leakage makes the dual-engine
  agreement non-independent evidence.

Both AGREE on the physics: the channel is real O(εη), and the **reconcile's "seal-5 = harmless O(η²)" is WRONG**.

## Verification (rule 13) — the pivot is whether PY dropped the channel IRRECOVERABLY

The disagreement reduces to one factual question: is PY's `rho_br_bg_rho4_constant` carried **opaquely**
(0 derivatives → c2 recovers the channel by re-binding it to the live field) or **frozen irrecoverably**?
Inspected the real export (this session):

```
# PY s11c_c1_face_response_coeffs, (LAB_HELD,1,RHO4_CONSTANT), PRESSURE:
#   Lambda_A_0*omega/(rho_br_bg_rho4_constant * s11cc1_q_out_output * ...)   → density is an OPAQUE 1/ρ symbol
#   exprs combining rho_br_bg_rho4_constant WITH w1_profile: 0   (0 derivatives; no η·w₁ baked in)
# background_density_map (in the fold): Eq(rho_br_bg_rho4_constant, W_bg*rho_br/W_0) = ρ_br(1+η·w₁)   [live]
# Re-binding the opaque symbol to the live field, first order in η:
python3: mu_theta/rho_br_bg_rho4_constant  --subs-->  mu_theta/(rho_br*(1+eta_bg*w1_profile))
         series(eta,0,2) = mu_theta/rho_br  -  eta_bg*mu_theta*w1_profile/rho_br
         O(eta) coefficient RECOVERED: -mu_theta*w1_profile/rho_br   (== the dropped channel)
```

⇒ PY carries the density **opaquely, 0 derivatives**; re-binding recovers **exactly** the O(εη) channel. The
channel is **NOT irrecoverably dropped** — c2 **is specified to** recover it (c2 v2 §3d.1 mandates binding
`rho_br_bg_rho4_constant` to `background_density_map`; ⚠ c2 is not yet built, so this is prospective). **Codex's
reopen severity is over-stated; Grok's disposition (no engine defect that reopens c1) is right.**
✅ **A fresh Codex-sol verify pass (this session) independently CONFIRMED the override** — re-inspected the export
(opaque `1/ρ`, 0 derivatives, no `w1_profile` baked in), re-ran the recovery, and substituted the live relation
into the *actual* exported pressure coefficient (residual 0): **NO REOPEN, override correct**
(`_measurements/S11c_c1_adjudication_verify` evidence; verify prompt `_legs/S11c_c1_adjudication_verify_prompt.md`).

## Resolution — c1's computed objects STAND (no engine reopen); record corrections OWED (3 MUST-level + carry-forward)

⭐ **c1's engines and exports are SOUND and STAND** — the two-momentum kernel + δp_s/J_s response are independently
confirmed (both legs re-derived the kernel from first principles). ⛔ **Do NOT reopen the c1 engines/exports.**

But the retro-review earned its keep — record corrections to the c1 RECORDS (⛔ not the engines). **3 MUST-level:**
1. **seal-5 (reconcile `_measurements/S11c_c1_comparator_reconcile.md` §3.5 + step-record method note):** replace
   "harmless O(η²) / 0 derivatives" with: **O(εη), carried opaquely by both engines, recoverable in c2 by
   re-binding `rho_br_bg_rho4_constant` to `background_density_map`** — a field-vs-constant *representational*
   question (still UNDECIDED, c2 re-adjudication mandatory), NOT a dropped channel.
2. **grazing (step record):** mark **exact grazing `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`** — the inverse
   `[I+(Λ_A/ρ_m²)Z]⁻¹` is nonanalytic as `q_out→0` (~1/η); the first-shape `Z₁` is a valid non-grazing asymptotic
   coefficient only, with `‖N₀⁻¹N₁‖≪1` on both legs.
3. **independence scoping (step record):** the c1 spec supplied the composition recipe + expected values (rigid
   cancellation, flat `Z₀`, zero-jet) — so part of the cross-engine "agreement" is **fidelity to the supplied
   structure**, not independent discovery. The **kernel is independently confirmed** (both legs' first-principles
   solve); scope the step record's independence claim accordingly.

**Lower-severity claim/control items to disposition or carry forward** (both retro reports; per review-until-clear,
anything changing what may be CLAIMED counts): energy-balance residual orientation/convention; `h_s`
graph-vs-outward-displacement + DtN (`N`) vs impedance (`Z`) terminology; naming the live density as a
multiplication operator `M_{1/ρ_br,bg}`; `K_a=(Z−Z†)/(2i)` is **Hermitian** (not anti-Hermitian); the evanescent
caveat should cover all second-shape grades **η²/ησ_W/σ_W²** (not only η²); the drain-projection `O(σ_W²)` wording;
the flat `Z₀=ρ_m ω/q_out` / rigid-shift "must cancel" expected-value leakage (partly under #3).

⚠ All of the above are corrections to committed c1 RECORDS (⛔ not the engines); per the (now-explicit) rule, a
physics-bearing record is review-until-clear → the corrections get their own 2 legs. c1's computed content is
untouched.
