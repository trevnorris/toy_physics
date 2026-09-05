# Independent physics review — S11c-c2 SHARED_PHYSICS **v2**, DECISION-LIST gate (round 2)

You are one of two independent legs re-reviewing an **orchestrator-written physics spec** *before any builder is
launched* (CLAUDE.md rule 7). The spec settles the structure, physics, provenance, and scope that the S11c-c2 build
directive and its two blind engines (SymPy + blind Wolfram) will honour. A defect here makes **both** engines agree
on the same wrong fold, so it is the one artifact the builder trusts. This is **v2** — folded once from a prior
two-leg gate. Your job: (a) verify each prior finding is **genuinely fixed and bites** — not reworded or relocated;
(b) independently hunt for **new** defects, especially any the fold introduced; (c) re-settle the load-bearing
computational claims on v2's actual text. ⛔ Do not rubber-stamp; ⛔ do not soften to reach agreement.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`  (this is v2)

## What S11c-c2 is
c1 exported the **closed permeable face response** `(δp_s,J_s,t_s)(V_s,μ_θ)`, the nonlocal DtN operator, its flat
symbol, and its two-momentum kernel. c2 **folds** the closed response into the S11c-b variable-coefficient slab
operator `slab_operator` (whose θ-row and mechanical rows still carry `δp_s` **symbolically**) and **re-extracts**
the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the **closed** full operator → the nonlocal **self-energy**.
CODE build; this is the pre-build decision gate, round 2.

## What you are handed (read the SOURCES OF TRUTH first, form your own view, THEN read v2)
Read in this order (a method instruction, not a blindness control):
1. Parent decision list `directives/S11c_decisions.md` (N-series scope).
2. c1 step record `steps/S11c_c1_curved_bulk_closure.md` + reconcile `directives/_measurements/S11c_c1_comparator_reconcile.md`.
3. S11c-b step record `steps/S11c_b_variable_coefficient_operator.md` (the DEFERRED cross-engine residual; the two
   whole-row sign conventions + #90's two flags that are cross-engine-UNVALIDATED, lines 112–115).
4. c1 spec `directives/S11c_c1_SHARED_PHYSICS.md` §1d (Λ-channel placement; "the routing … is S11c-c2's") + §3b
   (operator inverse; the energy note that the traction-vs-slab pairing "is S11c-c2's, after the fold").
5. S11c-b spec `directives/S11c_b_SHARED_PHYSICS.md` §3b/§3c; the gated bind-closure design
   `directives/export_ledger_bind_closure_design.md` (esp. lines 148–153 on `coupling_kernel` and 186–188 on c1 keys).
6. The two real export files `scripts/S11c_b_exports.py`, `scripts/S11c_c1_exports.py`, and `scripts/ledger_fold.py`.
   ⭐ You may fold them in Python and inspect the **real** `slab_operator` / `closure_shape_deriv` /
   `background_density_map` / `s11c_c1_face_response` rows (`load_model(base, c1_delta)`).
7. **The prior gate record** `directives/_measurements/S11c_c2_spec_review.md` (the 9 folded findings) and the prior
   leg reports `directives/_measurements/S11c_c2_spec_review_{grok,codex_sol}.txt`.
8. **Then** v2 of the spec.

## Required method (DOCUMENT branch + computation where a claim is computational)
Form your own independent view from the sources first, then read v2. Quote both sides for every finding.

⛔ **A prose derivation is worth nothing.** Where a claim is computational, **write your own SymPy script and save
both the script and its literal stdout to named absolute paths** (under `/tmp` or the `_legs/` dir); without them
your derivation claims are discarded. Do NOT spawn Mathematica kernels. Copy anything you run to `/tmp`; never modify
the working tree.

## A. Verify each prior finding is genuinely fixed AND bites (not reworded/relocated)
For each, quote the v2 text that resolves it and confirm the fix is *sound*, or flag it as unfixed/half-fixed:
1. **Isolation claim (the central one).** v2 §3c now says the increment drops only the bulk/kinetic base and that the
   face-force + #90 closure-fold sign conventions **carry through** and are surfaced/adjudicated. ⚙ Re-derive: does
   the increment really cancel the bulk base but retain the two carrier conventions? Is the new **§3d.4 traction-vs-
   slab mechanical-power pairing** actually able to adjudicate the face-force sign (i.e. does a one-sided `t_s`-sign
   corruption move that residual)? Or is §3d.4 a tautology / non-independent route?
2. **The fold operation.** v2 §1c/§3a now substitute **closed `δp_s` + w-jets** into `delta_p_±`/`d_w_delta_p_±`
   slots (not a closed `J_s`), and put `Λ_X` in traction only. ⚙ Check against the **real** `closure_shape_deriv`
   row and c1 §1d: is that the correct, complete, non-double-counting fold? Are the pressure **jets**
   (`d_w_delta_p_±`) and c1's `V_s`/opaque-`μ_θ` symbols all accounted for?
3. **Two-face assembly.** v2 §3a emits the assembled two-face operator per `(α,ρ)`. Correct against the `Σ_s` in the
   mass/virtual-work laws?
4. **Six re-adjudications** (density, `t_s`, DtN kernel/whole-form, traction-slab pairing, flat-symbol usage,
   `μ_R,bg` form). Is the set now complete, or does the fold make yet another c1-UNDECIDED/DEFERRED item load-bearing
   that v2 still omits (e.g. ENERGY, the giants)?
5. **Three named operators** (closed-coupling-kernel / ordering-commutator / self-energy-increment) — unambiguous now?
6. **Export wiring** — positional `load_model`; consume-set uses `s11c_c1_*` write-keys (not bare `face_response`);
   `coupling_kernel` regression-only, not a §3c construction operand. ⚙ Fold the real files and confirm the intended
   consume-set binds the **closed** c1 rows, and that the guard's existence-only pass is correctly flagged.
7. **`Z→0 ≠ Λ_A⁰=0 ≠ impermeable`** — v2 §5e split into three limits. ⚙ Confirm the three solves differ.
8. **N6 two-route control** — v2 §5c uses the two anchorings as the representation-invariance pair, then one-sided
   corruption. Is that a valid N6 independent-route construction, or still under-defined?
9. **Density provenance** — v2 binds `rho_br_bg_rho4_constant` to `background_density_map`; the two density
   representatives are the field-vs-constant pair; no `∇ρ` jet; not `ρ_m`. ⚙ Confirm the map really carries the live
   relation and the two representatives differ as claimed.

## B. Independent fresh hunt (rule 15 — did the fold introduce anything?)
Beyond the prior findings: any NEW under-specification, manufactured recipe (rule 3), tautological residual (rule 2
corollary 3), leaked expected value (rule 5), rule-17 freeze, F9 false-equal, or scope error introduced by the v2
edits. In particular stress: does the new §3d.4 pairing or §5c rep-invariance residual smuggle a tautology? Does the
"two-face assembled operator" decision conflict with the per-`(α,ρ,s)` provenance anywhere? Is the symbol-map
deferral to the build directive safe, or does it hide a load-bearing decision?

## Physics filter
Report a finding only if it catches a way the c2 **decisions** would make the physics wrong, under-specified, or
un-reviewable. Rank most-load-bearing first; mark each MUST-FIX / SHOULD-FIX / NIT with the exact v2 line and the
source line it contradicts.

## Output
Part A: per prior finding, FIXED / HALF-FIXED / UNFIXED + the quote or computation. Part B: new findings, ranked.
End with an explicit verdict: **ALL-CLEAR (no MUST-FIX remaining, sound to build against)** or **BLOCK** with the
remaining must-fix list. Two legs finding real defects is the expected outcome — but if v2 is genuinely sound, say so
plainly.
