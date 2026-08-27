# Measurements — S11c-a PY projection current fix directive

Grounding for the directive's factual claims (rule 2). Repo root `/var/projects/toy_physics`.

## PY code facts (the defect)
```
$ sed -n '176,181p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
176: j_bulk = (delta_j_bulk_1 … delta_j_bulk_4)                      # bare face-value symbols
177: dw_delta_j_bulk = { s: (d_w_delta_j_bulk_{plus,minus}_{i}) }    # NORMAL jets (declared)
181: grad_j_bulk = ((delta_j_bulk_{i}_d{j}))                          # in-plane jets
```
- `projection_terms` (lines 1114–1171) builds the projection from `j_bulk[i]` (constant),
  `current_divergence = Σ_i grad_j_bulk[i][i]` (in-plane jets only), and
  `WINDOW_NORMAL_CURRENT = -ε·∫ j_bulk[3]·window_normal dw` (δj_w a `w`-constant). It does **not** reference
  `dw_delta_j_bulk`.
- The normal jets ARE used in the *trace* construction:
  ```
  $ sed -n '650,652p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
  current_perturbation = tuple(affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face) for i in range(4))
  ```
  The trace's affine mechanism is FACE-keyed (`affine_bulk_perturbation(...,face)`, expands about ±W₀/2)
  and the projection is FACE-LESS — so rev-1's "reuse line 651" recipe was WRONG (2-leg review, §21). rev-2
  names the OBJECT only (∂_wδj_w must enter the projected ∇₄·δj), forbids the face-trace recipe, forbids a
  double-count against the existing post-IBP WINDOW_NORMAL_CURRENT, and scopes both branches + the §5c twin.

## Why this is the fix (finding provenance — NOT a build target; the builder never sees this file)
- CAS-verify (both engines, runnable SymPy + stdout, `~/.s11_build/proj_cas_*`): freezing the current at the
  face vs. keeping its `w`-profile gives a genuinely different projection (Q1 = DIFFERENT, witnesses
  `40000*log(10)/9999` and `sqrt(2π)·W0·c·(-A+B)·exp(-W0²/2)/2`); the window shape-derivative form is benign
  (Q2 = SAME).
- From-spec adjudication (Codex + Grok, `~/.s11_build/proj_wdep_adj_*`): unanimous — §1b requires the full
  `w`-dependent bulk current (WL faithful, PY defect), and it is an implementation defect, not a spec
  ambiguity. Orchestrator rule-13 self-verified against §1b/§3c and the raw engine source.
- Corroboration (§20): PY's UNIFORM_LIMIT_RESIDUAL is nonzero only in PROJECTION dynamic objects
  (`delta_j_bulk_i_di` + ∫O_window survive), while WL's is 72/72 zero — PY's frozen-current projection fails
  the §5c reduce-to-S11b regression that WL passes.

## §1b spec text (the governing law the fix implements)
```
j = ρ_4D v_bulk ,   ∂_tρ_4D + ∇₄·j = 0 .
Ω is a smooth slab window … integrating this conservation law against Ω and integrating by parts in w.
```
`∇₄·j = ∂_i j_i + ∂_w j_w` includes the normal-derivative term. §3c's zero-background/trace language is scoped
to TRACED fields at the face `h_s⁰`, not the bulk current under the `w`-integral.

## Rule-5 note
The directive states the OBLIGATION (integrate the full ∇₄·δj including ∂_wδj_w) and explicitly withholds
which normal-jet terms survive ("computed by the engine, not asserted here"). No expected residual value, no
"reduce-to-S11b", and no cross-engine match target is supplied to the builder.
