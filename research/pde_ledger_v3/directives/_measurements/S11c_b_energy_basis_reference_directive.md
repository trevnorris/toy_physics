# Measurements — backing the #86 energy-basis REFERENCE directive (rule 2)

Every factual claim the directive makes about an artifact, with the command that produced it. The **corrected**
quotient dimension is the withheld computed result (rule 5) and is NOT stated here; it is recorded separately
in the #86 result record after the reference build. Only the FROZEN (public) values appear below.

## 1 · Spec authority quotes (§1a / §1d / §3a) are verbatim
```
$ sed -n '150,171p' directives/S11c_b_SHARED_PHYSICS.md      # §1d
  → "the uniform quotient does not lift trivially … c∇·F ≡ −(∇c)·F modulo a boundary term … a first-jet
     invariant that is physics in the operator/kernel, not a representational identity."
$ sed -n '242,271p' directives/S11c_b_SHARED_PHYSICS.md      # §3a
  → "together with the supplied background first jets {∂W_bg, ∂μ_R,bg} treated as symmetry-breaking spurion
     data … a second spatial derivative of W_bg is still first order in background amplitude … and is not
     dropped." ; "impose the exact map e_W,bg=(W_0/W_bg)e_W before any independence/rank test."
$ sed -n '44,92p'   directives/S11c_b_SHARED_PHYSICS.md      # §1a
  → "the non-uniform background breaks in-plane translation invariance, so u (either part) may enter
     undifferentiated when contracted with a background gradient — these are the N15 spurion couplings §3a
     constructs"; "whether u enters undifferentiated, is computed in §3a/§3c, not stated here."
```

## 2 · The frozen-spurion mechanism is code-verified in the engine
```
$ sed -n '936,970p' scripts/S11c_b_brane_operator_sympy_audit.py    # basis_euler_signatures
  → derivative_maps built by `for field, first, second in fields:` only; basis_dx differentiates atoms in that
    map. The background spurion is NOT among `fields`.
$ sed -n '1025,1029p' scripts/S11c_b_brane_operator_sympy_audit.py  # basis_fields
  → basis_fields = the DOF fields only: (bu[a],bG[a],bU2), (btheta,bq,bTheta2), (be,br,bE2). No W_bg/∂W_bg row.
$ sed -n '1270,1273p' scripts/S11c_b_brane_operator_sympy_audit.py  # §3a call site
  → candidates = enumerate_new_candidates(g_vector); signatures = basis_euler_signatures(expressions,
    basis_fields).  g_vector (=∂W_bg) is in the candidates but absent from basis_fields ⇒ basis_dx treats it
    as constant ⇒ the ∂²W_bg term §3a retains is never generated. This is the frozen (uniform) quotient §1d
    forbids at variable coefficients.
```

## 3 · The frozen quotient is UNDERCOMPLETE — the merges it makes are spurious (own CAS, rule 13)
`~/.s11_build/S11c_b_86_probes/orchestrator_crux_check.py` (own construction; imports nothing). Codex's decision
leg emitted the frozen null-relation `c01 + c02 ≡ 0 (mod div)`. With W_bg a real field (`∂_iW=g_i`,
`∂_i g_j=H_ij` symmetric), the check computes the total-divergence current `J_j = (g·u) u_j`:
```
c01 = (∂_k u_k)(g·u) ,  c02 = g_i (∂_j u_i) u_j
∂_j J_j − (c01+c02)  =  H_ij u_i u_j        (residual vs H·uu = 0)
⇒ c01 + c02 ≡ −H_ij u_i u_j  (mod in-plane divergence)          # NONZERO with the Hessian retained
frozen view (H→0): c01 + c02 ≡ 0                                 # the spurious merge
```
So a pair the frozen quotient merges is genuinely independent once the background carries its Hessian; the
dropped invariant carries real physics (`H·uu`). Direction: the frozen §3a basis is undercomplete.

## 4 · The FROZEN basis count (the public regression anchor; NOT the corrected answer)
```
$ grep -a 'S11CB_ENERGY_BASIS_COUNT' scripts/out/S11c_b_brane_operator_sympy_audit.out
  → PY_S11CB_ENERGY_BASIS_COUNT: … Integer(26) …   (both anchorings LAB_HELD, MATERIAL_ADVECTED)
```
Control 1 of the reference regresses against this emitted frozen selection only.

## 5 · The two decision legs' machinery is validated against ground truth
Durable probes: `~/.s11_build/S11c_b_86_probes/{codex_leg_probes,grok_leg_probes}/`. Both independently
reproduce the engine's emitted FROZEN per-source selection `[1,4,5,6,7,9,10,13]` (rank 8) — i.e. each leg's
reconstruction matches the committed engine output, validating its quotient machinery before its corrected
computation. Leg reports: `~/.s11_build/S11c_b_86_reference_directive_review_{codex,grok}.txt`. Six findings
folded into rev 2 of the directive (finite seed; grade projection; drop the false family-separation control;
strip the acceptance-outcome leaks; add μ_R Hessian + map-ordering diagnostics + divergence certificates; emit
rank+span+certificates, not a non-invariant index selection).
