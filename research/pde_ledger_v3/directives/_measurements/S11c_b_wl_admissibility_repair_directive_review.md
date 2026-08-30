# Decision-list review of the STEP 3 WL-repair directive (v1) — two legs, both independent CAS

Orchestrator-written directive → legs = **Codex + Grok** (rule 7). Launched 2026-08-30 on
`directives/S11c_b_wl_admissibility_repair_directive.md` (v1). Both legs derived the spec-question in SymPy and saved
script + literal stdout; the gate did **not** clear.

- Grok leg: `~/.s11_build/S11c_b_wl_repair_review_grok.txt`; derivations
  `/tmp/S11c_b_itemA_theta_body_force_derive.{py,stdout}`,
  `/tmp/S11c_b_itemB_thickness_inertia_derive.{py,stdout}`.
- Codex leg: `~/.s11_build/S11c_b_wl_repair_review_codex.txt`; derivation
  `/tmp/S11c_b_step3_independent_spec_derivation.{py,stdout}`.

## JOB-1 verdicts (the spec-question-first) — the two items go OPPOSITE ways

### ITEM B (thickness kinetic inertia) = OVER-PRODUCTION — WL is correct; DROP the repair
Both legs, independent CAS, and the source-of-truth definition agree:
- §1a L61: `e_W ≡ δW/W₀` (constant `W₀`). ⇒ `δW = W₀·e_W` identically.
- §1c `T = ½ μ_W (∂_t δW)²` ⇒ in the named `e_W` coordinate `T = ½ μ_W W₀² (∂_t e_W)²`.
- Codex stdout (literal): `coefficient of d_t^2 e_W: W_0**2*mu_W`; `derivative … w.r.t local W_bg: 0`;
  assert `inertia_coefficient == mu_W * W_0**2` PASSED.
- Grok stdout (literal): `inertia_coeff_of_eW_tt = W0**2*mu_W`; `inertia_coeff_depends_on_W_bg = False`;
  STEP 6: `δW := W_bg·e_W` conflates `e_W` with `e_W,bg` and violates the exact map `e_W,bg=(W₀/W_bg)e_W`
  unless `W_bg≡W₀`.
- §2a scope: the promotion `W₀→W_bg(y)` is scoped to explicit `W₀` factors in `U` and the listed bindings;
  `T` is not `U`, and `δW=W₀e_W` is the §1a field definition, not a listed binding; `μ_W` is a modulus that
  does not vary in-plane.

⇒ WL `kineticEw = muW WZero^2 D[eWField,{time,2}]` (L838/L923) is the faithful emit. The sibling's `W_bg²`
is an over-production (or an `e_W`/`e_W,bg` normalization conflation). **Do not repair.** This REVERSES the
§8 KINETIC verdict (see the correction folded into `S11c_b_step2_adjudication.md`).

### ITEM A (admissibility θ body force) = GENUINE UNDER-RETENTION — but the mechanism was mis-diagnosed
Both legs, independent CAS:
- The §3a-independent mixed invariant `∇θ·∇e_W`, with its thickness-gradient factor lifted to the §3d full
  width `∇(W_bg+δW)`, produces a nonzero θ Euler covector at 𝔅⁰.
- Codex stdout (literal): `theta Euler covector at first background-jet order: -K_thetaW*Derivative(W_bg(y),
  (y,2))/W_bg(y)`; nonzero in both density representatives (RHO4 and RHOBR specializations identical);
  the analogous e_W covector `-W_0*kappa_W*W_bg^2*W_bg''` is the (already-agreed) E_W body force.
- Grok: same object; coefficient `κ_θW (−W_bg W_bg'' + (W_bg')²)/W_bg²`, nonzero on a test profile.
- This independent derivation reproduces the sibling's θ body force structure (`κ_θ_W·∇²W_bg`) ⇒ the
  sibling's Item-A content is spec-correct and WL under-retains it (rule-1 cross-check: independent route
  matches PY here, opposite to Item B).
- Both legs: the mandated mechanism is the **full-field thickness-gradient lift of the mixed invariant**, NOT
  the `ρ_4D=ρ_bg⁰(1+θ)` density-weighting the v1 directive diagnosed. §1c/§3a do not name `ρ_4D` as an energy
  multiplier; density-weighting is an extra, unmandated mechanism. (Codex noted it *also* yields θ-linear
  terms, but it is not the named path.)
- Confirmed by orchestrator code read (rule 13): `constructFullFieldBackgroundEnergy` lifts the pure-thickness
  invariant [[7]] to `fullWidth` (L554) but leaves the mixed invariant [[8]] `Dot[gradTheta, gradLocalEw]`
  (L555) on the wave-only `gradLocalEw` (L543) — the asymmetry that zeroes the θ body force.

## JOB-2 directive defects (both legs), folded into v2
1. **Rule-5 leak (negation-bound + current values).** "Do not freeze the physical thickness at `W₀`" bounds
   the withheld object (and is now known-wrong); stating `Integer(0)` and `muW WZero²` and "insert the missing
   θ channel" leak values/direction. → v2 removes all current values, directional prose, and the (dropped)
   Item B.
2. **Rule-3 over-specification / wrong mechanism.** v1 pinned Item A on a `ρ_4D` density-multiplier recipe.
   → v2 names the OBJECT (the §3d full-field first variation with full-field thickness-gradient content
   `∇(W_bg+δW)` in every invariant carrying a thickness gradient) and does not prescribe a construction.
3. **Wrong WL line cites.** `fullWidth` is L541 (not L565), the full-field gradient invariant [[7]] L554
   (not L560), `thetaVariation` L538 (not L562); the mixed invariant [[8]] is L555, `gradLocalEw` L543,
   `localEw` L540. L1334–1340 and L838/L923 were accurate. → v2 corrects the cites.
4. **Leak-gate claim was incomplete.** The v1 grep gate missed the negation-bounds and current-value
   citations. → replaced by semantic review; v2 re-gated.
5. **Blindness clean** (both legs): no sibling import/copy requested. Preserved in v2.

## Net
- v2 directive = **Item A only** (admissibility θ full-field mixed-invariant lift), leaks/over-spec/cites
  fixed. Re-run 2 decision-list legs (Codex + Grok) on v2 before any builder.
- Item B **withdrawn**; the §8 KINETIC verdict corrected to: WL correct, sibling over-produces `W_bg²`
  (re-adjudicate as representational `e_W`/`e_W,bg` conflation vs a genuine sibling bug at step 4).
