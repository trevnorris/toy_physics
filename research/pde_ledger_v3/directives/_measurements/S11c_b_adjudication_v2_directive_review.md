# S11c-b adjudication v2 directive — decision-list review, round 1 (two legs, NOT SOUND → folded)

**Legs (orchestrator-written ⇒ Codex + Grok):** prompt `_legs/S11c_b_adjudication_v2_directive_review.md`.
Raw: Codex `~/.s11_build/S11c_b_adjud_v2_directive_codex.txt`; Grok `~/.s11_build/S11c_b_adjud_v2_directive_grok.txt`.
**Both legs: NOT SOUND.** The v2 design (Bridge D + generic divergence classification) had two load-bearing
errors — a wrong jet expansion and an over-broad divergence classifier — both of which would manufacture false
agreement. Caught before any build.

## Orchestrator verification (rule 13) — Bridge D jet map (Grok B1 / Codex f1)
My directive wrote the naive continuum chain rule `W_bg_d{i} ↦ W_0·eta_bg·w1_profile_d{i}`. WRONG. PY's committed
`PROFILE_GRADE_SUBS` (`scripts/S11c_b_brane_operator_sympy_audit.py:648-668,1862-1868`) is:
```
W_bg   -> W0*(1 + eta_bg*w1_profile)          mu_R_bg -> mu_R*(1 + eta_bg*m1_profile)
grad_W[i]  (W_bg_d{i})   -> sigma_W*w1_grad[i]                    # first jet: sigma_W, NOT W_0*eta_bg
grad_W[i]  (2nd)         -> sigma_W*w1_profile_second/L_W         # second jet: /L_W
grad_mu[i] (mu_R_bg_d{i})-> mu_R*sigma_W*m1_grad[i]/W0
grad_mu[i] (2nd)         -> mu_R*sigma_W*m1_profile_second/(W0*L_W)
rho4_rho4 -> rho_br/W0 ;  rhobr_rho4 -> rho_br*(1 + eta_bg*w1_profile)
```
`sigma_W` (DERIVED "first-background-jet bookkeeper", L156) and `eta_bg` (KNOB "zero-jet contrast bookkeeper",
L157) are INDEPENDENT symbols. The zero-jet `W_bg = W0(1+eta_bg·w1)` is correct; the JET map is sigma_W/L_W.
⇒ **Bridge D must REUSE the engine's own `PROFILE_GRADE_SUBS`** (covers W_bg AND mu_R_bg AND jets AND density
bookkeepers — resolves Grok B2), applied to the renamed operands. Do not hand-write a chain rule.

## Blocking findings (folded into v2-corrected)
- **B1/f1 Bridge D jet map wrong** → reuse `PROFILE_GRADE_SUBS` (above); forbid identifying `sigma_W ≡ W_0·eta_bg`.
- **B2 mu_R_bg omitted** → `PROFILE_GRADE_SUBS` includes it; fold the full map.
- **Codex f2 — the divergence classifier is UNSOUND for STRONG operators** (the deepest finding). A strong
  operator `O=∇·V` has `∫ψO=−∫∇ψ·V≠0`; only a formed density difference `∇·(ψV)` is boundary-equivalent.
  SLAB_OPERATOR, SLAB_OPERATOR_TERM_ORIGINS, MU_THETA_OPERATOR, ADMISSIBILITY_OPERATOR_OPERAND, kinetic/container
  are STRONG (WL `variationalSource = ∂L/∂q − Inactive[Div](...)`, WL:274-277,793-795); dropping their HeldDiv
  (comparator:547-552) erases bulk derivatives. ⇒ the divergence classification applies ONLY to weak-pairing
  DENSITIES — currently only `COUPLING_KERNEL` (WL:1075-1124). Strong operators are compared EXACTLY (Bridge A +
  Bridge D only); their nonzero survivors are genuine differences.
- **Codex f3 — V reconstruction mandatory + verified.** `REPRESENTATIONAL_DIVERGENCE` requires raw residual R +
  explicit 3-component V + printed verification `R − Σ∂_iV_i = 0` in the jet algebra. An Euler-signature zero
  only screens; it is not the verdict and not "the bulk part". Unsupported → `DIVERGENCE_INCOMPLETE`/op-error.
- **Codex f4 — fixtures too weak.** Add the variable-coefficient product-rule pair: `a·φ_,1 → RESIDUAL_BULK`;
  `a_,1·φ + a·φ_,1 = D_1(a·φ) → REPRESENTATIONAL_DIVERGENCE, V=(a·φ,0,0)`; a route fixture that a divergence in a
  strong-operator family is NOT eligible; a fixture-local field registry (the comparator base registry 494-525
  doesn't know placeholder f,g,h).
- **Codex f5 — protected reps hide inside operator/coupling residuals.** Gate on exact protected ATOMS
  (`gamma_s11cb_{w_bg,mu_r_bg}_{07,10}`, `gamma{Width,Modulus}DivGrad{Theta,Ew}`), not just the ENERGY_BASIS
  family → a `PROTECTED_UNREDUCED` route keeps them raw; Bridge D also barred from protected families.
- **B4/f6 rule-5 leak** (MY recurring defect, again): the preamble "decompose into missed-bridge/protected/
  divergence" + "AGREE on ADMISSIBILITY (38)" leaks the expected taxonomy AND misstates (38 is global MATCH; the
  4 ADMISSIBILITY BODY_FORCE/THETA are FLAG — PY nonzero, WL 0 — the candidate finding). Strip to the neutral
  measured fact (v1: 38 MATCH, 40 FLAG).

## Cleared
Case-multiset accounting + drop protection adequate; `--drop-divergence`/`--drop-bridge-d` are genuine ablations
after the above; fixture placeholders are placeholders (good).

## Disposition
Fold into v2-corrected: Bridge D = `PROFILE_GRADE_SUBS`; strong operators compared EXACTLY (divergence ONLY on
COUPLING_KERNEL density, with mandatory verified V + protected-atom gating); strengthened fixtures; strip the
rule-5 preamble. Re-leg. If round 2 breeds new blockers, delegate the directive to Codex (rule 15).

---

# Round 2 — v2 directive review (two legs, NOT SOUND → final fold; build legs gate)

Raw: Codex `~/.s11_build/S11c_b_adjud_v2_directive_r2_codex.txt`; Grok `~/.s11_build/S11c_b_adjud_v2_directive_r2_grok.txt`.
**Both legs: round-1 fixes VERIFIED clean** (Bridge D=PROFILE_GRADE_SUBS, strong/weak split, V-mandatory,
rule-5). Remaining blockers are subtle physics; folded into the final directive, then build legs gate (2 leg
rounds done — rule 7 fold-and-go; the subtleties are not cleanly pre-enumerable in prose,
[[feedback_delegate_build_when_prose_directive_repeats]]).

## Blocking findings (folded)
- **Codex1/GrokN2 — import the FULL `PROFILE_GRADE_SUBS`** (all four density bookkeepers `rho4_rho4,
  rhobr_rho4, rho4_rhobr, rhobr_rhobr`, PY:667; NOT `e_W_bg`); build check compares the complete key/value set.
- **Codex2 — consume RAW `case.value`, not `case.compare_value`** (adjudication:524): the latter is replaced by
  `modulo_total_divergence(...)` (comparator:617) which DROPS every HeldDiv (547); adjointness is marked
  `reduce_divergence=True` (comparator:872,909). v2 must use raw values, disable the legacy prepass, and assert
  NO `DIVERGENCE_REDUCED` context reaches the classifier.
- **Codex3 — Bridge D and ∂ DON'T COMMUTE** (sigma_W ⟂ W_0·eta_bg): `B_D(∂_1(W_bg·φ)) = sigma_W·w1_d1·φ + …`
  but `∂_1(B_D(W_bg·φ)) = W_0·eta·w1_d1·φ + …` — DIFFERENT. ⇒ expand HeldDiv/divergences in the ANCHORED
  (pre-Bridge-D) jet algebra with the engine's formal derivative; verify a weak witness by
  `BridgeD(R_preD − Σ_i ∂_i V_i) == 0`, NEVER `∂_i(BridgeD(V_i))`; apply Bridge D LAST. Add an anchored-profile
  non-commutation fixture (the generic a,φ fixtures miss it).
- **Codex4 — eligibility = source-proven scalar-DENSITY path, not `family==COUPLING_KERNEL`.** Exclude
  `ADJOINTNESS_RELATION` (relational). `COUPLING_KERNEL_TERM_ORIGINS` are ALSO weak densities (WL extractCoupling
  per origin WL:1193; PY weak-block per origin PY:1996) — either total-bijection-adapt those origin-density
  leaves (eligible) or leave STRUCTURE_INCOMPLETE/COVERAGE; NEVER exact-FLAG them as strong. Add
  `ADMISSIBILITY_SUPPORT_OPERAND` to the exact non-density list.
- **Codex5 — `ENERGY_BASIS_COUNT` compared EXACTLY** (it supplies 2 of the 38 MATCH); VARIABLE/NEW_INVARIANTS/
  OMISSIONS structurally protected; none enter Bridge D/divergence.
- **GrokB1 — explicit route ORDER:** renames + Bridge A + Bridge D → if residual free-syms ∩ protected atoms →
  `PROTECTED_UNREDUCED` (no divergence; Bridge D already applied, safe since PROFILE_GRADE_SUBS doesn't touch
  protected symbols) → else COUPLING scalar density → V-classify → else strong exact `MATCH`/`FLAG`.
- **GrokB2 — HeldDiv expand = `HeldDiv(V) ↦ Σ_i formal_∂_i(V_i)`** (the anchored formal derivative), NEVER
  `_drop_held_divergences`/`reduce_divergence=True` on strong families.
