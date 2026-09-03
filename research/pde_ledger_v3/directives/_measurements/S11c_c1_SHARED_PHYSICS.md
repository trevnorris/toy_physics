# S11c-c1 SHARED_PHYSICS spec — decision-leg review record (2 legs, rule 7 TRIGGER)

Artifact `directives/S11c_c1_SHARED_PHYSICS.md` (v1, the first unit of the 2-way split of S11c-c, authored
2026-09-03). Orchestrator-written physics-bearing spec ⇒ 2 decision legs (Codex + Grok, document review). Prompt
`directives/_legs/S11c_c1_spec_review.md` (12135 b). Logs `~/.s11_build/S11c_c1_spec_review/{codex,grok}.log`.
Both exit 0; each saved independent SymPy derivations + literal stdout:
- Codex: `directives/_legs/S11c_c1_independent_physics_derivations.{py,out}`.
- Grok: `directives/_legs/S11c_c1_grok_leg_derivations.{py,out}` + `..._grok_leg_energy_tautology.{py,out}`.

## Both verdicts: SAFE AFTER THE LISTED FOLDS — NOT A RE-AUTHOR. Both confirm the object + the split.
The core is DERIVED-CORRECT by both (with scripts): the two-momentum DtN kernel
`Z_1(k,k′)=iρ_mω·Ŵ(k−k′)·[q(k)q(k′)+k·k′−ω²/c_s0²]/(q(k)q(k′))` (depends on BOTH legs, on-shell rigid-shift
`K(k,k)=0`); face parity (outward in-plane tilt `−½∇W` on BOTH faces, even in s); Λ placement (`J`=Λ_A/Λ_V, `t`=Λ_X,
T-i has no Λ_X) + the operator inverse `[I+(Λ_A/ρ_m²)Z]^{-1}`; the Fredholm-vs-algebraic-locus boundary; the
radiation-preserving-vs-secular-global-scaling distinction (`δZ` first variation `∝w′` secular); the grazing
domain (`|q v/ω|→0` as `q→0` while `|δZ/Z|=|ωv|/(c²|q|)→∞`); the split seam (c1 exports `(δp,J,t)(V,μ_θ)`+DtN,
does not fold); N1/N8/N12/N14; the two inherited S11b caveats are legitimate setup, not leaks.

## Union of folds (all localized; some converge across legs). c1 v1→v2, ONE pass.

### §3a — the DtN operator form
- **[Grok F1 BLOCKER] Delete the "physical-space product … OR equivalently the two-momentum kernel" fork.** A
  left-quantized `a(x,k)=h(x)σ(k)` IS "a product of W_bg with a pseudodifferential factor" and still has live
  `W_bg(x)` (so it is NOT the forbidden locally-constant `Z(k;∇W_bg)`), yet it carries only `q(k)` and drops
  `q(k′)` → re-licenses the one-leg freeze (Grok energy_tautology.out: two-sided `q(k)q(k′)` vs left-quantized
  `q(k)²` differ). ⇒ require BOTH: the physical-space **composition** `N_0∘M_h∘N_0` + the `Div(h∇)` + `κ²h`
  terms, AND the two-momentum kernel with BOTH `q_out(k)`, `q_out(k′)` explicit; ⛔ forbid a left-quantized
  single-`q(k)` object; and give the rigid-shift / `Ŵ(0)` cancellation a NAMED emit tag `S11CC1_DTN_RIGID_SHIFT_*`
  (it was required in prose but had no tag → skippable).
- **[Codex F5 MODERATE] Reactive/imaginary language.** §3a forbids entrywise `Re/Im` then asks for "the sign of
  every imaginary part" + "the coefficient of `∂_t²`" — but off-diagonal kernel entries change phase under
  profile/basis translation (Codex out: entrywise `Re(A)=0` while the Hermitian quadratic form ≠0), and a
  dispersive mode-mixing impedance need not be a local `m_add·∂_t²`. ⇒ on real ω emit the Hermitian
  `H_a=(Z+Z†ₐ)/2` AND the reactive `K_a=(Z−Z†ₐ)/2i`; define per-face added mass ONLY on a named purely-reactive
  block where `p=m_add∂_t²ζ` actually exists; keep branch/sheet signs as separate `q_out` data; ⛔ remove "sign of
  every imaginary part" of individual kernel entries.

### §3b — the dissipation/energy package (the two blockers; both legs converge that it is toothless)
- **[Codex F1 + Grok F2 BLOCKER/MAJOR] The "independent energy route" cannot see a `t_s` sign error and
  wrongly cites a slab/EOM object.** `½Re(δp_s V_s*)` at the face IS the bulk Poynting flux (an acoustic
  identity for any radiation-satisfying φ with `V=v_n`) → the residual is a structural 0 that never examines `t_s`
  (Codex: `P_slab=−pV`, `P_out=+pV`; Grok: `spec pairing − bulk flux = 0`, `residual contains s_t? False`). And
  the S11b `:463-474` teeth come from pairing the **slab EOM** with velocities — c1 has no slab rows (that fold
  is c2), so citing it invites c1 to import `S11CB_SLAB_OPERATOR` and breaks the split. ⇒ face operand = the
  true-area TRACTION pairing `½Re Σ_s∫ a_s^{0} t_s·v_face,s*` (= `−½Re Σ_s∫ a_s^0 p_s V_s*` in the impermeable
  `Λ_X⁰=0` slice) from the §3b `t_s` object; bulk operand = the outgoing Poynting flux from φ on a control
  surface **in the half-space / far field**, ⛔ NOT `δp` at the face; emit both + residual; a one-sided
  traction-sign corruption must move ONLY the face operand; ⛔ drop the "`t_s` fold" language and the S11b EOM
  citation from c1 (that check is c2's, after the fold). Keep the Hermitian route.
- **[Codex F2 BLOCKER] The dissipation operator is misidentified as the BARE bulk `Z`.** `Z_diss=(Z+Z†ₐ)/2` on
  the bulk DtN audits only bulk RADIATION — the bare `Z` depends on none of the `Λ_I` (Codex: `p,J` depend on
  `Λ_A,Λ_V`, `t` on `Λ_X`, bare `Z` on none). Permeable-response dissipation (Λ_I / ωτ_I dependence) needs the
  **two-port power form** `P_out=½Re[(p+Λ_X𝒜)V*+μ_s J*]` (S11b `:705-717`) with its block Hermitian matrix on
  input `(V,μ_s)`. ⇒ emit TWO distinct real-ω objects: (a) bare bulk radiation `H_a[Z]`; (b) the power-conjugate
  closed-port response `(V,μ_s)→(p+Λ_X𝒜, J)` with its block-Hermitian form under the true-area pairing; use (b)
  for the `Λ_I`/`ωτ_I` dissipativity; specify the background measure `a_s^{0,α}`; restrict to real ω.
- **[Codex F4 MAJOR] Sign-definiteness is not inferable from a first-shape-order Hermitian operator on the flat
  evanescent nullspace.** Flat evanescent channels have zero zeroth-order Hermitian/radiative part; the O(η)
  mixing's leakage is O(η²), exactly the omitted completion (Codex: PSD completion `[[a,ηb],[ηb,η²b²/a]]`,
  first-order truncation `det=−η²b²` — an apparent negative eigenvalue that is itself O(η²)). An ordered
  propagating↔evanescent regime block is not itself Hermitian either. ⇒ assemble adjoint-related regime pairs
  BEFORE forming the Hermitian part; test the power identity through O(η); restrict sign-definiteness claims to
  subspaces where the zeroth-order Hermitian form is nondegenerate; on its nullspace emit
  `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` (the O(η²) leakage is S11c-e's).

### §5 — controls
- **[Codex F3 + Grok clean MAJOR] Remove the `Σ_E` advection mutation (§5a) and the `μ_R,bg` ablation (§5b) from
  c1 — structurally absent from c1's `(V,μ_θ)` response.** With `μ_θ` symbolic, density advection does not enter
  `μ_s` at linear wave order (Codex: `εμ₁/(ρ₀−εσδρ)=εμ₁/ρ₀+O(ε²σ)`), and `μ_R,bg` lives INSIDE the supplied
  `μ_θ` operator, not c1's response coefficients. Emitting them would be a forbidden `A−A` control OR force c1 to
  expand `μ_θ` into slab DOFs (prematurely doing c2's composition → double-count). ⇒ keep the MANDATORY `W_bg`
  tilt slope-flip on the DtN; MARK `Σ_E` advection and `μ_R,bg` structurally absent from c1; reserve their
  end-to-end mutation for c2. (If an N4 test on the response is wanted in c1, mutate `ρ_br,bg` in
  `μ_s=μ_θ/ρ_br,bg`, ⛔ not `Σ_E`.)
- **[Grok F4 MAJOR] §5a route 2 mislabels the Hanzawa construction with S11c-a's MATERIAL/N4 tags.** The route is
  correctly Hanzawa/layer-potential (not the secular global scaling), but "mapped back to Eulerian (N4:
  Δρ=δρ_E+u·∇ρ⁰)" + the `..._MATERIAL_OPERAND` tag is S11c-a's slab-density map, which does not apply to the
  BULK DtN (`ρ_m` constant; `Σ_E` absent from `Z`) → a builder implements the forbidden secular map + a spurious
  `u·∇ρ` in the bulk problem → wrong radiation branch. ⇒ rename the operands to EULERIAN vs HANZAWA/LAYER-POTENTIAL
  (drop MATERIAL); ⛔ drop the `N4:Δρ` citation from the DtN route; keep the anchoring (LAB_HELD/MATERIAL_ADVECTED
  face geometry) fixed across routes.
- **[Grok F3 MAJOR] §5d zero-jet regression target can cancel the cavity error it catches.** "the S11b flat `Z`
  at that same thickness `±W̄₀(1+η)/2`" tells the builder to re-solve B0b with gap `W̄₀(1+η)` — which IS the
  forbidden finite-gap cavity → both operands become the same cavity → residual 0 (the vacuous limit §5d exists
  to stop). Flat `Z_0=ωρ_m/q_out` is thickness-independent (Grok OBJ8; a uniform thickness change is a rigid
  shift, on-shell `K(k,k)=0`). ⇒ operand A = c1's DtN at `σ_W→0`, `η` retained, `w_1` a constant symbol; operand
  B = the UNMODIFIED S11b half-space `Z` (`S11B_Z_IMPERMEABLE`), ⛔ NO `W_0→W̄₀(1+η)` substitution, ⛔ NO two-face
  re-solve; ⛔ do not state the residual is zero.
- **[Codex F6 + Grok F1 MODERATE] §5e momentum-freeze convention.** §3a writes `Z(k,k′)` with `Ŵ(k−k′)` (k
  output, k′ input) but §5e says "replace `q(k′)` by `q(k)` on the output leg" — k′ is the INPUT leg → two engines
  corrupt differently → false comparator disagreement; and the freeze is ill-posed on a single-leg object (fixed
  by the §3a two-leg requirement). ⇒ with the two-momentum kernel now mandatory, define the convention once and
  run BOTH one-leg freezes (input-leg and output-leg) for strongest coverage.

## Convergence / divergence notes (rule 13)
- CONVERGENT (both, highest confidence): the energy route is toothless (Codex F1 / Grok F2) — both derived the
  `½Re(δpV*)`-is-the-flux identity; the DtN must carry both momentum legs (Codex confirmed / Grok sharpened the
  left-quantization hole); `Σ_E`/`μ_R,bg` structurally absent from c1 (Codex F3 / Grok clean).
- Codex-unique, orchestrator-verified: the bare-`Z`-vs-two-port dissipation misidentification (F2), the
  evanescent-nullspace sign-definiteness (F4), the reactive-part `K_a` (F5). Grok CONFIRMED the Hermitian-part
  DEFINITION (its clean axis) — compatible: definition right, applied operator wrong (Codex).
- Grok-unique, orchestrator-verified: the §5d cavity-cancellation target (F3), the §5a Hanzawa/MATERIAL/N4
  mislabel (F4), the §3a left-quantization / rigid-shift-tag hole (F1).
- Both explicitly: DO NOT launch a builder until these fold; the 2-way split is physics-sound and respected.

## Disposition
Fold all above into c1 v2 in ONE pass (rule 7: fold and go, ⛔ do not re-leg to green). The build's 2 legs +
the cross-engine comparator are the next review layer; a from-scratch WL engine that disagrees surfaces any
residual spec ambiguity. Then dual-engine build of c1 (SymPy + blind Wolfram, each 2 build legs).
