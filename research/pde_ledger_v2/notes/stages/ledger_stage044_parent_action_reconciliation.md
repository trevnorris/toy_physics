# ledger_stage044_parent_action_reconciliation

## Status

**Part VII — Integration. VII-2a (build-order 044; the SPINE stage) — the
conservative parent action ASSEMBLED at the completeness floor, with the
`r_B`↔`χ_B` wall reconciliation, `PARENT_ACTION_ASSEMBLED_AT_COMPLETENESS_FLOOR`.**
This is a genuinely-NEW synthesis, not a reshape: it assembles Parts I–VI
(medium/gravity, wall, light, charge, magnetism) into ONE **conservative parent
action** over ONE field set, RECONCILES the two historical wall descriptions
(G0-card `r_B` ↔ Part-I `χ_B`) at the STATIC WALL-ENERGY level, and NAMES the two
candidate (non-variational) drain interfaces — **at the structure /
completeness-floor level only.** The card uses **`\StatusOpen`**, NOT
`\StatusExactClosure`. Scope class: **SYNTHESIS**, dual-engine.

**⭐ KEY FRAMING — carried verbatim in intent, NEVER softened, NEVER inflated.**
044 assembles the STRUCTURE of the conservative action and verifies its internal
consistency; it does **NOT** close the action, does **NOT** select or equate the
drain, and does **NOT** perform any two-body / far-field / throat solve. Honest
framing:

- `closed_action = false` — the closed unified action awaits the shared throat
  solve (`throat_solve = SIM_DEFERRED_CENTRAL_DEBT`).
- `g0_card = DRAFT_V0` — the G0 closure card
  (`software/em_charge_attribute/g0_closure_card_v0.md`) is DRAFT-v0 / postulated,
  not a ratified stage; many gates are DEFERRED.
- `drain_equivalence = UNRESOLVED`, `drain_selection = DEFERRED_045_AND_USER_GATE`
  — **044 makes NO drain decision.** The genuine drain-placement mini-gate is the
  USER's, at 045.
- `Z_χ` is an **unearned inertial postulate** (`[POSTULATE]`, no reduction route) —
  the honest inertial-completion cost of building a conservative wall from Part-I's
  dissipative one; NOT a "no new knob" stage.

The verdict is NOT "one medium proven" and NOT "second substance found": it is the
assembly of a conservative parent-action STRUCTURE that is internally consistent,
with the drain equivalence UNRESOLVED and the closure SIM-deferred.

## The assembled two-object action

Two NESTED objects, kept DISTINCT (never conflated):

- **`S_cons^G0`** = the G0 card's six conservative summands (minus-sign Legendre
  convention):
  ```text
  S_cons^G0 = S_bulk + S_χ + S_scalar + S_mouth + S_hold + S_geon,const
  ```
- **`S_assembled`** = the FULL parent action over all sectors:
  ```text
  S_assembled = S_cons^G0 + S_brane[light] + S_move[magnetism moving coupling]
  ```
  where `S_move`'s ONLY new conservative term is counted **ONCE**, and `S_scalar`
  is EITHER the parent triple `S_H + S_u + S_mix` OR the reduced `S_Lh`, **never
  both** (`scalar_branch = REDUCED_S_Lh | PARENT_S_H_plus_S_u_plus_S_mix`,
  exclusive).

Summand roster (verdict `S_cons_G0_summands` + `S_assembled_extra`):

- `S_bulk` — medium + gravity-flow (Madelung; stages 002/004; fields `ρ` [`L⁻⁴`],
  `θ` [dimensionless Madelung FLOW phase, `v=(ℏ/m)∇₄θ`]); `U=(K/4)ρ⁵`, `P=Kρ⁵`.
- `S_χ` — the ONE wall, INERTIAL/conservative form (stages 006/007):
  `½Z_χ(∂_t r_B)² − ½κ_χ|∇₄r_B|² − (λ_χ/4)r_B²(1−r_B)²`, `λ_χ=2κ_χ/ℓ²`, `r_B∈[0,1]`.
- `S_scalar` — the reduced `(u_L,h)` charge block (stage030) `S_Lh`, XOR its parent
  triple `S_H+S_u+S_mix` (`S_mix=−C_hu∇u_L·∇h` is INSIDE the block).
- `S_mouth` — the fixed-source `χ↔h` Robin sheet (stage031); the odd charge kernel
  `Q_χ[r_Σ,s_i]=s_i` is a FIXED normalized source/orientation relation. What is
  R63-deferred is the physical electric BC / pair-force SIGN, NOT the `s_i=±1`
  label.
- `S_hold` — the holonomic wall-freeze multiplier (stage032), CONSERVATIVE, IN
  `S_cons^G0`: `−∫dt∫_{Γ_Σ} d³A λ_Σ(r_B−½)`, pinning `r_B=½` on the mid-surface.
  `λ_Σ` is a **non-propagating auxiliary Lagrange multiplier** (in the variational
  inventory, EXCLUDED from the physical DOF count). ⚠ Notation firewall: write the
  ambient background piece `Γ_Σ^amb` — the G0 glyph COLLIDES with Part-I's `Γ_B`
  (the `T⁻¹` order-conversion RATE); DISTINCT objects.
- `S_geon,const` — additive separation-independent rest constant (no field
  content; in G0-v0 every field-dependent geon coupling is ZERO).
- `S_brane[light]` — the transverse light seam (stages 003/007): the `g_ℓ`-LOCALIZED
  `∫dt d⁴X g_ℓ(w)[L_Mac+L_uw]`, MacCullagh transverse shear `Ω_u=∇_∥×u` + the gapped
  `u_w` DOF; light speed `c_γ²=μ_R/ρ_br`; two transverse photons. `u^a` (light
  transverse) ≡ `u_T` (magnetism), projected divergence-free (`∇·u_T=0`) so its
  longitudinal inertia is NOT double-counted with `u_L`.
- `S_move[magnetism moving coupling]` — the ONE new conservative magnetism term
  (stages 034–039): `∫dt d³x q_T Σ_i s_i η_a V_i·u_T` (`∇·u_T=0`, `q_T=λ_T τ_d`),
  counted ONCE; the kinetic/gradient rows `½ρ_br|u̇_T|²−½μ_R|∇×u_T|²` are the
  IMPORTED `L_Mac` (cited, not re-counted). `BOOST_STRUCTURAL_RELATION_HOLDS`
  (magnetism structurally IS the boost of electricity; only `r_BA=q_T²/(ρ_br A_E)`
  open); `b_T=∇×u_T` is T-EVEN (`B_TIME_REVERSAL_EVEN`, 039).

**Field-set union (branch-dependent; verdict `field_set_*`).** `θ_B` (brane phase)
is ELIMINATED by Schur into `A_eff=ρ_br+C_J²/κ_phase` — DISTINCT from the RETAINED
Madelung flow phase `θ`:

- PARENT: `{ρ, θ, r_B≡χ_B, H, u_L, u_T(divfree), u_w, λ_Σ(aux)}` (`u_L`
  independent, `h=P₀H` derived).
- REDUCED: `{ρ, θ, r_B≡χ_B, u_L, h, u_T(divfree), u_w, λ_Σ(aux)}`.

## The wall static-energy reconciliation

`r_B` (G0) and `χ_B` (Part I) are the SAME `[0,1]` wall order parameter. Under the map

```text
r_B ≡ χ_B ,   κ_χ ≡ κ_B ,   λ_χ = 4 a_B
```

the build CONFIRMS (computed, able-to-fail): `ℓ = √(2κ_χ/λ_χ) = √(κ_χ/2a_B) = δ`;
`σ_χ = κ_χ/(6ℓ) = √(2a_Bκ_B)/6 = σ_wall`; coincident minima {0,1}; and equality of
the **STATIC WALL-ENERGY densities** (potential + gradient + resulting
kink/width/tension) under the map. This is `wall_dedup =
STATIC_WALL_ENERGY_IDENTITY`.

⚠ **Scoped explicitly to the static wall energy + order parameter — NOT a
full-dynamical-functional identity.** The FULL functionals DIFFER
(`FULL_FUNCTIONALS_DIFFER`): G0's `S_χ` is INERTIAL (`½Z_χ ṙ_B²`); Part I supplies a
free energy + DISSIPATIVE adjunct dynamics `D_t χ_B = −M_χ μ_χ (+ Γ_B)` with NO
`Z_χ χ̇²` term. Do NOT claim a full-Lagrangian identity. Able-to-fail: perturb the map
(`λ_χ→2a_B`) ⇒ width/tension mismatch.

**Precision correction (NOT "three descriptions collapse to one").** It is (a) an
EXACT static-energy dedup `r_B≡χ_B`; PLUS (b) a **SCOPE-SPLIT for `g_ℓ`** — the
Gaussian localization envelope `g_ℓ(w)=exp(−(w/ℓ_g)²)/(√π ℓ_g)` is NOT merged/removed.
`χ_B` supersedes only `g_ℓ`'s MATERIAL-STATE-CLOSURE role (register edge **R21**, a
scope split, NOT a reduction), while `g_ℓ` PERSISTS as the light-brane localization
envelope in `S_brane`. The `g_ℓ` scope-split CITES the existing R21 — it does NOT
merge, and no new edge is minted for it (`g_ell = SCOPE_SPLIT_R21_NOT_MERGED`).

**Light gating honesty.** The FROZEN light action uses the `g_ℓ` localization
envelope; its dynamic `χ_B`-projection is the DISTINCT stage006 term
`f_shear = χ_B·½μ_R⁽⁴⁾(curl₄ u_d)²` with `μ_R=∫χ_B μ_R⁽⁴⁾ dw` (edge **R17**,
PENDING; firewall `[μ_R]≠[μ_R⁽⁴⁾]`, R22). The two are NOT the same object and must
NOT be double-counted (`light_gating = g_ell_LOCALIZED; chi_B_projection R17
PENDING`).

## The drain interfaces — NAMED BOTH, adjudicate NEITHER

044 does **NOT** assert the two candidate drains are the same, and does **NOT**
select one. Both are **non-variational** (neither is a term in a conservative
action), so **`S_cons^G0` / `S_assembled` contain NEITHER** — that is the computed
anchor (tooth G: `δS_assembled/δθ ⇒ ∂_tρ+∇·(ρv)=0`, `conservative_continuity =
SOURCE_FREE`). The two candidate interfaces are DIFFERENT balance-law objects:

- **Part-I `Γ_B`** — an INTERNAL order-conversion:
  `∂_t(χ_B n)+∇₄·(χ_B n u+J_χ)=n Γ_B`, `Γ_B=Γ_return−Γ_drain`, dim `T⁻¹`; converts
  ORDERED↔DISORDERED material at the throat; **total `n` conserved exactly**;
  NON-variational (an added RHS source, NOT in `F`; the `−M_χ μ_χ` gradient-flow
  piece IS variational, `Γ_B` is not).
  (`part_I_Gamma_B = internal_order_conversion_total_n_conserved`.)
- **G0-card `S_drain`** — a LOCAL SOURCE in TOTAL-`ρ` continuity with a remote
  return: `S_drain=−Σ_i Γ_0 D_i`, `∂_tρ+∇·(ρv)=S_drain+S_leakage`, charge-EVEN;
  "deliberately non-variational."
  (`g0_S_drain = total_rho_sink_plus_remote_return`.)

No holonomic reduction mapping one to the other is supplied, and `S_hold` does NOT
toggle between them — the earlier "option (c) regime toggle" claim is WITHDRAWN.
044's deliverable: (i) PROVE `S_assembled` is source-free (tooth G); (ii) NAME both
non-variational interfaces + their distinct balance-law type (cited provenance,
tooth H); (iii) record their EQUIVALENCE as **UNRESOLVED** and their SELECTION as
**DEFERRED to the 045 non-variational block + the standing user mini-gate**:

```text
drain_interfaces  = TWO_NONVARIATIONAL_NAMED
drain_equivalence = UNRESOLVED
drain_selection   = DEFERRED_045_AND_USER_GATE
```

## P-retirement — four category-clean manifests

Decision 16 retired the brane polar field `P` under
`INSTABILITY_CONFIRMED_STRUCTURAL` (helical Lifshitz instability
`det=k²(2A_P μ_R k−λ_Pu²)<0`, no gyroscopic rescue); the wall is now UNAMBIGUOUSLY
the scalar `χ_B` (pathA_24 T1 falsified the polar-vector wall). The bookkeeping is
FOUR SEPARATE, CATEGORY-CLEAN, COMPUTED manifests (fields, parameters, and
cross-stage terms are NOT mixed in one set; verdict `P_retirement`):

1. **stage007 ACTION-summand partition** — operative `{S_GNLS, gL_Mac, gL_uw}` vs
   retired `{L_pol, gL_Pu}`, DISJOINT. (These tokens are STAGE007's — cited to
   stage007, NOT to decision-16, which records the DECISION, not the token set.)
2. **stage007 FIELD/DOF removal** — `P` removed ⇒ DOF `8→4`.
3. **stage007 DRIFT reduction** — `11→7` (`−4`, from `λ_Pu` / polar dynamics).
4. **stage006 DRIFT reduction** — `6→5` (`−1`, from `α_aniso` — this belongs to
   STAGE006, NOT the stage007 action summands).

Net route-less reduction `−5 = −4 (stage007 drift) − 1 (stage006 drift)` — **in the
DRIFT counts, NOT the DOF count.** A symbolic summand-set PARTITION over the immutable
stage007-T0 record (banked PENDING; no built v2 stage is touched now), NOT byte
surgery. The polar route is retired-but-NOT-foreclosed (re-entry needs a NEW T0-level
freeze).

## `Z_χ` — the new DRAFT-v0 knob + count consequence

The conservative (inertial) wall `S_χ` needs the inertial normalization `Z_χ`, which
the G0 card declares `[POSTULATE]` and which has NO derivation/reduction route in the
register (the register has `a_B, κ_B` and the DISSIPATIVE `M_χ`, but no `Z_χ`). It
CANNOT be removed by the static wall map: `r_B∈[0,1]` with fixed well minima {0,1}
pins the field normalization, so `Z_χ` is NOT absorbable by rescaling `r_B`. 044
records `Z_χ` as a **NEW DRAFT-v0 action input, `[POSTULATE]` status**, with NO
reduction route — the honest inertial-completion cost of building a conservative wall
from Part-I's dissipative one (`new_draft_knob = {Z_χ:
POSTULATE_inertial_normalization_no_reduction_route (G0 DRAFT-v0)}`).

**Count consequence (do NOT claim "unaffected").** By 043's OWN counting rule (a
route-less continuous ACTION input is counted unless a reduction is
DERIVED-and-EXECUTED), an independent `Z_χ` shifts the continuous codimension range

```text
[40, 49]  →  [41, 50]      (both endpoints +1; spread 9 unchanged)
```

044 records this as an explicit **count-reconciliation / sensitivity edge (R97)**,
NOT as "unaffected" (`Z_chi_count_consequence = SHIFTS_CONTINUOUS_[40,49]->[41,50]
(043-rule; both endpoints +1; spread 9 unchanged); committed re-count DEFERRED_046
(may substitute dissipative M_chi)`).

⚠ **The 043-R92 parallel.** 043's own recorded sensitivity R92 already flags
`K_θ/κ_phase`-as-2-DOF as an INDEPENDENT `+1` source that would also carry the
continuous range to `[41,50]`. These are **TWO INDEPENDENT `+1` sources**, each
DEFERRED to 046. `Z_χ`'s shift is recorded as ITS OWN sensitivity; 046 (the
count-consumer stage) reconciles all such shifts into the committed re-count — do NOT
double-add and do NOT contradict 043's headline. Whether `Z_χ` is a genuine net `+1`
or SUBSTITUTES for the dissipative `M_χ` (which lives in the non-variational sector /
045, not the conservative action) is a 046 calibration-map reconciliation. **The
`[40,49]` 043 headline STANDS as the committed Parts-I–VI count until 046 ratifies the
parent-action re-count.**

## Wave speeds + factorizations

Three DISTINCTLY-attributed wave speeds (`wave_speeds`; the scalar `c_±²` is NOT a
bulk speed, `c_γ` is a light import):

```text
bulk    c_s0²  = 5 K ρ0⁴ / m                     (from S_bulk)
scalar  c_±²   = (3 ± √2)/2 = {0.7928932188, 2.2071067812}   (units c_E²; D* = 7/4)
light   c_γ²   = μ_R / ρ_br                       (edge R4 import)
```

The scalar pair are the roots of `det(K−zM)=0` with `M=diag(1,1)`,
`K=[[2,½],[½,1]]`; positivity margin `D* = B_eff K_h − C_hu² = 7/4`. Stability is
`POSITIVE_DEFINITE` via the **Sylvester** criterion (leading minor `B_eff>0` AND
`D*=7/4>0` — det-positivity alone would permit negative-definite); wall double-well
minima (not maxima) at {0,1} (`V_χ''>0`); transverse block positive-semidefinite with
gap `Ω_w²≥0`.

Two operator factorizations, each with its annihilated zero mode
(`wall_factorization` + the charge block):

```text
L_χ^(2)/κ_χ = A_χ†A_χ = −∂_x² + (1 − (3/2)sech²(x/2ℓ))/ℓ²,   A_χ = ∂_x + tanh(x/2ℓ)/ℓ,   A_χ r₀' = 0
O_⊥       = A†A = −∂_w² + [4 − 6 sech²(w/ℓ)]/ℓ²,             A   = ∂_w + 2 tanh(w/ℓ)/ℓ,   O_⊥ f₀ = 0
                                                              f₀ = 1/[ℓ cosh²(w/ℓ)],       N₀ = 8/(3ℓ)
```

`L_χ^(2)` is the second-variation / Hessian (NOT the full Lagrangian `L_χ`); `A_χ`
annihilates the kink slope `r₀'` (the translational zero mode). Dimensional
homogeneity (tooth F): EVERY term of `S_assembled` carries the SAME action dimension
once all field + coefficient dims are restored ([L,T,M], no natural-unit hiding of a
dropped Jacobian / `a⁵`), free-carrier-independent.

## Gravitomagnetism — prose pointer + honest asymmetry

⚠ **PROSE / pointer item ONLY in 044's note + card — NO register row (that is 047's
job).** The gravitomagnetic (velocity-dependent flow) terms in the gravity sector —
frame-dragging / Lense–Thirring — are the velocity/boost part of gravity, reproduced
by the audited GR-matched **1PN→4PN ladder** (`research/4d_*pn*`; `g_{0i}` / `v`-terms).
The unification highlight: the **"B-analog = boost of the E-analog" structure holds on
BOTH forces** — magnetism = boost of the electric ±w throat (Part V);
gravitomagnetism = boost of the gravitoelectric drain/return flow
(`gravitomagnetism = BOOST_OF_GRAVITOELECTRIC_FLOW`).

⚠ **The HONEST asymmetry (never imply perfect force-symmetry):** gravity matches GR
CLEANLY (incl. gravitomagnetism, through 4PN); EM genuinely **DEPARTS** from exact
Maxwell (`NATIVE_P_NO_EMERGENT_GAUSS`, 033; `b_T` T-even, 039) — the gravitomagnetic
side is BETTER-covered.

## DRAFT-v0 / completeness-floor / throat-debt caveats

Recorded honestly, NEVER softened:

- `g0_card = DRAFT_V0` — read from the card header flag (not a literal). The card is
  postulated, not ratified; many gates DEFERRED.
- `closed_action = false` — DERIVED as `closed_action = (no DEFERRED gate remains)`;
  the DEFERRED gate set is non-empty, so the boolean is FALSE. 044 assembles the
  structure at the completeness FLOOR; it does not close the action.
- `throat_solve = SIM_DEFERRED_CENTRAL_DEBT` — the shared throat solve is the central
  sim-deferred reduction debt (the SAME solve carrying the gravity/electric/magnetism
  R1s + the `{ρ_br, μ_R, c_E}` reduction debt). 044 NAMES it, does not perform it.
- No two-body / far-field solve is performed. The reconciliation is at the STRUCTURE /
  completeness-floor level — NOT a proof that the closed unified action exists.

## The tooth roster

Dual-engine, both exit 0, **20 executable teeth each** (SymPy `.py` + materially
independent Mathematica `.wl`), **43 mutations** in the per-tooth ablation. The
tooth order:

```text
A_WALL_STATIC_ENERGY_DEDUP          — wall static-energy dedup under the map (§ wall reconciliation)
B_WALL_HESSIAN_ZERO_MODE            — L_χ^(2)/κ_χ = A_χ†A_χ, A_χ r₀'=0
C1_CHARGE_FACTORIZATION_ZERO_MODE   — O_⊥=A†A, O_⊥f₀=0, f₀=1/[ℓcosh²(w/ℓ)]
C2_CHARGE_ZERO_MODE_NORM            — norm N₀=8/(3ℓ)
D1_BULK_SOUND_SPEED                 — c_s0²=5Kρ0⁴/m from S_bulk
D2_SCALAR_WAVE_SPEEDS               — c_±²=(3±√2)/2 as roots of det(K−zM)=0
D3_LIGHT_WAVE_SPEED                 — c_γ²=μ_R/ρ_br
E1_SCALAR_SYLVESTER                 — B_eff>0 AND D*=7/4>0 (leading minor + det)
E2_WALL_MINIMA                      — double-well minima {0,1}, V_χ''>0
E3_TRANSVERSE_GAP                   — Ω_w²≥0 positive-semidef gap
F_DIMENSIONAL_HOMOGENEITY           — units-restored [L,T,M], whole S_assembled
G_SOURCE_FREE_CONTINUITY            — δS/δθ ⇒ ∂_tρ+∇·(ρv)=0 (the drain anchor)
H_DRAIN_PROVENANCE                  — TWO_NONVARIATIONAL_NAMED (cited source-facts)
I_P_RETIREMENT                      — the four disjoint P-retirement manifests
J_FIELD_UNION_SCHUR                 — θ≠θ_B derived (incidence + Schur), field sets
K1_DEFERRED_GATE_COMPLETENESS       — extracted deferred-gate set == expected
K2_DRAFT_OPEN_ACTION                — closed_action=(no deferred gate)=False; g0_card=DRAFT_V0
L_ACTION_INCIDENCE                  — keyed summand→field-incidence map == expected
M_PROVENANCE_STATUS_BINDING         — R21/R17/Z_χ/count/gravitomag/throat provenance binds
REPRODUCTION                        — first-match verdict + field binding
```

## Verification

- **Dual-engine, both exit 0, byte-identical verdict.** SymPy 20/20 + Mathematica
  20/20; the canonical-ordered-JSON verdict object is byte-identical across engines
  (**2196 bytes**).
  `scripts/ledger_stage044_parent_action_reconciliation_sympy_audit.py` — **SymPy 20
  teeth**.
  `mathematica/ledger_stage044_parent_action_reconciliation_mathematica_audit.wl` —
  **Mathematica 20 teeth**. The verdict is COMPUTED via a documented first-match over
  the per-check results, never a stored literal; the runner performs the cross-engine
  byte comparison.
- **Per-tooth ablation: `43/43` fired across BOTH engines** (after the can't-fail
  fix). Every able-to-fail check is individually ablated via `LEDGER_STAGE044_MUTATION`
  to a non-zero exit AT that check's assert; the guards
  `FIRED_AT_OWN_ASSERT`/`MUTATION_DID_NOT_FIRE`/`UNKNOWN_MUTATION` are genuine
  (adversarially probe-verified). Tooth K uses the resolve-ALL (not drop-one) mutation
  for `closed_action`; tooth L fires on a DELETED summand; every §4 JSON field maps to
  a producing tooth (A–M / REPRODUCTION).
- **Directive bookend Codex→Grok→Codex GREEN** — Grok reproduced all 6 numeric
  identities independently.
- **Tri-review: FIDELITY CLEAN.** The TRANSLITERATION + ADVERSARIAL legs substantiated
  3 can't-fail SUB-conjuncts (tooth F `carrier_normalization_dependencies`, tooth F
  `target_action_dimension`, tooth E1 `B_eff`) — all **FIXED** (removed 2 redundant
  echoes; added `E_BEFF_PERTURB` making the Sylvester `B_eff>0` leading-minor genuinely
  ablated). **044 ships with ZERO can't-fail conjuncts.**
- **The ONE honest non-blocking robustness note.** A subset of ELEMENTARY teeth
  (B / D1 / E2 / G) + the data manifests (H–M) are `.wl`↔`.py` transliterations
  (CAS-engine diversity still cross-checks them; the HARD falsifiable teeth
  C1/C2/D2/E1 ARE method-independent via DSolve/Integrate/Eigenvalues/Reduce;
  method-diversity is inherent-impossible for provenance manifests) — an accepted
  robustness limitation, disclosed not hidden.
- ⛔ ~~**GLM-5.2 tertiary: NOT yet run** — HELD for the user's explicit confirm; the one deferred check.~~
  **CLOSED, NOT OWED (user decision, 2026-07-29/30):** the tertiary pass is retired, so nothing is
  outstanding on this stage. 044's fresh-agent fidelity + adversarial legs were run and were CLEAN.

## Ledger accounting — what 044 does NOT do

- **Does NOT close the action** (`closed_action=false`); the closed unified action
  awaits the shared sim-deferred throat solve.
- **Makes NO drain decision** — records `drain_equivalence=UNRESOLVED` +
  `drain_selection=DEFERRED_045_AND_USER_GATE`; it does NOT select, equate, or build
  the drain (that is 045 + the standing user mini-gate).
- **Does NOT claim "no new knobs"** — `Z_χ` is disclosed as a NEW DRAFT `[POSTULATE]`
  input with a recorded `[40,49]→[41,50]` count sensitivity.
- **Does NOT re-count** — the committed count re-reconciliation is 046's; the `[40,49]`
  043 headline STANDS until then.
- **Mints NO gravitomagnetism register row** (047's job) — the gravitomagnetism item
  is PROSE only.
- **Does NOT build the non-variational block** (`S_drain`/`S_leakage`, the controllers,
  the BCs, the zero-mode quotients, the `F_var+F_flux+F_𝔅+F_rad` partition) — that is
  stage 045.

## Downstream consumers

- **Stage 045 (VII-2b, the non-variational block)** — CONSUMES both named drain
  interfaces (`TWO_NONVARIATIONAL_NAMED`) + the source-free conservative anchor, and
  carries the drain-placement crux to the **standing USER mini-gate** (equivalence
  UNRESOLVED, selection DEFERRED to here).
- **Stage 046 (the calibration map / count-consumer)** — CONSUMES the `Z_χ`
  count-reconciliation (`[40,49]→[41,50]`, may substitute the dissipative `M_χ`) and
  reconciles it with 043's R92 `K_θ/κ_phase` `+1` sensitivity into the committed
  re-count; also consumes the star-witness calibration map.
- **Stage 047 (the permanent registers)** — mints the gravitomagnetism one-line edge
  (044 supplies the PROSE pointer + the honest EM-departs asymmetry only).

## Provenance

- **Source:** there is NO `pathA_44` source gate — 044 is a genuinely-NEW synthesis.
  Its ground truth is the G0 closure card
  (`software/em_charge_attribute/g0_closure_card_v0.md`, DRAFT-v0) + the assembled
  Parts-I–VI ledger (`parameter_register.md`, R1–R92) + the recon-extracted action
  terms. The card summands + identities are TRANSCRIBED into in-engine data structures
  as the stage's self-contained input facts (no runtime file-reads, no source-grep).
- **Consumes:** cites the assembled register (R1–R92) + the G0 card; the shared
  sim-deferred throat solve (`SIM_DEFERRED_CENTRAL_DEBT`) that would close the action is
  NAMED, not performed.
- **Governing:** `notes/part7_integration_atomic_split.md` (the RATIFIED 7-stage split;
  the spine `044→044+045`); `notes/stage044_parent_action_prep.md` (the prep note);
  `_scratch/stage044/stage044_synthesis_directive.md` (the build directive, rev 2 — Codex
  `ISSUES_FOUND(12)` folded) and `OUT_stage044_ablation.txt` — ⛔ **neither is retained**; both lived in
  gitignored `_scratch/` and no copy survives, so these names record that the artifacts existed and are
  not auditable citations. The verdict pair `verdict_py.json`/`verdict_wl.json` **was** rescued, to
  `notes/stage044_evidence/`; `notes/parameter_register.md` (R1–R92 + the new R93–R97
  + the `Z_χ` master row); `research/pde_ledger/notes/MATHEMATICA_MIRROR_POLICY.md`;
  `notes/ledger_v2_blueprint.md` §5 (standalone engine spec) + §6 (per-tooth ablation).
