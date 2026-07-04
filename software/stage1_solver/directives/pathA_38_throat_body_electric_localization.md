# pathA_38 — Charge as a 4D throat-body interaction: does it give a brane-localized `1/r²` (Coulomb)?

**Status:** ⭐ §2 MEDIATING FIELD **RESOLVED** (Claude↔Codex rounds 1–2) + **§2A executable primitives folded, Codex design-review
GREEN** (r1→r2→r3) + **GLM tertiary pass folded** (`SOUND_WITH_CONCERNS`, 5 concerns + 4 notes all folded) + **Codex re-confirm r4 =
`SOUND_AS_IS`** (after the one stale-label cleanup, 2026-07-03). ✅ **DIRECTIVE DONE — READY FOR DUAL-ENGINE EXECUTION.** NEXT = commission
SymPy+Mathematica scripts (§2A) → ENGINE_AGREE → tri-review. Logs in `_scratch/pathA_38_*`. This is the **carry-over spec** written while full context is fresh (2026-07-03, pre-`/compact`). **Type:** throat/bulk
sector (4D), linearized inter-throat interaction. **Supersedes the mechanism of** `pathA_37` (the flow/counterflow gate — RETIRED;
kept as a documented waypoint). **Author:** orchestrator (scaffolding).

---

## §0. Why this gate — the crux in one paragraph

The v7 conceptual pivot (`docs/conceptual_foundation.md` §3/§4 ⭐ v7; memory `feedback-native-em-mechanisms`): **electric charge is not
a brane flow or deformation** (all three flow/deformation routes were explored and set aside — deformation → wrong energy `1/r⁶`;
one-fluid flow → charge=mass; counterflow → needs a 2nd component the `T=0` condensate lacks in flat 3D). Instead, **charge is the
interaction of two throats' 4D BODIES beyond the mouth** — a geometric interaction of the extended throat structures in the bulk, with
`±w` = the sign. **THE MAKE-OR-BREAK:** an electric field is long-range (`1/r²` force), so the throat-body interaction must come out
**`1/r²` in the 3D brane.** A source acting freely through the **4D bulk** gives the WRONG falloff (4-space → `1/R²` potential →
`1/R³` force); it is correct 3D Coulomb **only if the throat's influence is brane-localized.** This gate tests exactly that.

---

## §1. The leverage — pathA_29 already did this for GRAVITY (reuse the machinery)

**Do not reinvent the localization analysis — pathA_29 is the template.** `pathA_29` (`RETURN_RESIDUAL_PREDICTION`, tri-reviewed;
`reports/pathA_29_*`) proved that **gravity's field survives the finite brane slab as `1/r²`, not `1/r³`**, because the mediating field
has a **normalizable transverse (`w`-direction) zero mode localized to the slab** → the in-brane field solves a real 3D-radial
equation → `p=2` (`1/r²`). A *delocalizing* warp gives `p=3` (`1/r³`) and was the counterfactual that got rejected. **This is the
Randall–Sundrum-style brane-localization mechanism**, already demonstrated in-repo for the gravity sector.

**The pathA_38 question is the same PRINCIPLE applied to the CHARGE (throat-body) sector:** does the field that mediates the inter-throat
electric interaction inherit a **normalizable localized transverse zero mode → `1/r²`**, or does a warped/nonlocalizing/operator/source
failure **restore bulk `1/r³` (FAIL)**? **The shared item is the PRINCIPLE, not the mechanism (N1).** Gravity (pathA_29) localized via a
FINITE SLAB and the scalar operator `−d²/dw²` with constant zero mode `1/√d` (which is NOT normalizable as `d→∞`;
`pathA_29_results.yaml:199-211`); charge localizes via the Z₂ WALL PROFILE, `f₀∼sech²(w/ℓ)`, normalizable on `ℝ` with a flat measure. So
the mechanisms are DISTINCT and **charge's is stronger** (no finite slab needed). A PASS is therefore "the electric sector realizes the
same normalizable-zero-mode principle, by a distinct stronger mechanism"; a FAIL (delocalized/gapped/no-monopole/ghost) is a first-class
no-go (charge can't be long-range Coulomb in this medium).

---

## §2. ⭐ THE MEDIATING FIELD — RESOLVED (Claude↔Codex rounds 1–2, 2026-07-03)

The load-bearing open piece is now decided (discussion logs: `_scratch/pathA_38_field_mechanism_codex_round{1,2}.md`). The mediator is
**NOT** the free internal `±w` director — that is generically **gapped** by the stable Z₂ ordering (tilting the arrows off-axis costs
energy → `M_or²>0`) → the radial problem becomes `g''+(2/r)g'−μ²g=0` → **Yukawa `e^{−μr}/r` → short-range → `FAIL`.** That is the
SAME wrong-sign / stabilizing-stiffness wall as `pathA_36` (C5) and `pathA_25`; a continuous-director escape (to make it massless)
would reconnect the `±w`/Z₂ vacua and reopen the connected-vacuum wall problem. So a free director cloud is generically `FAIL_YUKAWA`.

**The mediator is the gapless transverse-embedding / orientation-lock Goldstone `h(x)`** — the wall's transverse displacement into
`w`, with the `±w` arrows **locked to the wall normal**. The gaplessness is protected because it is the **Goldstone of spontaneously
broken `w`-translation**: a localized wall in a translation-invariant bulk has the exact zero mode `f₀(w) ∝ ∂_w[background wall
profile]`, with finite norm equal to a wall-tension-like integral. Role split (supersedes the round-0 "3 candidates" list):
- **the native `±w` orientation order (A) is the mediator ONLY as this locked Goldstone**, not as a free director cloud;
- **the throat 4D funnel/body (C) is the SOURCE geometry** — it sets a signed boundary datum on `h`, it is NOT itself the mediator
  (this keeps it off the retired in-plane-deformation `u_L` route: a transverse SCALAR bending mode `h~1/r` gives `(∇h)²~1/r⁴`,
  self-energy `~1/a` — the RIGHT scaling — whereas the dead route was the in-plane vector `∇(∇·u)=0 → u_L~1/r² → (∂u)²~1/r⁶`);
- **the scalar Coulomb potential (B) is the DERIVED OUTPUT** after the localized zero-mode projection, never a postulate (postulating
  B bakes in `1/r` and is circular).

**The operator is a COUPLED fluctuation matrix, not a scalar Sturm–Liouville problem with `M_or²→0` by declaration.** Fluctuation
vector `Ψ = (wall/order transverse displacement, director-normal component, relative director tilt, …)` on the finite-slab + bulk
geometry. One eigenvector is the **exact Goldstone** (`f₀ ∝ ∂_w[profile]`: "translate/bend the whole wall and let the arrows follow the
normal"); the orthogonal **relative-tilt** combination is **massive** and must be projected OUT. If the arrows are NOT locked to the
normal, the mediator collapses to a gapped Frank/director cloud → `FAIL_YUKAWA`.

The three sub-questions are now concrete **DERIVE-not-assume** requirements (each able-to-fail — see §3):

1. **Gaplessness — the branon must stay massless.** `f₀ ∝ ∂_w(profile)` is the wall-translation Goldstone. It is massless ONLY if the
   wall is **not pinned**. If the medium/slab/stack pins the wall, or the earlier `u_w` gap (`pathA_35` Gate L) applies to this same
   mode, `h` is massive → `FAIL_PINNED_BRANON`. Check, do not assume.
2. **The throat-body geometry + monopole overlap.** Model the throat as a funnel (area > mouth; extrinsic curvature into `w`), reduced
   à la pathA_29 (NO nonlinear interior — sim-deferred), leveraging `[[project-light-4d-throat-hypothesis]]` (the 4D-throat/geon
   structure). The effective charge is the overlap `q_eff ∼ ⟨f₀, S_throat⟩`. A symmetric bend/dimple with zero net signed normal flux
   sources a **dipole**, not a monopole → no Coulomb → `FAIL_NO_MONOPOLE`. The signed monopole overlap must be nonzero.
3. **The sign — DERIVED from the throat boundary condition, not assumed.** `±w` (puncture direction) = the sign of the signed source.
   Like-repel/unlike-attract must fall out of the actual throat BC + positive field energy, **NOT** from "tensioned-membrane"
   intuition (a bare normal force on a membrane can give scalar-like *attraction*). Wrong sign via an unstable / indefinite mediator
   sector → `FAIL_GHOST_INSTABILITY`; a bare `σ_mouth` flip is only a non-physical classifier-sensitivity control (the like-repel sign
   is FORCED on a positive-definite PASS branch — §2A). (The earlier GLM Green's-identity + Bernoulli sign result carries forward as
   corroboration, but the sign must be re-derived from the localized mode + BC here.)

---

## §2A. Reduced-model EXECUTABLE PRIMITIVES (Claude↔Codex authored 2026-07-03 — folds the r1 design-review `NOT_SOUND`)

**Provenance:** authored by Codex (model/route), reviewed by Claude for anti-circularity (Goldstone norm + `w→−w` parity checked by
hand). Concreteness matches pathA_29.

**⭐ SIGN/COUPLING PROVENANCE — RESOLVED (Claude↔Codex confirm-r2, 2026-07-03).** The main-branch sign is NOT a free knob:
- **`σ_mouth ≡ s`** — the mouth-work sign IS the `±w` puncture direction, not an independent branch label.
- **The like-repel/unlike-attract STRUCTURE is DERIVED, not chosen:** on the PASS branch (normalizable gapless zero mode + positive-
  definite Hessian) the static zero-mode Green function is `G₀(R)=f₀(w)f₀(w')/N₀ · 1/(4πR) > 0`, so the two-throat cross term is
  `U_int = s₁s₂·(q_eff²/N₀)·1/(4πR)` → `s₁s₂=+1` repels, `s₁s₂=−1` attracts. The script MUST EMIT this from the SOLVED `G`
  (`U_{++}>0`, `U_{+−}<0` as outputs), never from choosing `σ_mouth`; the `σ_mouth=−1` row in §2A(3) is a NON-PHYSICAL able-to-fail
  control.
- **What stays CALIBRATED (first-principles sim-deferred):** the charge MAGNITUDE `Q_E` (universal — one `e`) and the modeling that a
  throat sources a NONZERO odd monopole (`q_eff≠0`). `u_E` is the projection of the (sim-deferred, nonlinear) mouth BC into `(η,n,τ)`;
  here it is a declared calibrated datum, and `FAIL_NO_MONOPOLE` (`S_orth`) is the guard that the overlap is a REAL computed quantity,
  not assumed.
- **⭐ HONEST HEADLINE / PROVENANCE SPLIT (GLM fold, 2026-07-03):**
  - **PREMISE (conditional Goldstone theorem):** `h` is gapless ONLY if the Z₂ wall is an unpinned self-localizing wall in an exactly
    `w`-translation-invariant bulk. `Ô f₀=0` is a consistency identity (differentiate the background EOM), NOT a discovery; it can fail
    only when explicit pinning/confinement is added (`FAIL_PINNED_BRANON`).
  - **DERIVED (given that premise):** the zero-mode NORMALIZABILITY + the extracted power `p` are the real localization teeth; the
    `w→−w` charge≠mass parity split is derived at the linear/source level.
  - **GENERAL THEOREM (not model-specific):** like-repel/unlike-attract follows from ANY stable positive-definite mediator with linear
    source coupling (`G₀>0`, `q_eff²>0` ⇒ `U_int=s₁s₂(q_eff²/N₀)/(4πR)`). The model-specific content is only {that the interaction is
    mediated by this Goldstone at all; that the mediator sector is stable}.
  - **CALIBRATED (first-principles sim-deferred):** charge magnitude `Q_E` (universal — one `e`); that the nonlinear throat mouth
    sources a nonzero odd monopole (`q_eff≠0`); the full throat-derived positive-definite Hessian/stability.
  - This is the project's calibrate-predict landing ([[feedback-calibrate-predict-methodology]]).

### (1) Fields, background, operator, BCs, parity
- **Fields** on the doubled slab `I_d=[−d,d]` (implemented as `0≤w≤d` with parity BCs at `w=0`; wall at `w=0`):
  `Ψ(x,w)=(η,n,τ)ᵀ` — `η` = wall/order-profile (translation) fluctuation [physical branon `h(x)` via `Ψ_h=−h(x)f₀(w)`]; `n` =
  director-normal component; `τ` = relative director tilt (after subtracting the geometric wall-normal rotation). Nonlinear throat
  interior NOT modeled (sim-deferred).
- **Background** (translation-invariant Z₂ wall family): `Φ₀(w)=(χ₀,N₀,0)ᵀ`, `χ₀=χ_*T(w/ℓ)`, `N₀=N_*T(w/ℓ)`, `T` odd, `T(±∞)=±1`
  (e.g. `tanh`); no tuned profile beyond solving the background EOM.
- **Static transverse Hessian:** `Ô=−∂_w(K_w∂_w)+M²(w)`, `K_w=diag(κ_χ,κ_N,κ_τ)>0`, `K∥(w)=diag(Z_χ,Z_N,Z_τ)>0`,
  `M²(w)=Hess_U[χ,N,τ]|_{Φ₀(w)}+diag(0,0,Λ_τ(w))`. Eigenproblem `Ô f_n=m_n² K∥ f_n`, inner product `⟨f,g⟩=∫_{I_d} f†K∥ g dw`.
- **Goldstone candidate:** `f₀(w)=∂_wΦ₀=(χ₀',N₀',0)ᵀ`. Acceptance requires COMPUTED `Ô f₀=0` (incl. BC residuals) and finite
  `N₀_norm=⟨f₀,f₀⟩<∞`. (Claude-checked, `T=tanh`/const weights: `N₀_norm ∝ 2[tanh(d/ℓ)−tanh³(d/ℓ)/3]/ℓ` → finite.)
- **BCs at `w=d`** (DC completions, à la pathA_29): `destructuring_absorbing_DC`: `(K_w∂_w f_n)|_d=0` (drain accounting in the source,
  not a tuned impedance); `bloch_stack_q0`: `f_n(d)=f_n(−d)`, `(K_w∂_w f_n)(d)=(K_w∂_w f_n)(−d)`; `sommerfeld_AC_only`: finite-`ω`
  only, EXCLUDED from the DC falloff verdict.
- **Parity:** `P_wΨ(w)=ΣΨ(−w)`, `Σ=diag(−1,−1,+1)`. Electric Goldstone `P_w f₀=−f₀` (parity-ODD); electric sources odd, gravity/drain
  sources even. (Claude-checked: `χ₀` odd → `χ₀'` even; `Σ=−1` on `η,n` → `P_w f₀=−f₀`. ✓)
- **⭐ pathA_35 reconciliation / explicit overturn (GLM concern 1).** On the PASS branch, `pathA_38` **OVERTURNS** the computed
  `pathA_35` Gate-L `u_w` gap (`Ω_w²>0`; `pathA_35_gateL_light.md:69`, `:93`, `:111-112`, `:132`): the transverse brane displacement
  that `pathA_35` required GAPPED is this gate's gapless electric mediator `h`. The distinction is **physical**: `pathA_35` used a
  flat/imposed-`V_conf` brane, where `w`-translation is EXPLICITLY broken, so no Goldstone protects `u_w`; `pathA_38` assumes a genuine
  self-localizing Z₂ material-state wall (rung-W-passed) in a translation-invariant bulk, so `w`-translation is SPONTANEOUSLY broken and
  `h=∂_wΦ₀` is gapless by the Goldstone theorem. The `pathA_35` fifth-force objection (`FAIL_BENDING_MASSLESS_FIFTH_FORCE`) is resolved
  HERE by **parity, not by a gap**: `h` is odd, mass/gravity sources are even, so `⟨f₀,S_grav_even⟩=0`. This overturn is **CONDITIONAL
  on exact, unpinned bulk `w`-translation symmetry** (the open premise, tested by `FAIL_PINNED_BRANON`).

### (2) Throat source law `S_throat` (compact VOLUME source — NO hand-written radial potential)
- `S_throat^s(x,w;X_i)=s·Q_E·ρ_a(|x−X_i|)·ρ_b(w)·u_E`, `s=±1`, `u_E=(u_η,u_n,u_τ)ᵀ` fixed before solving; `ρ_a,ρ_b` compact
  normalized (`∫d³x ρ_a=1`, `∫_{I_d}dw ρ_b=1`, `supp ρ_a≤a`, `supp ρ_b≤b≪d`). Dims `[ρ_a]=L⁻³`, `[ρ_b]=L⁻¹`, `[S_A]=E/(L⁴[Ψ_A])`.
- `q_eff=⟨f₀,S_w⟩=∫_{I_d} f₀†K∥(Q_E ρ_b u_E)dw`, computed AFTER solving. `U_int(R)=½∫S_i†(x,w)G(x−x',w,w')S_j(x',w')…` from the
  SOLVED slab Green function `G` of `Ô`; **NO `R⁻ⁿ`/`1/R`/force law inserted by hand.** `p` extracted from the solved zero-mode radial
  dependence (static AND dynamic routes).
- Controls: `FAIL_NO_MONOPOLE` via `S_orth=S−f₀⟨f₀,S⟩/⟨f₀,f₀⟩` (→ computed `q_eff=0`). `σ_mouth=−1` is retained ONLY as a NON-PHYSICAL
  classifier-sensitivity control — it is NOT a physical FAIL mode on the PASS branch (the sign is forced once the Hessian is positive-
  definite). The physical wrong-sign FAIL is **`FAIL_GHOST_INSTABILITY`** (a non-positive-definite Hessian → `G₀<0`), §2A(3).

### (3) `fail_witnesses`
| fail mode | ablation knob | expected COMPUTED classifier input |
|---|---|---|
| `FAIL_YUKAWA` | source only the relative tilt: `u_E=e_τ`, `Λ_τ>0` | lowest source-overlapping eigenvalue `m_src²>0`; Green tail `e^{−m_src R}/R`; `gapped_yukawa` |
| `FAIL_PINNED_BRANON` | two-stage `u_w`-descendant witness: first compute the unpinned Z₂-wall descendant `m_desc,Z2²=⟨f₀,Ô_Z2 f₀⟩/⟨f₀,K∥f₀⟩`; then ablate by adding the pathA_35-style explicit confinement curvature `ΔÔ_conf` (flat/imposed `V_conf` branch), normalized so its projected `u_w` curvature is the computed `Ω_w,35²` | unpinned: `m_desc,Z2²=0`, `Ô_Z2 f₀=0`, `pathA35_gap_overturned=true`; confinement ablation: `m_desc,conf²=Ω_w,35²>0`, Yukawa far field, `FAIL_PINNED_BRANON` |
| `FAIL_DELOCALIZED_BULK_1_OVER_R3` | **replace the compact slab with a NONCOMPACT anti-localizing control** (`w∈[0,∞)` or `ℝ`, same near-wall source), anti-localizing measure `e^{2k_w w}K∥`, with `k_w` **exceeding the zero-mode transverse decay rate** (for the tanh wall `f₀∼sech²(w/ℓ)∼e^{−2w/ℓ}` ⇒ integrand `∼e^{(2k_w−4/ℓ)w}` ⇒ require `k_w>2/ℓ`) — pathA_29's half-line warp. (A finite `[−d,d]` weight, or `k_w` below threshold, keeps the norm finite and does NOT delocalize; the domain MUST be noncompact AND `k_w` above threshold.) | `∫f₀†K∥f₀dw=∞`; continuum Green; computed `p=3` |
| `FAIL_NO_MONOPOLE` | `S_orth`, or a compact dipole with zero `f₀`-overlap | `q_eff=0`; no `1/R` term |
| `SIGN_CLASSIFIER_SENSITIVITY_NONPHYSICAL` | flip `σ_mouth=−1`, or collapse labels by hand | classifier response changes — but this is NOT a physical FAIL branch (retained only as a sensitivity control) |
| `FAIL_GHOST_INSTABILITY` | flip the sign of a diagonal source-overlapping `K_w` (or `K∥`) entry | Hessian non-positive-definite / negative-norm sector; computed zero-mode Green component `G₀<0`; `U_{++}<0` or sign matrix wrong; `ghost_instability_wrong_sign` |
| `FAIL_GRAVITY_MIXING` | contaminate mass source: `S_mass=S_grav_even+ε_mix S_E_odd` | `⟨f₀,S_mass⟩≠0`; neutral mass excites `h`; `gravity_mixes_with_h` |
| `FAIL_OPERATOR_PARITY_MIXING` | **SIM-DEFERRED:** nonlinear throat modifies the operator, `Ô→Ô+δÔ_throat` with `[δÔ_throat,P_w]≠0` | even/odd block split fails; the even (mass) sector overlaps the odd Goldstone; `⟨f₀,S_grav_even⟩≠0`; reopens the pathA_35 fifth-force concern (distinct from source-level `FAIL_GRAVITY_MIXING`) |
| `FAIL_SOURCE_NOT_COMPACT` | **SIM-DEFERRED:** replace compact `ρ_b` by the actual throat-funnel body if its `w`-support is noncompact | continuum spectral weight survives; a `1/R²` potential / `1/R³` force component can remain even with a normalizable zero mode; `source_not_compact_sim_deferred` |

**Note (N3):** because `f₀∼sech²(w/ℓ)` is normalizable on flat `ℝ` with a flat measure (unlike pathA_29's constant slab mode, which
needs the finite slab), `FAIL_DELOCALIZED_BULK_1_OVER_R3` requires a WARPED/anti-localizing noncompact control — scripts must NOT expect
a flat noncompact domain alone to fire it. (Charge localization is thus *stronger* than gravity's — §1.)

### (4) Goldstone / locking checks
- **Locking energy:** `E_lock=½∫[Λ_N(n−(N₀'/χ₀')η)²+Λ_τ τ²]`, `Λ_N,Λ_τ>0`. "Arrows locked to the normal" = the computed locked
  residuals VANISH on the Goldstone: `n−(N₀'/χ₀')η=0` and `τ=0` for `f₀=(χ₀',N₀',0)` (⇒ `f₀` is a flat direction of `E_lock` ⇒ stays
  gapless — Claude-checked).
- Required outputs: `goldstone_residual: Ô_Z2 f₀=0` (+ BC residuals); `m_desc,Z2²=0` for the `u_w` descendant (the pathA_35 overturn);
  `N₀_norm<∞` (explicit `2[tanh(d/ℓ)−tanh³(d/ℓ)/3]/ℓ` for tanh/const); `relative_tilt_eigenvalue>0` (bounded below by `min(Λ_τ/Z_τ)`);
  the pathA_35-style `V_conf` ablation gives `m_desc,conf²=Ω_w,35²>0` and fires `FAIL_PINNED_BRANON`.

### (5) Parity / distinctness computation
- `S_charge±=±Q_E ρ_a ρ_b u_E`; `S_mass/drain=M ρ_a ρ_b u_G+ε_mix S_charge+`; `P_w S_charge±=−S_charge±`, `P_w u_G=+u_G`.
- `q_h(+)=⟨f₀,S_charge+⟩`, `q_h(−)=−q_h(+)` ⇒ `q_h(+)+q_h(−)=0` (neutral `+/−` composite carries no `h` charge). `FAIL_GRAVITY_MIXING`
  fires iff `⟨f₀,S_mass/drain⟩=ε_mix⟨f₀,S_charge+⟩≠0`. **Parity-protected:** pure-even mass has `⟨f₀,S_grav_even⟩=0` by odd·even ⇒
  **charge≠mass is a `w→−w` selection rule** (at the linear/source level; the operator-level nonlinear check is `FAIL_OPERATOR_PARITY_MIXING`,
  sim-deferred). (N.B. the pathA_29 gravity zero mode `1/√d` belongs to a DIFFERENT scalar operator `−d²/dw²`, not to `Ô`; the argument
  needs only odd `f₀` and even mass/gravity source — no cross-operator `f_g` overlap is invoked.)

### (6) `ENGINE_AGREE` + dimensional firewall
- SymPy & Mathematica independently agree on: `Ô`, `K_w`, `K∥`, all BC residuals, parity eigenvalues, `Ô f₀`, `N₀_norm`,
  relative-tilt eigenvalue, pinning-gap shift, `q_eff`, source dims, static Green, dynamic Green, `p_static`, `p_dynamic`, `U_int`
  signs, all `fail_witnesses`, neutral-composite separation.
- **Dim firewall:** units-restored for every term in `E_quad`, `E_lock`, `E_pin`, `SΨ`, `G`, `U_int`; ≥2 able-to-fail dim ablations
  that MUST fire (omit `ρ_b` → source `L⁻³` not `L⁻⁴`; incompatible `K∥` weight → `m_n²` not `L⁻²`; corrupt locking by dropping
  `N₀'/χ₀'`).

---

## §3. The able-to-fail test (FIVE reachable physical FAIL modes + the distinctness guard + 2 sim-deferred risks)

Compute the **inter-throat interaction energy `U_int(R)`** vs 3D brane separation `R`, from the §2 coupled operator + throat source on
the finite-slab + bulk geometry.

- **PASS (`THROAT_ELECTRIC_LOCALIZED_COULOMB`):** `U_int ∼ 1/R` (Coulomb potential → `1/R²` force), arising from the **normalizable,
  gapless embedding Goldstone** with **nonzero throat monopole overlap** `q_eff=⟨f₀,S_throat⟩≠0`, and `±w` giving
  like-repel/unlike-attract.

- **FAIL modes (each physical, each reachable — the gate must be able to land on any of them):**
  1. **`FAIL_YUKAWA`** — the mediator is a **gapped** internal director (the relative-tilt combination), not the locked Goldstone →
     `U_int ∼ e^{−μR}/R`, short-range.
  2. **`FAIL_PINNED_BRANON`** — the wall-translation zero mode is itself **gapped** (pinned wall / inherited `u_w` gap from
     `pathA_35` Gate L) → not massless → Yukawa.
  3. **`FAIL_DELOCALIZED_BULK_1_OVER_R3`** — the zero mode is **non-normalizable** (delocalizes into the 4D bulk) → `U_int ∼ 1/R²`
     (→ `1/R³` force). This is the pathA_29 `p=3` counterfactual, now a reachable verdict.
  4. **`FAIL_NO_MONOPOLE`** — `q_eff=⟨f₀,S_throat⟩=0` (the throat sources only a **dipole**) → no Coulomb term.
  5. **`FAIL_GHOST_INSTABILITY`** — the reduced mediator sector is not stable/positive-definite; a negative-norm/ghost component makes
     `G₀<0` → like-attract or otherwise wrong sign. (The like-repel sign is FORCED on a positive-definite PASS branch — see §2A(3): a
     bare `σ_mouth=−1` flip is a non-physical classifier control, NOT a physical FAIL; the physical wrong-sign route is a ghost.)
  - **Plus 2 SIM-DEFERRED risks (named, not executed here):** `FAIL_OPERATOR_PARITY_MIXING` (the nonlinear throat modifies `Ô` and
    breaks the `w→−w` even/odd split at the operator level → reopens the pathA_35 fifth-force; distinct from the source-level
    `FAIL_GRAVITY_MIXING`) and `FAIL_SOURCE_NOT_COMPACT` (a noncompact throat-funnel body samples continuum modes → a `1/R³` component
    can survive a normalizable zero mode). See §2A(3).

- **⭐ The distinctness guard (the hardest — Codex round-2).** Electric `h` must be **ODD under `w→−w`** (sourced by the puncture
  orientation `±w`); **gravity** is the **EVEN** drain/compression channel (sourced by mass/throughput); **light** is in-plane
  MacCullagh shear at `c_γ`. If ordinary **mass/throat energy** also sources this same massless `h`, the mode is merely a **scalar
  branon / fifth-force partner of gravity** (the `pathA_35` Gate-L `u_w` concern) → `FAIL_GRAVITY_MIXING`. **Clean separation test: a
  neutral `+/−` throat composite must carry ZERO `h` charge but NONZERO gravity drain.** (This is the parity discriminator that keeps
  charge ≠ mass at the field level.)

- **Counterfactual guard (à la pathA_29):** a non-localizing warp / non-normalizable mode MUST flip the exponent to `1/R³` — the test
  genuinely fails on the delocalizing family.

- **Each FAIL mode's exact ablation knob + expected COMPUTED classifier input is tabulated in §2A(3) `fail_witnesses`** — every verdict
  is an OUTPUT of a solve, every FAIL genuinely reachable.
- **Dual-engine (SymPy + Mathematica), `ENGINE_AGREE`** on the full §2A(6) checklist (operator/BC residuals, parity eigenvalues, the
  gap, zero-mode normalizability `q_eff`, static+dynamic exponent `p`, `U_int` sign, all `fail_witnesses`, neutral-composite
  separation) + the §2A(6) dimensional firewall. Each engine assembles its own headline (no `x−x` tautologies); ≥1 able-to-fail
  ablation per decisive claim.

**If PASS:** the electric sector realizes pathA_29's **normalizable-transverse-zero-mode PRINCIPLE** — by a *distinct and stronger*
wall-profile mechanism (§1, N1) — yielding `1/r²` Coulomb with `±w` sign from the throat's 4D body. The forces are then **connected
through one localization principle** (not one identical mechanism), and the charge sector becomes a calibrated entry in the central
`pde_ledger` alongside pathA_29's gravity localization ([[project-calibrated-pde-goal]]).

---

## §4. Honest scope (the algebra/sim boundary)

- **In-scope for algebra:** the *linearized* inter-throat interaction + the *transverse-mode normalizability* analysis. pathA_29 did
  exactly this for gravity with reduced (algebra-tractable) machinery — the charge version is the same class of calculation.
- **Sim-deferred:** the full nonlinear throat interior, the self-consistent throat-body shape, dynamics. Those stay posed as
  sim-dependent open items ([[project-simulation-deferred-complete-pde-strategy]]) — **completing this gate completes the reduced
  far-field localization/sign gate conditional on the declared throat source law; it does NOT complete the nonlinear throat-interior
  theory or prove the full charge sector.**
- **Harder than the flat-brane far-field:** this is throat-interior / 4D / bulk territory (Gate-T neighborhood), so expect the setup
  (§2) to take real work before the gate is executable.

---

## §5. References (read first, next session)
- `docs/conceptual_foundation.md` §3 + §4 ⭐ v7 (the charge = 4D-throat-body mechanism + this `1/r²` crux).
- `reports/pathA_29_*` + `directives/pathA_29_brane_bulk_return.md` (the gravity brane-localization template — the machinery to reuse).
- `directives/pathA_37_c5_throat_electrostatics.md` (the RETIRED flow-gate — the documented exploration that led here; the GLM reviews
  `_scratch/pathA_37_v5_glm_review.md` established the sign result + the flat-3D counterflow no-go).
- Memories: `[[project-brane-existence-defect-structure]]` (UPDATE 2026-07-03b), `[[feedback-native-em-mechanisms]]` (v7),
  `[[project-light-4d-throat-hypothesis]]` (4D throat structure), `[[project-pn-gravity-ladder]]` / `[[project-calibrated-pde-goal]]`.

---

## §6. Changelog
- v0 DRAFT (2026-07-03) — carry-over spec written pre-`/compact` while context is fresh. Nails the crux (`1/r²` via brane-localization,
  the pathA_29 template), lists the §2 setup sub-questions to resolve first, and the §3 able-to-fail test. NOT yet gauntleted.
- v1 (2026-07-03, post-`/compact`) — **§2 mediating field RESOLVED via Claude↔Codex rounds 1–2** (logs in
  `_scratch/pathA_38_field_mechanism_codex_round{1,2}.md`). Decision: the mediator is the **gapless transverse-embedding /
  orientation-lock Goldstone `h`** (wall displacement into `w`, arrows locked to the normal), **not** the free internal `±w` director
  (generically gapped → `FAIL_YUKAWA`, the pathA_36/pathA_25 wall). Roles: A(orientation)=mediator-only-as-locked-Goldstone,
  C(throat funnel)=signed source geometry, B(scalar Coulomb)=derived output after zero-mode projection. Operator = a COUPLED
  fluctuation matrix (Goldstone eigenvector `f₀∝∂_w[profile]`, relative-tilt partner massive → projected out). §3 rewritten to the
  five reachable FAIL modes (`FAIL_YUKAWA`, `FAIL_PINNED_BRANON`, `FAIL_DELOCALIZED_BULK_1_OVER_R3`, `FAIL_NO_MONOPOLE`,
  `FAIL_SIGN_STRUCTURE`) + the `w→−w` parity **distinctness guard** (`FAIL_GRAVITY_MIXING`; neutral `+/−` composite = zero `h` charge,
  nonzero gravity drain).
- v2 (2026-07-03) — **Codex directive design-review r1 = `NOT_SOUND`** (6 blockers: placeholder operator, undefined `S_throat`,
  asserted-not-witnessed FAIL reachability, uncheckable Goldstone/locking, non-computable parity guard, thin dual-engine plan). Folded
  a Codex-authored **§2A reduced-model executable primitives** block (fields `Ψ=(η,n,τ)`, Z₂-wall background, coupled Hessian `Ô`,
  Goldstone `f₀=∂_wΦ₀`, DC-completion BCs, parity `Σ=diag(−1,−1,+1)`; compact-volume `S_throat`; `fail_witnesses` table; locking-energy
  checks; `q_h(+)+q_h(−)=0` parity computation; full `ENGINE_AGREE` + dim firewall). Claude-reviewed: Goldstone norm + `P_w f₀=−f₀`
  hand-checked; flagged CONFIRM-REVIEW item = `σ_mouth`/`u_E` sign + coupling must be DERIVED from the throat mouth BC, not chosen.
- v3 (2026-07-03) — **Codex confirm design-review r2 = `NOT_SOUND`** (2 remaining blockers, down from 6). Folded both: (1)
  **sign/coupling provenance RESOLVED** (Claude↔Codex) — `σ_mouth≡s` (puncture direction, not a free knob); the like-repel/unlike-
  attract SIGN STRUCTURE is DERIVED from `G₀>0`+`s=±1` (`U_int=s₁s₂·q_eff²/N₀·1/4πR`, emitted from the solved `G`); HONEST HEADLINE =
  DERIVED{`p`, sign-structure, charge≠mass parity} + CALIBRATED{charge magnitude `Q_E`, nonzero-monopole modeling; first-principles
  mouth BC sim-deferred}. (2) **delocalizing control fixed** — must use a NONCOMPACT half-line (`w∈[0,∞)`/`ℝ`) anti-localizing measure
  (a finite `[−d,d]` weight keeps norm finite → wouldn't give `p=3`). NEXT: fresh Codex confirm r3.
- v4 (2026-07-03) — **Codex confirm design-review r3 = `SOUND_WITH_CONCERNS`, NO blockers** (both r2 blockers cleared; Codex
  independently reproduced the `U_int=s₁s₂·q_eff²/N₀·1/4πR` sign derivation). Folded the one non-blocking nicety: the
  `FAIL_DELOCALIZED` control now states the explicit threshold `k_w>2/ℓ` (the anti-localizing growth must beat the `sech²` zero-mode
  decay, else norm stays finite). **Codex GREEN.**
- v5 (2026-07-03) — **GLM tertiary pass = `SOUND_WITH_CONCERNS`** (could not break the core; independently reproduced the sign
  derivation + confirmed NOT pathA_25/36 in disguise). Folded all 5 concerns + 4 notes (Codex-authored, Claude-integrated): (1)
  explicit **pathA_35 `u_w`-gap overturn** block + `FAIL_PINNED_BRANON` upgraded to compute the `u_w`-descendant eigenvalue
  (`m_desc,Z2²=0` vs the `V_conf` ablation `Ω_w,35²>0`); (2) `FAIL_SIGN_STRUCTURE`→**`FAIL_GHOST_INSTABILITY`** (physical wrong-sign =
  indefinite Hessian; `σ_mouth=−1` reclassified non-physical); (3) gaplessness reframed as **PREMISE** (Goldstone theorem) not derived;
  (4) sign structure reframed as **GENERAL THEOREM** in the four-way headline (PREMISE/DERIVED/THEOREM/CALIBRATED); (5) added
  **`FAIL_OPERATOR_PARITY_MIXING`** + **`FAIL_SOURCE_NOT_COMPACT`** (sim-deferred); notes: qualified "connected"→shared PRINCIPLE
  (distinct/stronger `sech²`-on-ℝ mechanism), the ℝ-normalizability note, dropped the ambiguous `f_g`. NEXT = final Codex re-confirm.
