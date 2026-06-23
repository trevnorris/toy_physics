# Decision 15 — The EM sector is medium-native: light as BRANE SHEAR (physical picture)

**Date:** 2026-06-22.
**Status:** CONCEPTUAL FOUNDATION (physical reasoning, user + Claude, private-chat unification mode). This is the "what we
believe is physically happening" record for the EM re-founding — a HYPOTHESIS to be tested by derivation, not a result.
Pairs with directive `pathA_23_em_medium_native.md` (which now needs REVISION to match this picture — see §8).
**Falsification stance:** this is a test ([[feedback-falsification-is-the-goal]]); the brane-shear picture can break (see
§7), and a clean break is a valid, valuable result — never to be rescued.

---

## 0. What this doc is
The conceptual/physical foundation behind the EM re-founding. The math has not been done; this records the physical picture
we reasoned out so the next phase (formalization + Codex derivations) is aimed correctly. The journey: the GR-quadrupole
verdict needed `λγ = c_γ/c_s`; chasing it exposed that the canonical EM sector had **drifted** away from the model's
founding concept; correcting the drift led, through an honest falsification check, to a specific new physical picture for
what light *is* in this model.

## 1. The chain that forced the re-founding
1. The verdict needs `λγ`. β-status: genuine gap. EM-anchor scoping: only a **speed-cone** observable isolates `λγ`.
2. `c_γ(ρ0)` derivation (`reports/pathA_cgamma_of_rho_derivation.md`): the canonical action is **Type-4 DECOUPLED** — its
   EM is a *fundamental unified gauge field on a flat metric*, separate from the medium, so `c_γ` is `ρ0`-independent while
   `c_s ∝ ρ0²`; the cones can't lock.
3. That decoupling is the fingerprint of a **conceptual drift** ([[project-single-medium-concept-vs-math-drift]]): standard
   EM machinery (unified `F_MN` on its own metric) crept in and made the model *dualistic* (medium + a separate field),
   contradicting the single-medium concept the model was built on.
4. Decision: re-found EM **medium-native** (directive `pathA_23`; design-reviewed `SOUND-AS-IS`; GLM tertiary).
5. GLM tertiary (`_scratch/glm_pathA_23_tertiary_review.md`) found a **near-theorem**: a single-component SCALAR superfluid
   cannot carry transverse light, period. User accepted this for a *3D scalar* medium and asked to re-examine the medium.

## 2. Why a (scalar) superfluid cannot carry light — precisely
- **Scalar order parameter → spin-0 only.** `Ψ = √ρ e^{iθ}` is one complex scalar. Its excitations are `δρ` and `δθ`, both
  spin-0/longitudinal; one broken `U(1)` → one Goldstone phonon (Goldstone theorem / Son–Wingate / Watanabe–Brauner).
  A scalar's spectrum cannot contain a spin-1 photon.
- **Fluid → no shear modulus → no transverse waves.** `c_transverse = √(μ_shear/ρ)`, and a fluid has `μ_shear = 0`.
- **Vortices are topological/non-radiative** (Kelvin waves: `ω ∝ k²`, not `c_s k`; need a pre-existing core).
- The model's OWN paper admits this: `em_fields.tex:909` — irrotational bulk gives `B=0` outside cores and "does not yet
  support ... transverse electromagnetic waves in the vacuum, as these require non-zero vorticity or shear degrees of
  freedom that are absent in a simple irrotational bulk."
- **The GNLS has NONE of the needed properties** (it *is* a single complex scalar). Confirmed.
- This is dimension-independent for a pure scalar: KK reduction of a scalar gives a tower of scalars, never a vector.

## 3. The load-bearing insight: "no sideways push" is REQUIRED, not a defect
- In the model, magnetism = **Magnus force** = a circulation interacting with the bulk flow, giving a push perpendicular to
  motion (`qv×B`). It comes out right **because the bulk is an ideal shear-free fluid.**
- If the bulk had genuine shear rigidity, an extra transverse *elastic* stress would ride on top of the Magnus force and the
  calibrated magnetic force would no longer match. **So the absence of bulk shear is a precondition for the magnetism that
  already works.** (User's finding — repeatedly confirmed by pushing the math; promoted here to a hard constraint.)
- Consequence — separate two things that were tangled:
  - **(a) EM fields & forces** (Coulomb, Gauss, Biot–Savart, Lorentz/Magnus): WORK, and work *because* it's a pure fluid.
  - **(b) Radiation / light:** the only place the scalar obstruction bites. A pure scalar fluid radiates **longitudinally**
    (sound); real light is **transverse**. That is the irreducible mismatch.

## 4. The resolution direction: substructure → tension/elasticity, located on SURFACES not in the bulk
- The vacuum has **substructure** (its own constituents held by their own cohesive forces — user's hypothesis).
- Cohesion → **tension** (at surfaces) + **transient/elastic stiffness** (viscoelasticity: fluid at low frequency, elastic
  at high frequency). Either way it provides a **transverse restoring force WITHOUT bulk shear** — exactly the loophole §3 demands.
- Everyday anchor (user): the **oscillating jet** (a stream's cross-section oscillating at 90°, held coherent by surface
  tension) — surface tension is a transverse restoring force and a coherence (anti-spreading) mechanism, with zero bulk shear.
- Key refinement: a free-surface *capillary* wave is dispersive (`ω²∝k³`, wrong for light); a **taut elastic sheet/membrane**
  carries *non-dispersive* transverse waves (`c=√(T/ρ)` or `√(μ/ρ)`) — the right kind.

## 5. THE PICTURE: the brane as an elastic drumhead; LIGHT = brane in-plane SHEAR (user's gut, unpacked)
- Our 3D space = **the brane** = a taut **elastic** sheet embedded in the 4D bulk.
- A real elastic membrane has THREE kinds of motion: **bending/flexural** (off-brane displacement into `w` — one *scalar*
  mode), **in-plane longitudinal**, and **in-plane SHEAR**.
- **The bending mode alone is a scalar (1 DOF) — NOT the photon.** The photon needs 2 transverse polarizations.
- **LIGHT = the brane's in-plane SHEAR waves.** In a 3D elastic medium, shear waves carry exactly **two transverse
  polarizations** — the photon's polarization count — propagating *within* our 3D space, **non-dispersive at a single speed**
  `c = √(μ_brane/ρ_brane)` (ideal elastic sheet). The "sideways sliding within the skin," not the up-down of the drumhead.
- **Critical:** the shear rigidity lives **ON THE BRANE, not in the bulk.** The bulk stays a pure shear-free fluid →
  **magnetism (§3) is untouched.** Physically normal: a soap film / lipid membrane has real surface elasticity while the
  fluid it sits in stays perfect (an interface may be stiff in ways the bulk is not). Brane tension (coherence) and brane
  shear (light) come from the **same** cohesive substructure — one story.

## 6. The three-sector division (the payoff of this picture)
| sector | what it is | character | status |
|---|---|---|---|
| **Gravity** | the bulk fluid's flow `v_r` (medium draining into defects) | longitudinal / scalar (`l=2` quadrupole) | works |
| **Magnetism** | the bulk's swirl (Magnus force), **no bulk shear** | bulk, shear-free | works (and protected) |
| **Light** | the **brane's in-plane shear** waves | transverse, 2 polarizations, on the brane | the new hypothesis |
Each lives in its proper place; the bulk is a pure fluid throughout, so gravity and magnetism never feel the brane's stiffness.

## 7. The open checks / make-or-break (this is where it can die)
1. **No leak into the bulk (hit this FIRST).** Prove brane elasticity does NOT feed into the bulk and perturb the
   Magnus/magnetic force. This is the most plausible failure mode; if it leaks, magnetism breaks and the picture is falsified.
   **STAGE-1 RESULT (2026-06-23, `LEAK_CONDITIONS_DEFERRED`, tri-reviewed):** the kinematic audit derived the bulk-stress
   projection onto the *bending* brane and found the interface traction `T_na = T_wa + (T_ww δ_ab − T_ab)∂_b u_w` carries a
   **generically transverse part**. No-leak holds ONLY if, at the brane, `P_T T_wa = 0` AND `T_ab` is isotropic — i.e. no
   normal-tangent momentum flux and no anisotropic in-plane stress. **Physical read: these conditions are NOT generic near a
   draining defect** (`T_wa = mρ v_w v_a` is exactly fluid turning from in-plane flow into the throat; radial in-plane flow makes
   `T_ab` anisotropic). So the leak is **expected to be present**; survival rides on the Stage-3 **throat solve** — either the
   leak *magnitude* is too small to spoil the calibrated Magnus force, or a **non-fine-tuned** profile makes the transverse
   *projection* (not `T_wa` itself) cancel. The `v_n → bulk-vorticity` feedback (feeds the Magnus reservoir directly) is the top
   Stage-3 priority. **Do not bank on no-leak.** (Report: `reports/pathA_23_stage1_kinematic_leak_audit.md`.)
2. **It must be MAXWELL, not just "a transverse wave."** The brane-shear field must reproduce EM *structure* — gauge freedom,
   `E`/`B`, the field equations, coupling to charges — not merely two modes at `c`.
3. **Universal, non-dispersive `c`.** An ideal elastic sheet gives it; verify the brane's actual constitutive law does too,
   and that all "colors" share one speed.
4. **Polarization count, honestly.** Confirm the in-plane shear genuinely yields 2 transverse polarizations on the 3-brane
   (the relocated make-or-break; better-placed than "can a scalar carry light at all").

### 7a. The "leftovers" (a feature, per precedent — user's curiosity)
The elastic brane has extra modes beyond light: the **flexural/bending scalar** `u_w` and the **in-brane longitudinal** mode
(and their relation to bulk sound). Do NOT discard them — catalog them. Precedent: `S_leak` was an accidental leftover that
turned out to help later. Treat unexpected leftover modes as candidate features, not nuisances, and see what they map to.
**BUT (GLM, FIX 3 — a hard constraint, not a curiosity):** a leftover is only a "feature" if it is **gapped or confined**. A
**massless, unsuppressed `u_w` propagating in our 3D space is phenomenologically EXCLUDED** — it would mediate a long-range
scalar **fifth force**. So `u_w` (and any residual in-brane longitudinal mode) must be shown gapped / confined-to-`w` /
constrained; a massless one is a `FAIL` (`FAIL_BENDING_MASSLESS_FIFTH_FORCE`), tested in Stage 4, not deferred to the leftovers
catalog.

## 8. Relationship to directive pathA_23 (it needs REVISION)
`pathA_23` v2 was framed to derive medium-native `E`/`B` from the **GNLS bulk flow** (E = compression, B = vorticity of the
scalar). The reasoning above concludes that route **cannot** give transverse light (its Stage-2 make-or-break would FAIL as
written — which is consistent, not contradictory: the bulk scalar genuinely can't do it). The new picture **relocates light
to BRANE elasticity.** So the directive must be re-pointed: from "derive `E`/`B` from the scalar bulk" to **"model the brane
as an elastic membrane; derive light = brane in-plane shear; prove it reproduces Maxwell AND leaves the bulk magnetism
intact."** The staged scaffold (action-fork + DOF audit → map → polarization/dispersion → Lorentz/E↔B → sources/charge →
payoff) largely carries over, re-aimed at the **brane-elastic action**.

## 9. Deferred questions (parked — do NOT lose)
- **Why does fluid flow INTO the mouth but then LEAK BACK into the brane?** (the inflow-vs-brane-leak puzzle.) Explicitly
  deferred by the user (2026-06-22). Park for after the EM re-founding.

## 10. NEXT (post-/compact)
Identify where all the math goes and prep the Codex derivation runs:
1. Write down a **brane-elastic action** (the brane as a taut elastic membrane in the 4D bulk; its in-plane shear sector).
2. Derive the **shear-wave spectrum** → confirm 2 transverse polarizations, non-dispersive, at a universal `c`.
3. **Maxwell-structure check** — does brane shear reproduce gauge freedom + `E`/`B` + the field equations + charge coupling?
4. **No-leak proof** — brane elasticity must not perturb the bulk Magnus/magnetic force (§7.1).
5. **Catalog the leftover modes** (§7a).
Then: revise `pathA_23` to this brane-elastic framing, design-review, and hand the derivations to Codex.

---

## 11. The known template: MacCullagh ROTATIONAL elasticity (the constitutive crux) — added 2026-06-22
The brane-shear hypothesis is not starting from zero. The exact structure we need was built ~190 years ago and is a real,
citable result — this both sharpens the make-or-break to a *single* question and resolves the scariest open check for free.

- **MacCullagh (1839)** constructed an elastic ether whose **potential energy depends only on the *rotation* (curl) of the
  displacement**, `U = ½ μ (∇×u)²` — *not* on compression `∇·u` or general strain. It reproduced reflection, refraction, and
  polarization across crystal optics. **FitzGerald (1878/80), Kelvin (1892), Larmor (1894)** then showed it maps directly onto
  Maxwell: identify the displacement with one EM field and its curl with the dual field and the kinetic/potential energies
  coincide with Maxwell's. *An elastic medium of the rotational type maps onto electromagnetism — in that historical template.*
  (This is a template to TEST, not a proof for our brane action: curl-only potential energy alone does not by itself prove EM
  gauge redundancy — the kinetic term, source coupling, constraints, and boundary conditions can still leave physical
  longitudinal/zero modes. Whether the same map holds for OUR brane action+coupling+sources+constraints is exactly what the
  derivation must establish, not assume.)
  (Refs: Darrigol, "James MacCullagh's ether: An optical route to Maxwell's equations?"; Whittaker, *History of the Theories
  of Aether and Electricity*, ch. 5, en.wikisource.org.)
- **Kelvin's "gyrostatic adynamic ether" (1890)** is the MECHANICAL realization: a lattice of **spinning gyrostats** has
  *rotational* rigidity but **no translational rigidity**. That is precisely the user's substructure hypothesis — **spinning
  constituents → rotational elasticity → light** — anticipated by 130 years. (Ref: Kelvin, "On a Gyrostatic Adynamic
  Constitution for Ether," 1890.)
- **Modern corroboration:** **Cosserat / micropolar** elasticity (constituents carrying independent micro-rotation DOF)
  admits a vector-potential representation of the stress with EM-like gauge redundancy `δA^i_μ = ∂_μ α^i`, i.e. a
  Maxwell-type structure. (Refs: arXiv:1908.06984 "duality between Cosserat elasticity and fractons"; arXiv:physics/0202006
  "Elasticity and electromagnetism. Part 1"; arXiv:2206.02473 micropolar survey.)

**THE CRUX, sharpened.** A *generic* (Cauchy/Navier) elastic solid, energy `½λ(∇·u)² + μ(strain)²`, carries a **longitudinal
wave** too (speed `√((λ+2μ)/ρ)`, *faster* than the transverse) — an unwanted "longitudinal photon," which would be fatal. The
**rotational (MacCullagh) form** has potential energy `∝ (∇×u)²` only, so the longitudinal part of `u` (`u → u + ∇χ`) costs no
*potential* energy. In the historical template this makes the longitudinal sector gauge and delivers, together: (1) two
transverse polarizations, (2) no longitudinal photon, (3) gauge invariance, (4) Maxwell. **CAUTION (do not assume):** that the
*same* collapse holds for our brane theory is a HYPOTHESIS — the kinetic term, sources, constraints, and boundary conditions
must be checked to leave no physical longitudinal/zero mode (this is what Stages 3–5 of the directive test, with reachable
`FAIL_CURL_ONLY_NOT_GAUGE` / `FAIL_RESIDUAL_LONGITUDINAL_ZERO_MODE` outcomes). Modulo that proof, §7's four checks collapse
toward ONE:

> **Does the brane's substructure produce *rotational* (curl-only / MacCullagh) elasticity rather than ordinary (Cauchy)
> elasticity?** (Candidate mechanism: Cosserat/micropolar — oriented/spinning constituents — i.e. Kelvin's gyrostats.)

**Honest caution (do not assume away).** MacCullagh's ether was long dismissed as "a mathematical invention without a
mechanical analog" — its stress tensor is *antisymmetric*/unusual — until Kelvin's gyrostatic model. That unusual constitutive
form is a genuine subtlety to **confront and justify from the substructure**, not to assume. Also: which field maps to which
(displacement↔`A` vs displacement↔`B`-type; `E↔∂_t u` vs `curl u`) is a **dual convention to be pinned by the derivation**, not
asserted here.

**The SPECIFIC MacCullagh gauge obstruction (GLM, C5 — a named, likely failure, not a vague risk).** The curl-only *potential*
energy `½μ(∇×u)²` is invariant under `u→u+∇χ`, but the **kinetic** term `½ρ(∂_t u)²` is NOT invariant for time-dependent `χ`.
So `u→u+∇χ` is **not** a symmetry of the *full* action — only of the potential energy. The EOM then gives `∂_t²(∇·u)=0`: the
longitudinal mode is a **constrained zero mode, not a gauge artifact.** In Maxwell, by contrast, `A→A+∇χ` IS a full-action
symmetry because the **scalar potential** `φ→φ-∂_tχ` compensates in the kinetic term (`E=-∂_tA-∇φ`). The MacCullagh ether has
**no `φ`**, so its gauge structure does *not* match Maxwell as-is. **Expected resolution to TEST (Stage 5):** either (a) the
brane theory supplies a scalar-potential analog `φ` (recovering full Maxwell), or (b) a constraint eliminates the longitudinal
mode. The directive must state which it expects and test it — `FAIL_CURL_ONLY_NOT_GAUGE` is a *specific likely* outcome here,
not just one of many.

## 12. Where the math goes (routing → directive stages) — added 2026-06-22
**Engine:** heavily symbolic ⇒ **Mathematica leads** (constitutive tensors, principal symbols, polarization eigenvectors,
dispersion relations, the MacCullagh↔Maxwell correspondence), with **SymPy cross-checks** where a second engine can verify
([[feedback-dual-engine-required]]).

**Fields.**
- *Bulk* (KEEP — works): scalar `Ψ → ρ, θ`; flow `v` (gravity = `v_r`) + swirl (magnetism, Magnus). Shear-free.
- *Brane* (NEW): displacement `u(t,x,y,z)` = in-plane `u_∥` (3 comps; the **2 transverse = light**) + off-brane `u_w`
  (bending/flexural scalar — a leftover, §7a). Plus the **brane↔bulk coupling** term (the no-leak object).

**Ordered derivations (→ `pathA_23` stages):**
1. Brane-elastic action + full field content + DOF audit + **honest classification** (NEW_PARENT_ACTION) + candidate
   constitutive forms (Cauchy vs rotational/MacCullagh vs Cosserat).
2. **NO-LEAK (hit FIRST — the most likely death):** the brane sector must not source bulk shear / perturb the Magnus
   (magnetic) or gravity flow. FAIL here = concept-fatal.
3. **Constitutive form (the crux, §11):** does substructure → rotational (curl-only) elasticity, or Cauchy? Cosserat mechanism.
4. **Spectrum / polarization:** 2 transverse, non-dispersive, universal `c`; longitudinal = gauge (rotational) vs stray
   wave (Cauchy).
5. **Maxwell dictionary + gauge + Lorentz / E↔B mixing.**
6. **Sources / Gauss / charge coupling:** how a bulk drain/defect couples to the brane displacement to act as a charge.
7. **Leftovers catalog (§7a) + payoff (λγ, see §13) + bookkeeping + honest scope.**

## 13. The λγ subtlety the brane picture INTRODUCES (important — do not lose) — added 2026-06-22
Relocating light to the brane changes the `λγ = c_γ/c_s` story. `c_γ = √(μ_brane/ρ_brane)` is now a **brane** property, while
`c_s` (bulk sound / gravity-ripple) is a **bulk** property. So `λγ = 1` is **no longer automatic "same medium, same ripple"** —
it requires the brane shear modulus and the bulk compressibility to **both descend from the same cohesive substructure and
lock**. The empirical GW170817 cone-lock (`c_GW ≈ c_light`) thus becomes, in this model, a **derived relation between two
distinct elastic sectors**, not a triviality. This may make `λγ = 1` **harder** to obtain than the naive single-medium hope —
which is exactly the kind of thing worth surfacing now (the user's reason for chasing `λγ` early: "hammering this out now will
reveal something we'll need to know"). The payoff stage must **derive — or honestly fail to derive — this brane↔bulk speed
lock**, and label its provenance. **(GLM, C7 — manage expectations):** in ordinary elasticity shear and compression speeds
differ generically (`√(μ/ρ)` vs `√((K+4μ/3)/ρ)`), so **`CONE_LOCK_ABSENT` / `CONE_SCALING_ONLY` is the EXPECTED outcome**;
`CONE_LOCK_DERIVED` would require a non-obvious mechanism tying brane shear to bulk compression — if it appears, scrutinize it
HARDER (a suspiciously clean lock), do not celebrate it.

## 14. Honest conceptual costs surfaced by the GLM tertiary (2026-06-22 — do not lose)
The brane-shear re-founding is the most promising EM route, but the tertiary review made four costs explicit. They are not
blockers (the user's ethos: let the math falsify; surface costs, never hide them — [[feedback-falsification-is-the-goal]]), but
they must be carried honestly into every verdict.

1. **Brane elasticity is NOT a consequence of the GP/NLS mean-field (C1).** A GP/NLS superfluid is a *fluid* — zero shear
   modulus. The brane's shear rigidity therefore requires the **substructure** (constituents/cohesion *beneath* the mean-field),
   which the GP/NLS equation does not contain. Honest framing: **"GP/NLS as the effective/coarse-grained medium + a deeper
   substructure that supplies the brane elasticity."** This is still the user's single-substance picture (the substructure is
   the one stuff; GP/NLS is its coarse behavior), but it is a real cost to the simplest "one mean-field equation" reading.
   `FAIL_NO_TRANSVERSE_STIFFNESS` is the *most likely* Stage-2 failure ([[project-single-medium-concept-vs-math-drift]]).
2. **A postulated constitutive form ⇒ a CONDITIONAL verdict (FIX 1, the deepest point).** If the curl-only (rotational) form is
   *postulated* rather than *derived from an independently-motivated microstructure*, then a full stage-pass means only
   **"Maxwell structure follows from a POSTULATED rotational elasticity"** — it does **NOT** establish that *the medium*
   produces rotational elasticity. That weaker result is **not, by itself, sufficient** to reclassify `λγ` in `decisions/14` or
   to integrate into the papers without the user's explicit acceptance of the conditional status. Moreover, a microstructure
   that is reverse-engineered *to give* the MacCullagh form (e.g. "assume spinning gyrostats") is **circular** — the derivation
   must start from the medium's *independently-motivated* properties.
3. **The construction is bulk + brane + coupling (C9).** "Our 3D space = an elastic brane" is *assumed* (`NEW_PARENT_ACTION`),
   not derived from the medium; and the model now has three pieces (bulk GP/NLS, brane elasticity, brane↔bulk coupling). Whether
   that is honestly "one medium" is a conceptual cost the final verdict must acknowledge.
4. **The treatment is CLASSICAL; quantization is deferred (C10).** The brane elastic action, constitutive laws, and wave
   spectra are classical, while the bulk is a quantum superfluid. The classical scope must be stated and quantization flagged
   as a deferred question.

## 15. Stage-2 constitutive-crux RESULT (2026-06-23) — the medium does NOT derive the brane shear; user elected the CONDITIONAL path
The constitutive crux (`pathA_23` Stage 2, track B) was executed, caught tautological, reworked, and tri-reviewed. **Verified
result: `FAIL_UNSPECIFIED_SUBSTRUCTURE`.** (`reports/pathA_23_stage2_constitutive_form.md`.)

- **The honest finding.** With the **full symmetric first-gradient energy** kept — `U_∥ = ½K_br θ² + μ_br e⟨ab⟩e⟨ab⟩`, the
  deviatoric shear invariant `e⟨ab⟩` PRESENT with a free coefficient — and μ_br routed through a genuine able-to-fail
  substructure classifier, the medium's record **does not determine μ_br**. Simple-fluid facts give μ_br=0; coherent-network
  facts give μ_br>0; the actual record lands UNDETERMINED because *persistent in-plane neighbor memory* and an *affine/network
  free energy* are never independently specified (consistent with §14 C1: the GP/NLS mean-field is a fluid, shear must come from
  the substructure; and §14 C9: the substructure is not fixed by GP/NLS). **So brane-shear EM is NOT derivable from the current
  single-medium specification.**
- **The tautology that was caught (lesson).** The first attempt returned `FAIL_NO_TRANSVERSE_STIFFNESS` by selecting `W=W(det F)
  =W(J)` (a fluid EOS) *at the input* — `μ_shear=0` was then a definitional identity, a can't-fail gate. The fidelity audit
  confirmed the code was clean *as an implementation*; the adversarial review caught that it implemented a pre-ordained setup and
  tested the wrong object (the mean-field, not the substructure). The rework de-tautologized it (deviatoric invariant present;
  classifier demonstrably lands ZERO/NONZERO/UNDETERMINED on different inputs). Same discipline that caught Stage 1.
- **The trilemma / pincer (verified algebraically — sharpens §11's CRUX).** Spectrum of the symmetric law is
  `(λ−μ_br k²)²(λ−(K_br+4μ_br/3)k²)`:
  - **μ_br=0** (fluid) → eigenvalues `0,0,K_br k²` → **no transverse light.**
  - **μ_br>0** (ordinary/Cauchy symmetric shear; needs NO directors/spin → *not* reverse-engineering) → two transverse modes
    `μ_br k²` **but also** a longitudinal mode `(K_br+4μ_br/3)k²` = a stray **"second photon"** → `FAIL_CAUCHY_STRAY_LONGITUDINAL`
    at Stage 4.
  - **MacCullagh curl-only** `½μ_R(∇×u)²` (the *only* form with clean transverse-only light, charpoly `λ(λ−μ_R k²)²`) → first-
    gradient stress is **antisymmetric** (no spin-couple closure without a gyrostat sector = reverse-engineering) and carries the
    **C5 gauge obstruction** (the longitudinal null is a constrained zero mode, not gauge, until a φ-analog/constraint is supplied).
  ⇒ the medium cannot produce clean transverse-only light without an **extra postulated ingredient**.
- **Consequence for the value map.** The EM sector / `λγ` (`c_γ=√(μ_br/ρ)` or `√(μ_R/ρ)`) is **not derivable medium-natively** at
  the current specification — `λγ` remains a **genuine free input** (consistent with `BETA_GENUINE_GAP`). (Do NOT yet edit
  `decisions/14`: the gate rule holds — a postulated form ⇒ CONDITIONAL, no value-map/paper updates without explicit acceptance.)
- **USER DECISION (2026-06-23) = Path 1 (proceed CONDITIONALLY).** Adopt the **rotational/MacCullagh form as an explicit
  POSTULATE**, flag the verdict CONDITIONAL throughout, and run Stages 3–6 to test whether the rest holds (no-leak; the C5 gauge
  test; charge firewall; the λγ cone payoff). Still able-to-fail at every stage. The gyrostat substructure stays an acknowledged
  gap, not a claim. **NEXT = Stage 3** (constitutive no-leak closure — decides the Stage-1 leak).

## 16. Stage-3 leak → the CURVATURE-LOCALIZED LEAK + GEON-THROAT hypothesis (2026-06-23, user)
Stage 3 (constitutive no-leak closure) returned `LEAK_BOUNDED_CONDITIONAL` (tri-reviewed: re-run 36/36 + `FIDELITY_CLEAN`; the
adversarial review flagged the verdict as too soft). The adversarial concern: a **defect-independent, intrinsic-to-light** channel —
the brane's own rotational stress `D_b σ^R_ba = −μ_R k² u_a` — would drive bulk vorticity for *every* photon, making D1
(no-leak/Magnus) effectively fatal unless a 4th postulate (brane→bulk impedance) is added. Adjudicating that produced a sharper,
falsifiable picture (user's insight + the ideal-fluid free-slip property).

**(A) The leak is CURVATURE-LOCALIZED, and the adversarial "intrinsic fatal" reading is an over-count.**
- An **ideal (shear-free) bulk transmits tangential traction ONLY convectively**: `T_wa = ρ v_w v_a` (the pressure part `P δ_wa` is
  purely normal ⇒ **free slip** — the very "no sideways push in the bulk" that makes Magnus work, §3). So the bulk can receive a
  sideways force only by a normal (w) flow *advecting* in-plane momentum.
- **Light is in-plane shear with `v_w = 0`** ⇒ on a **flat brane it transmits zero tangential traction ⇒ NO leak ⇒ Magnus/EM clean
  in the far field.** The brane's `σ^R` is its *internal* restoring force; by the action's variational structure it sources the
  *brane* EOM, not the *bulk* EOM (vary bulk fields ρ,θ,v → only S_ψ + S_cpl appear) — so folding `D_b σ^R_ba` into the bulk source
  is the over-count.
- The light↔`u_w` **mixing rides on the brane slope `s_a = D_a u_w`** (the `(T_ww δ−T_ab) D_b u_w` term in `T_na`): ≈0 between
  defects, **O(1) at a throat** where the brane bends hard into 4D and the small-slope expansion breaks. ⇒ `ε_leak ∝ (curvature)²`
  — a **derived, far-field-vanishing** suppression, NOT the unmotivated "assume ε_leak≪1" the adversarial flagged. Bonus: residual
  curvature everywhere (= the gravity field) ⇒ tiny light–bulk coupling ∝ local curvature = light bending in gravity.

**(B) The GEON-THROAT hypothesis (user — the bigger structural bet).** The throat of a defect = a **trapped brane-shear standing
wave** ("photons standing between 3D and 4D — beyond the mouth, before the bulk"), held by the throat geometry. Ancestor: Wheeler's
**geon** (1955). Trapped wave → quantized k → rest energy = the defect's **mass**; circulation → **charge/spin**. Ties the EM sector
(light) to the gravity sector (drain/throat) AT the defect — a real unification step, and a candidate answer to the parked §9
inflow-vs-back-leak puzzle. See memory `[[project-geon-throat-hypothesis]]`.

**Status / falsifiable ([[feedback-falsification-is-the-goal]]):** the free-space half (σ^R not a bulk source; convective-only
tangential traction; `v_w(light)=0 ⇒ T_wa(light)=0`; leak ∝ slope) is being **verified now** (focused Stage-3b re-derivation,
dual-engine + tri-review, able-to-fail in BOTH directions). If it holds, the fatal reading is retired and Stage 3 = a
curvature-localized (defect-only, throat-deferred) `LEAK_BOUNDED_CONDITIONAL`. The **make-or-break** for the geon throat is
**Stage-4 `u_w` confinement** (a massless unconfined `u_w` = fifth force = death) + showing the throat-localized strong mixing does
NOT reach the far-field Magnus. Do NOT bank it; let the math falsify it.

**Stage-3b OUTCOME (2026-06-23, tri-reviewed — re-run 30/30·29/29, `FIDELITY_CLEAN`, adversarial `RESTS_ON_MODELING_CHOICE_NEEDS_GLM`
+ `STAGE3B_OVERCLAIMS`).** Verdict `OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED` (`SIGMA_R_NOT_A_BULK_SOURCE`, `LIGHT_FREE_SLIPS_NO_LEAK`).
The over-count is **real and verified GIVEN the declared (separate-sector) action** — `σ^R` is a brane-internal force, not a bulk
source; flat density-preserving transverse light free-slips. **BUT** (adversarial, fair): the result is **model-contingent** — it
rests on the brane field `u` being a separate sector from the bulk fluid `v` (no in-plane `u̇^a↔v^a` link). Membrane-in-fluid → no
leak; strict "brane *is* the fluid" → `u̇^a` is a bulk `δv^a` → leak returns. AND the real (defect/throat) leak is **relocated, not
retired** (still open, deferred to a throat solve; the `|K|L_mix` scaling admittedly breaks at the throat). Claude's physics read:
the separation is **probably legitimate** — inviscid (shear-free) bulk ⇒ **free slip** ⇒ no tangential `u̇^a=v^a` matching (no-slip
holds only for a *viscous* fluid), so brane shear (`u`, neighbor-preserving elastic) and bulk flow (`v`, mass transport) are two
**distinct emergent DOF of one substructure**, correctly decoupled tangentially — the same free-slip that makes Magnus work. **§17
below upgrades this from "probably" to a mechanism** (the brane is a domain wall = a genuine emergent structure ⇒ separation is
single-medium-legitimate by construction). Report framing to be downgraded ("over-count corrected, defect/throat leak open"); GLM Q1
likely subsumed by §17.

## 17. WHY does the brane exist at all? Brane = a DOMAIN WALL / phase interface (2026-06-23, user — foundational)
The most-imposed, least-derived part of the whole construction. Currently the brane is **put in by hand**: the confinement potential
`V_conf(X;Σ)`, the `w`-profiles `Z(w)`/`W(w)`/`B_ℓ(w)`, and the `k_w u_w²` restoring term. None explains *why* there is a sheet at
`w=0` or what stops us flowing into the bulk.

- **User's mechanism (surface tension):** the brane = a **domain wall / phase interface** in the medium. A field with two degenerate
  vacua forms a codim-1 surface carrying **energy/area = surface tension** (the water-surface intuition); it is **topologically
  stable** (interpolates two vacua, can't unwind); we are **confined as bound zero modes** in the wall's potential well (leaving the
  wall costs energy). Grounding: **Rubakov–Shaposhnikov 1983** "Do we live inside a domain wall?"; **Volovik ³He** A–B phase interface.
- **The honest obstacle:** the current potential `U(ρ)=(K/4)ρ⁵` is **single-well** ⇒ no domain wall from the scalar alone ⇒ deriving
  the brane needs **more medium structure** (degenerate vacua / a 2nd component / self-trapping by the collective drain network).
  This is the **same "substructure"** §14 C1 said we need for shear, and the same gap Stage 2 hit.
- **THE DOOR — the convergence:** brane existence, the Stage-2 shear crux, the §14 substructure, and the Stage-3b separate-sector
  justification are **ALL one question: what is the wall's internal structure?**
  - **Structureless wall** (simple scalar) = a *fluid membrane*: tension + bending but **NO in-plane shear** = exactly the Stage-2
    failure (`U_∥=½K(∇·u)²`, μ_shear=0). Gives confinement but not light.
  - **Structured wall** (internal order: a director / broken symmetry across the wall / textured interface) **CAN carry in-plane
    shear** ⇒ that internal order **IS** the substructure §14 C1 needed ⇒ could **derive** the shear law and **retire the postulated
    MacCullagh (Stage-2 CONDITIONAL → derivation)**.
  - A wall as a genuine emergent topological object is the **single-medium-legitimate justification for the separate-sector no-leak**
    (answers Stage-3b/GLM Q1): one substance *organized into* a stable wall with its own collective surface DOF — not a 2nd medium,
    not "just the fluid."
  So deriving "why the brane exists" could yield **confinement + shear (light) + the legitimate separation in one move.**

## 18. The DEFECT = a PUNCTURE through the brane: CHARGE / MASS / THROAT unification (2026-06-23, user — the prize)
Extends §17 to tie the **electric** and **mass** sectors together at the defect. A defect/particle = a **puncture (throat) through the
brane** into the bulk; the brane's **surface tension wants to restore (close)** the puncture → stored deformation energy.

- **CHARGE = the topological PUNCTURE itself** (that it punctures at all + its winding/orientation). **Quantized and mass-independent
  because topological** → explains **charge universality** (same `|e|` despite wildly different mass). This is a **mechanism for the
  EXISTING `η_Q e_*` ontology** (`pde.tex:279-312`) and respects the **charge firewall** (`η_Q`, not circulation) — NOT a new free knob.
- **MASS = the geon standing-wave content** holding the throat **open** (trapped brane-shear; wave energy = rest mass, §16B).
  Continuous/variable.
- **DECOUPLING (user, load-bearing assumption):** charge ⊥ mass *at this level* — the standing-wave photons (mass) do **not** directly
  interact with the electric (charge/puncture) part. ("Maybe they interact at some point, but not here.") ⇒ the EM-charge sector and
  the mass/geon sector are **separable**; this is what makes charge mass-independent.
- **ELECTRIC FIELD ENERGY = the brane-tension deformation energy around the puncture.** The tension trying to close the hole **↔ the
  electric charge**; the stored field energy is the brane's restoring energy.
- **THROAT SIZE/STRUCTURE = the energy BALANCE:** brane-tension-closing (∝ charge) **vs** standing-wave/drain holding-open (∝
  mass/energy). Sets a **finite throat radius** → candidate **resolution of the classical point-charge self-energy divergence**
  (finite `a` ⇒ finite Coulomb self-energy `~q²/a`). Necessary to model the particle defect at a fundamental level.
- **THE PRIZE:** a candidate **fundamental model of the charged massive particle** — a punctured, tense domain-wall brane with a
  trapped standing wave: **charge from topology, mass from wave energy, size from the tension-vs-holding-open balance** — tying EM
  (charge=tension, light=shear) to gravity/mass (throat=drain, mass=geon) through the brane's existence.
- **Falsifiable tests (do NOT bank — [[feedback-falsification-is-the-goal]]):** (1) does the puncture's tension-energy reproduce the
  Coulomb `1/r²` field + `~q²/a` self-energy? (2) does the balance give a finite, sensible throat size? (3) does charge come out
  quantized/universal (mass-independent)? (4) does the charge⊥mass decoupling hold or break? (5) does the tie to brane moduli feed
  `λγ` (`c_γ²=μ_brane/ρ_brane`) and the verdict? Any of these can fail and falsify (this construction of) the picture.

- **CALIBRATION RUNWAY (user, 2026-06-23 — why this matters).** Because charge is universal (one `e`) and the brane tension `σ` is one
  universal number, `tension↔charge` is **one equation** → **calibrate `σ` on the measured `e`/`α`**, then predict the **held-out
  surplus** against experiment (classical electron radius `r_e`; `λγ=c_γ/c_s≈1` cone-lock → closes the verdict gap; lepton mass ratios;
  `g−2`; stretch = `α` from brane/bulk moduli). **NON-NEGOTIABLE first gate (the Gate-4 lesson):** the `tension→charge` map must derive
  **constant-free** (fewer free constants than downstream predictions) — else calibrating `σ` on `e` merely ABSORBS `e` (zero surplus =
  the `g_G`-absorbs-`G` trap), which we report plainly. The cascade (one `σ` → many predictions) further requires §17's wall structure to
  **relate the moduli** (`σ`, `μ_brane`→`c_γ`, `ρ_brane`). This is calibrate-predict ([[feedback-calibrate-predict-methodology]]) at the
  EM coupling — the deepest open gap. Folded into `pathA_24` Phase B (B6).

**NEXT PHASE = directive `pathA_24` (brane existence + defect structure) — the next many steps.** See `directives/pathA_24_brane_existence_defect_structure.md` (DRAFT → design-review → GLM → execute). Memory: `[[project-brane-existence-defect-structure]]`.
