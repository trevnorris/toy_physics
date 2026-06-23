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
