# S9b v10 — mechanical authoring lookups

Commands run from `/var/projects/toy_physics`; literal stdout and exit codes follow.
These are text retrieval/difference and version lookups, with no CAS or physical computation.

```bash
git diff --no-ext-diff c2f1cf2b -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
diff --git a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
index 0d59269d..ccaa2fe5 100644
--- a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
+++ b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
@@ -1,7 +1,25 @@
-# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v8)
+# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)
 
-**Authors:** Codex (v5–v6), with v7–v8 edits by Claude (orchestrator). **Status:** v8, cleared for the build
-on 2026-10-06.
+**Authors:** Codex (v5–v6, v10), with v7–v8 edits by Claude (orchestrator). **Status:** v10 authored
+2026-10-08; no v10 review clearance or computed result is claimed.
+
+**Deliverable:** specify the v8 optical objects in Parts A–C and the neutral, linked steady in-plane
+brane balance and its conditional optical requirements in Part D, retaining every unsupplied O2 input.
+
+**Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions
+are `directives/S9b_repair_decision_list.md` at `f069c37c` (D1–D5). v9's Part D is not an input. Paths
+below are relative to `research/pde_ledger_v3/`. The additional source abbreviations are:
+
+- **O2-R:** `steps/O2_steady_brane_balance.md`, accepted at `72866fcf`, especially §§4 and 8.
+- **O2-S:** `directives/O2_SHARED_PHYSICS.md` at `4680e251`.
+- **O2-C:** `directives/O2_input_contract.md` at `217a92e9`.
+- **O2-P:** `directives/O2_premise_decision_list.md` at `77d2c39a`.
+
+`CLAUDE.md` M1–M3 and E1–E2 govern the artifact. Equations labelled **supplied** are conditional
+inputs that this build cannot test. Adopted premises retain that status and their provenance. Every
+dependent result is flagged with the supplied identification or adopted premise it uses. Reference
+objects are comparison inputs, rather than premises for deriving the brane's observables. The spec
+names computed objects and solution conditions; it supplies no outcome or expected-value acceptance test.
 
 ## Why
 
@@ -30,28 +48,56 @@ Each piece is a supplied identification. Flag any result that depends on one.
   Here `μ_⊥` is the basis-invariant transverse stiffness (`steps/S11b_interface_coupling_law.md:74–87`);
   S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
   remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
-  for both polarizations, measured relative to the shear-carrying material.
-- **Advection.** `V(x)` is the in-plane velocity of the material that `u` displaces. Its kinetic symbol is
-  `(ω − V^i k_i)²`: the coefficient of `ω²` is `1`, the coefficient of `ω` is `−2V^i k_i`, and the `kk`
-  slot contains `V^iV^j k_i k_j`.
+  for both polarizations, measured relative to the shear-carrying material; the supplied polarization
+  identification is `c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`.
+- **Advection (supplied).** `V(x)` is the in-plane velocity of the material that `u` displaces. Its
+  supplied kinetic symbol is
+
+  ```
+  K_adv(ω,k;x) ≡ (ω − V^i(x) k_i)² .
+  ```
 - **Steady brane mass balance.** The supplied balance is
 
   ```
   ∇·(ρ_br(x)V(x)) = −j_n(x) .
   ```
 
+  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
+  coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as
+  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. For an induced-measure
+  interpretation, keep `∂_r[(∂ξ_w)²]` live. D4 supplies the optional additional scale condition
+  `∂_r[(∂ξ_w)²] = O(ε/r)`; it is not imposed here. No order for the relative correction to `j_n`
+  is supplied from the slope-amplitude grade alone.
+
   `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
   owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object
   (`directives/S11b_SHARED_PHYSICS.md:99–111`); the relation between `j_n` and `v_dr` is open.
-- **Embedding.** `ξ_w(x)` is the brane's displacement into the bulk direction, as a length. The ledger's `h`
-  is dimensionless, with `ξ_w = ℓh` (`directives/S11b_SHARED_PHYSICS.md:95–96`). `g^{ij}` is the inverse of
-  the induced spatial metric `g_ij` above.
-- **Anchor.** With `V = 0`, `ξ_w = 0` and constant coefficients, the relation reduces to S11b's uniform
-  transverse branch, `ρ_br⁰ω² = μ_⊥k²` (`steps/S11b_interface_coupling_law.md:74–87`).
-- **Background anchoring of `c_γ`.** The steady object in this step is the LAB_HELD profile `c_γ(x)`
-  (`directives/S11c_a_SHARED_PHYSICS.md:232–242`). The S11c-a MATERIAL_ADVECTED anchoring is not used:
-  for a steady radial profile, `D_t c_γ = 0` and `∂_t c_γ = 0` give `V_r ∂_r c_γ = 0`, inconsistent with
-  keeping both a radial drain and a radially varying `c_γ` live.
+- **Embedding (supplied).** `ξ_w(x)` is the brane's displacement into the bulk direction, as a length.
+  The ledger's `h` is dimensionless. The supplied equations are
+
+  ```
+  w = ξ_w(x) ,      ξ_w = ℓh ,      g^{ij} = (g_ij)⁻¹ .
+  ```
+
+  `ℓ` is the fixed reduction scale (`directives/S11b_SHARED_PHYSICS.md:95–96`; O2-C §5), rather than
+  a selected slab width.
+- **Uniform anchor (supplied on its historical domain).** S11b's uniform transverse input is
+
+  ```
+  ρ_br⁰ω² = μ_⊥k² ,      V = 0 ,      ξ_w = 0 ,      ∂_iμ_⊥ = ∂_iρ_br = 0 .
+  ```
+
+  Source: `steps/S11b_interface_coupling_law.md:74–87`. It supplies no live flowing material law.
+- **Background anchoring of `c_γ` (supplied).** This step selects LAB_HELD from the two distinct
+  physical anchorings (`directives/S11c_a_SHARED_PHYSICS.md:232–242`; O2-C §6):
+
+  ```
+  Q_bg^L(x,t) ≡ Q_bg(x) ,      Q_bg^M(x,t) ≡ Q_bg(χ(x,t)) ,
+  c_γ^L(x,t) ≡ c_γ(x) .
+  ```
+
+  `χ` is the inverse material map. MATERIAL_ADVECTED is not selected. LAB_HELD supplies neither a
+  material-reference evolution law nor a physical holder.
 
 **Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
 - direction-dependent (radial versus tangential) stiffness;
@@ -59,22 +105,23 @@ Each piece is a supplied identification. Flag any result that depends on one.
 - polarization-dependent propagation;
 - mixed `ωk` content other than advection.
 
-**Reference speed.** `c₀` is the asymptotic value of `c_γ`, identified with the measured light speed. The bulk
+**Reference speed (supplied).** `c₀ ≡ lim_{r→∞} c_γ(r)`, identified with the measured light speed. The bulk
 sound speed `c_s` is separate and is not identified with `c₀`. Local light speed is not an observable here,
 because rulers and clocks are made of the same medium. Only the far-field observables below are compared.
 
-**Order.**
+**Order (supplied optical retained set and counting).**
 - **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
   coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
   terms of the underlying operator, which affect amplitude and polarization transport.
-- **Smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained set is every monomial
+- **Optical smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained optical set is every monomial
 
   ```
   δ^a (V/c₀)^b ((∂ξ_w)²)^c,      0 ≤ a ≤ 1,  0 ≤ b ≤ 2,  0 ≤ c ≤ 1,
   ```
 
-  with nonnegative integer `a`, `b` and `c`. Each retained monomial is printed separately, and nothing outside
-  this set is computed.
+  with nonnegative integer `a`, `b` and `c`. Each retained optical monomial is printed separately;
+  optical terms outside this set are not computed. This is not a truncation of Part D's mechanical
+  or inherited energy accounting.
 - **Order counting** (supplied; the build cannot test it; owned by the gravity sector or S12). In the far zone,
   let `ε(r) ≡ GM/(c₀² r)`, with `GM` as supplied below. The supplied counting is
 
@@ -82,13 +129,15 @@ because rulers and clocks are made of the same medium. Only the far-field observ
   δ = O(ε),      (∂ξ_w)² = O(ε),      V/c₀ = O(ε^{1/2})   ⇒   (V/c₀)² = O(ε) .
   ```
 
-  Under it, every other monomial in the retained box is `o(ε)`. Within these bounds, `δ`, `V` and `ξ_w` stay
-  free radial profiles. Part B compares at first order in `ε`. For each Part B observable, it sums the
-  contributions of the grades `δ`, `V/c₀`, `(V/c₀)²` and `(∂ξ_w)²` and subtracts the reference. The `V/c₀`
-  contribution is computed and printed at its own order, `ε^{1/2}`, not assumed. Every other retained monomial
-  is printed as its own higher-order object and is not compared with the first-order references.
-- **Profiles.** Prefer general radial functions. If an engine restricts itself to a family, it says so and keeps
-  every exponent and coefficient symbolic.
+  Within these bounds, `δ`, `V` and `ξ_w` stay free radial profiles in Parts A–C. Part B's comparison
+  object is the sum of contributions of the grades `δ`, `V/c₀`, `(V/c₀)²` and `(∂ξ_w)²`, minus the
+  reference. The `V/c₀` contribution is computed and printed at its supplied grade `ε^{1/2}`, with
+  no outcome assumed. Every other retained monomial is printed separately with its order under the
+  supplied counting, and is not compared with the first-order references.
+- **Profiles.** Prefer general radial functions in Parts A–C. If an engine restricts those objects to a
+  family, it says so and keeps every exponent and coefficient symbolic. Part D's unsupplied responses
+  remain general live unknowns, including admissible gradients and material history; no engine-chosen
+  response family or derivative/history cutoff is permitted (D3; O2-R §8).
 - **Inherited limits and freezes.** S9 took the sharp zero-width-sheet limit, no dissipation,
   frequency-independent moduli, the continuum limit and amplitude `→ 0`
   (`steps/S9_light_requires_shear.md:349–350`). It also took two distinct background limits
@@ -104,10 +153,15 @@ because rulers and clocks are made of the same medium. Only the far-field observ
   the live normal-exchange profile `j_n`. The effect of lifting every other inherited limit remains outside
   this step's claim.
 
+  Adopted P1 below revisits S9's no-dissipation and frequency-independent-moduli limits for the
+  steady-load material/reference response (O2-P P1; O2-C §4). Those historical limits do not constrain
+  Part D's relaxation response. Parts A–C retain their supplied optical dispersion and idealizations;
+  compatibility of that optical regime with the material's relaxation remains a later question.
+
 ## The observables
 
-**Quantifier.** Each Part B comparison is required to hold for **every** far-zone `b`, not only at selected
-`b`.
+**Quantifier.** The Part B conditions to compute are the profile conditions for matching the reference
+for **every** far-zone `b`. Part D uses that same quantifier. No matching outcome is supplied.
 
 **Setting.** One isolated, spherically symmetric mass at rest, with its drain flowing: steady, not frozen,
 with `V` live. Far field; linear waves. Time is the lab time of the brane's far-field rest frame.
@@ -126,17 +180,25 @@ trip in the required directions. Outside either condition, report the branch typ
 or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
 them.
 
-## Supplied references (inputs the build cannot test)
+## Reference and bulk input objects
 
 - **`GM`** (owned by the gravity sector). The far-field Newtonian mass parameter, as measured by the orbits of
   slowly moving test matter. It is an independent symbol, not identified with any profile amplitude. Every
   Part B and Part C condition is stated relative to it.
-- **GR reference, in PPN form with `c₀`:**
-  - full-flyby deflection: `Δθ = (1+γ)·2GM/(b c₀²)`;
-  - for `Z_E, Z_R ≫ b`, the round-trip logarithmic term is
-    `2(1+γ)(GM/c₀³)·ln(4 r_E r_R/b²)`. Part B uses only its coefficient of `ln(1/b²)`;
-  - `γ = 1` in GR.
-- **Part C bulk:**
+- **GR reference objects (supplied comparison inputs), in PPN form with `c₀`:**
+
+  ```
+  Δθ_ref(b;γ) ≡ (1+γ)·2GM/(b c₀²) ,
+  Δt_RT,log,ref(b;γ) ≡ 2(1+γ)(GM/c₀³)·ln(4 r_E r_R/b²) ,
+  γ_GR ≡ 1 .
+  ```
+
+  The logarithmic reference has domain `Z_E, Z_R ≫ b`; Part B uses only its coefficient of
+  `ln(1/b²)`. These equations define the oracle objects to compare against. They do not specify
+  either observable computed from the supplied brane dispersion, or an acceptance test.
+  The reference family defines each effective `γ`; GR residuals and matching conditions use its
+  supplied `γ_GR` member.
+- **Part C bulk (supplied):**
   - `P = Kρ^n`, with `n` symbolic, `ρ` the bulk number density and `m` the particle mass (v3 convention). So
     `c_s² = nKρ^(n−1)/m`.
   - Define the fractional bulk-density change at the brane by
@@ -155,9 +217,10 @@ them.
   - An effective `γ` from `Δθ`, and one from the coefficient of `ln(1/b²)` in the round-trip excess time for
     `Z_E, Z_R ≫ b`.
   - Their residuals against the references.
-  - For each observable, the condition on the profiles, relative to `GM`, that sets the first-order sum (see
-    "Order counting") minus its reference to zero. Higher-order monomials are printed but not set against the
-    references.
+  - For each observable, the solution condition on the profiles, relative to `GM`, for the comparison
+    residual of the selected grade sum (see "Order counting") at every far-zone `b`. This is an object
+    to compute, with its domain and unresolved inputs; it is not a demanded residual value. Other
+    retained monomials are printed but not set against the references.
   - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
   - The difference between the two `γ`s, printed as an object.
 - **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
@@ -170,11 +233,244 @@ them.
   print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
   if it enters at all.
 
+- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
+  momentum balance below in force. Print the constructed balance pieces and conditional objects,
+  their premises and remaining OPEN inputs, which of `δ`, `V`, `ρ_br` and `j_n` each condition
+  determines relative to the independent orbital `GM`, which stay free, and the implied `j_n` for
+  each condition. Preserve the every-far-zone-`b` quantifier and print the domain of each condition.
+  No profile, relation among these quantities, or success in matching is supplied.
+
+## Part D inputs: linked steady in-plane brane
+
+**Object and domain (D3; O2-R §8; O2-S §§1, 4–7).** The object is the brane's far-field steady
+in-plane momentum balance: the in-plane component object of O2's `ℬ_hold^live`, conditional on the
+adopted inputs below, with O2's named accounting and force/power qualifications retained. It is a
+neutral-sector restriction of the supplied optical
+setting. The supplied coordinate/radial restrictions are
+
+```
+r ≡ |x| > 0 ,      V^i(x) = V_r(r) x^i/r ,
+Q(x,t) = Q(r)      (steady scalar profiles in this setting).
+```
+
+`V_r`, `ρ_br`, `j_n` and all unsupplied responses remain general live profiles/actions. Eulerian
+steadiness supplies no material constancy or history cutoff. Coordinate momentum storage, transport,
+exchange, force, load, energy and power densities use the same coordinate `d³x` measure as the mass
+law. Native face-area factors and reductions remain explicit through O6. Coordinate `w` content and
+graph-normal content are distinct; the centre graph fixes no finite-thickness face response.
+
+### Adopted equations and their scope
+
+All entries in this subsection are **supplied adopted substrate inputs to a conditional model**,
+rather than derived substrate laws. P1–P4 retain the label and date
+**adopted substrate input to a conditional model (2026-10-06)** from O2-P. The new inputs carry the
+dates and qualifications of D2. Each printed result carries the labels of every entry it uses.
+
+**P1 — material reference (2026-10-06; O2-P P1; O2-C §4; O2-S §§2, 3.2, 6).** The carrier responds
+elastically in the optical shear regime, represented here by the supplied dispersion and
+`c_γ² ≡ μ_⊥/ρ_br`. It relaxes under steady load. The named equations identifying its unsupplied
+relaxation and power inputs are
+
+```
+ℛ_ref/strain^live ≡ OPEN reference/strain evolution and relaxation response ,
+𝒫_ref/relax^live ≡ OPEN power of that response .
+```
+
+These operand identifications supply no kernel, time scale, rate, sign or material-reference law.
+The response keeps general live fields, gradients and material history. Its power is carried
+explicitly, with any net supply accompanied by its physical supplier and budget. Optical
+compatibility with that relaxation response remains a later light question.
+
+**P2 — drain drive (2026-10-06; O2-P P2; O2-S §§2, 4–5).** Adopt O2's supplied representation:
+
+```
+F_drive,separate ≡ 0 .
+```
+
+The drive is dynamical order conversion in conserved material, represented through the material
+stress, face/support loading and O3 exchange entries. Local source/controller functions and
+mouth/collar/return/IR/bulk-boundary data are separate OPEN inventories. Their reaction and supply
+systems stay named. The historical frozen-wall total-mass sink is not substituted for this drive.
+The source-side S14a bridge and response-side S16 interface to `GM` remain OPEN; no `F_drive(GM)`
+law is supplied.
+
+**P3 — exchange momentum (2026-10-06; O2-P P3; O2-S §5).** Use the signed outward-loss convention
+of the coordinate mass law. The supplied local carried-material equation is
+
+```
+(Π_n^carry)^i ≡ j_n V^i ,      i = 1,2,3 .
+```
+
+This fixes carried in-plane momentum only. Additional non-variational S12 partners and the system
+carrying their reaction remain OPEN. Bulk-direction carry retains `𝒥_map`, `𝒩_br^live`, native
+geometry and unsupplied face-to-material velocity identifications. P3 supplies no carried-total-energy
+formula. `j_n`, `V_s`, native bulk velocity and `v_dr` are not identified with one another.
+
+**P4 — bulk loading (2026-10-06; O2-P P4; O2-C §7; O2-S §§2, 3.3, 5).** Retain the postulated
+shear-free scalar bulk and its supplied native face-normal traction restriction:
+
+```
+t_bulk,s^live ≡ 𝒯_bulk,n,s^live n̂_s .
+```
+
+`𝒯_bulk,n,s^live` is a general signed OPEN amplitude, and `n̂_s` is the native face normal.
+There is no independent tangential bulk mechanical stress. Full `T_hold,s`, its support partition,
+native geometry and reduction remain OPEN; the bulk amplitude qualifies its bulk part, rather
+than being a second load beside it. A centre-graph restriction supplies no native-face reduction.
+Mechanical loading, carried momentum and any external support remain separately identifiable.
+
+**P5 — momentum density (2026-10-07; D2).** The supplied brane in-plane momentum-density map is
+
+```
+𝒫_br,inplane^live ≡ ρ_br V .
+```
+
+This supplies the previously OPEN in-plane momentum-density identification. If stressed brane
+material carried additional momentum from its stress, P5 would change. Whether that enters at a
+retained grade is OPEN and owned by S8. P5 supplies no total-energy law, normal-response map or
+independent closure of every O2 transport/history action.
+
+**P6 — steady in-plane pressure (2026-10-08; D2).** In Part D, the supplied steady in-plane stress
+is isotropic pressure, as the P1 material relaxes under steady load. In the Cartesian Cauchy-stress
+convention where stress contracted with a unit normal gives traction, the premise is
+
+```
+𝒯_br,inplane^live,ij ≡ −p_br(ρ_br) δ^{ij} ,
+c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
+```
+
+`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
+This supplies the steady in-plane part of `𝒯_br^live` only. It supplies no normal stress or normal
+material law, relaxation evolution/power law, or optical stiffness law. In the optical regime light
+continues to see `μ_⊥`. No order or value is assigned to `c₀/c_comp`.
+
+**Density link (2026-10-06; D2).** Part D alone adopts the supplied one-exponent stiffness response.
+Writing the proportionality relative to the live asymptotic density and stiffness gives its input
+equations, together with the supplied asymptotic optical identification:
+
+```
+ρ_br⁰ ≡ lim_{r→∞} ρ_br(r) ,      μ_⊥⁰ ≡ lim_{r→∞} μ_⊥(r) ,
+μ_⊥(r)/μ_⊥⁰ ≡ [ρ_br(r)/ρ_br⁰]^α ,
+c₀² ≡ μ_⊥⁰/ρ_br⁰ ,      c_γ(r)² ≡ μ_⊥(r)/ρ_br(r) .
+```
+
+`α` is one live symbol and `ρ_br⁰` stays live. The domain of the symbolic response is printed.
+This supplies `μ_⊥` as a function of `ρ_br` alone in Part D, setting aside its other O1 `ℳ_⊥`
+dependences there; every result using this restriction is flagged. Parts A–C keep the two local
+profiles independent subject to their supplied optical ratio, and Part C's bulk-density route stands.
+The link supplies no `ρ_br(f)` law or separate smallness grade for the brane density or stiffness.
+The engines compute its implication for `δ`; no such implication is stated here.
+
+**w-parity — neutral sector (2026-10-07; D2, D4).** For an electrically neutral mass the adopted
+far-field `w → −w` symmetry restricts the graph displacement and material bulk-direction velocity:
+
+```
+ξ_w(r) ≡ 0 ,      U_material^w(r) ≡ 0 .
+```
+
+Print **neutral-sector restriction** with Part D and each dependent result, including the reduction
+of the supplied `g_ij`. Parts A–C keep `ξ_w` live, including the charged case. This restriction adds
+no force term and supplies no native-face, thickness, exchange-map or normal constitutive law.
+
+### O2 balance pieces that remain live
+
+The engines compose the balance from these physical roles and the supplied equations above.
+No assembled momentum balance, chosen transport tensor, cancellation or profile solution is supplied
+here (D3; O2-S §§3–6; O2-R §8). The author identifies the supplied content as follows:
+
+| O2 content | Supplied content in Part D | Content retained as OPEN |
+|---|---|---|
+| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density. | Remaining flowing/embedded inertia, normal and transport/history actions and their unresolved relations. The conservative antecedent is named, rather than added as another momentum species. |
+| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the steady in-plane pressure part of the full stress. | Normal material response, conservative antecedents, reference evolution and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; they are not independently added forces. |
+| Optical/material constitutive inputs; O1 `ℳ_⊥`, O7 `ℛ_br` | The density link supplies Part D's optical stiffness response only. | Brane-density response and every other unsupplied material or reduction identification. Parts A–C retain the v8 local profiles. |
+| Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric and neutral-sector restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, normal response and material/order/projection map. No second normal equation is imposed. |
+| Mechanical face/support loading; `T_hold,s`, `𝒯_bulk,n,s^live` | P4 supplies the bulk part's native direction only. | Full load/support partition, native normal-load amplitude, projections/maps and application-point velocities. No external support is selected by LAB_HELD. |
+| Exchange; O3 `Π_n`, S12 partners/reactions | P3 supplies the coordinate in-plane carry `j_n V^i`. | Additional momentum partners/reactions and bulk-direction carry with O6/normal material content. Native convective transfer and the same O3 carry are not counted twice. |
+| Boundary/source/core inputs; O5 `ℋ_core`, `𝔅_A13` | P2 supplies the drain-drive representation only. | A13 branch; distinct local conversion/controller and boundary/domain inventories; physical core holder/mouth data and their response. Q2/S22 retain O5; S12 retains conversion and return partners. |
+
+All OPEN entries are general unknown actions, with admissible live fields, gradients, entire material
+history and native dependences (O2-S §3.2; O2-C §1; O2-R §8). A name supplies no closed argument list,
+locality, finite internal-variable set, derivative order, stress split or constitutive family.
+Retain every named operand and complete live-object dependence printed by either O2 engine in the
+content no adopted premise supplies, including PY-only native/chart and core/material-compatibility
+content. Neither engine's finite inventory, nor their union, exhausts admissible dependences.
+Historical homogeneous kinetic, uniform quadratic, static embedding and frozen-profile/face laws
+retain their original domains (O2-S §§7–8; O2-C §§2–7); they supply no additional live law here.
+
+**Paired energy content (P1; O2-S §6; O2-C §8; O2-R §8).** Carry the inherited force/power and
+energy-accounting qualifications with each Part D condition wherever its material/load/exchange
+content requires them. The named OPEN object is `ℬ_E^steady`, retaining `ℰ_br^live`, `𝒥_E^live`,
+`𝒫_ref/relax^live`, `𝒫_convert/exchange^live`, `𝒫_boundary^live`, `𝒮_E,net` and `𝒫_E,supply`.
+P5, P6, the density link and P3 do not supply the total energy, relaxation power, exchanged energy,
+or net supplier/budget. Preserve the unresolved `C_ref` energy-reference/improvement convention
+where applicable. Pair material stress, face loads and generalized/couple actions with the actual
+application-point velocities or generalized rates on compatible measures/maps. Identify shared
+work occurrences without adding them twice; the channel names prescribe no additive split.
+If a conditional closure requires net supply, it carries the physical supplier and budget as OPEN
+inputs. A selected non-passive interface response would separately require its reservoir and stated
+budget; none is selected. Part D closes no energy law or physical supplier.
+
+**Counting (D3; O2-R §8; O2-S §7; O2-C §9).** The mechanical and inherited energy accounting is
+untruncated. The optical grade box removes no material, force, exchange, boundary or power term.
+Individual density/stiffness, inertia/stress/normal response, exchange/source/load, relaxation/power,
+holder/mouth/embedding and derivative grades remain OPEN wherever unsupplied. Only the optical
+stiffness/density ratio inherits the speed-change grade. First order in bulk `f` supplies no `f`–`ε`
+relation. `j_n`, `p_br(ρ_br)`, `c_comp`, `α` and `ρ_br⁰` remain live. No new scale is used to remove
+an OPEN operand or select a profile.
+
+### Returned O2 interpretation and dependence locations
+
+**D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7
+obligation returned to the orchestrator under M1**. S9b does not discharge it; the S9b record carries
+that status. It is distinct from ownership of a physical input. Preserve the accompanying seven
+WL-only first-derivative keys and the complete O2 difference ledger wherever it bears on a condition.
+
+The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
+`mu_perp`, `o2_rho_br_live` at order 1 and `xi_w` at order 2, each at the stored argument
+`sqrt(x1²+x2²+x3²)` (O2-R §4). The names denote `V_r`, `δ`, `f`, `h`, `j_n`, `μ_⊥`, `ρ_br` and
+`ξ_w` respectively. Retain all eight at every following O2 role location in content no adopted
+premise supplies; the location table is provenance, rather than a prescribed CAS representation:
+
+| O2 balance location | Roles carrying the eight keys |
+|---|---|
+| `hold_inplane[0]` | `OPEN_MomentumFlux_0_0`, `OPEN_MomentumFlux_0_1`, `OPEN_MomentumFlux_0_2` |
+| `hold_inplane[1]` | `OPEN_MomentumFlux_1_0`, `OPEN_MomentumFlux_1_1`, `OPEN_MomentumFlux_1_2` |
+| `hold_inplane[2]` | `OPEN_MomentumFlux_2_0`, `OPEN_MomentumFlux_2_1`, `OPEN_MomentumFlux_2_2` |
+| `hold_bulk` | `OPEN_MomentumFlux_3_0`, `OPEN_MomentumFlux_3_1`, `OPEN_MomentumFlux_3_2` |
+| `hold_normal` | All twelve `OPEN_MomentumFlux_i_j` roles, `i=0,1,2,3`, `j=0,1,2` |
+| `energy_balance` | `OPEN_MaterialEnergyFlux_0`, `OPEN_MaterialEnergyFlux_1`, `OPEN_MaterialEnergyFlux_2` |
+
+P5 supplies the in-plane momentum density, P6 the steady in-plane pressure stress, and the density
+link the optical stiffness response, each with the dependence stated in its adopted equation.
+Their use in a constructed flux/action is flagged; they are not blanket replacements of every O2
+momentum-flux or energy-flux action. Keep the general OPEN objects, their full dependence declarations,
+and the labelled neutral-sector restrictions visible with their restricted component objects.
+Neutral centre-graph geometry does not remove unsupplied native or material-history content.
+
+Carry O2-R §4's other differences alongside the conditions they affect: named/native/geometry/map,
+generalized-work and material-compatibility content; the one-sided `OPEN_MaterialCompatibility`
+actions in `energy_storage`, `energy_transport`, `energy_power` and `energy_balance`; one-sided
+material density/flux actions in `energy_power`; the four density/flux orientation differences in
+`energy_balance`; and the density/flux occurrences inside `OPEN_JointPowerAccounting` in
+`energy_power` and `ℬ_E^steady` alongside their separate `energy_balance` occurrences. Retain the
+force/power-compatibility and duplicate-power questions with these differences. No limited inventory
+match supplies full OPEN equality, compatibility, or a duplicate-power conclusion (O2-R §§4–5, 8).
+
+**Register handoff (D5; O2-R §8).** The later record carries the adopted-premise provenance and the
+undischarged interpretation. A register entry follows only from an established sourced requirement
+on which the conditional object rests: material/phase identifications or expressly carried
+accounting/admissibility obligations. Unsupplied closure forms remain OPEN handoffs. An unresolved
+engine difference alone creates no requirement. S21 owns the later integration/sort; this spec
+performs no register edit or physical reconciliation.
+
 ## Engines, review, scope
 
 - **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).
+- **Spec review.** A fresh non-author Claude agent and Grok review v10 until clear, after this
+  authoring stop (D7). The orchestrator runs those reviews.
 - **Build review.** Codex-written, so a fresh Claude agent and Grok, each with a mandatory FORM ablation.
-- **Model point.** As in "Setting": leading eikonal with the retained multigraded set above. Results do not
+- **Optical model point.** As in "Setting": leading eikonal with the retained multigraded set above.
+  Part D's mechanical and inherited energy content retains its stated untruncated domain. Results do not
   transfer to:
   - the strong field;
   - the throat mouth or interior;
@@ -184,10 +480,14 @@ them.
   - anisotropic or coupled branches.
 - **Deferred to the build** (implementation, not new physics; the build directive owns it):
   - the symbolic handling of the every-`b` requirement;
-  - the mass balance on the induced metric. Under the supplied counting, its difference from the flat form is a
-    relative `O(ε)` correction to the implied `j_n`.
+  - component/measure calculus on the supplied coordinate mass law, with the D4 gradient-scale
+    qualification above for any induced-measure interpretation;
+  - representation of general OPEN actions and executable controls, including FORM ablation (E2).
 - **Stop and report**, without choosing, when any of these happens:
   - a second method failure;
   - a sub-problem this spec does not name;
   - a premise this spec does not supply.
 - **The step record** interprets the results.
+
+**Authoring STOP.** Write v10 and report D1–D5 changes against v8, source conflicts and missing
+sourced pieces. No CAS, build, review launch, commit, push or spawned agent is part of this task.
```

Exit code: `0`.

```bash
rg -n '^# S9b|^\*\*Authority|coordinate `d³x`|^  interpretation|^  `∂_r|^  Adopted P1|^\*\*Quantifier|^  `ln|^  The reference family|^  - For each observable|^- \*\*Part [ABCD]|^\*\*Object and domain|^\*\*P[1-6]|^\*\*Density link|^\*\*w-parity|^### O2 balance|^All OPEN|^\*\*Paired energy|^\*\*Counting|^### Returned|^\*\*D5|^The full inherited|^\| `hold_|^\| `energy_balance|^P5 supplies|^Carry O2-R|^\*\*Register|^- \*\*Spec review|^\*\*Authoring STOP' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
1:# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)
9:**Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions
65:  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
68:  interpretation, keep `∂_r[(∂ξ_w)²]` live. D4 supplies the optional additional scale condition
69:  `∂_r[(∂ξ_w)²] = O(ε/r)`; it is not imposed here. No order for the relative correction to `j_n`
156:  Adopted P1 below revisits S9's no-dissipation and frequency-independent-moduli limits for the
163:**Quantifier.** The Part B conditions to compute are the profile conditions for matching the reference
197:  `ln(1/b²)`. These equations define the oracle objects to compare against. They do not specify
199:  The reference family defines each effective `γ`; GR residuals and matching conditions use its
201:- **Part C bulk (supplied):**
211:- **Part A.** The branch-existence and path-traversal conditions. Then, for the LAB_HELD profile, as expressions
216:- **Part B.**
220:  - For each observable, the solution condition on the profiles, relative to `GM`, for the comparison
226:- **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
236:- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
245:**Object and domain (D3; O2-R §8; O2-S §§1, 4–7).** The object is the brane's far-field steady
258:exchange, force, load, energy and power densities use the same coordinate `d³x` measure as the mass
269:**P1 — material reference (2026-10-06; O2-P P1; O2-C §4; O2-S §§2, 3.2, 6).** The carrier responds
284:**P2 — drain drive (2026-10-06; O2-P P2; O2-S §§2, 4–5).** Adopt O2's supplied representation:
297:**P3 — exchange momentum (2026-10-06; O2-P P3; O2-S §5).** Use the signed outward-loss convention
309:**P4 — bulk loading (2026-10-06; O2-P P4; O2-C §7; O2-S §§2, 3.3, 5).** Retain the postulated
322:**P5 — momentum density (2026-10-07; D2).** The supplied brane in-plane momentum-density map is
333:**P6 — steady in-plane pressure (2026-10-08; D2).** In Part D, the supplied steady in-plane stress
347:**Density link (2026-10-06; D2).** Part D alone adopts the supplied one-exponent stiffness response.
364:**w-parity — neutral sector (2026-10-07; D2, D4).** For an electrically neutral mass the adopted
375:### O2 balance pieces that remain live
391:All OPEN entries are general unknown actions, with admissible live fields, gradients, entire material
400:**Paired energy content (P1; O2-S §6; O2-C §8; O2-R §8).** Carry the inherited force/power and
413:**Counting (D3; O2-R §8; O2-S §7; O2-C §9).** The mechanical and inherited energy accounting is
421:### Returned O2 interpretation and dependence locations
423:**D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7
428:The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
436:| `hold_inplane[0]` | `OPEN_MomentumFlux_0_0`, `OPEN_MomentumFlux_0_1`, `OPEN_MomentumFlux_0_2` |
437:| `hold_inplane[1]` | `OPEN_MomentumFlux_1_0`, `OPEN_MomentumFlux_1_1`, `OPEN_MomentumFlux_1_2` |
438:| `hold_inplane[2]` | `OPEN_MomentumFlux_2_0`, `OPEN_MomentumFlux_2_1`, `OPEN_MomentumFlux_2_2` |
439:| `hold_bulk` | `OPEN_MomentumFlux_3_0`, `OPEN_MomentumFlux_3_1`, `OPEN_MomentumFlux_3_2` |
440:| `hold_normal` | All twelve `OPEN_MomentumFlux_i_j` roles, `i=0,1,2,3`, `j=0,1,2` |
441:| `energy_balance` | `OPEN_MaterialEnergyFlux_0`, `OPEN_MaterialEnergyFlux_1`, `OPEN_MaterialEnergyFlux_2` |
443:P5 supplies the in-plane momentum density, P6 the steady in-plane pressure stress, and the density
450:Carry O2-R §4's other differences alongside the conditions they affect: named/native/geometry/map,
459:**Register handoff (D5; O2-R §8).** The later record carries the adopted-premise provenance and the
469:- **Spec review.** A fresh non-author Claude agent and Grok review v10 until clear, after this
492:**Authoring STOP.** Write v10 and report D1–D5 changes against v8, source conflicts and missing
```

Exit code: `0`.

```bash
bash -c 'git show c2f1cf2b:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | nl -ba | sed -n '"'"'91,97p;195,201p'"'"''
```

```text
    91	  every exponent and coefficient symbolic.
    92	- **Inherited limits and freezes.** S9 took the sharp zero-width-sheet limit, no dissipation,
    93	  frequency-independent moduli, the continuum limit and amplitude `→ 0`
    94	  (`steps/S9_light_requires_shear.md:349–350`). It also took two distinct background limits
    95	  (`directives/S9_wl_rebuild_directive.md:359`; `steps/S10_two_transverse_photons.md:816–819`):
    96	  - S9's in-plane background flow `v₀ → 0` removed all convective terms. This step lifts that freeze by
    97	    keeping the in-plane material velocity `V` live.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '"'"'27,41p;63,94p'"'"''
```

```text
    27	| **P5** (2026-10-07) | The brane's momentum density is `ρ_br V`. | If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
    28	| **P6** (2026-10-08) | In steady flow the brane's in-plane stress is an isotropic pressure `p_br(ρ_br)` that depends on `ρ_br` only, because shear relaxes under steady load (P1). `p_br` stays a general function. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. |
    29	| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |
    30	| **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |
    31	
    32	## D3. Part D: the linked brane (new)
    33	**The object.** Part D is the brane's far-field steady in-plane momentum balance. Its ingredients:
    34	- O2's conditional hold balance, from the record's §8 handoff and its sources;
    35	- P3, P5, P6 and w-parity, as supplied;
    36	- v8's mass balance `∇·(ρ_br V) = −j_n`;
    37	- the density link.
    38	
    39	**The deliverable.** Part D re-expresses the Part B conditions with that balance in force. It prints:
    40	- which of `δ`, `V`, `ρ_br` and `j_n` it determines relative to `GM`;
    41	- which of them stay free;
    63	The O2 record returned one interpretation to the orchestrator: the WL-only `ξ_w''` dependence, together with the
    64	seven WL-only first-derivative keys at the same locations.
    65	
    66	**Where they sit.** The record places all eight in the momentum-flux and energy-flux roles (its §4 location table).
    67	
    68	**Routing:**
    69	- **The obligation stays the orchestrator's** (M1; record §8). S9b does not discharge it. The S9b record carries it
    70	  as undischarged. It is distinct from any step's ownership of a physical input.
    71	- **Retained keys.** Part D retains all eight keys, at every §4 location, for every piece of O2 content that no
    72	  adopted premise supplies (record §8).
    73	- **Supplied content.** Where P5, P6 or the density link supplies an object, the premise states that object's
    74	  dependence, and every result that uses it is flagged.
    75	  - The spec author determines, from the O2 sources, which content each premise supplies.
    76	  - The author does not compose the balance (D3).
    77	  - The user chose P6 over keeping the stress OPEN (2026-10-08).
    78	- **Register.** A register entry follows only if a sourced requirement is established. An unresolved difference
    79	  alone does not create one (record's rest-on criterion).
    80	
    81	## D6. Owed build findings: routed to the repair build directive
    82	These come from the build review preserved at `bb94b885`, as listed in that commit's message. They concern implementation, not spec physics. They go to the
    83	repair build directive, which gets its own G2 pass, and the build legs verify them (E2):
    84	- **B1.** The exports bind `w`, `q` and `L` onto unrelated upstream rows (bulk normal coordinate, wave-norm
    85	  coordinate, half-interval size) by bare-symbol equality (`F9B_EQUAL`). The rows are retagged as corroborated S9b
    86	  KNOBs. Both the binding and the corroboration/status promotion are unsupported.
    87	- **B2.** The every-`b` requirement is restated, not reduced. No condition is solved relative to `GM` per stratum,
    88	  and no implied `j_n` is printed (with the `V` sign symbolic).
    89	- **B3.** The radar `ln(1/b²)` coefficient comes from a hand rule, not from Part A's round trip.
    90	- **B4.** Path independence is computed on a placeholder.
    91	- **B5.** The demo refusal is a typed literal.
    92	- **B6.** The SymPy Part C rows conjoin domain predicates built from the pre-response amplitude `d` (seven
    93	  occurrences in the `δ = 0` row).
    94	- **B7.** The tag sets are not parallel across engines.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/O2_premise_decision_list.md | sed -n '"'"'16,27p;41,48p'"'"''
```

```text
    16	1. **Material reference: viscoelastic (the "glacier" choice).** The shear carrier responds elastically in the
    17	   optical shear regime and relaxes under steady load. The user chose this on 2026-10-06 over a separately offered
    18	   option: first resuming the paused retained-reference/relaxed-reference comparison (inventory §2). This
    19	   supersedes the order of that earlier steer.
    20	   - The choice revisits S9's no-dissipation and frequency-independent-moduli limits
    21	     (`steps/S9_light_requires_shear.md:349–350`).
    22	   - The reference and strain evolution law is the relaxation response. It is a general unknown, not an
    23	     engine-chosen family. The input contract (sub-step 2) supplies its form, or keeps it as a named operand.
    24	   - The input contract carries the steady state's energy balance explicitly, including the relaxation
    25	     response's power. Its sign is left to the computation, and any net power names its supplier.
    26	   - The relaxation response's consequences for light in the optical regime are a later light-compatibility
    27	     question. They are not an O2 deliverable.
    41	## Retained OPEN operands
    42	Each stays named, as a general unknown:
    43	- the material branch (A13);
    44	- the live brane stress, inertia and normal material response. Premise 1 fixes its character, not its form.
    45	- the stiffness and density responses `ℳ_⊥` and `ℛ_br`;
    46	- the live embedding relation `ℰ_h^live` (O4). L3's field identity `ξ_w = ℓh` is retained, and O4's identity with
    47	  O2's normal part stays unsettled.
    48	- `T_hold,s`. Its bulk part is face-normal under premise 4. A declared support would be a held input.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/O2_input_contract.md | sed -n '"'"'94,102p;160,177p'"'"''
```

```text
    94	## 3. Live stress, conservative antecedents, inertia and normal response
    95	
    96	**Status of live-stress character: adopted premise 1 (user, 2026-10-06); Label: adopted substrate input
    97	to a conditional model (2026-10-06).** The decision list retains “the live brane stress, inertia and
    98	normal material response. Premise 1 fixes its character, not its form.” (`O2_premise_decision_list.md:44`).
    99	**Form: OPEN `𝒯_br^live`.** The live in-plane and normal stress is that of the premise-1 material,
   100	with the evolving reference/strain state of `ℛ_ref/strain^live` and its material history kept live and
   101	general. Its form, including any non-conservative content, is OPEN; no conservative/dissipative split is
   102	adopted.
   160	## 4. Material-reference and strain evolution
   161	
   162	**Status: adopted premise 1 (user, 2026-10-06); Label: adopted substrate input to a conditional model
   163	(2026-10-06).** The optical shear response is elastic while the material relaxes under steady load.
   164	**Status of its form: OPEN `ℛ_ref/strain^live`.** This single name denotes the unsupplied reference/strain
   165	evolution and relaxation response, including its physical reference carrier, transport, formation or
   166	renewal through conversion/return, and work content where applicable. It has no prescribed argument
   167	list, tensor realization, relaxation kernel, rate, time scale or zero-frequency form. Premise 1 fixes
   168	character, not an equation (`O2_premise_decision_list.md:16–25,44`). §3's OPEN `𝒯_br^live` keeps its
   169	dependence on this evolving reference/strain state and its history.
   170	
   171	**Recorded limits:** S9 took no dissipation and frequency-independent moduli, as well as a sharp sheet,
   172	rest background, continuum and vanishing wave amplitude (`steps/S9_light_requires_shear.md:349–350`).
   173	v9:93–106 lifts in-plane background flow and isotropic background-strain freezes while retaining its
   174	other stated scope boundaries. Premise 1 expressly revisits the no-dissipation/frequency-independent
   175	limits; those limits cannot be imposed on its steady-load relaxation response. The consequences for
   176	optical propagation remain a later light-compatibility question, as required by the decision list
   177	`:20–27`. LAB_HELD speed anchoring is not a reference-evolution law.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/steps/O2_steady_brane_balance.md | sed -n '"'"'283,318p;478,481p;509,515p'"'"''
```

```text
   283	Every nonempty balance `live_arguments` difference consists of the following eight complete keys:
   284	`ProfileDerivative` of **`V_r`, `delta`, `f`, `h`, `j_n`, `mu_perp`, `o2_rho_br_live` at order 1**, and
   285	**`xi_w` at order 2**, all evaluated at the stored argument `sqrt(x1²+x2²+x3²)` and all with signed
   286	PY-minus-WL value **−1** (WL-only). M6 retains the full keys, including derivative order and argument.
   287	Each of those eight keys occurs in **every** role in the following exhaustive location table:
   288	
   289	| Balance component | Every entry key carrying all eight differences |
   290	|---|---|
   291	| `hold_inplane[0]` | `["OPEN_MomentumFlux_0_0"]`, `["OPEN_MomentumFlux_0_1"]`, `["OPEN_MomentumFlux_0_2"]` |
   292	| `hold_inplane[1]` | `["OPEN_MomentumFlux_1_0"]`, `["OPEN_MomentumFlux_1_1"]`, `["OPEN_MomentumFlux_1_2"]` |
   293	| `hold_inplane[2]` | `["OPEN_MomentumFlux_2_0"]`, `["OPEN_MomentumFlux_2_1"]`, `["OPEN_MomentumFlux_2_2"]` |
   294	| `hold_bulk` | `["OPEN_MomentumFlux_3_0"]`, `["OPEN_MomentumFlux_3_1"]`, `["OPEN_MomentumFlux_3_2"]` |
   295	| `hold_normal` | `["OPEN_MomentumFlux_0_0"]`, `["OPEN_MomentumFlux_0_1"]`, `["OPEN_MomentumFlux_0_2"]`, `["OPEN_MomentumFlux_1_0"]`, `["OPEN_MomentumFlux_1_1"]`, `["OPEN_MomentumFlux_1_2"]`, `["OPEN_MomentumFlux_2_0"]`, `["OPEN_MomentumFlux_2_1"]`, `["OPEN_MomentumFlux_2_2"]`, `["OPEN_MomentumFlux_3_0"]`, `["OPEN_MomentumFlux_3_1"]`, `["OPEN_MomentumFlux_3_2"]` |
   296	| `energy_balance` | `["OPEN_MaterialEnergyFlux_0"]`, `["OPEN_MaterialEnergyFlux_1"]`, `["OPEN_MaterialEnergyFlux_2"]` |
   297	
   298	All other comparator-printed balance live-difference arrays are empty **at the unioned balance-role
   299	live-key inventory level**; their
   300	named-operand differences remain (M6). These are OPEN cross-engine content differences. The comparator
   301	acceptance specifically calls the WL-only `xi_w''(r)` a **“cross-engine difference in OPEN content”**,
   302	with filed reviewer measurements `{{1, 43}, {2, 12}}` for WL XiW derivatives versus `{1: 105}` for
   303	PY `OPEN_MomentumFlux_*` occurrences (M1, acceptance's measurement section). These are quotations
   304	of that evidence, not a new count or computation. That source assigns the interpretation:
   305	
   306	> Interpreting it belongs to the O2 record (sub-step 7), under M1: it is preserved, never designed away.
   307	
   308	The record directive bounds this assignment (M1, directive L54–56):
   309	
   310	> Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named stored object.
   311	> A question that needs computation beyond retrieval is listed as open, not computed.
   312	
   313	Retrieval establishes the stored order-2 key, radial argument, signed WL-only difference and all
   314	six balance-component locations above; it does not establish physical equivalence, a transport
   315	interpretation or a force/power consequence. Those questions require more than retrieval.
   316	**The interpretation remains an undischarged sub-step-7 obligation, returned to Claude (orchestrator)
   317	under M1.** No answer is computed or reconciliation attempted. The named operands and other
   318	live-key differences remain OPEN, with no adjudication owner named by their sources.
   478	The mass-law qualification travels explicitly: `∇·(ρ_br V)=−j_n` is a **coordinate-`d³x` measure**
   479	input. For a claim reading `j_n` or `ρ_br` per induced measure, or comparing an induced-metric mass
   480	law, the recorded qualification is a **relative `O(ε)` correction to `j_n`**, not a live O6 law;
   481	O2 supplies and derives **no induced-metric mass balance** (spec L330–336). Fixed-`ℓ` counting,
   509	In particular, O2 supplies **no momentum-density map `ρ_br V`**. Its open inertia/momentum mapping
   510	cannot be cancelled or converted into a chosen advective kinetic law using the mass relation.
   511	The directive's folded disposition keeps the separate user premise for S9b repair at that later
   512	decision-list gate (M1). Part D must source any additional kinetic, constitutive, projection,
   513	profile, support, order or energy input through its own authorized premise/specification work.
   514	It inherits no `GM` normalization, drain-to-mouth identification, static flowing extension or
   515	light-observable compatibility from O2.
```

Exit code: `0`.

```bash
git rev-parse f069c37c:research/pde_ledger_v3/directives/S9b_repair_decision_list.md HEAD:research/pde_ledger_v3/directives/S9b_repair_decision_list.md 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md HEAD:research/pde_ledger_v3/steps/O2_steady_brane_balance.md 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md HEAD:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md 217a92e9:research/pde_ledger_v3/directives/O2_input_contract.md HEAD:research/pde_ledger_v3/directives/O2_input_contract.md 77d2c39a:research/pde_ledger_v3/directives/O2_premise_decision_list.md HEAD:research/pde_ledger_v3/directives/O2_premise_decision_list.md
```

```text
ab3047b6544580c9554523406873965e3b462272
ab3047b6544580c9554523406873965e3b462272
84745153f47d4fe2c07871af708be5aa4dd2658e
84745153f47d4fe2c07871af708be5aa4dd2658e
4ba41e1fc548c61d9f043fd37fbc6673df58c531
4ba41e1fc548c61d9f043fd37fbc6673df58c531
2b5f663e22fa25335e01c4ac66b60049366d74e3
2b5f663e22fa25335e01c4ac66b60049366d74e3
7d8a2885c2b3a500924974680231719896e48eec
7d8a2885c2b3a500924974680231719896e48eec
```

Exit code: `0`.

```bash
git diff --check -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
```

Exit code: `0`.
