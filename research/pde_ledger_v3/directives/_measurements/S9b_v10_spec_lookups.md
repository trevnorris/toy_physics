# S9b v10 — mechanical authoring lookups

## Original authoring: preserved baseline `af1674e5`

The following lookups predate amendment 1. The repair 1 lookups below use the amended authority.

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

## Repair 1: baseline `af1674e5`, amendment 1 at `79055918`

The current user instruction routes historical O2-emission retention checks to the comparator and
record; engine-facing obligations use the explicit inputs in the repaired spec. The commands below
retrieve the amended inputs, unchanged O2 sources, repaired text and literal text/version differences.
No CAS or physical computation is performed.

```bash
git diff --no-ext-diff af1674e5 -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
diff --git a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
index ccaa2fe5..015d369c 100644
--- a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
+++ b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
@@ -1,13 +1,14 @@
 # S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)
 
-**Authors:** Codex (v5–v6, v10), with v7–v8 edits by Claude (orchestrator). **Status:** v10 authored
-2026-10-08; no v10 review clearance or computed result is claimed.
+**Authors:** Codex (v5–v6, v10), with v7–v8 edits by Claude (orchestrator). **Status:** v10 repair 1,
+2026-10-08, against the preserved, not accepted baseline `af1674e5`; no repair clearance or computed
+result is claimed.
 
 **Deliverable:** specify the v8 optical objects in Parts A–C and the neutral, linked steady in-plane
 brane balance and its conditional optical requirements in Part D, retaining every unsupplied O2 input.
 
 **Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions
-are `directives/S9b_repair_decision_list.md` at `f069c37c` (D1–D5). v9's Part D is not an input. Paths
+are `directives/S9b_repair_decision_list.md` at `79055918` (amendment 1). v9's Part D is not an input. Paths
 below are relative to `research/pde_ledger_v3/`. The additional source abbreviations are:
 
 - **O2-R:** `steps/O2_steady_brane_balance.md`, accepted at `72866fcf`, especially §§4 and 8.
@@ -64,10 +65,12 @@ Each piece is a supplied identification. Flag any result that depends on one.
 
   The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
   coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as
-  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. For an induced-measure
-  interpretation, keep `∂_r[(∂ξ_w)²]` live. D4 supplies the optional additional scale condition
-  `∂_r[(∂ξ_w)²] = O(ε/r)`; it is not imposed here. No order for the relative correction to `j_n`
-  is supplied from the slope-amplitude grade alone.
+  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. Any induced-measure claim
+  names its reading: the same densities re-expressed per induced volume, or a mass law imposed on
+  the induced measure. It keeps `∂_r[(∂ξ_w)²]` live and attaches no relative order to the correction
+  to `j_n` unless it states a condition bounding that correction relative to `j_n` itself (amended
+  D4). No such bound is supplied here. Neither the supplied slope counting nor the historical
+  gradient-scale condition `∂_r[(∂ξ_w)²] = O(ε/r)` supplies that bound. The latter is not imposed.
 
   `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
   owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object
@@ -237,7 +240,9 @@ them.
   momentum balance below in force. Print the constructed balance pieces and conditional objects,
   their premises and remaining OPEN inputs, which of `δ`, `V`, `ρ_br` and `j_n` each condition
   determines relative to the independent orbital `GM`, which stay free, and the implied `j_n` for
-  each condition. Preserve the every-far-zone-`b` quantifier and print the domain of each condition.
+  each condition. For each condition, also print what the OPEN pieces left in the in-plane balance
+  must supply, each with its general live dependence (amended D3). Preserve the every-far-zone-`b`
+  quantifier and print the domain of each condition.
   No profile, relation among these quantities, or success in matching is supplied.
 
 ## Part D inputs: linked steady in-plane brane
@@ -319,30 +324,44 @@ native geometry and reduction remain OPEN; the bulk amplitude qualifies its bulk
 than being a second load beside it. A centre-graph restriction supplies no native-face reduction.
 Mechanical loading, carried momentum and any external support remain separately identifiable.
 
-**P5 — momentum density (2026-10-07; D2).** The supplied brane in-plane momentum-density map is
+**P5 — momentum density and carriage (2026-10-07; carriage 2026-10-08; amended D2).** Scoped to
+Part D, this **adopted substrate input to a conditional model** supplies the brane's in-plane
+momentum density and its carriage with the brane material at `V`:
 
 ```
-𝒫_br,inplane^live ≡ ρ_br V .
+𝒫_br,inplane^live ≡ ρ_br V ,
+𝒥_br,carry^ij ≡ ρ_br V^i V^j .
 ```
 
-This supplies the previously OPEN in-plane momentum-density identification. If stressed brane
-material carried additional momentum from its stress, P5 would change. Whether that enters at a
-retained grade is OPEN and owned by S8. P5 supplies no total-energy law, normal-response map or
-independent closure of every O2 transport/history action.
+`𝒥_br,carry^ij` names P5's contribution to the in-plane momentum current. Any other in-plane
+momentum current stays a general OPEN action, counted beside this supplied current and P6's stress;
+it contains neither again. Each Part D condition prints what the remaining OPEN pieces must supply.
+The material current is distinct from P3's outward exchange momentum. If stressed brane material
+carried additional momentum from its stress, P5 would change. Whether that enters at a retained
+grade is OPEN and owned by S8. P5 supplies no total-energy law or normal-response map.
 
-**P6 — steady in-plane pressure (2026-10-08; D2).** In Part D, the supplied steady in-plane stress
-is isotropic pressure, as the P1 material relaxes under steady load. In the Cartesian Cauchy-stress
-convention where stress contracted with a unit normal gives traction, the premise is
+**P6 — adopted steady in-plane stress (2026-10-08, revised the same day; amended D2).** Scoped to
+Part D, this **adopted substrate input to a conditional model** supplies isotropic pressure plus
+linear viscous stress in the Cartesian Cauchy convention:
 
 ```
-𝒯_br,inplane^live,ij ≡ −p_br(ρ_br) δ^{ij} ,
+𝒯_br,inplane^live,ij ≡ T^{ij} ,
+T^{ij} ≡ −p_br(ρ_br) δ^{ij}
+         + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k)
+         + ζ δ^{ij} ∂_kV^k ,
+t_br,inplane^i ≡ T^{ij} n_j ,      F_br,inplane^i ≡ ∂_j T^{ij} ,
 c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
 ```
 
 `p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
-This supplies the steady in-plane part of `𝒯_br^live` only. It supplies no normal stress or normal
-material law, relaxation evolution/power law, or optical stiffness law. In the optical regime light
-continues to see `μ_⊥`. No order or value is assigned to `c₀/c_comp`.
+The shear viscosity `η` and bulk viscosity `ζ` are live general profiles, including their gradients.
+P6 supplies the steady in-plane part of `𝒯_br^live` only. It is an adopted steady form, not a
+consequence of P1; P1's relaxation response remains OPEN outside this form. The linear viscous
+form excludes power-law creep. P6 supplies no normal stress or normal material law, relaxation
+evolution law, or optical stiffness law. In the optical regime light continues to see `μ_⊥`.
+The power this stress dissipates has no supplied identification with an O2 energy operand. It
+stays in the OPEN energy accounting, with its physical supplier and budget OPEN. No order or value
+is assigned to `c₀/c_comp`, `η` or `ζ`.
 
 **Density link (2026-10-06; D2).** Part D alone adopts the supplied one-exponent stiffness response.
 Writing the proportionality relative to the live asymptotic density and stiffness gives its input
@@ -375,13 +394,13 @@ no force term and supplies no native-face, thickness, exchange-map or normal con
 ### O2 balance pieces that remain live
 
 The engines compose the balance from these physical roles and the supplied equations above.
-No assembled momentum balance, chosen transport tensor, cancellation or profile solution is supplied
+No assembled momentum balance, total transport tensor, cancellation or profile solution is supplied
 here (D3; O2-S §§3–6; O2-R §8). The author identifies the supplied content as follows:
 
 | O2 content | Supplied content in Part D | Content retained as OPEN |
 |---|---|---|
-| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density. | Remaining flowing/embedded inertia, normal and transport/history actions and their unresolved relations. The conservative antecedent is named, rather than added as another momentum species. |
-| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the steady in-plane pressure part of the full stress. | Normal material response, conservative antecedents, reference evolution and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; they are not independently added forces. |
+| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density and the material-carried contribution `ρ_br V^i V^j` to the in-plane current. | Any other in-plane momentum current, with general live dependence, counted beside the supplied current and P6 stress and containing neither again; normal/embedded responses and unresolved relations wherever unsupplied. The conservative antecedent is named, rather than added as another momentum species. |
+| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the full steady in-plane stress in its adopted pressure-plus-linear-viscous form, with the stated Cauchy traction/force convention. | Normal material response, conservative antecedents, reference evolution outside the adopted form, and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; no additional in-plane stress duplicates P6. |
 | Optical/material constitutive inputs; O1 `ℳ_⊥`, O7 `ℛ_br` | The density link supplies Part D's optical stiffness response only. | Brane-density response and every other unsupplied material or reduction identification. Parts A–C retain the v8 local profiles. |
 | Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric and neutral-sector restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, normal response and material/order/projection map. No second normal equation is imposed. |
 | Mechanical face/support loading; `T_hold,s`, `𝒯_bulk,n,s^live` | P4 supplies the bulk part's native direction only. | Full load/support partition, native normal-load amplitude, projections/maps and application-point velocities. No external support is selected by LAB_HELD. |
@@ -391,9 +410,11 @@ here (D3; O2-S §§3–6; O2-R §8). The author identifies the supplied content
 All OPEN entries are general unknown actions, with admissible live fields, gradients, entire material
 history and native dependences (O2-S §3.2; O2-C §1; O2-R §8). A name supplies no closed argument list,
 locality, finite internal-variable set, derivative order, stress split or constitutive family.
-Retain every named operand and complete live-object dependence printed by either O2 engine in the
-content no adopted premise supplies, including PY-only native/chart and core/material-compatibility
-content. Neither engine's finite inventory, nor their union, exhausts admissible dependences.
+Each engine retains the named operands and complete live-object dependences supplied in this spec
+in content no adopted premise supplies. This includes native/chart, core/material compatibility,
+the derivative keys and component locations below, and general admissible dependences outside those
+explicit keys. No finite dependence inventory exhausts the general OPEN actions. Checking retention
+against either O2 engine's historical emissions belongs to the comparator and record, as stated below.
 Historical homogeneous kinetic, uniform quadratic, static embedding and frozen-profile/face laws
 retain their original domains (O2-S §§7–8; O2-C §§2–7); they supply no additional live law here.
 
@@ -402,7 +423,9 @@ energy-accounting qualifications with each Part D condition wherever its materia
 content requires them. The named OPEN object is `ℬ_E^steady`, retaining `ℰ_br^live`, `𝒥_E^live`,
 `𝒫_ref/relax^live`, `𝒫_convert/exchange^live`, `𝒫_boundary^live`, `𝒮_E,net` and `𝒫_E,supply`.
 P5, P6, the density link and P3 do not supply the total energy, relaxation power, exchanged energy,
-or net supplier/budget. Preserve the unresolved `C_ref` energy-reference/improvement convention
+or net supplier/budget. The power dissipated by P6's adopted stress remains in this OPEN accounting,
+with no supplied identification with an O2 energy operand or physical supplier/budget. Preserve
+the unresolved `C_ref` energy-reference/improvement convention
 where applicable. Pair material stress, face loads and generalized/couple actions with the actual
 application-point velocities or generalized rates on compatible measures/maps. Identify shared
 work occurrences without adding them twice; the channel names prescribe no additive split.
@@ -415,21 +438,24 @@ untruncated. The optical grade box removes no material, force, exchange, boundar
 Individual density/stiffness, inertia/stress/normal response, exchange/source/load, relaxation/power,
 holder/mouth/embedding and derivative grades remain OPEN wherever unsupplied. Only the optical
 stiffness/density ratio inherits the speed-change grade. First order in bulk `f` supplies no `f`–`ε`
-relation. `j_n`, `p_br(ρ_br)`, `c_comp`, `α` and `ρ_br⁰` remain live. No new scale is used to remove
-an OPEN operand or select a profile.
+relation. `j_n`, `p_br(ρ_br)`, `c_comp`, `η`, `ζ`, `α` and `ρ_br⁰` remain live. No order or value is
+assigned to `c₀/c_comp`, `η` or `ζ`. No new scale is used to remove an OPEN operand or select a profile.
 
 ### Returned O2 interpretation and dependence locations
 
 **D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7
 obligation returned to the orchestrator under M1**. S9b does not discharge it; the S9b record carries
-that status. It is distinct from ownership of a physical input. Preserve the accompanying seven
-WL-only first-derivative keys and the complete O2 difference ledger wherever it bears on a condition.
+that status. It is distinct from ownership of a physical input. The engine-facing dependence inputs
+are stated explicitly below; preservation of the historical difference ledger is the comparator/record
+handoff at the end of this subsection.
 
 The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
 `mu_perp`, `o2_rho_br_live` at order 1 and `xi_w` at order 2, each at the stored argument
 `sqrt(x1²+x2²+x3²)` (O2-R §4). The names denote `V_r`, `δ`, `f`, `h`, `j_n`, `μ_⊥`, `ρ_br` and
-`ξ_w` respectively. Retain all eight at every following O2 role location in content no adopted
-premise supplies; the location table is provenance, rather than a prescribed CAS representation:
+`ξ_w` respectively. Each engine retains all eight as admissible live dependences at every following
+component/role location in content no adopted premise supplies. These are complete dependence
+inputs in this spec, rather than a request to inspect an O2 emission or reproduce its serialization.
+Equivalent component/action notation is permitted, with the roles identifiable:
 
 | O2 balance location | Roles carrying the eight keys |
 |---|---|
@@ -440,21 +466,39 @@ premise supplies; the location table is provenance, rather than a prescribed CAS
 | `hold_normal` | All twelve `OPEN_MomentumFlux_i_j` roles, `i=0,1,2,3`, `j=0,1,2` |
 | `energy_balance` | `OPEN_MaterialEnergyFlux_0`, `OPEN_MaterialEnergyFlux_1`, `OPEN_MaterialEnergyFlux_2` |
 
-P5 supplies the in-plane momentum density, P6 the steady in-plane pressure stress, and the density
-link the optical stiffness response, each with the dependence stated in its adopted equation.
-Their use in a constructed flux/action is flagged; they are not blanket replacements of every O2
-momentum-flux or energy-flux action. Keep the general OPEN objects, their full dependence declarations,
-and the labelled neutral-sector restrictions visible with their restricted component objects.
+P5 supplies the in-plane momentum density and its material-carried current contribution; P6 supplies
+the steady in-plane pressure-plus-linear-viscous stress; the density link supplies the optical
+stiffness response, each with the dependence stated in its adopted equations. Other in-plane
+momentum-current content stays OPEN beside the supplied current and P6 stress, containing neither
+again. Their use in a constructed action is flagged; they do not replace unsupplied momentum-current
+or energy-flux content. Keep the general OPEN objects, their full dependence declarations and the
+labelled neutral-sector restrictions visible with their restricted component objects.
 Neutral centre-graph geometry does not remove unsupplied native or material-history content.
 
-Carry O2-R §4's other differences alongside the conditions they affect: named/native/geometry/map,
-generalized-work and material-compatibility content; the one-sided `OPEN_MaterialCompatibility`
-actions in `energy_storage`, `energy_transport`, `energy_power` and `energy_balance`; one-sided
-material density/flux actions in `energy_power`; the four density/flux orientation differences in
-`energy_balance`; and the density/flux occurrences inside `OPEN_JointPowerAccounting` in
-`energy_power` and `ℬ_E^steady` alongside their separate `energy_balance` occurrences. Retain the
-force/power-compatibility and duplicate-power questions with these differences. No limited inventory
-match supplies full OPEN equality, compatibility, or a duplicate-power conclusion (O2-R §§4–5, 8).
+**Additional engine-facing OPEN content.** Retain named/native geometry, chart/measure/map,
+core/material compatibility and generalized/rotational work as general OPEN actions wherever they
+bear on a condition. The explicit material-compatibility role is `OPEN_MaterialCompatibility` in
+`energy_storage`, `energy_transport`, `energy_power` and `energy_balance`. Material energy-density
+and energy-flux roles are `OPEN_MaterialEnergyDensity` and `OPEN_MaterialEnergyFlux_i`, `i=0,1,2`;
+their content remains OPEN in energy storage/transport and power accounting. The named joint action
+`OPEN_JointPowerAccounting` keeps force/power compatibility and the unresolved identification of
+shared work/energy content visible, without choosing a decomposition or counting the same work twice.
+These role names denote the existing material/energy operands above, not additional independently
+additive energy species or a prescribed nesting/occurrence count.
+
+**Comparator and record handoff (O2-R §§4–5, 8).** The comparator and record check that the Part D
+content no adopted premise supplies retains every named operand and complete live-object dependence
+printed by either O2 engine, including all eight derivative keys at every listed location and the
+PY-only native/chart and core/material-compatibility content. They carry the complete O2 difference
+ledger wherever it bears on a condition. Each engine constructs from the explicit inputs in this
+spec alone; it is not charged with checking historical O2 emissions it was not given.
+The handoff includes the one-sided material-compatibility actions in all four named energy rows;
+one-sided material density/flux actions in `energy_power`; four density/flux orientation differences
+in `energy_balance`; and density/flux occurrences inside `OPEN_JointPowerAccounting` in
+`energy_power` and `ℬ_E^steady` alongside separate `energy_balance` occurrences. Those differences
+travel with the force/power-compatibility and duplicate-power questions. No limited inventory match
+supplies full OPEN equality, compatibility or a duplicate-power conclusion. This check discharges
+neither the returned interpretation nor any unsupplied physical input.
 
 **Register handoff (D5; O2-R §8).** The later record carries the adopted-premise provenance and the
 undischarged interpretation. A register entry follows only from an established sourced requirement
@@ -470,7 +514,10 @@ performs no register edit or physical reconciliation.
   authoring stop (D7). The orchestrator runs those reviews.
 - **Build review.** Codex-written, so a fresh Claude agent and Grok, each with a mandatory FORM ablation.
 - **Optical model point.** As in "Setting": leading eikonal with the retained multigraded set above.
-  Part D's mechanical and inherited energy content retains its stated untruncated domain. Results do not
+  Part D's mechanical and inherited energy content retains its stated untruncated domain. Its results
+  are conditional on P5's adopted momentum carriage and P6's adopted steady pressure-plus-linear-viscous
+  form, with general live `η` and `ζ`. P6 is not derived from P1; power-law creep is outside this model
+  point, and relaxation outside the adopted form and energy accounting remain OPEN. Results do not
   transfer to:
   - the strong field;
   - the throat mouth or interior;
@@ -480,8 +527,8 @@ performs no register edit or physical reconciliation.
   - anisotropic or coupled branches.
 - **Deferred to the build** (implementation, not new physics; the build directive owns it):
   - the symbolic handling of the every-`b` requirement;
-  - component/measure calculus on the supplied coordinate mass law, with the D4 gradient-scale
-    qualification above for any induced-measure interpretation;
+  - component/measure calculus on the supplied coordinate mass law, with amended D4's reading,
+    live-gradient and relative-to-`j_n` bound requirements above for any induced-measure claim;
   - representation of general OPEN actions and executable controls, including FORM ablation (E2).
 - **Stop and report**, without choosing, when any of these happens:
   - a second method failure;
@@ -489,5 +536,5 @@ performs no register edit or physical reconciliation.
   - a premise this spec does not supply.
 - **The step record** interprets the results.
 
-**Authoring STOP.** Write v10 and report D1–D5 changes against v8, source conflicts and missing
-sourced pieces. No CAS, build, review launch, commit, push or spawned agent is part of this task.
+**Authoring STOP.** Repair v10 in place and report changes for repair items 1–4, source conflicts
+and missing sourced pieces. No CAS, build, review launch, commit, push or spawned agent is part of this task.
```

Exit code: `0`.

```bash
git diff --no-ext-diff f069c37c 79055918 -- research/pde_ledger_v3/directives/S9b_repair_decision_list.md
```

```text
diff --git a/research/pde_ledger_v3/directives/S9b_repair_decision_list.md b/research/pde_ledger_v3/directives/S9b_repair_decision_list.md
index ab3047b6..dd0cc090 100644
--- a/research/pde_ledger_v3/directives/S9b_repair_decision_list.md
+++ b/research/pde_ledger_v3/directives/S9b_repair_decision_list.md
@@ -2,6 +2,10 @@
 
 **Author:** Claude (orchestrator), 2026-10-08. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
 both NEEDS REVISION; dispositions in `directives/_measurements/S9b_repair_dl_review_disposition.md`).
+**Amendment 1** (2026-10-08) rewrites D2's P5 and P6 rows, D3's live symbols and D4 in place, after spec v10 review
+round 1 (`af1674e5`; dispositions in `directives/_measurements/S9b_spec_v10_r1_review_disposition.md`). The user
+re-selected P5 and P6 the same day. The amendment had its own single Codex + Grok pass (both NEEDS REVISION) and was folded once (dispositions in
+`directives/_measurements/S9b_repair_dl_amend1_review_disposition.md`).
 
 This list sets what the next S9b spec version (v10) must contain, and routes the build findings still owed. It names
 objects and sources. It states no expected value, sign or relation between symbols.
@@ -24,8 +28,8 @@ Every result that depends on one is flagged with it.
 
 | Premise | Content | Qualification carried with it |
 |---|---|---|
-| **P5** (2026-10-07) | The brane's momentum density is `ρ_br V`. | If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
-| **P6** (2026-10-08) | In steady flow the brane's in-plane stress is an isotropic pressure `p_br(ρ_br)` that depends on `ρ_br` only, because shear relaxes under steady load (P1). `p_br` stays a general function. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. |
+| **P5** (2026-10-07; carriage 2026-10-08) | The brane's in-plane momentum density is `ρ_br V`, and that momentum is carried with the brane material at `V`. Its contribution to the in-plane momentum current is `ρ_br V^i V^j`. | Scoped to Part D. Any other in-plane momentum current stays OPEN. It is counted beside the supplied current and P6's stress, which it does not contain again. Part D states what the OPEN pieces must supply. If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
+| **P6** (2026-10-08, revised the same day) | In steady flow the brane's in-plane stress is an isotropic pressure plus a linear viscous stress: `T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`, in the Cauchy convention (traction `T^{ij} n_j`; force density `∂_j T^{ij}`). `p_br` is a general function of `ρ_br` only. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. The shear viscosity `η` and bulk viscosity `ζ` are live general profiles. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. It is an adopted steady form, not a consequence of P1: P1's relaxation response stays OPEN outside it. The linear viscous form excludes power-law creep. The power this stress dissipates is not identified with any O2 energy operand. It stays in the OPEN energy accounting, with its physical supplier and budget OPEN. The first 2026-10-08 row (pressure only, justified by relaxation under steady load) is withdrawn: steady flow keeps straining, so relaxation does not remove the viscous stress (spec v10 review C2). |
 | **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |
 | **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |
 
@@ -39,23 +43,27 @@ Every result that depends on one is flagged with it.
 **The deliverable.** Part D re-expresses the Part B conditions with that balance in force. It prints:
 - which of `δ`, `V`, `ρ_br` and `j_n` it determines relative to `GM`;
 - which of them stay free;
-- for each condition, the `j_n` it implies.
+- for each condition, the `j_n` it implies;
+- for each condition, what the OPEN pieces left in the in-plane balance must supply, each with its general live
+  dependence.
 
 **What must be true:**
 - **Supplied pieces.** The spec supplies O2's balance pieces and the premises as equations, from the O2 sources. It
   does not compose them in advance; the engines compose them.
 - **OPEN operands.** Every O2 OPEN operand that P3–P6 and w-parity do not supply stays a general live unknown, not
   an engine-chosen family (v9 review). This includes its admissible gradient and history dependence (O2 record §8).
-- **Live symbols.** `j_n`, `p_br(ρ_br)` (and with it `c_comp`), `α` and the asymptotic `ρ_br⁰` stay live. No order or
-  value is assigned to `c₀/c_comp`.
+- **Live symbols.** `j_n`, `p_br(ρ_br)` (and with it `c_comp`), `η`, `ζ`, `α` and the asymptotic `ρ_br⁰` stay live.
+  No order or value is assigned to `c₀/c_comp`, `η` or `ζ`.
 - **Bulk-density route.** v8's Part C is unchanged.
 
 ## D4. Induced metric and order
 - **v8's sentence.** v8 says the induced-metric mass balance differs from the flat form by "a relative `O(ε)`
-  correction to the implied `j_n`". That holds only if `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok).
-- **One rule for v10.** The supplied mass balance is on the coordinate `d³x` measure (O2 record §8). Any reading on
-  the induced measure carries the scale `∂_r[(∂ξ_w)²]` live, or states the condition above. No order is attached
-  without one or the other.
+  correction to the implied `j_n`". No source supplies that order. The gradient-scale condition
+  `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok) does not supply it either (spec v10 review C4).
+- **One rule for v10.** The supplied mass balance is on the coordinate `d³x` measure (O2 record §8). An induced-measure
+  claim names which reading it uses: the same densities re-expressed per induced volume, or a mass law imposed on the
+  induced measure. Such a claim keeps `∂_r[(∂ξ_w)²]` live. It attaches no relative order to the correction to `j_n`
+  unless it states a condition that bounds that correction relative to `j_n` itself.
 - **Neutral Part D.** With `ξ_w = 0` (w-parity), the supplied `g_ij` reduces to `δ_ij`, and the engines print that
   reduction. This limits the claim; it adds no term.
 
@@ -103,6 +111,7 @@ repair build directive, which gets its own G2 pass, and the build legs verify th
 4. **Afterwards.** Build legs review until clear, with a FORM ablation each. Then the comparator, then the record.
 
 ## Not decided here
-- Any expected value, sign or limit. This includes the relation of `c_comp` to `c₀`, and of `j_n` to zero.
+- Any expected value, sign or limit. This includes the relation of `c_comp` to `c₀`, of `j_n` to zero, and the size
+  or sign of `η` and `ζ`.
 - The laws for S8 inertia and stress, and O1/O3–O7. The exception is Part D, where P5, P6 and the density link
   supply content as scoped above.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/S9b_repair_decision_list.md | sed -n '"'"'5,8p;27,68p;70,90p'"'"''
```

```text
     5	**Amendment 1** (2026-10-08) rewrites D2's P5 and P6 rows, D3's live symbols and D4 in place, after spec v10 review
     6	round 1 (`af1674e5`; dispositions in `directives/_measurements/S9b_spec_v10_r1_review_disposition.md`). The user
     7	re-selected P5 and P6 the same day. The amendment had its own single Codex + Grok pass (both NEEDS REVISION) and was folded once (dispositions in
     8	`directives/_measurements/S9b_repair_dl_amend1_review_disposition.md`).
    27	**New for S9b:**
    28	
    29	| Premise | Content | Qualification carried with it |
    30	|---|---|---|
    31	| **P5** (2026-10-07; carriage 2026-10-08) | The brane's in-plane momentum density is `ρ_br V`, and that momentum is carried with the brane material at `V`. Its contribution to the in-plane momentum current is `ρ_br V^i V^j`. | Scoped to Part D. Any other in-plane momentum current stays OPEN. It is counted beside the supplied current and P6's stress, which it does not contain again. Part D states what the OPEN pieces must supply. If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
    32	| **P6** (2026-10-08, revised the same day) | In steady flow the brane's in-plane stress is an isotropic pressure plus a linear viscous stress: `T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`, in the Cauchy convention (traction `T^{ij} n_j`; force density `∂_j T^{ij}`). `p_br` is a general function of `ρ_br` only. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. The shear viscosity `η` and bulk viscosity `ζ` are live general profiles. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. It is an adopted steady form, not a consequence of P1: P1's relaxation response stays OPEN outside it. The linear viscous form excludes power-law creep. The power this stress dissipates is not identified with any O2 energy operand. It stays in the OPEN energy accounting, with its physical supplier and budget OPEN. The first 2026-10-08 row (pressure only, justified by relaxation under steady load) is withdrawn: steady flow keeps straining, so relaxation does not remove the viscous stress (spec v10 review C2). |
    33	| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |
    34	| **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |
    35	
    36	## D3. Part D: the linked brane (new)
    37	**The object.** Part D is the brane's far-field steady in-plane momentum balance. Its ingredients:
    38	- O2's conditional hold balance, from the record's §8 handoff and its sources;
    39	- P3, P5, P6 and w-parity, as supplied;
    40	- v8's mass balance `∇·(ρ_br V) = −j_n`;
    41	- the density link.
    42	
    43	**The deliverable.** Part D re-expresses the Part B conditions with that balance in force. It prints:
    44	- which of `δ`, `V`, `ρ_br` and `j_n` it determines relative to `GM`;
    45	- which of them stay free;
    46	- for each condition, the `j_n` it implies;
    47	- for each condition, what the OPEN pieces left in the in-plane balance must supply, each with its general live
    48	  dependence.
    49	
    50	**What must be true:**
    51	- **Supplied pieces.** The spec supplies O2's balance pieces and the premises as equations, from the O2 sources. It
    52	  does not compose them in advance; the engines compose them.
    53	- **OPEN operands.** Every O2 OPEN operand that P3–P6 and w-parity do not supply stays a general live unknown, not
    54	  an engine-chosen family (v9 review). This includes its admissible gradient and history dependence (O2 record §8).
    55	- **Live symbols.** `j_n`, `p_br(ρ_br)` (and with it `c_comp`), `η`, `ζ`, `α` and the asymptotic `ρ_br⁰` stay live.
    56	  No order or value is assigned to `c₀/c_comp`, `η` or `ζ`.
    57	- **Bulk-density route.** v8's Part C is unchanged.
    58	
    59	## D4. Induced metric and order
    60	- **v8's sentence.** v8 says the induced-metric mass balance differs from the flat form by "a relative `O(ε)`
    61	  correction to the implied `j_n`". No source supplies that order. The gradient-scale condition
    62	  `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok) does not supply it either (spec v10 review C4).
    63	- **One rule for v10.** The supplied mass balance is on the coordinate `d³x` measure (O2 record §8). An induced-measure
    64	  claim names which reading it uses: the same densities re-expressed per induced volume, or a mass law imposed on the
    65	  induced measure. Such a claim keeps `∂_r[(∂ξ_w)²]` live. It attaches no relative order to the correction to `j_n`
    66	  unless it states a condition that bounds that correction relative to `j_n` itself.
    67	- **Neutral Part D.** With `ξ_w = 0` (w-parity), the supplied `g_ij` reduces to `δ_ij`, and the engines print that
    68	  reduction. This limits the claim; it adds no term.
    70	## D5. The O2 record's returned obligation
    71	The O2 record returned one interpretation to the orchestrator: the WL-only `ξ_w''` dependence, together with the
    72	seven WL-only first-derivative keys at the same locations.
    73	
    74	**Where they sit.** The record places all eight in the momentum-flux and energy-flux roles (its §4 location table).
    75	
    76	**Routing:**
    77	- **The obligation stays the orchestrator's** (M1; record §8). S9b does not discharge it. The S9b record carries it
    78	  as undischarged. It is distinct from any step's ownership of a physical input.
    79	- **Retained keys.** Part D retains all eight keys, at every §4 location, for every piece of O2 content that no
    80	  adopted premise supplies (record §8).
    81	- **Supplied content.** Where P5, P6 or the density link supplies an object, the premise states that object's
    82	  dependence, and every result that uses it is flagged.
    83	  - The spec author determines, from the O2 sources, which content each premise supplies.
    84	  - The author does not compose the balance (D3).
    85	  - The user chose P6 over keeping the stress OPEN (2026-10-08).
    86	- **Register.** A register entry follows only if a sourced requirement is established. An unresolved difference
    87	  alone does not create one (record's rest-on criterion).
    88	
    89	## D6. Owed build findings: routed to the repair build directive
    90	These come from the build review preserved at `bb94b885`, as listed in that commit's message. They concern implementation, not spec physics. They go to the
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | sed -n '"'"'121,147p;170,181p;272,305p'"'"''
```

```text
   121	### 3.2 OPEN material and constitutive inputs
   122	
   123	Every entry below is a **general OPEN operand**. Notation supplies no closed argument list, locality,
   124	instantaneous response, tensor realization, constitutive family or derivative cutoff. Live fields,
   125	gradients and material history remain admissible dependences (C §1).
   126	
   127	| Operand; status/source | Physical content and entry into the object | Owner retained |
   128	| --- | --- | --- |
   129	| `𝔅_A13`; **OPEN**, C §2 | Real/dissipative versus complex/inertial order-field branch, its action, degrees of freedom and conversion content. Branch dependence stays attached to the material and source responses; premise 1 does not decide A13. | S1's A13 gate, propagated through S5/S12. |
   130	| `𝒯_br^live`; **OPEN form**, with character fixed by **adopted premise 1**, C §3 | Full live in-plane and normal stress, including any relaxation content. It enters the internal material-force accounting once, with the live reference/strain state. No conservative/dissipative split is adopted. | O2 carries it; S1.5/S8 antecedents, S22 nonlinear completion; relaxation ownership unassigned. |
   131	| `𝒫_br^cons`; **OPEN conservative antecedent**, C §3 | Conservative material momentum input to the live momentum description. It is not an additional external force or a second momentum species. No live reduction into `ℐ_br^live` is supplied. | S1.5/S8. |
   132	| `𝒯_br^cons`; **OPEN conservative antecedent**, C §3 | Conservative in-plane/normal-stress content, retained as an antecedent to the material-force description. It is not added beside the full live stress as another stress. Its relation to that stress remains unspecified. | S1.5/S8; nonlinear completion S22. |
   133	| `ℐ_br^live`; **OPEN**, C §3 | Flowing/embedded inertial and kinetic response, entering momentum storage and transport. No identification with `ρ_br V`, a quadratic live kinetic energy, or a constant inertia is supplied. | S8, with S1.5 antecedents when available. |
   134	| `𝒩_br^live`; **OPEN**, C §3 | Normal material response, keeping embedding, centre and thickness content distinct. It qualifies the normal material accounting jointly with the stress and inertia inputs; their overlap/identification remains unresolved, rather than being treated as three additive normal forces. | S5–S8 ingredients, Q1/Q2 static sector, S8/S22 live completion. |
   135	| `𝒜_rot^live`; **OPEN**, C §3 | Whether internal angular momentum/couple stress is present and its physical rotational reference frame. Retain its effect on the admissible momentum/stress action and on the corresponding power accounting. No stress symmetry, vanishing couple content or chosen carrier is assumed. | S8 requirements; register assignments retain pass-2-review-pending status. |
   136	| `ℛ_ref/strain^live`; **OPEN form**, **adopted premise 1**, C §4 | Reference/strain evolution and relaxation, including carrier, transport, formation/renewal through conversion/return, and work content. It enters through the material's evolving state/history and its explicit energy partner, not as an independently postulated body force. | Relaxation/reference ownership **unassigned**; O2 carries it. S8/S22 links are inferences; conversion/return functions remain S12's. |
   137	| `ℳ_⊥` (O1); **OPEN**, C §6 | General stiffness response, including frequency regime and loading/material history, bulk/flow/embedding/thickness/projection dependence and gradients. It is a constitutive input to the linked optical/material description; no relation to `𝒯_br^live` is supplied. It is not a separate force. | S8, substrate reduction and S22 completion. |
   138	| `ℛ_br` (O7); **OPEN**, C §6 | General brane-density response with live bulk, flow, embedding, thickness/projection dependence and gradients. It is input to the density in mass/momentum/energy accounting; it is not inferred from the bulk EOS or slab factorization. | S8, substrate reduction and S22 completion. |
   139	| `ℰ_h^live` (O4); **OPEN**, C §5 | Live embedding/longitudinal relation with flow, exchange, variable coefficients and gradients. Keep it as a coupled input with `ξ_w=ℓh`. Its identity with or independence from O2's normal relation remains **unsettled**. Do not impose a second normal equation or use it as a duplicate normal force. | Q1/Q2 static ingredients; O4 live identification, S8/S22 completion. |
   140	
   141	**Accounting convention, not a constitutive decomposition:** the material entries describe one
   142	brane-material momentum/force accounting object. C's names supply no equation saying which response
   143	contains or determines another. Keep the unresolved relations visible within that object. Where the
   144	force action, momentum map, normal response or their compatibility is unsupplied, print the formal
   145	component action of the named OPEN inputs. Do not manufacture an explicit action by choosing a
   146	stress measure or by independently adding every named response. An auxiliary name for an unevaluated
   147	action is notation for these existing OPEN inputs, not a new closed physical response.
   170	- **Material momentum and inertia:** use the live flowing/embedded response `ℐ_br^live`, with the
   171	  conservative momentum antecedent `𝒫_br^cons` still named. Keep storage, spatial transport and
   172	  material-history effects where the response requires them. Eulerian steadiness does not remove
   173	  convective transport. The mass law identifies the mass-source convention; it does not close the
   174	  momentum-to-velocity relation. The material velocity on which this response acts is that of §1,
   175	  with in-plane components `V^i` and the bulk-direction component that the supplied graph and `V`
   176	  determine. `ℐ_br^live` still supplies no map from that velocity to momentum, and normal content that
   177	  the graph does not determine stays with `𝒩_br^live` and `𝒥_map`.
   178	- **Internal material force:** use the full `𝒯_br^live` once. Retain the conservative antecedent
   179	  `𝒯_br^cons`, evolving `ℛ_ref/strain^live` and admissibility/frame content `𝒜_rot^live` as its
   180	  unresolved material inputs, without asserting a split or containment relation. Do not add the
   181	  conservative antecedent or a separately invented relaxation stress to the full stress.
   272	**`ℬ_E^steady` is OPEN**, with the accounting requirement fixed by **adopted premise 1** and its S21
   273	label in §2 (C §8). It denotes the steady energy relation to be constructed, not a precomputed
   274	residual. Its required inputs are:
   275	
   276	| OPEN operand (C §8) | Required content and entry | Owner |
   277	| --- | --- | --- |
   278	| `ℰ_br^live`, `𝒥_E^live` | Material energy storage and transport compatible with the live stress, inertia, normal response, reference evolution and rotational content. Retain live transport in the steady setting. | S1.5 antecedents, S8 material, S22 completion. |
   279	| `𝒫_ref/relax^live` | Explicit power associated with `ℛ_ref/strain^live`, with no formula, sign or vanishing assumed. | O2 requirement; relaxation/reference ownership unassigned. |
   280	| `𝒫_convert/exchange^live` | Order-conversion work, energy carried with exchanged material, and additional non-variational energy partners, with their reaction/supply systems. Premise 3 fixes transported momentum, not the carried total energy. | S12. |
   281	| `𝒫_boundary^live` | Mechanical face/support work and other energy transfer through the live bulk/boundary data; any explicitly supplied hold is identified. | O2 accounting; S12 native data; Q2/S22 core response. |
   282	| `𝒮_E,net`, `𝒫_E,supply` | Identity of the physical supplier of any net power and its stated power budget. Both stay general OPEN inputs. Naming the drain does not specify available energy or a budget. | O2 requirement; source/holder/supplier forms retain their owners. |
   283	
   284	Pair each mechanical contribution with the velocity or generalized rate of the material point or
   285	degree of freedom on which it acts, on the same measure (§1) and with the same geometric map as its
   286	momentum occurrence. For a brane material point on the supplied graph, that velocity is the material
   287	velocity of §1, the same velocity as in its momentum occurrence. In particular:
   288	
   289	- The stress/normal-response contribution has its corresponding material stress work and energy
   290	  transport. Any couple/frame content retains the matching rotational/generalized-rate work as
   291	  OPEN `𝒜_rot^live` content. No ordinary elastic stored-energy functional is substituted for the
   292	  full viscoelastic response.
   293	- A face/support traction is paired with the actual velocity of its load application point.
   294	  Where the live face-to-material identification is unavailable, retain it through the OPEN
   295	  normal/material response and map. Do not use the light perturbation velocity or silently equate
   296	  `V`, `V_s` and `v_dr`.
   297	- Carried exchange energy and additional conversion/source power are accounted for separately from
   298	  mechanical face work. Transported momentum at the local brane material velocity alone does not
   299	  authorize a formula for total energy per converted mass, or an energy-free change of material
   300	  reference.
   301	- Keep `𝒫_ref/relax^live` explicit even when its work is represented within the material stress and
   302	  internal-energy accounting. Identify its occurrence there instead of adding the same work again
   303	  as an independent external power. The same rule applies when boundary or source power is also
   304	  described within another energy operand. These names are accounting obligations, not an adopted
   305	  additive decomposition into independent channels.
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/directives/O2_input_contract.md | sed -n '"'"'94,102p'"'"''
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
```

Exit code: `0`.

```bash
bash -c 'nl -ba research/pde_ledger_v3/steps/O2_steady_brane_balance.md | sed -n '"'"'248,278p;283,296p;469,515p'"'"''
```

```text
   248	| Joined row / structure stream line | Declared roles with no WL occurrence |
   249	|---|---|
   250	| `energy_storage` / 195 | `OPEN_MaterialCompatibility` |
   251	| `energy_transport` / 199 | `OPEN_MaterialCompatibility` |
   252	| `energy_power` / 203 | `OPEN_MaterialCompatibility`, `OPEN_MaterialEnergyDensity`, `OPEN_MaterialEnergyFlux_0`, `OPEN_MaterialEnergyFlux_1`, `OPEN_MaterialEnergyFlux_2` |
   253	| `energy_balance` / 207 | `OPEN_MaterialCompatibility` |
   254	
   255	These are differences in the action inventories beneath the energy balance entries. They remain
   256	OPEN with no owner named; the record's literal stored classification tuple count match in §3 (M7) does
   257	not establish a match of the actions beneath those entries. No occurrence-level pairing, complete
   258	work equality or duplicate-power conclusion is formed from either engine's one-sided content.
   259	
   260	At the **row-unioned paired-action inventory level**, the comparator prints `[]` for `head` and
   261	`role` differences in all 234 paired role groups; M4's literal counts give 234 of each empty field,
   262	230 empty `orientation` fields and four `[["argument",1]]` fields (M6 exhaustive retrieval).
   263	These are counts of printed grouped comparisons, not occurrence comparisons. Every nonempty paired head/role/
   264	orientation difference is therefore in the following table. These are **OPEN cross-engine
   265	differences; no owner named for adjudication**:
   266	
   267	| Joined row / structure stream line | Paired role | Printed orientation difference (PY minus WL) |
   268	|---|---|---|
   269	| `energy_balance` / 207 | `OPEN_MaterialEnergyDensity` | `[["argument",1]]` |
   270	| `energy_balance` / 207 | `OPEN_MaterialEnergyFlux_0` | `[["argument",1]]` |
   271	| `energy_balance` / 207 | `OPEN_MaterialEnergyFlux_1` | `[["argument",1]]` |
   272	| `energy_balance` / 207 | `OPEN_MaterialEnergyFlux_2` | `[["argument",1]]` |
   273	
   274	For each of these four roles, the stored PY occurrence collection has length **2**, with orientation
   275	fields `{"argument":1}` and `{"1":1}`; WL has length **1**, with `{"1":1}` (M12, line 207).
   276	The raw PY function-name nodes and their paths in `energy_power` (line 202) and `energy_balance`
   277	(line 206) show the density and all three flux actions inside `OPEN_JointPowerAccounting`.
   278	In `energy_balance`, each also has an occurrence outside that joint-power action (M12).
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
   469	Part D must retain **every named operand and complete live-object dependence printed by either engine**,
   470	including all eight WL-only derivative keys at every location in §4 and all PY-only native/chart,
   471	core/material-compatibility content. It must also preserve the spec's general admissible live fields,
   472	gradients, entire material history and other native dependences (spec L124–125). Neither engine's
   473	finite OPEN inventory is the object's dependence list, and even taking both lists cannot impose a
   474	closed argument list or derivative/history cutoff. The full inventories/operands remain in the
   475	filed production stream; their differences are M6. No surrounding held calculus, coefficient algebra,
   476	position or multiplicity outside §5's inventory scope becomes a comparison result.
   477	
   478	The mass-law qualification travels explicitly: `∇·(ρ_br V)=−j_n` is a **coordinate-`d³x` measure**
   479	input. For a claim reading `j_n` or `ρ_br` per induced measure, or comparing an induced-metric mass
   480	law, the recorded qualification is a **relative `O(ε)` correction to `j_n`**, not a live O6 law;
   481	O2 supplies and derives **no induced-metric mass balance** (spec L330–336). Fixed-`ℓ` counting,
   482	independent orbital `GM`, only the stiffness/density ratio inheriting the speed-change grade, and
   483	first order in `f` with no `f`–`ε` relation also travel unchanged. Historical static/homogeneous,
   484	frozen/uniform or supplied-profile relations retain exactly §2's restricted domains.
   485	
   486	Part D must keep coordinate versus graph-normal content distinct, use the same material and measure
   487	with compatible native maps, retain untruncated unknown terms/grades, avoid duplicate stress/exchange/
   488	power content, and name the supplier/budget if a closure requires net power. Residual `0` at the
   489	balance closed-part level (printed, M8), unioned balance-role inventory matches (printed, M6),
   490	and joint stored classification tuple count matches (record literal counts, M7) supply neither
   491	a solved material state nor full OPEN equality. Each condition inherits the whole difference
   492	ledger in §4 wherever a difference bears on its content, including native-map/measure, material
   493	compatibility, named/live dependences and generalized-work scope (M6/M9/M12).
   494	In particular, the force/power-compatibility and duplicate-power conditions carry the one-sided
   495	`OPEN_MaterialCompatibility` actions in all four energy rows, the one-sided density/flux actions in
   496	`energy_power`, **and** the four paired density/flux orientation differences `[["argument",1]]`
   497	in `energy_balance` (line 207). The stored PY density/flux actions inside `OPEN_JointPowerAccounting`
   498	travel with those conditions in both `energy_power` and `ℬ_E^steady` itself, alongside the separate
   499	occurrences in `energy_balance` (lines 202/206; §4; M12). The record's literal stored classification tuple
   500	count match (M7) and the comparator's unioned balance-role inventory matches (M6)
   501	do not settle compatibility or whether the complete work/energy content is counted once.
   502	Part D must carry these differences alongside those conditions; it may not discard or reconcile
   503	them, or judge either engine to duplicate power from either limited match. The `xi_w''` difference
   504	also travels with the admissible-dependence and force/power conditions as **OPEN interpretation,
   505	undischarged at sub-step 7 and returned to the orchestrator under M1** (§4). Any physical resolution
   506	requires work beyond retrieval.
   507	Source: spec §§1–10, contract §§1–10 and M6–M9/M12.
   508	
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
rg -n '^\*\*Authors|^are `directives|^  The divergence|^  `ρ_br`\. This input|^  names its reading|^  the induced measure|^  to `j_n`|^  D4|^  gradient-scale|^- \*\*Part D|^  each condition\. For|^  must supply|^\*\*P5|^\*\*P6|^𝒫_br,inplane|^𝒥_br,carry|^𝒯_br,inplane|^T\^|^         [+=]|^t_br,inplane|^c_comp|^`𝒥_br,carry|^momentum current|^The shear viscosity|^P6 supplies|^The power this stress|^is assigned|^\| Material momentum|^\| Internal material force|^Each engine retains|^against either O2|^relation\. `j_n`|^assigned to|^The full inherited|^`ξ_w` respectively|^P5 supplies|^\*\*Additional engine|^\*\*Comparator|^  are conditional|^  form, with general|^  point, and relaxation|^\*\*Authoring STOP' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
3:**Authors:** Codex (v5–v6, v10), with v7–v8 edits by Claude (orchestrator). **Status:** v10 repair 1,
11:are `directives/S9b_repair_decision_list.md` at `79055918` (amendment 1). v9's Part D is not an input. Paths
66:  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
68:  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. Any induced-measure claim
69:  names its reading: the same densities re-expressed per induced volume, or a mass law imposed on
70:  the induced measure. It keeps `∂_r[(∂ξ_w)²]` live and attaches no relative order to the correction
71:  to `j_n` unless it states a condition bounding that correction relative to `j_n` itself (amended
72:  D4). No such bound is supplied here. Neither the supplied slope counting nor the historical
73:  gradient-scale condition `∂_r[(∂ξ_w)²] = O(ε/r)` supplies that bound. The latter is not imposed.
239:- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
243:  each condition. For each condition, also print what the OPEN pieces left in the in-plane balance
244:  must supply, each with its general live dependence (amended D3). Preserve the every-far-zone-`b`
327:**P5 — momentum density and carriage (2026-10-07; carriage 2026-10-08; amended D2).** Scoped to
332:𝒫_br,inplane^live ≡ ρ_br V ,
333:𝒥_br,carry^ij ≡ ρ_br V^i V^j .
336:`𝒥_br,carry^ij` names P5's contribution to the in-plane momentum current. Any other in-plane
337:momentum current stays a general OPEN action, counted beside this supplied current and P6's stress;
343:**P6 — adopted steady in-plane stress (2026-10-08, revised the same day; amended D2).** Scoped to
348:𝒯_br,inplane^live,ij ≡ T^{ij} ,
349:T^{ij} ≡ −p_br(ρ_br) δ^{ij}
350:         + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k)
351:         + ζ δ^{ij} ∂_kV^k ,
352:t_br,inplane^i ≡ T^{ij} n_j ,      F_br,inplane^i ≡ ∂_j T^{ij} ,
353:c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
357:The shear viscosity `η` and bulk viscosity `ζ` are live general profiles, including their gradients.
358:P6 supplies the steady in-plane part of `𝒯_br^live` only. It is an adopted steady form, not a
362:The power this stress dissipates has no supplied identification with an O2 energy operand. It
364:is assigned to `c₀/c_comp`, `η` or `ζ`.
402:| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density and the material-carried contribution `ρ_br V^i V^j` to the in-plane current. | Any other in-plane momentum current, with general live dependence, counted beside the supplied current and P6 stress and containing neither again; normal/embedded responses and unresolved relations wherever unsupplied. The conservative antecedent is named, rather than added as another momentum species. |
403:| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the full steady in-plane stress in its adopted pressure-plus-linear-viscous form, with the stated Cauchy traction/force convention. | Normal material response, conservative antecedents, reference evolution outside the adopted form, and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; no additional in-plane stress duplicates P6. |
413:Each engine retains the named operands and complete live-object dependences supplied in this spec
417:against either O2 engine's historical emissions belongs to the comparator and record, as stated below.
441:relation. `j_n`, `p_br(ρ_br)`, `c_comp`, `η`, `ζ`, `α` and `ρ_br⁰` remain live. No order or value is
442:assigned to `c₀/c_comp`, `η` or `ζ`. No new scale is used to remove an OPEN operand or select a profile.
452:The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
455:`ξ_w` respectively. Each engine retains all eight as admissible live dependences at every following
469:P5 supplies the in-plane momentum density and its material-carried current contribution; P6 supplies
478:**Additional engine-facing OPEN content.** Retain named/native geometry, chart/measure/map,
489:**Comparator and record handoff (O2-R §§4–5, 8).** The comparator and record check that the Part D
518:  are conditional on P5's adopted momentum carriage and P6's adopted steady pressure-plus-linear-viscous
519:  form, with general live `η` and `ζ`. P6 is not derived from P1; power-law creep is outside this model
520:  point, and relaxation outside the adopted form and energy accounting remain OPEN. Results do not
539:**Authoring STOP.** Repair v10 in place and report changes for repair items 1–4, source conflicts
```

Exit code: `0`.

```bash
git rev-parse af1674e5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md HEAD:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md HEAD:research/pde_ledger_v3/directives/S9b_repair_decision_list.md 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md HEAD:research/pde_ledger_v3/steps/O2_steady_brane_balance.md 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md HEAD:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md 217a92e9:research/pde_ledger_v3/directives/O2_input_contract.md HEAD:research/pde_ledger_v3/directives/O2_input_contract.md 77d2c39a:research/pde_ledger_v3/directives/O2_premise_decision_list.md HEAD:research/pde_ledger_v3/directives/O2_premise_decision_list.md
```

```text
ccaa2fe51b924fb3880593783d715443e8ec0fd8
ccaa2fe51b924fb3880593783d715443e8ec0fd8
dd0cc090c5cd6486909617dbfacfd6e6baffdf2c
dd0cc090c5cd6486909617dbfacfd6e6baffdf2c
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

## Repair 2: baseline `478b5e7b`, decisions unchanged at `79055918`

Author: Codex. These are mechanical source/version and text lookups for the three requested repairs.
The same-content density equations define D4's density-re-expression reading; they are not a new
mass law. The graph-velocity notation identifies the source-constructed material tangent lift;
unsupplied normal material/map content stays OPEN. Review and authoring workflow is outside the
engine-facing spec. No CAS or physical computation was run. `cat -E` marks line endings in the
new literal stdout blocks; previous lookup sections remain preserved as filed.

```bash
bash -o pipefail -c 'git diff --no-ext-diff 478b5e7b -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | cat -E'
```

```text
diff --git a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md$
index 015d369c..fa57ad80 100644$
--- a/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md$
+++ b/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md$
@@ -1,10 +1,8 @@$
 # S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)$
 $
-**Authors:** Codex (v5–v6, v10), with v7–v8 edits by Claude (orchestrator). **Status:** v10 repair 1,$
-2026-10-08, against the preserved, not accepted baseline `af1674e5`; no repair clearance or computed$
-result is claimed.$
+**Version:** v10, repair 2 (2026-10-08).$
 $
-**Deliverable:** specify the v8 optical objects in Parts A–C and the neutral, linked steady in-plane$
+**Objects:** the v8 optical objects in Parts A–C and the neutral, linked steady in-plane$
 brane balance and its conditional optical requirements in Part D, retaining every unsupplied O2 input.$
 $
 **Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions$
@@ -16,8 +14,8 @@ below are relative to `research/pde_ledger_v3/`. The additional source abbreviat$
 - **O2-C:** `directives/O2_input_contract.md` at `217a92e9`.$
 - **O2-P:** `directives/O2_premise_decision_list.md` at `77d2c39a`.$
 $
-`CLAUDE.md` M1–M3 and E1–E2 govern the artifact. Equations labelled **supplied** are conditional$
-inputs that this build cannot test. Adopted premises retain that status and their provenance. Every$
+Equations labelled **supplied** are conditional inputs that the computation cannot test. Adopted$
+premises retain that status and their provenance. Every$
 dependent result is flagged with the supplied identification or adopted premise it uses. Reference$
 objects are comparison inputs, rather than premises for deriving the brane's observables. The spec$
 names computed objects and solution conditions; it supplies no outcome or expected-value acceptance test.$
@@ -29,7 +27,7 @@ needs in order to be deflected and delayed by a mass. This step states that requ$
 local wave equation that brane light sees near a mass. It does not derive that equation from the substrate:$
 the nonuniform brane operator is unfinished S11c work, and the substrate comes at the knit.$
 $
-## Governing object (supplied; the build cannot test it)$
+## Governing object (supplied; the computation cannot test the input)$
 $
 Light is the transverse branch of the brane's in-plane displacement `u`, and `u` is the material displacement$
 of the stuff whose density is `ρ_br` (`steps/S11_stray_longitudinal.md:32–35`; S10). Near the mass, at$
@@ -65,12 +63,31 @@ Each piece is a supplied identification. Flag any result that depends on one.$
 $
   The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`$
   coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as$
-  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. Any induced-measure claim$
-  names its reading: the same densities re-expressed per induced volume, or a mass law imposed on$
-  the induced measure. It keeps `∂_r[(∂ξ_w)²]` live and attaches no relative order to the correction$
-  to `j_n` unless it states a condition bounding that correction relative to `j_n` itself (amended$
-  D4). No such bound is supplied here. Neither the supplied slope counting nor the historical$
-  gradient-scale condition `∂_r[(∂ξ_w)²] = O(ε/r)` supplies that bound. The latter is not imposed.$
+  `ρ_br`. This input supplies no induced-measure or finite-slab replacement mass law. Every$
+  induced-measure claim names which of the following readings it uses (D4):$
+$
+  - **Density re-expression.** The same mass and normal exchange are re-expressed as densities per$
+    induced volume. The measure and same-content identifications defining this reading are$
+$
+    ```$
+    dvol_g ≡ √det(g_ij) d³x ,$
+    ρ_br^(g) dvol_g ≡ ρ_br d³x ,      j_n^(g) dvol_g ≡ j_n d³x .$
+    ```$
+$
+    These are supplied identifications of the reading, using the supplied metric; they do not$
+    replace the mass law. Print the re-expressed densities and correction objects with their$
+    orders under the supplied metric and slope counting. No additional independent bound on the$
+    drain divergence is a premise of this density re-expression. The imposed-law restriction$
+    below does not apply to this reading.$
+  - **Imposed mass law.** An independently imposed mass law on the induced measure, without the$
+    same-content density re-expression above, is a separate physical input; none is supplied here.$
+    A comparison on this reading names that law and keeps$
+    `∂_r[(∂ξ_w)²]` live. A relative order for its correction to the implied `j_n` requires a$
+    condition bounding that correction relative to `j_n` itself. **For this imposed-law reading$
+    only**, the unqualified recorded relative-order statement in O2-R L478–480 is not carried.$
+    The supplied slope counting and the historical gradient-scale condition$
+    `∂_r[(∂ξ_w)²] = O(ε/r)` supply no such uniform drain-relative bound for an imposed-law$
+    comparison (D4). The historical derivative condition is not imposed here.$
 $
   `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is$
   owned by the gravity sector or S12. S11b's uniform background normal drain `v_dr` is a different object$
@@ -125,7 +142,7 @@ because rulers and clocks are made of the same medium. Only the far-field observ$
   with nonnegative integer `a`, `b` and `c`. Each retained optical monomial is printed separately;$
   optical terms outside this set are not computed. This is not a truncation of Part D's mechanical$
   or inherited energy accounting.$
-- **Order counting** (supplied; the build cannot test it; owned by the gravity sector or S12). In the far zone,$
+- **Order counting** (supplied; the computation cannot test it; owned by the gravity sector or S12). In the far zone,$
   let `ε(r) ≡ GM/(c₀² r)`, with `GM` as supplied below. The supplied counting is$
 $
   ```$
@@ -262,7 +279,8 @@ Q(x,t) = Q(r)      (steady scalar profiles in this setting).$
 steadiness supplies no material constancy or history cutoff. Coordinate momentum storage, transport,$
 exchange, force, load, energy and power densities use the same coordinate `d³x` measure as the mass$
 law. Native face-area factors and reductions remain explicit through O6. Coordinate `w` content and$
-graph-normal content are distinct; the centre graph fixes no finite-thickness face response.$
+graph-normal content are distinct. The centre graph fixes neither finite-thickness face nor interior$
+normal velocity content; any content it does not determine stays with `𝒩_br^live` and `𝒥_map`.$
 $
 ### Adopted equations and their scope$
 $
@@ -381,30 +399,43 @@ The link supplies no `ρ_br(f)` law or separate smallness grade for the brane de$
 The engines compute its implication for `δ`; no such implication is stated here.$
 $
 **w-parity — neutral sector (2026-10-07; D2, D4).** For an electrically neutral mass the adopted$
-far-field `w → −w` symmetry restricts the graph displacement and material bulk-direction velocity:$
+far-field `w → −w` symmetry restricts the graph displacement and the constructed centre-graph$
+material velocity. Define the supplied graph embedding and its material tangent lift by$
 $
 ```$
-ξ_w(r) ≡ 0 ,      U_material^w(r) ≡ 0 .$
+X_graph(x) ≡ (x¹,x²,x³,ξ_w(x)) ,      U_graph ≡ V^i ∂_i X_graph .$
+```$
+$
+`U_graph^w` denotes the bulk-coordinate component of this centre-graph material velocity, constructed$
+from the steady graph and `V` (O2-S L68–79; O2-R §§1–2). It is not an independent normal velocity$
+operand. The adopted neutral-sector equations are$
+$
+```$
+ξ_w(r) ≡ 0 ,      U_graph^w(r) ≡ 0 .$
 ```$
 $
 Print **neutral-sector restriction** with Part D and each dependent result, including the reduction$
 of the supplied `g_ij`. Parts A–C keep `ξ_w` live, including the charged case. This restriction adds$
-no force term and supplies no native-face, thickness, exchange-map or normal constitutive law.$
+no force term. Normal velocity content not determined by the centre graph remains OPEN in$
+`𝒩_br^live` and `𝒥_map`, including face and interior content and unsupplied face-to-material$
+identifications. The reduced bulk-direction carry `(Π_n^carry)^w` remains an OPEN action of those$
+operands with native geometry; it is not supplied as `j_n` times one centre-graph velocity$
+(O2-S §5). No native-face, thickness, exchange-map or normal constitutive law is supplied.$
 $
 ### O2 balance pieces that remain live$
 $
 The engines compose the balance from these physical roles and the supplied equations above.$
 No assembled momentum balance, total transport tensor, cancellation or profile solution is supplied$
-here (D3; O2-S §§3–6; O2-R §8). The author identifies the supplied content as follows:$
+here (D3; O2-S §§3–6; O2-R §8). The supplied content and remaining OPEN actions are:$
 $
 | O2 content | Supplied content in Part D | Content retained as OPEN |$
 |---|---|---|$
 | Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density and the material-carried contribution `ρ_br V^i V^j` to the in-plane current. | Any other in-plane momentum current, with general live dependence, counted beside the supplied current and P6 stress and containing neither again; normal/embedded responses and unresolved relations wherever unsupplied. The conservative antecedent is named, rather than added as another momentum species. |$
 | Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the full steady in-plane stress in its adopted pressure-plus-linear-viscous form, with the stated Cauchy traction/force convention. | Normal material response, conservative antecedents, reference evolution outside the adopted form, and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; no additional in-plane stress duplicates P6. |$
 | Optical/material constitutive inputs; O1 `ℳ_⊥`, O7 `ℛ_br` | The density link supplies Part D's optical stiffness response only. | Brane-density response and every other unsupplied material or reduction identification. Parts A–C retain the v8 local profiles. |$
-| Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric and neutral-sector restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, normal response and material/order/projection map. No second normal equation is imposed. |$
+| Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric, constructed `U_graph` and neutral centre-graph restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, face/interior normal velocity content not determined by the graph, normal response and material/order/projection map. No second normal equation is imposed. |$
 | Mechanical face/support loading; `T_hold,s`, `𝒯_bulk,n,s^live` | P4 supplies the bulk part's native direction only. | Full load/support partition, native normal-load amplitude, projections/maps and application-point velocities. No external support is selected by LAB_HELD. |$
-| Exchange; O3 `Π_n`, S12 partners/reactions | P3 supplies the coordinate in-plane carry `j_n V^i`. | Additional momentum partners/reactions and bulk-direction carry with O6/normal material content. Native convective transfer and the same O3 carry are not counted twice. |$
+| Exchange; O3 `Π_n`, S12 partners/reactions | P3 supplies the coordinate in-plane carry `j_n V^i`. | Additional momentum partners/reactions and bulk-direction carry with O6/normal material content, native geometry and unsupplied face-to-material velocities. Centre-graph parity supplies no reduced bulk-carry formula. Native convective transfer and the same O3 carry are not counted twice. |$
 | Boundary/source/core inputs; O5 `ℋ_core`, `𝔅_A13` | P2 supplies the drain-drive representation only. | A13 branch; distinct local conversion/controller and boundary/domain inventories; physical core holder/mouth data and their response. Q2/S22 retain O5; S12 retains conversion and return partners. |$
 $
 All OPEN entries are general unknown actions, with admissible live fields, gradients, entire material$
@@ -444,10 +475,9 @@ assigned to `c₀/c_comp`, `η` or `ζ`. No new scale is used to remove an OPEN$
 ### Returned O2 interpretation and dependence locations$
 $
 **D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7$
-obligation returned to the orchestrator under M1**. S9b does not discharge it; the S9b record carries$
-that status. It is distinct from ownership of a physical input. The engine-facing dependence inputs$
-are stated explicitly below; preservation of the historical difference ledger is the comparator/record$
-handoff at the end of this subsection.$
+obligation returned to the orchestrator under M1**. This is an unresolved input qualification,$
+distinct from ownership of a physical input. Each dependent Part D condition carries that status;$
+no physical interpretation or reconciliation is supplied. The dependence inputs are explicit below.$
 $
 The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,$
 `mu_perp`, `o2_rho_br_live` at order 1 and `xi_w` at order 2, each at the stored argument$
@@ -486,33 +516,27 @@ shared work/energy content visible, without choosing a decomposition or counting$
 These role names denote the existing material/energy operands above, not additional independently$
 additive energy species or a prescribed nesting/occurrence count.$
 $
-**Comparator and record handoff (O2-R §§4–5, 8).** The comparator and record check that the Part D$
-content no adopted premise supplies retains every named operand and complete live-object dependence$
-printed by either O2 engine, including all eight derivative keys at every listed location and the$
-PY-only native/chart and core/material-compatibility content. They carry the complete O2 difference$
-ledger wherever it bears on a condition. Each engine constructs from the explicit inputs in this$
-spec alone; it is not charged with checking historical O2 emissions it was not given.$
-The handoff includes the one-sided material-compatibility actions in all four named energy rows;$
+**O2 emission provenance and comparison limits (O2-R §§4–5, 8).** Retention against historical$
+O2 emissions is comparator/record scope, outside engine construction from this spec. The historical$
+content includes every named operand and complete live-object dependence from either O2 engine,$
+the eight keys at every listed location, and PY-only native/chart and core/material-compatibility$
+content. The explicit engine inputs above retain general admissible dependences independently of$
+that finite inventory. The O2 difference ledger is unresolved wherever it bears on a condition.$
+It includes one-sided material-compatibility actions in all four named energy rows;$
 one-sided material density/flux actions in `energy_power`; four density/flux orientation differences$
 in `energy_balance`; and density/flux occurrences inside `OPEN_JointPowerAccounting` in$
 `energy_power` and `ℬ_E^steady` alongside separate `energy_balance` occurrences. Those differences$
 travel with the force/power-compatibility and duplicate-power questions. No limited inventory match$
-supplies full OPEN equality, compatibility or a duplicate-power conclusion. This check discharges$
-neither the returned interpretation nor any unsupplied physical input.$
-$
-**Register handoff (D5; O2-R §8).** The later record carries the adopted-premise provenance and the$
-undischarged interpretation. A register entry follows only from an established sourced requirement$
-on which the conditional object rests: material/phase identifications or expressly carried$
-accounting/admissibility obligations. Unsupplied closure forms remain OPEN handoffs. An unresolved$
-engine difference alone creates no requirement. S21 owns the later integration/sort; this spec$
-performs no register edit or physical reconciliation.$
-$
-## Engines, review, scope$
-$
-- **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).$
-- **Spec review.** A fresh non-author Claude agent and Grok review v10 until clear, after this$
-  authoring stop (D7). The orchestrator runs those reviews.$
-- **Build review.** Codex-written, so a fresh Claude agent and Grok, each with a mandatory FORM ablation.$
+supplies full OPEN equality, compatibility or a duplicate-power conclusion. The returned interpretation$
+and unsupplied physical inputs remain OPEN.$
+$
+**Requirement status (D5; O2-R §8).** An established sourced requirement is a condition on which the$
+conditional object rests: a material/phase identification or an expressly carried accounting/$
+admissibility obligation. Unsupplied closure forms remain OPEN. An unresolved engine difference$
+alone supplies no requirement.$
+$
+## Model point and scope$
+$
 - **Optical model point.** As in "Setting": leading eikonal with the retained multigraded set above.$
   Part D's mechanical and inherited energy content retains its stated untruncated domain. Its results$
   are conditional on P5's adopted momentum carriage and P6's adopted steady pressure-plus-linear-viscous$
@@ -525,16 +549,3 @@ performs no register edit or physical reconciliation.$
   - a drain that changes while light crosses;$
   - polarization transport;$
   - anisotropic or coupled branches.$
-- **Deferred to the build** (implementation, not new physics; the build directive owns it):$
-  - the symbolic handling of the every-`b` requirement;$
-  - component/measure calculus on the supplied coordinate mass law, with amended D4's reading,$
-    live-gradient and relative-to-`j_n` bound requirements above for any induced-measure claim;$
-  - representation of general OPEN actions and executable controls, including FORM ablation (E2).$
-- **Stop and report**, without choosing, when any of these happens:$
-  - a second method failure;$
-  - a sub-problem this spec does not name;$
-  - a premise this spec does not supply.$
-- **The step record** interprets the results.$
-$
-**Authoring STOP.** Repair v10 in place and report changes for repair items 1–4, source conflicts$
-and missing sourced pieces. No CAS, build, review launch, commit, push or spawned agent is part of this task.$
```

Exit code: `0`.

```bash
bash -o pipefail -c 'git show 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md | nl -ba | sed -n '"'"'27,34p;59,68p;76,90p'"'"' | cat -E'
```

```text
    27	**New for S9b:**$
    28	$
    29	| Premise | Content | Qualification carried with it |$
    30	|---|---|---|$
    31	| **P5** (2026-10-07; carriage 2026-10-08) | The brane's in-plane momentum density is `ρ_br V`, and that momentum is carried with the brane material at `V`. Its contribution to the in-plane momentum current is `ρ_br V^i V^j`. | Scoped to Part D. Any other in-plane momentum current stays OPEN. It is counted beside the supplied current and P6's stress, which it does not contain again. Part D states what the OPEN pieces must supply. If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |$
    32	| **P6** (2026-10-08, revised the same day) | In steady flow the brane's in-plane stress is an isotropic pressure plus a linear viscous stress: `T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`, in the Cauchy convention (traction `T^{ij} n_j`; force density `∂_j T^{ij}`). `p_br` is a general function of `ρ_br` only. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. The shear viscosity `η` and bulk viscosity `ζ` are live general profiles. In the optical regime, light still sees the elastic transverse stiffness `μ_⊥`. | Scoped to Part D. It supplies the steady in-plane part of O2's stress `𝒯_br^live` only. It is an adopted steady form, not a consequence of P1: P1's relaxation response stays OPEN outside it. The linear viscous form excludes power-law creep. The power this stress dissipates is not identified with any O2 energy operand. It stays in the OPEN energy accounting, with its physical supplier and budget OPEN. The first 2026-10-08 row (pressure only, justified by relaxation under steady load) is withdrawn: steady flow keeps straining, so relaxation does not remove the viscous stress (spec v10 review C2). |$
    33	| **Density link** (2026-10-06) | v8's `c_γ² = μ_⊥/ρ_br` is kept, with `μ_⊥ ∝ ρ_br^α` and `α` one live symbol. This is the user's selected live exponent. | Scoped to Part D. There it supplies `μ_⊥` as a function of `ρ_br` alone. Its other O1 `ℳ_⊥` dependences are set aside by this premise, and results are flagged. Parts A–C keep `μ_⊥(x)` and `ρ_br(x)` independent, as in v8. |$
    34	| **w-parity** (2026-10-07) | For an electrically neutral mass, the far field is symmetric under `w → −w`. So `ξ_w` and the material `w`-velocity vanish there. | This is a **labelled neutral-sector restriction**. Parts A–C keep `ξ_w` live as in v8, which covers the charged case. Part D is stated for the neutral sector, and prints the label. |$
    59	## D4. Induced metric and order$
    60	- **v8's sentence.** v8 says the induced-metric mass balance differs from the flat form by "a relative `O(ε)`$
    61	  correction to the implied `j_n`". No source supplies that order. The gradient-scale condition$
    62	  `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok) does not supply it either (spec v10 review C4).$
    63	- **One rule for v10.** The supplied mass balance is on the coordinate `d³x` measure (O2 record §8). An induced-measure$
    64	  claim names which reading it uses: the same densities re-expressed per induced volume, or a mass law imposed on the$
    65	  induced measure. Such a claim keeps `∂_r[(∂ξ_w)²]` live. It attaches no relative order to the correction to `j_n`$
    66	  unless it states a condition that bounds that correction relative to `j_n` itself.$
    67	- **Neutral Part D.** With `ξ_w = 0` (w-parity), the supplied `g_ij` reduces to `δ_ij`, and the engines print that$
    68	  reduction. This limits the claim; it adds no term.$
    76	**Routing:**$
    77	- **The obligation stays the orchestrator's** (M1; record §8). S9b does not discharge it. The S9b record carries it$
    78	  as undischarged. It is distinct from any step's ownership of a physical input.$
    79	- **Retained keys.** Part D retains all eight keys, at every §4 location, for every piece of O2 content that no$
    80	  adopted premise supplies (record §8).$
    81	- **Supplied content.** Where P5, P6 or the density link supplies an object, the premise states that object's$
    82	  dependence, and every result that uses it is flagged.$
    83	  - The spec author determines, from the O2 sources, which content each premise supplies.$
    84	  - The author does not compose the balance (D3).$
    85	  - The user chose P6 over keeping the stress OPEN (2026-10-08).$
    86	- **Register.** A register entry follows only if a sourced requirement is established. An unresolved difference$
    87	  alone does not create one (record's rest-on criterion).$
    88	$
    89	## D6. Owed build findings: routed to the repair build directive$
    90	These come from the build review preserved at `bb94b885`, as listed in that commit's message. They concern implementation, not spec physics. They go to the$
```

Exit code: `0`.

```bash
bash -o pipefail -c 'git show 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md | nl -ba | sed -n '"'"'62,97p;115,123p;469,515p'"'"' | cat -E'
```

```text
    62	The engine acceptance reports the computed metric/inverse and `det g = 1+ξ′²`, graph normal from$
    63	tangents, graph material velocity with `U^w = V_r ξ′` and zero graph-normal component, `(V·∇)U`,$
    64	the unsolved coordinate-measure mass residual, and in-plane carried `j_n V^i`. It also explicitly$
    65	keeps momentum density/flux, internal force, native reductions, bulk carry and energy inputs OPEN.$
    66	These are the acceptance disposition's scoped statements, retrieved by M1, not new derivations.$
    67	The comparator's `metric_determinant`, geometry/velocity container leaves, mass rows, profile and$
    68	gradient rows print their own results (M8). The material-acceleration vector has no emitted WL$
    69	counterpart and is unjoined (M9); the engine acceptance does not make it a cross-engine comparator$
    70	agreement claim.$
    71	$
    72	## 2. Exact domain, adopted premises and supplied inputs$
    73	$
    74	The supplied setting is one isolated, spherically symmetric mass at rest with its drain flowing:$
    75	far field, Eulerian steady material profiles, lab time in the brane's far-field rest frame, linear$
    76	optical waves and leading eikonal. There are three in-plane coordinate directions `x^i` and bulk$
    77	coordinate `w`. Scalar profiles are general live functions of `r=|x|` on the far-field `r>0` domain, with radial$
    78	`V^i=V_r(r)x^i/r`. The material velocity is that of the material displaced by `u`, distinct from$
    79	wave velocity, outward face velocity `V_s`, bulk-normal drain `v_dr` and native bulk velocity.$
    80	The steady graph determines its centre graph velocity; it does not identify finite-thickness face$
    81	velocities or the reduced bulk-direction exchange. Material identity is **recorded** on the$
    82	homogeneous displacement anchor and a **supplied S9b identification**, not a live kinetic law.$
    83	The historical ordered/disordered split of one conserved material is **postulated**; A13 remains$
    84	OPEN. The static localized reduction is conditional on its **postulated** parent sector.$
    85	Sources: spec §§1, 5, 7–8; contract §§1–2, 5–7 (M1). WL separately emits its `r>0` inequality;$
    86	PY has no paired inequality emission (M9), so domain use is not an inequality-comparison result.$
    87	$
    88	Every density in this object is per coordinate `d³x`, including mass/source, material momentum,$
    89	exchange, loads, energy and power; optical `μ_⊥` uses the same measure as `ρ_br`. Induced metric and$
    90	native area factors remain explicit. The supplied geometric/optical/mass inputs are$
    91	`g_ij=δ_ij+∂_iξ_w∂_jξ_w`, its inverse, `ξ_w=ℓh`,$
    92	`c_γ²≡μ_⊥/ρ_br`, `c_γ=c₀(1+δ)` and `∇·(ρ_br V)=−j_n`.$
    93	LAB_HELD is supplied spatial speed anchoring, not a holder or reference-evolution law. The bulk$
    94	inputs `P=Kρ^n`, `c_s²=nKρ^(n−1)/m`, `f=ρ/ρ₀−1` use bulk number density, particle mass and symbolic$
    95	EOS exponent. They supply no bulk profile, brane-density response or live traction. These are inputs$
    96	on their recorded domains, not results earned by engine agreement. Sources: spec §§1, 3.1, 7–8;$
    97	contract §§5–7 (M1).$
   115	Supplied counting is `ε=GM/(c₀²r)`, `δ=O(ε)`, `(∂ξ_w)²=O(ε)`, `V/c₀=O(ε^{1/2})`,$
   116	`(V/c₀)²=O(ε)` and the optical monomial box `0≤a≤1, 0≤b≤2, 0≤c≤1` in$
   117	`δ^a(V/c₀)^b((∂ξ_w)²)^c`, with nonnegative integer indices. Only the stiffness/density ratio$
   118	inherits the speed-change grade. Individual density/modulus, stress/inertia/normal response,$
   119	exchange/source/load, relaxation/power, holder/embedding and derivative grades remain OPEN. O2 is$
   120	untruncated; this optical box removes no mechanical or energy term. Bulk first order in `f` is$
   121	separate, with no supplied `f`–`ε` relation. The mass law carries its recorded relative-`O(ε)`$
   122	qualification for claims transferring `j_n` or `ρ_br` to induced measure; no induced-measure mass$
   123	law is supplied or derived. Sources: spec §7; contract §9 (M1).$
   469	Part D must retain **every named operand and complete live-object dependence printed by either engine**,$
   470	including all eight WL-only derivative keys at every location in §4 and all PY-only native/chart,$
   471	core/material-compatibility content. It must also preserve the spec's general admissible live fields,$
   472	gradients, entire material history and other native dependences (spec L124–125). Neither engine's$
   473	finite OPEN inventory is the object's dependence list, and even taking both lists cannot impose a$
   474	closed argument list or derivative/history cutoff. The full inventories/operands remain in the$
   475	filed production stream; their differences are M6. No surrounding held calculus, coefficient algebra,$
   476	position or multiplicity outside §5's inventory scope becomes a comparison result.$
   477	$
   478	The mass-law qualification travels explicitly: `∇·(ρ_br V)=−j_n` is a **coordinate-`d³x` measure**$
   479	input. For a claim reading `j_n` or `ρ_br` per induced measure, or comparing an induced-metric mass$
   480	law, the recorded qualification is a **relative `O(ε)` correction to `j_n`**, not a live O6 law;$
   481	O2 supplies and derives **no induced-metric mass balance** (spec L330–336). Fixed-`ℓ` counting,$
   482	independent orbital `GM`, only the stiffness/density ratio inheriting the speed-change grade, and$
   483	first order in `f` with no `f`–`ε` relation also travel unchanged. Historical static/homogeneous,$
   484	frozen/uniform or supplied-profile relations retain exactly §2's restricted domains.$
   485	$
   486	Part D must keep coordinate versus graph-normal content distinct, use the same material and measure$
   487	with compatible native maps, retain untruncated unknown terms/grades, avoid duplicate stress/exchange/$
   488	power content, and name the supplier/budget if a closure requires net power. Residual `0` at the$
   489	balance closed-part level (printed, M8), unioned balance-role inventory matches (printed, M6),$
   490	and joint stored classification tuple count matches (record literal counts, M7) supply neither$
   491	a solved material state nor full OPEN equality. Each condition inherits the whole difference$
   492	ledger in §4 wherever a difference bears on its content, including native-map/measure, material$
   493	compatibility, named/live dependences and generalized-work scope (M6/M9/M12).$
   494	In particular, the force/power-compatibility and duplicate-power conditions carry the one-sided$
   495	`OPEN_MaterialCompatibility` actions in all four energy rows, the one-sided density/flux actions in$
   496	`energy_power`, **and** the four paired density/flux orientation differences `[["argument",1]]`$
   497	in `energy_balance` (line 207). The stored PY density/flux actions inside `OPEN_JointPowerAccounting`$
   498	travel with those conditions in both `energy_power` and `ℬ_E^steady` itself, alongside the separate$
   499	occurrences in `energy_balance` (lines 202/206; §4; M12). The record's literal stored classification tuple$
   500	count match (M7) and the comparator's unioned balance-role inventory matches (M6)$
   501	do not settle compatibility or whether the complete work/energy content is counted once.$
   502	Part D must carry these differences alongside those conditions; it may not discard or reconcile$
   503	them, or judge either engine to duplicate power from either limited match. The `xi_w''` difference$
   504	also travels with the admissible-dependence and force/power conditions as **OPEN interpretation,$
   505	undischarged at sub-step 7 and returned to the orchestrator under M1** (§4). Any physical resolution$
   506	requires work beyond retrieval.$
   507	Source: spec §§1–10, contract §§1–10 and M6–M9/M12.$
   508	$
   509	In particular, O2 supplies **no momentum-density map `ρ_br V`**. Its open inertia/momentum mapping$
   510	cannot be cancelled or converted into a chosen advective kinetic law using the mass relation.$
   511	The directive's folded disposition keeps the separate user premise for S9b repair at that later$
   512	decision-list gate (M1). Part D must source any additional kinetic, constitutive, projection,$
   513	profile, support, order or energy input through its own authorized premise/specification work.$
   514	It inherits no `GM` normalization, drain-to-mouth identification, static flowing extension or$
   515	light-observable compatibility from O2.$
```

Exit code: `0`.

```bash
bash -o pipefail -c 'git show 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | nl -ba | sed -n '"'"'49,79p;103,114p;234,248p;330,347p'"'"' | cat -E'
```

```text
    49	The component object has three in-plane directions and normal material content on the supplied$
    50	embedded brane. State the component basis and its relation to the far-field coordinates `x^i` and$
    51	bulk direction `w`. Distinguish coordinate projections from projections onto the embedded graph's$
    52	normal. Geometric component changes may be constructed from the supplied graph; a component basis$
    53	does not choose a sharp-sheet or finite-slab material reduction. Where a material or face quantity$
    54	cannot be projected without O6 or a live material response, its component action stays explicitly$
    55	unevaluated with that operand named.$
    56	$
    57	**Declared measure (a convention).** Every density in this object is a density per coordinate volume$
    58	`d³x` of the far-field coordinates `x^i`, the measure on which the supplied mass law's divergence is$
    59	written (§3.1). This covers `ρ_br`, `j_n`, material momentum storage and transport, the carried and$
    60	other exchange occurrences, and the force, load, energy and power densities. In the optical ratio$
    61	`c_γ² ≡ μ_⊥/ρ_br`, `μ_⊥` is referred to the same measure as `ρ_br`. Induced-metric factors, such as$
    62	`g_ij`, `g^{ij}`, `det g_ij` and the graph normal, enter explicitly in the geometric actions; no$
    63	occurrence is moved to another measure. A quantity defined per native face area enters through its$
    64	native geometric factor and `𝒥_map` (§5). This fixes the volume measure of densities only. It$
    65	selects no stress measure, supplies no induced-measure mass balance, and does not remove the supplied$
    66	law's recorded qualification (§7).$
    67	$
    68	`V` is the background in-plane velocity of the material displaced by light's `u`. It is distinct$
    69	from the perturbation velocity, outward face velocity `V_s`, bulk-normal drain `v_dr`, and a native$
    70	bulk velocity. No value for a normal material response or identification among these velocities$
    71	follows from those names (C §§2, 5–7). On the steady supplied graph `w = ξ_w` (§3.1), a brane$
    72	material point with in-plane coordinate velocity `V^i` has the bulk-direction coordinate velocity$
    73	that the graph and `V` determine. Construct that component from them; it stays live through `V` and$
    74	`ξ_w`, and there is no separate bulk-direction velocity operand. Normal velocity content that the$
    75	graph does not determine belongs to `𝒩_br^live` and `𝒥_map`. One example is different local$
    76	velocities at the faces of a finite-thickness realization, whose centre and thickness content `ξ_w`$
    77	does not fix. Such content stays an unevaluated action with those operands named. This material$
    78	velocity enters the material momentum entries (§4) and the force/power pairings (§6); the$
    79	bulk-direction carried exchange momentum is the OPEN reduction of §5.$
   103	## 3. Live input register$
   104	$
   105	### 3.1 Material identity, geometry, optical identifications and mass source$
   106	$
   107	| Status and source in C | Equation or operand | Domain and use |$
   108	| --- | --- | --- |$
   109	| **Recorded material identity; supplied S9b identification**, C §2 | `u` displaces the material carrying `ρ_br`; `V` is its live background in-plane velocity. | Fixes the material whose momentum is accounted for. The homogeneous kinetic/continuity anchor is reported in §8.1; it supplies no flowing kinetic law. |$
   110	| **Supplied geometric identification**, C §5 | `g_ij = δ_ij + ∂_iξ_w ∂_jξ_w`, `g^{ij} = (g_ij)⁻¹`, `ξ_w = ℓh`. | Induced spatial metric on the recorded graph and retained L3 field identity. `ξ_w` is the brane's displacement into the bulk direction `w`, so the supplied graph is `w = ξ_w` over the far-field coordinates `x^i`. `ξ_w` has length, `h` is dimensionless, and `ℓ` is the fixed reduction scale, not a selected slab width. Geometry and slopes remain live. |$
   111	| **Supplied optical-regime live identifications**, C §6 | `c_γ(r)² ≡ μ_⊥(r)/ρ_br(r)`, `c_γ(r) ≡ c₀[1+δ(r)]`. | `μ_⊥` is premise 1's optical elastic stiffness. The ratio supplies no full stress, steady-load stiffness or inertial law. `c₀` is the supplied asymptotic light speed. |$
   112	| **Supplied steady mass balance**, C §6 | `∇·(ρ_br V) = −j_n`. | Use it as written on the declared measure (§1), retaining the full density factor and derivatives in the v9 source convention. This is the live mass input, not a proof of O6. A finite-slab or induced-measure replacement is not supplied by this equation. Its recorded relative-`O(ε)` qualification is a limit on claims, not a term in the object (§7). |$
   113	| **Supplied speed-profile anchoring**, C §6 | `Q_bg^L(x,t) = Q_bg(x)`, `Q_bg^M(x,t) = Q_bg(χ(x,t))`; v9 selects LAB_HELD for `c_γ`. | These are distinct physical anchorings on the recorded supplied background. `χ(x,t)` is the inverse material map, not `χ_B`. LAB_HELD does not impose a material-reference law or supply a holder. |$
   114	| **Supplied bulk-density inputs**, C §6 | `P = Kρ^n`, `c_s² = nKρ^(n−1)/m`, `f(r) ≡ ρ(r)/ρ₀ − 1`. | `ρ` is bulk number density, `m` particle mass and `n` a symbolic EOS exponent. The profile is unsolved; the recorded Part C response is first order in `f`. Neither `ρ_br(f)` nor a live normal traction follows. |$
   234	On a native face description, transported momentum uses the native relative exchanged mass current$
   235	and the local brane material velocity of the converted material. A full normal/component reduction$
   236	retains `𝒥_map`, the live geometry and any unsupplied face-to-material velocity identification. Do$
   237	not choose an independent converted-material velocity, replace it with `v_dr`, or identify it with a$
   238	native bulk velocity. The premise-3 identification applies to the material after the specified$
   239	conversion/transfer; it does not fix how a native bulk current acquires that momentum. Any additional$
   240	non-variational conversion momentum partner and its reaction system remain OPEN with S12.$
   241	$
   242	**Bulk-direction carried component.** Premise 3 applies to each native transfer at its own location.$
   243	It does not supply the reduced bulk-direction component `(Π_n^carry)^w` as `j_n` times one$
   244	bulk-direction velocity. In a finite-thickness realization, the transfers at different faces need$
   245	not share one current or one local bulk-direction velocity, and premise 3 supplies no such equality.$
   246	The reduced component therefore depends on the OPEN `𝒥_map` and the normal material response$
   247	`𝒩_br^live`, with the live geometry and any unsupplied face-to-material velocity identification. Print$
   248	it as an unevaluated OPEN action with those operands named, in the same outward convention.$
   330	`GM` is the independent slow-test-matter orbital parameter, not a source, profile, mouth or drain$
   331	amplitude. Only `μ_⊥/ρ_br` inherits the speed-change grade. Keep the full density factor and$
   332	derivatives in the mass input. The recorded induced-metric qualification is a relative `O(ε)`$
   333	correction to `j_n`, not a live O6 law. O2 uses the law on the declared coordinate measure (§1) and$
   334	supplies no induced-metric mass balance. The qualification is carried as a limit on claims, not as a$
   335	term in the object: a claim that reads this `j_n` or `ρ_br` as a density per induced measure, or$
   336	compares it with an induced-metric mass balance, carries it. Fixed `ℓ` transfers the live `h` slope$
   337	to the geometric grade without supplying a mouth-amplitude grade. First order in bulk `f` is a$
   338	separate recorded response domain; no relation between `f` and `ε` is supplied.$
   339	$
   340	**OPEN grades and derivative scales (C §9):** individual density and stiffness; inertia; stress and$
   341	normal response; reference/relaxation and power; source and exchange momentum; force and$
   342	traction/support; holder and mouth data; embedding-sector coefficients/source and longitudinal$
   343	field, wherever no grade is recorded. Retain these as named unknown grades/scales attached to their$
   344	inputs. Do not assign them a convenient higher order or separate density/modulus `O(ε)` variations.$
   345	$
   346	O2 is **untruncated**. The optical monomial box records the model's optical counting; it is not an$
   347	order contract for removing material, exchange, boundary or power terms from O2.$
```

Exit code: `0`.

```bash
bash -o pipefail -c 'git show 217a92e9:research/pde_ledger_v3/directives/O2_input_contract.md | nl -ba | sed -n '"'"'107,114p;197,211p;247,253p;332,347p'"'"' | cat -E'
```

```text
   107	| Operand | Physical input whose form is missing | Owner / boundary |$
   108	| --- | --- | --- |$
   109	| `𝒯_br^live` | Full live brane stress, including any relaxation contribution to the in-plane and normal force content. | O2 carries the operand; conservative/material antecedents S1.5/S8 and nonlinear material completion S22. Relaxation/reference ownership remains unassigned (§4). |$
   110	| `𝒫_br^cons` | Conservative brane material momentum antecedents. | S1.5 supplies conservative substrate balance antecedents; S8 supplies brane material content. No live reduction is supplied. |$
   111	| `𝒯_br^cons` | Conservative in-plane and normal-stress antecedents from S1.5/S8. | S1.5/S8; nonlinear material completion S22. These do not supply the full `𝒯_br^live`; no stress measure, symmetry or constitutive form is selected here. |$
   112	| `ℐ_br^live` | Inertial/kinetic response on the flowing, embedded material background. | S8, using S1.5 antecedents when available. The optical density identification alone does not supply this response. |$
   113	| `𝒩_br^live` | Normal material response, retaining distinctions among embedding, centre and thickness content until a model connects them. | S5–S8 for wall/width/compression ingredients; Q1/Q2 for their recorded static embedding sector; live identification with O4 remains OPEN. |$
   114	| `𝒜_rot^live` | Required internal angular-momentum/couple-stress content and the physical rotational reference frame, including whether such content is present. | S8's OPEN requirements; no carrier or frame is chosen. |$
   197	## 5. Geometry, embedding and O4$
   198	$
   199	**Status: supplied geometric identification.**$
   200	$
   201	```text$
   202	g_ij = δ_ij + ∂_iξ_w ∂_jξ_w ,       g^{ij} = (g_ij)⁻¹ ,       ξ_w = ℓh .$
   203	```$
   204	$
   205	**Source/domain:** v9:19–21,46–48,216–228; the induced spatial metric on its supplied graph and the$
   206	retained L3 field identity. Underlying identity:$
   207	v2:`stages/ledger_stage031_puncture_deflection_field_identity_source.md:60–76`, within the postulated$
   208	G0 sector. `ξ_w` is a length, `h` dimensionless and `ℓ` that reduction's fixed scale; `ℓ` is not a newly$
   209	selected slab width. `ξ_w`, `h` and their spatial derivatives remain live. `ζ_c` and `W` are independent$
   210	face-centre/thickness variables in S11b, rather than replacements for `ξ_w`$
   211	(`directives/S11b_SHARED_PHYSICS.md:89–96`). **Owner:** static identity Q1/Q2; live applicability O4.$
   247	**Status: OPEN `ℰ_h^live` (O4).** It denotes the missing live embedding/longitudinal relation with flow,$
   248	exchange, variable coefficients and their gradients (v9:295–300). The recorded action supplies no map$
   249	from `u_L` to steady `V`, or from its mouth source to the mass's drain. **Relation to O2's normal content:$
   250	unsettled**, as fixed by `O2_premise_decision_list.md:46–47`. O4 remains separately named as a coupled$
   251	input; this contract chooses neither an identification with `𝒩_br^live` nor an independent equation$
   252	count. **Owner:** recorded static ingredients Q1/Q2; the live identification remains O4 for the later$
   253	O2 spec, with material completion at its S8/S22 owners.$
   332	**Status: OPEN `𝒥_map` (O6).** The material/order weighting, live projection/window and measure$
   333	identifications connecting sheet `j_n`, slab fluxes and bulk-normal `v_dr`, with bulk/return data, remain$
   334	general (v9:305–308; `O2_premise_decision_list.md:50`). The window in$
   335	`directives/S11c_a_SHARED_PHYSICS.md:394–402` is tied to its own face maps; it is not a newly selected$
   336	O2 projection. Stage 006's shear projection$
   337	$
   338	```text$
   339	μ_R = ∫ dw χ_B μ_R⁽⁴⁾$
   340	```$
   341	$
   342	is **postulated/PENDING**, with dimensional consistency asserted only$
   343	(v2:`stages/ledger_stage006_two_phase_chiB_ontology.md:98`); it is not a supplied `ℳ_⊥` law or a complete$
   344	live map. **Owner:** S12 for dynamical conversion and separate source/boundary inventories$
   345	(`V3_STEP_PLAN.md:579–617`; `steps/S11c_PARTIAL_CLOSEOUT.md:21,33`); S14a retains the distinct$
   346	projected order-loss/far-field flux bridge (`V3_STEP_PLAN.md:626–641`). Sub-step 3 must state its use of$
   347	O6 in the balance; this contract does not perform that reduction.$
```

Exit code: `0`.

```bash
bash -o pipefail -c 'nl -ba research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n '"'"'1,21p;58,90p;278,284p;400,439p;477,480p;519,552p'"'"' | cat -E'
```

```text
     1	# S9b — what brane light needs in order to bend and be delayed like GR (question spec, v10)$
     2	$
     3	**Version:** v10, repair 2 (2026-10-08).$
     4	$
     5	**Objects:** the v8 optical objects in Parts A–C and the neutral, linked steady in-plane$
     6	brane balance and its conditional optical requirements in Part D, retaining every unsupplied O2 input.$
     7	$
     8	**Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions$
     9	are `directives/S9b_repair_decision_list.md` at `79055918` (amendment 1). v9's Part D is not an input. Paths$
    10	below are relative to `research/pde_ledger_v3/`. The additional source abbreviations are:$
    11	$
    12	- **O2-R:** `steps/O2_steady_brane_balance.md`, accepted at `72866fcf`, especially §§4 and 8.$
    13	- **O2-S:** `directives/O2_SHARED_PHYSICS.md` at `4680e251`.$
    14	- **O2-C:** `directives/O2_input_contract.md` at `217a92e9`.$
    15	- **O2-P:** `directives/O2_premise_decision_list.md` at `77d2c39a`.$
    16	$
    17	Equations labelled **supplied** are conditional inputs that the computation cannot test. Adopted$
    18	premises retain that status and their provenance. Every$
    19	dependent result is flagged with the supplied identification or adopted premise it uses. Reference$
    20	objects are comparison inputs, rather than premises for deriving the brane's observables. The spec$
    21	names computed objects and solution conditions; it supplies no outcome or expected-value acceptance test.$
    58	- **Steady brane mass balance.** The supplied balance is$
    59	$
    60	  ```$
    61	  ∇·(ρ_br(x)V(x)) = −j_n(x) .$
    62	  ```$
    63	$
    64	  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`$
    65	  coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as$
    66	  `ρ_br`. This input supplies no induced-measure or finite-slab replacement mass law. Every$
    67	  induced-measure claim names which of the following readings it uses (D4):$
    68	$
    69	  - **Density re-expression.** The same mass and normal exchange are re-expressed as densities per$
    70	    induced volume. The measure and same-content identifications defining this reading are$
    71	$
    72	    ```$
    73	    dvol_g ≡ √det(g_ij) d³x ,$
    74	    ρ_br^(g) dvol_g ≡ ρ_br d³x ,      j_n^(g) dvol_g ≡ j_n d³x .$
    75	    ```$
    76	$
    77	    These are supplied identifications of the reading, using the supplied metric; they do not$
    78	    replace the mass law. Print the re-expressed densities and correction objects with their$
    79	    orders under the supplied metric and slope counting. No additional independent bound on the$
    80	    drain divergence is a premise of this density re-expression. The imposed-law restriction$
    81	    below does not apply to this reading.$
    82	  - **Imposed mass law.** An independently imposed mass law on the induced measure, without the$
    83	    same-content density re-expression above, is a separate physical input; none is supplied here.$
    84	    A comparison on this reading names that law and keeps$
    85	    `∂_r[(∂ξ_w)²]` live. A relative order for its correction to the implied `j_n` requires a$
    86	    condition bounding that correction relative to `j_n` itself. **For this imposed-law reading$
    87	    only**, the unqualified recorded relative-order statement in O2-R L478–480 is not carried.$
    88	    The supplied slope counting and the historical gradient-scale condition$
    89	    `∂_r[(∂ξ_w)²] = O(ε/r)` supply no such uniform drain-relative bound for an imposed-law$
    90	    comparison (D4). The historical derivative condition is not imposed here.$
   278	`V_r`, `ρ_br`, `j_n` and all unsupplied responses remain general live profiles/actions. Eulerian$
   279	steadiness supplies no material constancy or history cutoff. Coordinate momentum storage, transport,$
   280	exchange, force, load, energy and power densities use the same coordinate `d³x` measure as the mass$
   281	law. Native face-area factors and reductions remain explicit through O6. Coordinate `w` content and$
   282	graph-normal content are distinct. The centre graph fixes neither finite-thickness face nor interior$
   283	normal velocity content; any content it does not determine stays with `𝒩_br^live` and `𝒥_map`.$
   284	$
   400	$
   401	**w-parity — neutral sector (2026-10-07; D2, D4).** For an electrically neutral mass the adopted$
   402	far-field `w → −w` symmetry restricts the graph displacement and the constructed centre-graph$
   403	material velocity. Define the supplied graph embedding and its material tangent lift by$
   404	$
   405	```$
   406	X_graph(x) ≡ (x¹,x²,x³,ξ_w(x)) ,      U_graph ≡ V^i ∂_i X_graph .$
   407	```$
   408	$
   409	`U_graph^w` denotes the bulk-coordinate component of this centre-graph material velocity, constructed$
   410	from the steady graph and `V` (O2-S L68–79; O2-R §§1–2). It is not an independent normal velocity$
   411	operand. The adopted neutral-sector equations are$
   412	$
   413	```$
   414	ξ_w(r) ≡ 0 ,      U_graph^w(r) ≡ 0 .$
   415	```$
   416	$
   417	Print **neutral-sector restriction** with Part D and each dependent result, including the reduction$
   418	of the supplied `g_ij`. Parts A–C keep `ξ_w` live, including the charged case. This restriction adds$
   419	no force term. Normal velocity content not determined by the centre graph remains OPEN in$
   420	`𝒩_br^live` and `𝒥_map`, including face and interior content and unsupplied face-to-material$
   421	identifications. The reduced bulk-direction carry `(Π_n^carry)^w` remains an OPEN action of those$
   422	operands with native geometry; it is not supplied as `j_n` times one centre-graph velocity$
   423	(O2-S §5). No native-face, thickness, exchange-map or normal constitutive law is supplied.$
   424	$
   425	### O2 balance pieces that remain live$
   426	$
   427	The engines compose the balance from these physical roles and the supplied equations above.$
   428	No assembled momentum balance, total transport tensor, cancellation or profile solution is supplied$
   429	here (D3; O2-S §§3–6; O2-R §8). The supplied content and remaining OPEN actions are:$
   430	$
   431	| O2 content | Supplied content in Part D | Content retained as OPEN |$
   432	|---|---|---|$
   433	| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density and the material-carried contribution `ρ_br V^i V^j` to the in-plane current. | Any other in-plane momentum current, with general live dependence, counted beside the supplied current and P6 stress and containing neither again; normal/embedded responses and unresolved relations wherever unsupplied. The conservative antecedent is named, rather than added as another momentum species. |$
   434	| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the full steady in-plane stress in its adopted pressure-plus-linear-viscous form, with the stated Cauchy traction/force convention. | Normal material response, conservative antecedents, reference evolution outside the adopted form, and rotational/couple/frame content wherever unsupplied. Their overlapping descriptions remain explicit within one material accounting object; no additional in-plane stress duplicates P6. |$
   435	| Optical/material constitutive inputs; O1 `ℳ_⊥`, O7 `ℛ_br` | The density link supplies Part D's optical stiffness response only. | Brane-density response and every other unsupplied material or reduction identification. Parts A–C retain the v8 local profiles. |$
   436	| Geometry, O4 `ℰ_h^live`, O6 `𝒥_map` | Supplied centre graph, `ξ_w=ℓh`, metric, constructed `U_graph` and neutral centre-graph restriction. | Live embedding/longitudinal relation and its unsettled identity with O2 normal content; native face geometry/measure, finite-thickness reduction, face/interior normal velocity content not determined by the graph, normal response and material/order/projection map. No second normal equation is imposed. |$
   437	| Mechanical face/support loading; `T_hold,s`, `𝒯_bulk,n,s^live` | P4 supplies the bulk part's native direction only. | Full load/support partition, native normal-load amplitude, projections/maps and application-point velocities. No external support is selected by LAB_HELD. |$
   438	| Exchange; O3 `Π_n`, S12 partners/reactions | P3 supplies the coordinate in-plane carry `j_n V^i`. | Additional momentum partners/reactions and bulk-direction carry with O6/normal material content, native geometry and unsupplied face-to-material velocities. Centre-graph parity supplies no reduced bulk-carry formula. Native convective transfer and the same O3 carry are not counted twice. |$
   439	| Boundary/source/core inputs; O5 `ℋ_core`, `𝔅_A13` | P2 supplies the drain-drive representation only. | A13 branch; distinct local conversion/controller and boundary/domain inventories; physical core holder/mouth data and their response. Q2/S22 retain O5; S12 retains conversion and return partners. |$
   477	**D5; O2-R §§4, 8.** The WL-only `ξ_w''` interpretation remains an **undischarged sub-step-7$
   478	obligation returned to the orchestrator under M1**. This is an unresolved input qualification,$
   479	distinct from ownership of a physical input. Each dependent Part D condition carries that status;$
   480	no physical interpretation or reconciliation is supplied. The dependence inputs are explicit below.$
   519	**O2 emission provenance and comparison limits (O2-R §§4–5, 8).** Retention against historical$
   520	O2 emissions is comparator/record scope, outside engine construction from this spec. The historical$
   521	content includes every named operand and complete live-object dependence from either O2 engine,$
   522	the eight keys at every listed location, and PY-only native/chart and core/material-compatibility$
   523	content. The explicit engine inputs above retain general admissible dependences independently of$
   524	that finite inventory. The O2 difference ledger is unresolved wherever it bears on a condition.$
   525	It includes one-sided material-compatibility actions in all four named energy rows;$
   526	one-sided material density/flux actions in `energy_power`; four density/flux orientation differences$
   527	in `energy_balance`; and density/flux occurrences inside `OPEN_JointPowerAccounting` in$
   528	`energy_power` and `ℬ_E^steady` alongside separate `energy_balance` occurrences. Those differences$
   529	travel with the force/power-compatibility and duplicate-power questions. No limited inventory match$
   530	supplies full OPEN equality, compatibility or a duplicate-power conclusion. The returned interpretation$
   531	and unsupplied physical inputs remain OPEN.$
   532	$
   533	**Requirement status (D5; O2-R §8).** An established sourced requirement is a condition on which the$
   534	conditional object rests: a material/phase identification or an expressly carried accounting/$
   535	admissibility obligation. Unsupplied closure forms remain OPEN. An unresolved engine difference$
   536	alone supplies no requirement.$
   537	$
   538	## Model point and scope$
   539	$
   540	- **Optical model point.** As in "Setting": leading eikonal with the retained multigraded set above.$
   541	  Part D's mechanical and inherited energy content retains its stated untruncated domain. Its results$
   542	  are conditional on P5's adopted momentum carriage and P6's adopted steady pressure-plus-linear-viscous$
   543	  form, with general live `η` and `ζ`. P6 is not derived from P1; power-law creep is outside this model$
   544	  point, and relaxation outside the adopted form and energy accounting remain OPEN. Results do not$
   545	  transfer to:$
   546	  - the strong field;$
   547	  - the throat mouth or interior;$
   548	  - a moving or rotating mass;$
   549	  - a drain that changes while light crosses;$
   550	  - polarization transport;$
   551	  - anisotropic or coupled branches.$
```

Exit code: `0`.

```bash
rg -n -i 'author|review|stop|launch|push|spawn|orchestrator|claude|build' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
8:**Authority and sources.** The base is v8 at `c2f1cf2b`, with its cited sources. The governing decisions
164:  (`directives/S9_wl_rebuild_directive.md:359`; `steps/S10_two_transverse_photons.md:816–819`):
478:obligation returned to the orchestrator under M1**. This is an unresolved input qualification,
```

Exit code: `0`.

```bash
git rev-parse 478b5e7b:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md HEAD:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md 79055918:research/pde_ledger_v3/directives/S9b_repair_decision_list.md HEAD:research/pde_ledger_v3/directives/S9b_repair_decision_list.md 72866fcf:research/pde_ledger_v3/steps/O2_steady_brane_balance.md HEAD:research/pde_ledger_v3/steps/O2_steady_brane_balance.md 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md HEAD:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md 217a92e9:research/pde_ledger_v3/directives/O2_input_contract.md HEAD:research/pde_ledger_v3/directives/O2_input_contract.md
```

```text
015d369c7226c1e351e5320c78e003f3563d993b
015d369c7226c1e351e5320c78e003f3563d993b
dd0cc090c5cd6486909617dbfacfd6e6baffdf2c
dd0cc090c5cd6486909617dbfacfd6e6baffdf2c
84745153f47d4fe2c07871af708be5aa4dd2658e
84745153f47d4fe2c07871af708be5aa4dd2658e
4ba41e1fc548c61d9f043fd37fbc6673df58c531
4ba41e1fc548c61d9f043fd37fbc6673df58c531
2b5f663e22fa25335e01c4ac66b60049366d74e3
2b5f663e22fa25335e01c4ac66b60049366d74e3
```

Exit code: `0`.

```bash
git diff --check -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```text
```

Exit code: `0`.
