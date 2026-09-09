# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v1 (DOCUMENT branch)

VERDICT: NOT-SOUND (5 must-fix). Verbatim (a 1-line reasoning-trace preamble precedes the report).

---

I'll form an independent view from the governing sources first, then open the spec. Starting with the decision list, c2 record, house template, N6 disposition, and profile-class consult.The consult file I opened is only the prompt. Next I’ll pull the consult responses, the S11c-d row, and the inherited power-counting/profile premises from S11c-a.I have a working picture of the profile class, order counting, and c2 import honesty. Next I’ll pull the inherited N6 pattern, confinement language, and conversion-form notes before opening the spec.Independent picture is in place: localized interface, Born-in-η with live σ_W, honest c2 import, two kill channels, and F'(0) only for the strong edge. Opening the spec now.Order bookkeeping and the strong-edge counterexample check out; several other claims do not. I’m writing the item-by-item findings and verdict next.# Independent physics review — S11c-d SHARED PHYSICS (v1)

**Artifact:** `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`  
**Role:** orchestrator-written physics spec (review-until-clear)  
**Verdict: NOT-SOUND** — five must-fix findings change what an engine computes or what the spec may claim.

I formed the object from the decision list, the c2 record/spec, and the N6 disposition **before** opening this spec. What S11c-d must be: the **linear mixing** of the transverse and thickness sectors on a **named localized-interface class**, using c2’s closed diagonal resolvent and closed off-diagonal kernel, with honest c2 status (per-engine values, covariance-only cross-engine N6, operand DEBT live), two distinct photon-kill channels as **computed** objects, `(ε,η)` bookkeeping that keeps `σ_W` and `k L_W` live, and the strong edge named as a **new** construction, not a reduction.

---

## Independent view (before the spec)

- **Class.** N5 forbids `ω(k)` for generic `W₀(x)`. The slit-edge endgame is a localized gradient, not a lattice or an adiabatic ramp, so the matching object is a Born/scattering kernel on a **localized interface**: asymptotically constant, `∫ W' = ΔW ≠ 0`. A bump (`∫ W' = 0`) is a different low-`q` vertex (the DC piece of `Ŵ'(q)` vanishes).
- **Grades.** S11c-a already splits contrast `η` from first-jet sharpness `σ_W = η W̄₀/L_W`. `k L_W` is a third, independent kinematic grade (vary `k` at fixed `η,σ_W`). “Weak-gradient” collapses sharpness into contrast and invites `L_W→∞`, which is WKB, a different N5 class.
- **Orders.** In a **linear** theory the converted amplitude is one insertion of an `O(η)` integrated vertex, hence `O(εη)`; fluxes are quadratic, so `J_conv = O(ε²η²)`, `J_in = O(ε²)`, and `C = J_conv/J_in = O(η²)` with `ε²` cancelling. Absolute and fractional labels are different objects; neither is the N10 nonlinear program.
- **Strong edge.** N7: a slit edge is order-unity, not Born. A nonzero Born coefficient is `F'(0)` of an **amplitude**, not a lower bound on `C` at finite contrast. `C = sin²(ηG)` is a valid coupler counterexample.
- **c2 honesty.** Operator/kernel **values** are per-engine SOUND; cross-engine content is N6 covariance only (`(0)−(0)`, not operand AGREE); carrier 40 / source 76 / Φ 18 are UNADJUDICATED and **material**; leftover SHAPE uninspected; F/G withdrawn; preserve **both** `R_N6 = 18/288` and `R_cov` no-nonzero.
- **N6.** Independent Eulerian vs material-coordinate constructions of **this sub-step’s object**, then one-sided tilt/advection corruption; not `∇W→0`; no `A−A` on structurally absent `RHO4_CONSTANT` advection; not an anchoring swap.
- **N13.** Confinement = survival of the **transverse** channel. Continuum conversion and bound capture are distinct emitted objects. Existence of a bound pole is a computation, not a theorem to paste onto a step.

---

## Must-fix findings

### F1. `F'(0)` is the wrong weak object — slope of the conversion **fraction** vanishes

**Spec §3d (quote):**

> S11c-d honestly delivers `F′(0)` (the `O(η)` slope of the conversion); the lab bounds `F(1)`.

and the emitted object

> `⇒ S11CD_WEAK_MATCHING_COEFFICIENT (F′(0))`

**Why it is wrong.** §3c of the **same** spec (and N12/N5) sets the conversion **fraction** `C = J_conv/J_in = O(η²)`. For any such `C`,

\[
C(η) = c_2 η^2 + O(η^4) \implies C'(0) = 0.
\]

If an engine takes `F` to be “the conversion” in the N7 sense (the dimensionless quantity the lab bounds), it will emit `F'(0) = 0` as the weak matching coefficient. That is a structurally zero object, not the Born content.

The counterexample in the same subsection already shows the right split:

\[
A_H = -i\varepsilon\sin(η G),\qquad C = \sin^2(η G) = η^2 G^2 + O(η^4).
\]

- Weak object = leading amplitude coefficient \(A_H/(εη)\big|_{η\to 0} = -i G\), or equivalently \(C/η^2\big|_{η\to 0} = G^2\).
- Lab object = \(C\) at order-unity contrast, i.e. \(|A_H/ε|^2\) at \(η\sim O(1)\), **not** “`F(1)`” of an unnamed `F`.

N7: “the lab bounds a **dimensionless conversion** (photons lost at an edge)”. That is `C`, not an amplitude. Identifying “lab bounds `F(1)`” with “`F'(0)` is the `O(η)` slope of the conversion” mixes amplitude-slope language with fraction-at-unity language and will send the two engines to different objects (zero slope of `C` vs Born coefficient of `A_H`).

The rest of §3d (order-unity out of scope; new construction, not a reduction; no lower bound from a nonzero Born coefficient) is physically right. Only the named weak/lab pair is wrong.

**Minimal fix.** Define the weak object as the leading coefficient of the converted **amplitude**, e.g. \(A_H/(εη)\) as \(η\to 0\) (or of `C/η²`). State that the lab bounds the **fraction** `C` at order-unity contrast, not `F(1)` of that amplitude. Drop the phrase “`O(η)` slope of the conversion.”

---

### F2. Bound-mode channel asserted by the 1D-well theorem, which does not apply to the named **interface**

**Spec §3b (quote):**

> (ii) bound-mode capture — transverse → a **bound** thickness/breathing pole … (a weak attractive well binds a mode in 1D, so this is a real photon-kill channel at **small** `η`, and is **not** a Bloch band).

**Spec §1c (the named class):**

> a **localized interface** … `f(−∞) ≠ f(+∞)` ⇒ `∫ W₀′(x) dx = W₊ − W₋ ≠ 0` — ⛔ **not** a defect bump (`∫ W₀′ dx = 0`).

**Why it is wrong.** The 1D weak-binding theorem is for a **well** that returns to the same asymptotic value (`V(±∞)` equal, `∫ V < 0`): any attractive 1D well binds, at arbitrarily small depth. A monotone **interface** (step) interpolates two different asymptotes. For a second-order wave/Schrödinger operator, `V(x)` then lies above `min(V_-,V_+)`, so there is **no** state below both continua. A thickness step is that object, not a well.

The bump the spec correctly excludes (`∫ W' = 0`) is the profile to which the cited theorem **does** apply. Using the well theorem to guarantee a pole **on the interface class** is the wrong existence claim for the named object.

Two further ways this goes wrong:

1. **Manufacture.** “This is a real photon-kill channel at small `η`” can send an engine to produce a pole even when the thickness resolvent of `S11CC2_CLOSED_SLAB_OPERATOR` has none. N13 requires the two **channels as emitted objects**; it does **not** require that a bound pole exist for every named class. N10’s confinement question is whether survival is unconditional — that is a **computed** (possibly empty) pole set, which §3b otherwise asks for correctly.
2. **Internal collision with §2.** §2 tells the engine to “Exclude unresolved threshold/resonance enhancements”. A bound-mode capture **is** a resolvent pole. That exclusion, read onto §3a/§3b, can delete the very object §3b asserts is real.

A Jackiw–Rebbi kink mode would be an interface mechanism, but that is a first-order/Dirac mass kink, not “a weak attractive well,” and is not supplied as the thickness-sector structure (S11b’s breathing sector is a second-order wave operator). If that is the intended mechanism, the spec must name **that** operator structure; it currently names the well theorem.

**Minimal fix.** Emit poles/residues of the thickness-diagonal resolvent as a **computed** object; permit computed absence. Do not assert 1D-well existence on an interface. If a bound channel is expected only for a well/bump, say so and keep it off the interface class (or name a distinct well profile). Restrict the §2 “regular domain” clause so it cannot delete §3b(ii).

---

### F3. N6 re-enters at **face linearization**, which bypasses c2’s close-then-extract

**Spec §5a (quote):**

> route 1 (level-set / graph): derive the off-diagonal mixing by direct level-set / graph linearization of the localized-interface faces ;  
> route 2 (flattened material): derive it AGAIN after flattening the faces into material coordinates…

**Spec §5 intro:** “Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result.”

**Spec §2:** mixing = one insertion of `S11CC2_CLOSED_COUPLING_KERNEL` between resolvents of `S11CC2_CLOSED_SLAB_OPERATOR`.

**Why it is wrong.** c2 exists because extract and close do not commute (`S11c_c2_SHARED_PHYSICS.md` §2):

\[
\texttt{extract}(\texttt{close}(\mathrm{SLAB})) \;\neq\; \texttt{close}(\texttt{extract}(\mathrm{SLAB})).
\]

The object d consumes is `extract(close(·))`. Linearizing the faces and extracting is c2’s **open** kernel path (`close(extract)`), which c2 forbids as a construction route and uses only as a §5a ablation. An engine that follows §5a literally rebuilds mixing from geometry and **drops the self-energy fold**. That is a different physical object: photon conversion without the bulk DtN threaded through the diagonal and the vertex.

Parent pattern: S11c-a/b apply level-set vs material routes to **their** new operator/kernel; c2 applies Eulerian vs material routes to **the increment**, from the imported closed face response, not by re-deriving S11c-b from scratch. d’s new object is the **mixing response**. N6 must construct **that** two ways from the two representations of the **closed** operator/kernel.

N6 in the decision list (`:94–104`) is written for the kernel, at the family level. Copying its “linearize the faces” sentence into a **consumer** of the closed kernel inverts the fold.

**Also in §5a:** “require a **nonzero** residual” for the one-sided corruption. That contradicts the same subsection (“computed value is the finding, ⛔ no target value is supplied, … ⛔ never a builder exit condition”) and S11c-b §5a: “the resulting residual is the computed outcome, not an asserted value (an honest zero from symmetric cancellation is a valid finding)”. A “must be nonzero” exit is a builder target (N7/M2): junk will satisfy it; an honest zero will be rejected.

**Minimal fix.** Route 1 = mixing from the Eulerian closed operator/kernel; route 2 = mixing from the native material-coordinate closed construction, already in the common Eulerian face basis (no separate `T` on the differenced mixing, as c2). Re-enter at those imported operands, not at a new face linearization. Print corruption residuals; adjudicate nonzero/absence on our side; keep the `RHO4_CONSTANT` computed-absence / no `A−A` rule.

---

### F4. Flux `J` is unnamed, so `C` can hide bound capture (the N13 defect)

**Spec §3c (quote):**

> dimensionless conversion FRACTION `C = J_conv / J_in = O(ε⁰ η²) = O(η²)` (the incident `ε²` CANCELS).

**N13 (quote):**

> It means **survival of the transverse polarization channel**. Conversion into a **bound** breathing/thickness mode kills the photon exactly as bulk radiation does … ⛔ “Energy stays in the slab” is **not** confinement.

**Why it is wrong.** No bilinear form is named for `J`. Candidate fluxes already in the chain disagree on N13:

| Flux | Bound capture | Continuum / bulk escape |
|---|---|---|
| Total mechanical energy / slab stored + kinetic | **zero** (energy stays in the slab) | possibly zero if energy remains in-brane |
| Far-field bulk Poynting (c1 ENERGY) | **zero** (no radiation) | nonzero |
| **Transverse-channel** energy flux | nonzero (photon lost) | nonzero |

If either engine uses total energy or bulk Poynting, `S11CD_BOUND_MODE_CAPTURE` does not appear in `C`, and “energy stays in the slab” is reintroduced as the conversion observable. That is exactly the reading N13 forbids. N7 also notes that a bare coupling “still carries an undetermined normalization until ‘what couples to what’ is fixed” (`S11b_SHARED_PHYSICS.md:824–825`). Flux-normalization is that fix; it only works if `J` is the **transverse** flux.

This is a missing premise, not a style gap: two engines can agree on amplitudes and disagree on `C`, and one of those `C`s is physically the N13-wrong object.

**Minimal fix.** Name `J` as the inherited **transverse-polarization-channel** energy flux (incident = incoming transverse; converted = that channel’s lost flux, equal to continuum thickness/bulk plus bound-capture). Point at the S11b quadratic energy / c2 traction–slab pairing as the bilinear form; do not leave “flux” generic.

---

### F5. Class and vertex are stated on `W`/`∇w₁`; the load-bearing mixer is `∇μ_R` (`m₁` independent)

**Spec §1c (quote):**

> The named class is a LOCALIZED INTERFACE. The background thickness profile is  
> `W₀(x) = W̄₀·[ 1 + η·f(x/L_W) ]` …  
> The mixing vertex depends on `∇w₁`.

**Spec §2 / §0 / §5d:** mixing is “driven by `∇μ_R ≠ 0` (`∝ k·a`)”, “supported only where `∇μ_R ≠ 0`”.

**S11c-a §2a (inherited by pointer):** `w₁` and `m₁` are **independent** `O(1)` profiles,

\[
\partial_i W_\mathrm{bg} = σ_W\,\partial_i w_1,\qquad
\partial_i μ_{R,\mathrm{bg}} = (μ̄_R/W̄₀)\,σ_W\,\partial_i m_1.
\]

**S11c-b §5b (quote):** “`w₁` and `m₁` are independent with separate derivative maps … ablating both together cannot separate thickness-slope from modulus-gradient coupling (the very channels `N6` distinguishes).”

**V3_STEP_PLAN:1179:** “the families do **not** mix in a *homogeneous* brane (coupling `∝ k·a`); `∇μ_R ≠ 0` mixes them.”

**Why it is wrong.**

1. **Internal contradiction.** §1c says the vertex depends on `∇w₁`; §2 says it is supported only where `∇μ_R ≠ 0`. Those are S11c-b’s two distinct first-jet channels. Reducing the closed kernel to `∇w₁` drops the modulus-gradient mixer that the program names as the mixing source.
2. **Incomplete class.** The interface restriction (`f(−∞)≠f(+∞)`) is imposed only on `W`/`f`. Under the inherited ansatz `m₁` may be a bump, a constant, or a different interface. Then §5c’s edge-vs-bump discriminant on `∫ W'≠0` does **not** control the low-`q` content of the `∇μ_R` vertex (`∫ m₁'`). An engine that only shapes `w₁` can emit an “interface” mixing that is constitutively a bump.
3. **N14 name.** The varying field is written `W₀(x)` while `W_0` is a reserved **constant** key (`N14`; S11c-a’s varying field is `W_bg`). That is the F9 false-equal / freeze hazard N14 exists to stop. §7 forbids reusing `W_0`; the display equation in §1c still writes the class as `W₀(x)`.

**Minimal fix.** Identify the class as a restriction of the **inherited** `{W_bg, μ_{R,bg}, ρ…}` ansatz (`w₁`, `m₁`, density maps). Apply the localized-interface condition to every background that enters the kernel, or state an explicit constitutive tie `m₁ = m₁[w₁]`. State that the vertex is the **full** `S11CC2_CLOSED_COUPLING_KERNEL` (tilt `∇W` **and** `∇μ_R` **and** N4 advection), not `∇w₁` alone. Write the varying thickness as `W_bg`, never `W₀(x)`/`W_0`.

---

## Order bookkeeping (item 2) — derivation, holds

Linear coupled-mode, one insertion of the off-diagonal vertex (the §2 object):

- Incident transverse field \(ψ_T = O(ε)\).
- Local kernel density is first-jet, \(K \sim σ_W\). Integrated over an interface,
  \[
  \int K\,dx \sim σ_W L_W \int f'(ξ)\,dξ = η\,\overline{W}_0\cdot Δf = O(η),
  \]
  independent of `L_W`. The Born amplitude is that integrated vertex times the incident field:
  \[
  ψ_H = G_H K ψ_T = O(εη).
  \]
  (`σ_W` remains live in the form factor \(\hat f'(k L_W)\); it does not change the `η` power.)
- Energy flux is quadratic in amplitude:
  \[
  J_\mathrm{in} \sim |ψ_T|^2 = O(ε^2),\qquad
  J_\mathrm{conv} \sim |ψ_H|^2 = O(ε^2η^2).
  \]
- Fraction:
  \[
  C = \frac{J_\mathrm{conv}}{J_\mathrm{in}} = O(η^2).
  \]
  Doubling `ε` doubles both fluxes; `ε²` cancels **because the equations are linear**. An `O(ε²)` **vertex** would be N10 (intensity/harmonic). An `O(ε²η²)` **observable** of a linear mixing is not that program.

Emitting **both** `O(ε²η²)` (absolute) and `O(η²)` (fractional) is required: N12 names the first, N5/N7 the second. That part of §3c is right. It does not salvage F1 or F4.

---

## Strong-edge counterexample (item 3) — derivation, holds except F1

Lossless two-mode coupler with integrated coupling \(η G\):

\[
\partial_z A_T = -i η g\, A_H,\quad
\partial_z A_H = -i η g\, A_T,\quad
G=\int g\,dz.
\]

With \(A_T(0)=ε\), \(A_H(0)=0\):

\[
A_H = -iε\sin(η G),\qquad C=\sin^2(η G).
\]

At small `η`, \(C=η^2 G^2+O(η^4)\) (Born). At \(η G = nπ\), \(C=0\) exactly, while Born is \((nπ)^2 ≠ 0\). So a nonzero Born coefficient gives **no** positive lower bound on strong-edge conversion. Evaluating the Born coefficient at `η=1` is `F'(0)` vs `F(1)` in the invalid sense N7 names. Order-unity matching is a **new** construction (piecewise-uniform sides are `η`-exact; it is not a first-jet kernel). All of that in §0/§3d is honest.

---

## What holds (no finding after the filter)

**Profile class / grades (item 1), aside from F5.** Naming a localized interface with `∫W' ≠ 0` vs bump `= 0` is the right non-global N5 object for an edge. The three grades `{η, σ_W, k L_W}` are independent at fixed `W̄₀` (vary `η`; vary `L_W` against `η` to move `σ_W`; vary `k` to move `k L_W`). “Weak-gradient” is the wrong name: the Born form factor is

\[
M(q)=η\,\overline{W}_0\,\widehat{f'}(q L_W).
\]

Sending `L_W→∞` at fixed `η` (`σ_W→0`) concentrates \(\widehat{f'}\) at `q=0` and kills finite-`q` conversion (WKB, the class N5 says does not match an edge). Sending `η→0` at fixed `L_W` is weak **contrast** and keeps the form factor. First-shape-order operators at `η=O(1)` are an incomplete operator, not a non-perturbative theory; §1d is right to forbid that.

**c2 import (item 4).** §1b marks per-engine SOUND values; N6 as covariance / Reading B; matched zeros as `(0)−(0)` not operand AGREE; carrier 40 / source 76 / Φ 18 UNADJUDICATED and material; leftover SHAPE uninspected; “not known to be just thickness”; F/G withdrawn; **both** `R_N6=18/288` and `R_cov` no-nonzero preserved; c1 UNDECIDED carry-ins and the two S11c-b slot-sign conventions named. I did not find a quiet upgrade of any of those to “closed.”

**N6 probes (item 5), aside from F3.** Uniform limit is a secondary regression (§5b). Tilt (N3) and advection (N4) are the two same-order channels. `∇W₀→0`/`η→0` is rejected as a corruption. `RHO4_CONSTANT` is computed absence, not `A−A`. Anchoring swap is excluded. `Δρ` is not used to bridge `LAB_HELD ↔ MATERIAL_ADVECTED`.

**Two channels as objects (item 6), aside from F2.** Distinct tags; confinement emitted, not asserted; bound pole correctly distinguished from a Bloch band.

**Recipe / freeze (item 7), aside from F1.** One kernel insertion **is** the leading linear mixing (object, not a path among many); Born vs DWBA coinciding at leading `O(η²)` in `C` is correct (`G = G_0+O(η)` ⇒ amplitude difference `O(η²)` ⇒ `C` difference `O(η³)`). `O(1)`/grating reductio and the numeric lab bound are withheld. `σ_W` and `k L_W` are kept live.

**House / N11 (item 8), aside from F4–F5.** Chain is generate-over-frozen-base with c1+c2 deltas and `ledger_fold.py` in `BUILD_INPUT_DIGESTS`; T7 comparator surfaces the DEBT; denylist stays cut; blind WL re-derives; supplied vs computed is split; `N11a` rest-frame limit is present (`v_bulk_normal_0`, not `v_0`). N2 allows the spec to refine the a–e boundary, so placing the FORM in d is not by itself a defect.

---

## Nits (do not change the verdict by themselves)

1. **§0:** “a periodic profile is a later evaluation of the **same** kernel, not an S11c-d class.” Bloch bands are not Born scattering evaluated on a periodic `f` (N5: different object). Harmless because periodic is out of scope; drop “same kernel” so an engine does not emit a periodic evaluation of the interface S-matrix and call it Bloch.
2. **§2 formula** carries an explicit minus, `− G_H⁺ · K_HT · ψ_T^{in,D}`, against §3’s “states no … sign.” Overall sign is `i0`/Lippmann–Schwinger convention; `C` is insensitive, but the comparator’s load-bearing residual is the mixing **amplitude**. Omit the minus; let each engine’s resolvent convention fix it.
3. **N2 vs d:** FORM and confinement interpretation are S11c-e in the decision-list table. N2 allows spec-stage boundary refinement; just don’t let d’s FORM become e’s withheld `O(1)` target.
4. **“Frozen approximation”** (§1d) is ambiguous next to N11b’s frozen-wall-width freeze. Say “Born / static-profile.”
5. **House completeness:** a/b/c2 emit named homogeneity-corruption objects; d only points at dimensional able-to-fail. Optional, not load-bearing for mixing.
6. **N11** “large `k c_{s0}/|ω|` is necessary, ⛔ not sufficient” is not restated; “boundary-layer/subsonic” is close enough.
7. **`kL` vs `kL_W`** in the §1d lead-in vs the grade list.

---

## Verdict

**NOT-SOUND.**

Must-fix: **F1** (wrong weak object `F'(0)` / lab `F(1)`), **F2** (1D-well existence claimed on an interface; pole may be manufactured or deleted), **F3** (N6 re-enters at face linearization and asserts a nonzero target, bypassing the closed kernel), **F4** (`J` unnamed; total-energy flux hides bound capture), **F5** (class/vertex stated on `W`/`∇w₁` while `∇μ_R`/`m₁` is the independent mixer; `W₀(x)` N14 hazard).

Nothing in items 2 (orders), 3 (strong-edge physics aside from F1), or 4 (c2 honesty) failed the filter. F1–F5 each change what an engine computes or what the spec may claim; they are not nits.
