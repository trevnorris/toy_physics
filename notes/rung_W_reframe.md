# Rung W reframe — material-state brane (χ_B): adversarial pre-check

**Role:** independent adversarial reviewer. **Goal:** try to break the χ_B brane route BEFORE building the directive. Three prior routes died; a clean pass here is suspicious.

**Date:** 2026-07-03. **Sources read:** `conceptual_foundation.md` §0/§2(v5)/§3/§8; `pathA_35_shear_surface_brane_gates.md` §2.3–2.6/§8; `pathA_35_gateL_light.md`; `pathA_25_gateRC_cubic.md`; `pathA_24_T1_wall.md`; `pathA_35_G0_freeze.md`; `decisions/15` §11/§16/§17; `brane_bulk_handoff.md`; `medium_requirements_and_prior_art.md`.

---

## Rung W verdict: `GENUINE_COIN_FLIP` (leaning `LIKELY_PASS` on the narrow question, but the pass is cheap and the route still dies at Gate L)

**The narrow question:** is a double-well χ_B wall stable (escapes T1) AND light-permitting (escapes the density no-go)? **My answer: yes, on a load-bearing assumption — but it doesn't matter, because the C5 obstruction survives into Gate L and the θ rescue is structurally broken.** The χ_B wall is a better wall than any prior route, but it is not a better *light substrate*. The route dies one rung later than W, not at W.

**The single load-bearing assumption:** χ_B is an **independent scalar order parameter** with a postulated double-well `f_B(n,χ_B)` that has disconnected minima at fixed density n. If χ_B is derived from P^i (e.g. χ_B = |P_∥|²), the wall reverts to T1 (P's vacuum manifold is connected). If χ_B is truly independent, it is second-medium drift (conceded under the analog license, §0.6, but costly — see the ledger below).

---

## 1. Pressure-test: Escape-T1

**Claim:** a scalar χ_B double-well gives disconnected minima → π₀=Z₂ → topologically stable φ⁴-kink-like wall.

**Verdict: GENUINE — this is the one claim that holds cleanly, and it is confirmed by T1's own negative control.**

The T1 report (`pathA_24_T1_wall.md` lines 114–116) ran a φ⁴ kink as its able-to-fail negative control and it returned **stable** (no negative mode, `π₀=Z₂`). T1 died because P^i lives on a **connected** S³ vacuum (`π₀=0`, line 104) → the wall spreads to infinite width with zero barrier (`σ_L → 0` as `L_half → ∞`, lines 74–90). The diagnosis (line 141) explicitly said: "a stable antipodal wall would need disconnected vacua, for example an easy-axis or component **double-well** that makes `π₀=Z₂`." The χ_B route is exactly that structure. So the topological escape is real.

**But does coupling χ_B to the rest of the medium re-introduce a connected direction?**

Three sub-pressures, in order of severity:

**(a) The conservation-law pressure (the strongest).** The medium has fixed total `∫n d⁴X` (`brane_bulk_handoff.md` §6, line 196). The free energy `f_B(n,χ_B)` (`brane_bulk_handoff.md` §10, line 433) couples χ_B to n. If the two minima of `f_B` at fixed n prefer **different equilibrium densities** (n₁ for χ_B=1, n₀≠n₁ for χ_B=0), then at coexistence (common tangent) the wall becomes a **density interface** — a liquid-vapor boundary. This is exactly the "(a) two densities" alternative already considered and rejected in `conceptual_foundation.md` §2 (lines 211–213: "dense vs dilute are not mirror images → asymmetric charge, and a pure density wall is a fluid membrane with no in-plane shear → no light").

*Can this be avoided?* Yes, **if** `f_B(n,χ_B)` has disconnected minima in χ_B **at the same n** — a Landau double-well `½α(T−T_c)χ_B² + ¼βχ_B⁴` with α<0, where both minima sit at the same density. This is possible (it is the Ising/φ⁴ case) but it is **postulated** — the GNLS EOS `U(ρ)=(K/4)ρ⁵` is single-well and provides no such structure. The v5 block concedes this (`conceptual_foundation.md` line 129: "the double-well `f_B` is still postulated — the GNLS potential `U(ρ)∝ρ⁵` is single-well, so 'two coexisting phases' is the same degenerate-vacua obstacle the brane has always had"). It is real drift, honestly labeled.

**(b) The "is χ_B really independent of P?" pressure.** If χ_B is a *description* of P's state (e.g. χ_B=1 when P is in-plane ordered, χ_B=0 when P is disordered/along ŵ), then χ_B is not an independent field — it is a function of P. The wall is then a wall in P, and P's vacuum manifold is connected (S³, T1). The wall unwinds. This is the fatal version: **χ_B must be genuinely independent of P to escape T1**, but "genuinely independent" is exactly what makes it look like a second field.

The analog license (§0.6) permits postulating χ_B as an independent order parameter of the one medium (like the nematic director is an order parameter of a liquid crystal, not a second substance). But the cost must be counted (see ledger in §4 below).

**(c) The "does the χ_B–P coupling open a flat direction?" pressure.** The anisotropy needed to force P in-plane in the ordered phase (see §2 below) breaks P's O(4) symmetry to O(3). The vacuum manifold of P in the ordered phase is S² — also connected (`π₀(S²)=0`). But the wall is in χ_B (π₀=Z₂), not in P. P's manifold being connected doesn't matter because the wall doesn't live in P's field space. This is the structural difference from T1, and it is genuine.

**Escape-T1 bottom line:** The scalar double-well does give a stable wall. The conservation-law pressure (a) is real but avoidable with a fixed-n double-well (postulated). The independence pressure (b) is the fundamental cost — χ_B must be independent of P. This is the load-bearing assumption.

---

## 2. Pressure-test: Escape-density-no-go (the expected kill)

**Claim:** χ_B is a single domain wall in an abstract order field, not a periodic density stack → no layer normal to pin the arrows against → in-plane shear lives in a separate field that χ_B merely gates on.

**Verdict: GENUINE — the χ_B wall does NOT re-create the density no-go. But the escape costs an extra postulated anisotropy, and the reason it escapes is structurally illuminating about why the density route couldn't.**

### Why the density no-go fired (the mechanism to compare against)

The density no-go (`pathA_25_gateRC_cubic.md` lines 54–58) was precise: the Cpin coupling `χ_Cpin(P₀·k)²` **simultaneously** (i) opened the codim-1 density-smectic window and (ii) pinned P along the layer normal. At the static minimizer with no rotation drive (`Ω_u=0`), `P_parallel=0` — the in-plane arrows were **pinned to zero**. **One coupling did two jobs**, and the second job killed light.

### Why χ_B escapes (the structural difference)

For χ_B, the two jobs are done by **two independent couplings**:

1. **Wall creation:** the χ_B double-well `f_B(χ_B)` + gradient stiffness `½κ_B|∇χ_B|²`. This creates the codim-1 wall. It does NOT involve P.
2. **P orientation:** an anisotropy `α χ_B (P·ŵ)²` (or equivalent), which forces P in-plane when χ_B=1. This is a SEPARATE term.

Because the wall-creation coupling doesn't involve P, there is no mechanism for the wall's geometry to pin P. The gating is multiplicative and local: `E_shear = ½ χ_B μ_R (∇×u)²` — the shear modulus is just turned on/off by χ_B, with no preferred in-plane direction.

### The concrete static calculation (the proof that P_∥ is free)

Write the minimal static energy at the wall (P^i ∈ ℝ⁴, T0 O(4)-isotropic, plus the χ_B-gated anisotropy):

```
E/A₃ = ∫ dw [ ½κ_B(∂_w χ_B)² + f_B(χ_B)
              + ½κ_P(∂_w P^i)² + ¼r(|P|²−1)²
              + α χ_B (P^w)² ]
```

Static equations for P^i at the kink background χ_B(w):

- **w-component:** `κ_P ∂_w² P^w = r(|P|²−1) P^w + 2α χ_B(w) P^w`
  → at the wall centre (χ_B=1), effective mass of P^w is `2α > 0` → **P^w is gapped to zero**. Good.

- **in-plane components:** `κ_P ∂_w² P^∥ = r(|P|²−1) P^∥`
  → at the wall centre, with P^w=0, `|P|² = |P^∥|²`, and the vacuum is `|P^∥|=1` → **P^∥ is a free S² Goldstone** (2 transverse in-plane rotations, zero frequency at k=0).

**The χ_B gradient `∂_w χ_B` does NOT appear in the P^∥ equation.** The anisotropy `α χ_B (P^w)²` is local in χ_B — it doesn't couple gradients of χ_B to P. So there is **no gradient-direction pinning** of the in-plane arrows.

Compare to the density no-go: there, the coupling `χ_Cpin(P₀·k)²` explicitly coupled P to the layer normal k (the modulation direction), and the static minimizer gave `P_∥=0`. Here, the anisotropy penalizes P along ŵ (the wall normal), which forces P **into** the plane — the *opposite* of pinning. The in-plane components are Goldstone modes.

### The cost of the escape: the anisotropy is postulated

The anisotropy `α χ_B (P·ŵ)²` is **not derived from T0** (which is O(4)-isotropic, no easy axis — `pathA_24_T1_wall.md` line 25: the freeze asserts no `P·w` or `w_hat` term). It must be **postulated** as a new ingredient. The sign (penalize P^w, not P^∥) must be chosen to match the desired "shear-supporting" phase.

This is allowed under the analog license (§0.6: postulate the structure freely). But:
- It is an **additional ŵ-dependent coupling**, on top of the parity-repaired P–u coupling `-λ_Pu ŵ·(P_∥ × (∇×u))` already in pathA_35. Both break O(4) using the wall normal ŵ.
- The sign choice (penalize out-of-plane P) is **not innocent**: it is the choice that makes the ordered phase "shear-supporting" — i.e. it postulates the *outcome* (in-plane shear) via an *ingredient* (the anisotropy). The pathA_35 directive (§1, line 49) warns: "postulating an ingredient is allowed; postulating an outcome is not." The anisotropy is an ingredient, but it is an ingredient *chosen to produce the desired outcome*, with no independent mechanical provenance. It should be labeled `POSTULATED_ANISOTROPY` in the G0 ledger.

**Escape-density-no-go bottom line:** The χ_B wall genuinely does NOT re-create the density no-go. The structural reason is clean: the wall-creation and P-orientation are done by independent couplings, so the wall's geometry doesn't pin P. The in-plane P components are free Goldstone modes at the wall centre. The cost is a postulated anisotropy (one new coupling, ŵ-dependent).

---

## 3. The fresh no-go I see (NOT where you expected — but it bites at Gate L, not W)

**I do NOT find a fresh no-go between "codim-1 localization" and "in-plane rotational-elastic shear" for a scalar order parameter.** That specific pair is compatible. The χ_B wall is stable, and the in-plane P Goldstones are free. Rung W, narrowly read, **passes**.

**But I find a deeper structural trap that makes the pass Pyrrhic — the χ_B wall changes the confinement mechanism, not the light-sector dynamics.** The pathA_35 Gate L four-way no-go (`FAIL_COUPLE_STRESS_NOGO`) fired on the *light-sector constitutive package* (MacCullagh + P–u + couple-stress + C5), not on the confinement profile. The χ_B wall replaces the frozen `g(w)` profile with a dynamical χ_B kink, but the light-sector questions are **identical**:

| Gate L sub-hurdle | pathA_35 (frozen g(w)) | χ_B wall | Changed? |
|---|---|---|---|
| L(a-i) traction | PASS, ARROWS_SUPPLY_TRACTION | same P–u operator, same ŵ | **No** |
| L(a-ii) hidden modes | FAIL (live P) / PASS (slaved) | same P mass/gap decision | **No** |
| L(a-iii) C5 | FAIL (no φ) | **same — unless θ saves it (see §5)** | **Only if θ works** |
| L(b) closure | FAIL (live P) / PASS (slaved) | same couple-stress sector | **No** |
| L(c) leak | PASS | same free-slip traction | **No** |
| L(d) u_w gap | PASS | same gap term | **No** |

The χ_B wall is a better *wall* (dynamical, self-consistent, escapes T1) but it is the **same light substrate**. The four-way no-go will fire identically unless θ (the C5-φ sub-hypothesis) solves the C5 leg — and θ is structurally broken (§5 below).

**This is the most valuable negative finding: the χ_B route does not need a fresh no-go at rung W because it inherits the pathA_35 no-go at Gate L.** The wall-existence *and* light-permitting question (rung W) passes, but the light-*consistency* question (Gate L) fails for the same reasons it failed before. The χ_B field is a red herring for the light sector — it solves the wrong problem.

### The one genuinely new Gate-L concern the χ_B wall introduces

The χ_B wall has a **translation Goldstone** (the kink can be displaced in w). This is a gapless scalar on the brane — potentially a fifth force (the same hazard as u_w, L(d)).

- If the wall is "carried by the medium" (δw = u_w), the u_w gap term also gaps the wall Goldstone. This is the same "carried by the medium" condition needed for MacCullagh traction (`conceptual_foundation.md` §3, lines 300–305). It is not an additional constraint — but it must be verified.
- If the wall is NOT carried by the medium (it can move independently), the Goldstone is gapless → `FAIL_BENDING_MASSLESS_FIFTH_FORCE` at L(d), for a new mode.

This is not a rung-W no-go (it's a Gate-L question), but it is a **fresh** concern that pathA_35 didn't face (the frozen g(w) profile had no translation Goldstone).

---

## 4. The χ_B parameter ledger (second-medium drift count)

The v5 block (`conceptual_foundation.md` lines 129–134) honestly flags the drift. Here is the concrete count for rung W:

| Item | Count | Tag |
|---|---|---|
| χ_B field (scalar DOF) | 1 new field | postulated-ingredient |
| Double-well `f_B(χ_B)` (≥2 constants: α, β) | 2 | postulated-ingredient |
| Gradient stiffness `κ_B` | 1 | postulated-ingredient |
| χ_B–shear gating (`χ_B · μ_R`) | 0 (structural, uses existing μ_R) | structural |
| χ_B–P anisotropy `α χ_B (P·ŵ)²` | 1 new coupling + 1 new constant | postulated-ingredient |
| (if complex) Phase sector `½χ̃|χ_B|²(∂_t θ)²` + `½κ_θ|∇θ|²` | 2 | postulated-ingredient |
| (if θ = C5 φ) Maxwell-type electric term `½ε(∂_t u + ∇θ)²` | 1 new constant | postulated-ingredient |

**Minimum (real χ_B, no θ):** 1 field + 4 constants + 1 structural constraint = **6 new inputs**.
pathA_35 G0 cost 11 and was labeled `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`. The χ_B package costs **at least 6** (real) or **9** (complex with θ + Maxwell term). This is significant drift — less than pathA_35 (because the confinement profile is no longer a separate postulate — χ_B IS the confinement), but still well above the "zero new DOF" baseline the program wants.

---

## 5. C5-φ sub-hypothesis verdict: `FATAL_FLAW` (structural — do NOT carry into rung φ)

**Claim:** if χ_B is complex (amplitude + phase θ), the phase θ is a new on-brane scalar with genuine mechanical provenance, not u_w — a candidate for the C5 φ.

**Verdict: FATAL FLAW. The phase of an order parameter cannot serve as the MacCullagh scalar potential. The obstruction is structural, not numerical.**

### The structural argument

The C5 obstruction (`decisions/15` §11, lines 180–188; `pathA_35_shear_surface_brane_gates.md` §2.4) is precise: MacCullagh's curl-only *potential* energy `½μ_R(∇×u)²` is invariant under `u→u+∇χ`, but the *kinetic* term `½ρ_br(∂_t u)²` is NOT — the EOM gives `∂_t²(∇·u)=0`, a constrained physical zero mode. Maxwell escapes because the scalar potential transforms as `φ→φ−∂_tχ`, making the **full action** invariant through the gauge-invariant combination `E = −∂_t A − ∇φ`, with electric energy `½ε|E|²`.

For θ to serve as φ, the action must contain the Maxwell-type structure:

```
S = ∫ [ ½ε(∂_t u + ∇θ)²  −  ½μ_R(∇×u)² ]
```

so that under `u→u+∇χ, θ→θ−∂_tχ` the action is invariant. This requires:

**(i) A cross-coupling `ε (∂_t u)·(∇θ)` in the kinetic energy.** This is NOT present in the MacCullagh package (which has `½ρ_br(∂_t u)²` — pure inertia, no θ). It is NOT present in the Landau theory of a complex order parameter either: the phase θ has its own kinetic energy `½χ̃|χ_B|²(∂_t θ)²` (the Josephson term), which depends on `∂_t θ`, not `∇θ` coupled to `∂_t u`.

**(ii) The electric-type energy `½ε(∂_t u + ∇θ)²` must be postulated.** This is a **new constitutive term** — an "electric MacCullagh" stiffness — that is not part of the original MacCullagh package and has no independent mechanical provenance. It changes the nature of the light sector from **mechanical** (wave speed = `√(μ_R/ρ_br)`, stiffness/inertia) to **gauge-field** (wave speed = `1/√(μ_R ε)`, ratio of stiffnesses). The brane inertia `ρ_br` drops out of the wave speed entirely.

**(iii) θ and u have different physical natures.** θ is a **compact scalar phase** (periodic, the conjugate to the order-parameter density). u is a **non-compact vector displacement**. For the Maxwell structure, `(θ, u)` must form a 4-vector (scalar + vector potential), but they are fields of different tensor character with different mechanical origins. An order-parameter phase couples to its **conjugate density** (through the Josephson relation `∂_t θ = −μ`), NOT to the **displacement divergence** `∇·u`. There is no variational principle that produces the coupling `ε(∂_t u)·(∇θ)` from the medium's mechanics.

### Why this is not the same as the superfluid velocity structure

One might object: in a superfluid, the velocity is `v_s = (ℏ/m)∇θ`, and the kinetic energy `½ρ v_s²` involves `∇θ`. So phases DO couple to gradients. But:
- The superfluid phase θ is the phase of `ψ = √ρ e^{iθ}` — the **medium's** phase. It couples to the **density** ρ (its conjugate), not to a **displacement** u.
- The brane-order phase θ_B is the phase of `χ_B = |χ_B| e^{iθ_B}` — a **different** order parameter. It couples to the **brane-order density** (its conjugate), not to the displacement u.
- Neither phase naturally couples to `∇·u` through the Maxwell structure. The coupling would have to be **postulated** — and it is not a small postulate; it is a fundamental restructuring of the light sector from MacCullagh to Maxwell.

### The one escape (and why it fails)

Could θ_B be the **time component of a 4-vector whose spatial components are u**? If so, the Maxwell structure arises naturally. But:
- u is a **3-component** in-plane vector. θ_B is a **scalar**. A 4-vector would need 3+1 components: `(θ_B, u^a)`. But there is no reason the brane-order phase and the material displacement form a 4-vector — they are different fields with different transformation properties.
- If the ordered state is a **condensate of shear quanta** (u quanta), then χ_B is the amplitude and θ_B is the phase of the condensate. In that case, θ_B IS related to u. But u is a vector (3 components), while θ_B is a scalar — they cannot be amplitude and phase of the same field. At best, θ_B could be the phase of a **scalar** condensate of longitudinal shear quanta, while u is the transverse displacement — but then θ_B couples to `∇·u` (the longitudinal part), not to u as a whole, and the gauge structure doesn't close for the transverse modes.

### Bottom line on θ

The C5-φ sub-hypothesis is **structurally broken**. The phase of a complex order parameter has its own dynamics (Josephson relation) and couples to its conjugate density, not to the displacement divergence. The Maxwell-like coupling `½ε(∂_t u + ∇θ)²` would have to be postulated as a new constitutive term, changing the light sector from mechanical (MacCullagh) to gauge-field (Maxwell). This is not "θ supplies the missing φ" — it is "replace MacCullagh with Maxwell and add an electric-type energy," which is a bigger change than the sub-hypothesis claims, and one with no mechanical provenance.

**Do NOT carry θ into rung φ.** The C5 obstruction will persist at Gate L(a-iii) regardless of whether χ_B is real or complex. The χ_B route faces the same `FAIL_C5_LONGITUDINAL_ZERO_MODE` (or `FAIL_COUPLE_STRESS_NOGO` via the four-way chain) that killed pathA_35.

---

## 6. Cheapest decisive test for rung W

**The test:** static χ_B–P kink calculation — does the in-plane P survive at the wall centre?

**What to compute:**

1. Write the minimal static free energy (1D in w):
```
E/A₃ = ∫ dw [ ½κ_B(∂_w χ_B)² + f_B(χ_B)
              + ½κ_P(∂_w P^i)² + ¼r(|P|²−1)²
              + α χ_B (P^w)² ]
```
where `f_B(χ_B) = ½a(χ_B²−1)²` (φ⁴ double-well, minima at χ_B=±1).

2. Solve the static Euler–Lagrange equations for χ_B(w) and P^i(w) simultaneously.

3. **The kill/save criterion — check `P_∥(w=0)`:**

   - **KILL:** `P_∥(w=0) = 0` (the χ_B gradient or the χ_B–P coupling forces the in-plane P to vanish at the wall centre → in-plane shear starved → density no-go re-created by another route).
   - **SAVE:** `P_∥(w=0) ≠ 0` (or P_∥ is a free Goldstone — zero-frequency mode at k=0) → in-plane arrows live at the wall centre → shear is not starved → density no-go genuinely escaped.

**My prediction:** SAVE. The anisotropy `α χ_B (P^w)²` penalizes P along ŵ, not P_∥. At the wall centre (χ_B=1), P^w is gapped to zero and P_∥ is a free S² Goldstone. The χ_B gradient does not enter the P_∥ equation (the anisotropy is local in χ_B, not in ∂_w χ_B). There is no mechanism for gradient-direction pinning. The test should confirm P_∥ ≠ 0.

**Why run it anyway:** (a) it is cheap (static 1D ODE, analytical for φ⁴ kink); (b) it is the direct analog of the C.6 test that killed the density route (`pathA_25_gateRC_cubic.md` lines 54–58), so it is the fair comparison; (c) there could be an unexpected coupling through the simultaneous solution (the χ_B equation back-reacts to P, and the effective χ_B profile could shift in a way that affects P_∥ — unlikely, but able-to-fail); (d) the negative control (Frank-only, no anisotropy → P_∥ should be unconstrained, not necessarily zero) and the density-no-go control (Cpin-style coupling `χ_Cpin(P·k)²` → P_∥=0) should both fire.

**What this test does NOT decide:** whether the light sector is dynamically consistent (C5, closure, mode-count). That is Gate L, and the χ_B wall doesn't change those sub-hurdles. The rung-W test only checks whether the wall starves the shear — and I expect it doesn't.

---

## 7. Prior-art anchors

| System | Relevance | Key lesson |
|---|---|---|
| **³He A–B interface** (Volovik) | Closest physical analog: two phases of one medium, codim-1 interface, complex order parameter | Interface is **externally maintained** (pressure/temperature-tuned coexistence), not self-localizing (`medium_requirements_and_prior_art.md` line 66). The χ_B wall faces the same issue: the double-well must be tuned to coexistence (degenerate minima). The analog license concedes this. |
| **Smectic-A / Smectic-C** (de Gennes) | Director orientation relative to layer normal — the density no-go mechanism | In SmA, director IS along the layer normal (P_∥=0 — exactly the density no-go). SmC (tilted director, P_∥≠0) requires a **separate tilt instability**. Confirms: getting P in-plane requires a mechanism **independent** of the layering. The χ_B route provides this (the anisotropy is independent of the wall). |
| **Cosserat / micropolar continua** (1909–present) | Rotational elasticity with independent micro-rotation DOF | Cosserat continua CAN support MacCullagh-type shear, but require **independent micro-rotation DOF** (= the P^i arrows). The couple-stress closure is the hard part — exactly where pathA_35 Gate L failed. The χ_B wall doesn't change the Cosserat structure. |
| **Ferroelectric domain walls** (e.g. BiFeO₃) | π₀=Z₂ scalar domain wall, stable, carries bound modes | Confirms a scalar double-well gives a stable wall with bound states. BUT: in ferroelectrics, P IS the order parameter (no separate χ_B). The χ_B route has an extra layer of independence (χ_B gates P), which is the cost. |
| **Lifshitz points** | Boundary between homogeneous and modulated phases | The χ_B wall is NOT a Lifshitz phenomenon (single wall, not modulation). But the Lifshitz point is where the density-smectic route lived (finite-k instability → triad, `pathA_25_gateRC_cubic.md` lines 37–41). The χ_B route avoids the Lifshitz point entirely — a structural advantage. |
| **Superfluid vortices / interfaces** | Phase of complex order parameter, Goldstone modes | The phase of a superfluid order parameter couples to the **density** (Josephson relation), not to displacement fields. Supports the §5 finding that θ_B cannot serve as the C5 φ. |
| **Rubakov–Shaposhnikov** (1983) | Domain wall brane in higher dimensions | The original "do we live inside a domain wall?" paper. The wall is postulated (not derived from a single scalar). The χ_B route is in the same tradition — postulate the wall, test consistency. |

---

## 8. Summary — the three-death pattern and where χ_B lands

| Route | What died | Why | χ_B escape? |
|---|---|---|---|
| **T1** (little-arrows wall) | No stable wall | P on connected S³, π₀=0 → unwinds | **Yes** — χ_B double-well gives π₀=Z₂ |
| **Density-smectic** (pathA_25 R/C) | Codim-1 brane kills light | Cpin pins P along layer normal → P_∥=0 | **Yes** — independent couplings, no layer normal |
| **Shear-surface** (pathA_35 Gate L) | Light sector inconsistent | Four-way: closure ⊥ mode-count ⊥ C5/φ ⊥ u_w-gap | **No** — same light sector, same no-go |

The χ_B route escapes the first two deaths but **inherits the third**. The C5 obstruction is not touched by changing the confinement mechanism from a frozen profile to a dynamical order parameter. The θ sub-hypothesis — the proposed C5 fix — is structurally broken (§5).

**The honest assessment:** the χ_B material-state closure is the best brane-existence mechanism the program has produced. It is the first route that gives a genuinely stable, self-consistent wall that doesn't starve the in-plane shear. But it is a better *wall* than a *light substrate*. The light sector's problems are constitutive (MacCullagh + P–u + C5), not geometric (confinement profile). Changing the geometry doesn't fix the constitution.

**If the program proceeds:** run the cheapest test (§6) to confirm rung W passes, then go directly to the C5 question — NOT via θ (dead, §5), but by asking whether ANY on-brane scalar can serve as the Maxwell φ without being u_w or a new field. If the answer is no (and I believe it is), the χ_B route dies at Gate L(a-iii) for the same reason pathA_35 did, and the program should redirect to either (a) accepting MacCullagh without gauge closure (a toy light that is not Maxwell), or (b) the lattice/Wen route that was previously declined.

---

## 9. One-line verdicts (for quick reference)

- **Rung W (wall stable AND light-permitting):** `GENUINE_COIN_FLIP`, leaning `LIKELY_PASS` — the wall is stable and the shear is not starved, but it's the same light sector that already died.
- **Fresh no-go at rung W:** **No** — codim-1 localization and in-plane shear are compatible for χ_B. The no-go is inherited from Gate L, not fresh.
- **C5-φ sub-hypothesis (θ as Maxwell φ):** `FATAL_FLAW` — order-parameter phase couples to conjugate density, not to displacement divergence. Do NOT carry into rung φ.
- **Cheapest decisive test:** static χ_B–P kink, check P_∥(w=0). Expected SAVE. KILL if P_∥=0.
- **Load-bearing assumption:** χ_B is an independent scalar order parameter (not derived from P). If dependent → T1 returns.
- **Drift cost:** ≥6 new inputs (real χ_B) or ≥9 (complex with θ). Significant, conceded under analog license.
