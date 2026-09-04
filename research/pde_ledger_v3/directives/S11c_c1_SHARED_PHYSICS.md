# S11c-c1 — SHARED PHYSICS (curved bulk closure: the nonlocal DtN/impedance and the permeable face response)

**Orchestrator-authored, 2026-09-03, under the S11c decision list (`directives/S11c_decisions.md`, N-series)
and `CLAUDE.md`.** This is the physics authority read by the two blind S11c-c1 engines: the SymPy engine imports
the inherited model via the fold `load_model` (`scripts/ledger_fold.py`) over the atomic frozen base
`scripts/S11c_b_exports.py` — which already carries the whole F9-resolved S11b + S11c-a + S11c-b model — binding
only its declared `IMPORT_KEYS` from that fold (§7; the generate-over-a-frozen-base topology of
`directives/export_ledger_bind_closure_design.md`); the Wolfram engine imports nothing and reconstructs every
object from the supplied setup below and its cited sibling specs. This document is an **obligation-to-compute
specification**, not a script and not a record of results. There is no acceptance value to withhold beyond the
value-free packages named below; there is no terminal `VERDICT`, `PASS`, or `FAIL`.

**This is the FIRST unit of the two-way split of the ratified S11c-c row** (`directives/S11c_decisions.md`, N1/N2
table; the split refines that row's boundary per N2, decided with the two spec-review legs —
`directives/_measurements/S11c_c_spec_review.md`). S11c-c1 solves the curved bulk and closes the face; the fold
of that closed response into the S11c-b slab operator and the re-extraction of the coupling kernel are **S11c-c2**
(the self-energy operator). c1 exports the closed face response c2 consumes.

## 0 · Scope

S11c-c1 takes the non-uniform background, its two anchoring branches, and the tilted-face shape-derivative
substrate (S11c-a, T-a..T-i), and **closes the bulk at the face**: it solves the perturbed **curved two-face
outgoing bulk** acoustic problem to obtain the nonlocal **Dirichlet-to-Neumann / impedance** operator relating
the bulk pressure perturbation at each curved face to that face's bulk normal velocity, and composes it with the
inherited interfacial mass balance and face closure to obtain the permeable curved face response — the bulk
pressure `δp_s`, relative flux `J_s`, and traction `t_s` as functions of the face velocity `V_s` and the
constitutive operand `μ_θ`. Concretely it asks for: the first-shape-order curved-face impedance/DtN **operator**
`Z_s(ω;k,k′)` (a two-momentum nonlocal kernel — §3a), its bulk-regime structure, branch selection, face parity,
and its reactive/radiative split via the operator Hermitian part; the permeable curved face response (S11b's B0c
on curved faces); the response's degenerate/Fredholm loci; and its dissipation audit checked against an
independent outgoing-flux route. Every object is multigraded by `(ε,η,σ_W)` and carries restored physical
dimension.

The following remain outside this build:

- **S11c-c2** (the second unit of the split): the fold of the closed face response into the S11c-b slab operator,
  the re-extraction of the closed off-diagonal coupling kernel from the **closed full operator**, and the
  resulting coupled nonlocal **self-energy operator**. ⛔ S11c-c1 emits the closed face response `(δp_s, J_s,
  t_s)(V_s, μ_θ)` and the DtN operator; it does **not** substitute them into `S11CB_SLAB_OPERATOR` or
  `S11CB_COUPLING_KERNEL` and does **not** re-derive any slab row;
- S11c-b: the variable-coefficient stored-energy basis, the divergence-form slab operator with its symbolic
  `δp_s`/`Λ_I(ω)` face slots, the off-diagonal coupling kernel, `μ_θ`. These are **imported unchanged** (`N1`),
  per-engine-verified with the S11c-b cross-engine residual deferred
  (`steps/S11c_b_variable_coefficient_operator.md`); c1 consumes `S11CB_MU_THETA_OPERATOR` as the named operand
  `μ_θ` feeding `μ_s`, and re-derives nothing of S11c-b;
- S11c-d: any profile-conditioned scattering, Born, Bloch, WKB, resonance, or spectral object; ⛔ **no global
  dispersion relation** `ω(k)` for a generic `W_bg(x)` is requested (`N5`). c1 emits the DtN and response as
  **operators/kernels**, ⛔ not their spectrum, and the **profile-conditioned resolvent / singular set** of the
  operator loci is S11c-d's, not c1's (§3b);
- S11c-e: the flux-normalized dimensionless conversion form, the leakage observable, the confinement
  interpretation (`N13`), and the bench-optics comparison (`N7`);
- all nonlinear DC, harmonic, sideband, intensity, and soliton questions (`N10`/`G14`).

⛔ **The convective bulk operator is NOT an S11c task (`N11a`, settled in the decision list).** `v_bulk_normal_0`
is inherited as a standing rest-frame limit; the rest-frame bulk wave equation `∂_t²φ=c_s0²∇₄²φ` (§1b) stands
unmodified. Its exclusion is legitimate at first shape order (the drain-projection correction is `O(σ_W²)`; the
outward normal obeys `n̂_s·ŵ − s = O(σ_W²)`). ⚠ The stated validity domain is **non-uniform near grazing** and is
made precise in §2b; ⛔ `v_bulk_normal_0` is never aliased to `v_0`.

There is no acceptance value to withhold in S11c-c1. There is no terminal `VERDICT`, `PASS`, or `FAIL`.

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

### 1a · Inheritance, degrees of freedom, and the sector split

The supplied background ansatz, anchoring branches, face maps, level sets, orientation, measures, trace law,
dynamic window, and interface laws are **S11c-a §§1–3 imported unchanged** (`directives/S11c_a_SHARED_PHYSICS.md`);
the tilted-face first-order shape derivatives T-a..T-i are the S11c-a §4 exports. The SymPy engine imports
T-a..T-i and `S11CB_MU_THETA_OPERATOR` from the fold over the atomic frozen base `scripts/S11c_b_exports.py`
(§7); the Wolfram engine re-derives every consumed object from the sibling specs, importing nothing. This document amends
and extends that setup with the bulk-closure constructions of §§2–3; where it does not amend, the S11c-a text
governs.

The slab degrees of freedom are S11b's, with the two independent face variables `{ζ_+,ζ_-}` (S11c-a §3a):

```text
u = u_L + u_T ,   ∇×u_L = 0 , ∇·u_T = 0   (in-plane displacement; the sector labels of S11c-b §1a) ;
θ(x,t)     Eulerian densification, ρ_4D = ρ_4D⁰(1+θ) ;
ζ_+, ζ_-   the two independent face variables, δW ≡ ζ_+ − ζ_- , ζ_c ≡ (ζ_+ + ζ_-)/2 ,
           ζ_s = ζ_c + s δW/2 ,   e_W ≡ δW/W₀ ,   s ∈ {+1,−1} .
```

⛔ `ζ_c` is an independent face DOF; ⛔ no S11c-c1 computation may set `ζ_c=0` or replace the two face variables
by a thickness-only ansatz, except in an explicitly named centre-fixed uniform regression.

### 1b · Bulk acoustics, radiation condition, and the branch object — SUPPLIED

Inherited from S11b §1/§1b/§2 and S11c-a §1b, unchanged. Per the rest-frame limit (`N11a`),

```text
v_bulk = ∇₄φ ,        δp = −ρ_m ∂_tφ ,        ∂_t²φ = c_s0² ∇₄²φ ,
j = ρ_4D v_bulk ,     ∂_tρ_4D + ∇₄·j = 0 .
```

`ρ_m` (bulk mass density) and `c_s0` (bulk sound speed) are **constant** — they do **not** vary in-plane
(`directives/S11c_b_SHARED_PHYSICS.md:202`). Both half-spaces exterior to the slab are bulk: the upper beyond the
upper face, the lower beyond the lower face; the slab interior is **not** bulk, so the **two half-spaces are
disconnected** (they share no bulk region), and the impedance is a **per-face** operator — the curvature does not
create propagation between the two exterior regions. Only the perturbations carry harmonic time dependence.

The **impedance** is the S11b object (`directives/S11b_SHARED_PHYSICS.md:92`)

```text
Z ≡ (bulk pressure perturbation at a face) / (bulk material's OUTWARD normal velocity at that face) .
```

The **radiation condition** (`:93`) is supplied and binding: in each half-space retain **only** waves carrying
energy **away** from the slab (real normal wavenumber) or **decaying** away from it (imaginary normal
wavenumber); ⛔ **there are no incoming waves from infinity**. Boundedness alone does **not** select a branch
when the normal wavenumber is real; the sign of the selected normal wavenumber in each half-space is fixed under
the harmonic convention and is part of the computation, not stated here.

The **branch object** `q_out(ω,k)` (`:116`) is the root of the bulk dispersion selected by the radiation
condition — the bulk **normal** wavenumber for an in-plane wavenumber `k`; its branch points sit on the real axis
at the **sound cone `ω = ±c_s0|k|`**. A complex frequency is continued along the supplied path that keeps those
branch points fixed; ⛔ nothing in S11c-c1 may re-select the branch to delete, suppress, or reclassify a growing
or decaying object (`:145`). At `q_out = 0` the two bulk solutions coalesce — the **degenerate/grazing** point
(`:150`). ⚠ `q_out` is the radiation-selected root **inside** the bulk operator; ⛔ it is **not** an ansatz `k`
for the curved problem, and the three regimes (`q_out²` positive/negative/zero) are structure to keep **live**,
never a value to freeze (§3a).

### 1c · The flat-face reference (B0b/B0c) — SUPPLIED as the uniform-limit operand only

On the uniform (flat, `η=σ_W=0`) background S11b solved this bulk closure: the flat impermeable impedance `Z`,
its three regimes and per-face inertial loading (B0b), and the flat permeable face response (B0c). These are
supplied S11b objects, imported/re-derived **only** as the `(η,σ_W)→(0,0)` regression operand of §5c:

```text
S11B_Z_IMPERMEABLE , S11B_Z_BY_REGIME , S11B_Z_BY_PARITY , S11B_ADDED_MASS , S11B_GRAZING_BEHAVIOUR ,
S11B_FACE_RESPONSE , S11B_FACE_RESPONSE_COEFFS , S11B_DEGENERATE_LOCI_* ,
S11B_PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY .
```

S11c-c1 states nothing here about their values and does **not** obtain the curved-face objects by substituting a
profile into them; the curved objects are constructed from §§2–3 and the flat objects are the regression target
only. ⚠ Two S11b-measured properties are inherited and must be carried, not restated as answers: (1) a **real**
part of `Z` in the propagating regime is **bulk radiation to infinity**, present even with impermeable faces, so
⛔ `Re Z` is not by itself evidence of transfer **through** the interface (`:760`); (2) the per-face inertial
loading acts with the **same sign against each face's OUTWARD acceleration**, so the global-`w` signed pair
`(X,−X)` is a convention artifact and ⛔ is not consumed as an inertia (`:614`).

### 1d · The face closure laws the bulk solve closes — SUPPLIED (note the Λ-channel placement)

The flat-face closure, affinity, response kernels, traction, and kinematic balance are exactly S11c-a §1c/§3b
(inherited by S11c-b §1c), stated here as the laws c1 closes the bulk against, per anchoring `α∈{L,M}` and face `s`:

```text
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ / ρ_br,bg⁰ ,   μ_θ = S11CB_MU_THETA_OPERATOR (the held-fixed operand) ,
J_s = Λ_A(ω) 𝒜_s + Λ_V(ω) V_s ,           ⇐ the FLUX closure carries only Λ_A, Λ_V ,
Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,   I ∈ {A,V,X} ,   the three τ_I independent ,
t_s = −(δp_s + Λ_X(ω) 𝒜_s) n̂_s ,           ⇐ the TRACTION carries Λ_X (a reciprocal-traction channel) ,
n̂_s·v_bulk,s = V_s + J_s/ρ_m               (interfacial mass balance — kinematics, not a result) .
```

⚠ **Λ-channel placement (a fold from the spec review).** `Λ_X` appears **only** in the traction `t_s`, **not** in
the flux closure `J_s` and **not** in S11c-a's T-i closure shape-derivative (`S11CA_CLOSURE_SHAPE_DERIV`, which is
the shape derivative of `J_s − Λ_A𝒜_s − Λ_V V_s = 0` only). c1 computes `δp_s`, `J_s`, and `t_s` (the last
carrying `Λ_X`) as the face response; the routing of `J_s` into the mass row and `t_s` into the mechanical rows
is S11c-c2's, not c1's. The face objects `δp_s`, `V_s`, `J_s`, `𝒜_s`, `t_s`, `n̂_s`, and the face measure `a_s^α`
are the **curved-face** objects of the S11c-a substrate (T-a `S11CA_FACE_NORMAL`, T-a′ `S11CA_CONORMAL_DERIV`,
T-a″ `S11CA_FACE_MEASURE_SHAPE_DERIV`, T-b `S11CA_FACE_VELOCITY`, T-c `S11CA_RELATIVE_FLUX`, T-c′
`S11CA_KINEMATIC_BALANCE`, T-d `S11CA_TRACTION`, T-e `S11CA_FACE_SHIFT`, T-i `S11CA_CLOSURE_SHAPE_DERIV`), ⛔ not
flat Cartesian objects. In S11c-b the bulk pressure `δp_s` and the kernels `Λ_I(ω)` remained **symbolic**;
S11c-c1 supplies the bulk relation `δp_s = Z · v_bulk,s` (§3a) that eliminates `δp_s`, keeping `Λ_I(ω)` symbolic.

Equations of motion are obtained by S11b's method — balance laws, the binding virtual-displacement rule,
variational derivatives with held-fixed fields named, and prescribed external virtual work — **not** by putting
an irreversible response kernel in an ordinary action.

---

## 2 · The supplied curved-bulk problem and the standing rest-frame limit

All definitions in §§1–2 are inputs. The expansions, coefficients, kernels, orders, regimes, branch selections,
face-parity blocks, and cancellations in §§3–5 are **outputs** computed from them; ⛔ none is stated here.

### 2a · The two disconnected curved half-spaces and the driven boundary data — SUPPLIED

For each anchoring `α∈{L,M}`, the upper bulk half-space is `{(x,w) : w > h_+^α(x,t)}` and the lower is
`{(x,w) : w < h_-^α(x,t)}`, with the S11c-a Eulerian graph heights and level sets

```text
h_s^α(x,t) = s W_bg^α/2 + ζ_s(χ(x,t),t) ,   F_s^α ≡ w − h_s^α ,   (S11c-a §3a; W_bg^L=W_bg(x), W_bg^M=W_bg(χ)) ,
```

the oriented outward unit normal `n̂_s^α` (T-a, obeying the binding orientation law `s(n̂_s^α·ŵ)>0`, ⛔ **not**
`sign(∇₄F_s^α)`), the conormal derivative `n̂_s^α·∇₄` (T-a′), and the true-area measure `a_s^α =
√(1+|∇_x h_s^α|²)` (T-a″). The **lab-`w` graph slopes** are opposite at first shape order
(`∂_x h_s^α = (s/2)∂_x W_bg^α + …`), a supplied-map identity. ⛔ **Do not infer a parity of any face object from
that.** The oriented **outward** normal, on which `Z` is defined, is `n̂_s^α = (−½∇_x W_bg^α, s)+O(|∇W_bg|²)`, so
its in-plane (tilt) component is the **same** for both faces; whether the first-shape-order correction to any
outward face object is even, odd, or mixed under `s→−s`, and whether it couples the `{δW,ζ_c}` combinations, is a
**computed** result (§3) — ⛔ not asserted here, in either direction.

In each half-space the bulk field `φ` obeys §1b under the §1b radiation condition, with the boundary data set at
the curved face by the kinematic balance `n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m` (§1d, T-c′): the face normal
velocity drives the bulk, and the bulk returns the pressure `δp_s` at the face. ⛔ The face is a **level set
`F_s^α=0`**, evaluated by the shifted-trace law (T-e, `S11CA_FACE_SHIFT`) — do **not** evaluate the bulk field at
the flat reference `w=sW_0/2` and discard the face shift; the retained first-shape-order dependence through
`h_s^α` (carrying `η`,`σ_W` via `W_bg`) is the curvature that distinguishes c1 from S11b's flat B0b.

### 2b · The standing rest-frame limit and its validity domain (`N11a`) — SUPPLIED

`v_bulk_normal_0` is the steady bulk-normal drain, inherited as the **inert rest-frame scope limit** of §1b
(`directives/S11c_a_SHARED_PHYSICS.md:49`; never aliased to `v_0`). The bulk operator solved in §3 is the
rest-frame `∂_t²φ=c_s0²∇₄²φ`; ⛔ **the convective operator `(∂_t+v_bulk_normal_0∂_w)²φ=c_s0²∇₄²φ` is not
constructed** (`N11a`; legitimate at first shape order — §0). ⚠ **The validity domain is non-uniform near
grazing** (a fold from the spec review): the S11b smallness condition `|q·v_bulk_normal_0/ω| ≪ 1` becomes
vacuous as `q_out→0` (grazing) while the rest-frame correction it is meant to bound does not, so it does **not**
by itself license the grazing behaviour §3a requests. Therefore:

- the requested **grazing** behaviour (`q_out→0`, §3a) is the **strict `v_bulk_normal_0=0` result** and is stated
  as such in the step record;
- away from grazing, each result is conditional on the S11b-derived domain, made explicit as the pair
  `|q_out·v_bulk_normal_0/ω| ≪ 1` **and** `|ω·v_bulk_normal_0|/(c_s0²|q_out|) ≪ 1` (the boundary-layer exclusion),
  together with the independent **subsonic** speed condition inherited from S11b (`v_bulk_normal_0 < c_s0`); the
  order of the limits (`v_bulk_normal_0→0` vs `q_out→0`) is stated explicitly;
- ⛔ the engines still compute only the rest-frame object; this domain is recorded by the step record, not carried
  as a term in any operator.

### 2c · Multigrade bookkeeping — SUPPLIED

The background bookkeepers are S11c-a's: the zero-jet contrast `η` and the supplied first-jet bookkeeper
`σ_W ≡ η W̄₀/L_W`, varied independently (by `η` and `L_W`); ⛔ no engine may replace `σ_W` by `η` or assign them a
common order. Wave and bulk perturbations (`u`, `ζ_c`, `δW`, `θ`, `φ`, `δp`, `V_s`, `J_s`) carry the amplitude
bookkeeper `ε`. **Every computed object is multigraded by `(ε,η,σ_W)` from its actual data dependency**; the
requested truncation is **first order in wave amplitude `ε` and first shape order in each background bookkeeper
`η` and `σ_W`** (`N5`/`N12`). ⛔ No term is removed merely because it carries both a wave and a background
bookkeeper; ⛔ no object's grade is stated here. A second spatial derivative of `W_bg` remains first order in
background amplitude and is not dropped; `|∇_x h_s|²` is `O(σ_W²)` and is not a first-shape-order ingredient.

---

## 3 · The bulk-closure construction

### 3a · The curved-face impedance / DtN — a two-momentum nonlocal OPERATOR

Solve the perturbed **curved two-face outgoing bulk** problem to obtain the operator relating the bulk pressure
perturbation `δp_s` at each curved face to that face's bulk **outward** normal velocity `v_bulk,s ≡
n̂_s^α·v_bulk,s` — the curved-face generalization of S11b's flat B0b `Z`. The construction rule is S11b's B0b,
extended to the curved level-set face of §2a: solve §1b in each half-space under the §1b radiation condition,
evaluate `v_bulk,s` and `δp_s` at the level-set face by the T-a normal, the T-a′ conormal, and the T-e shifted
trace, and expand to first shape order.

⭐ **The object is an OPERATOR on the in-plane face fields carrying BOTH momentum legs, ⛔ not a single-`k`
multiplier and ⛔ not a one-leg left-quantized symbol.** Emit **all** of: (i) the zeroth-shape-order **flat**
piece — the unperturbed half-space Fourier **symbol** `Z_s^{(0)}(ω,k)`, built from the radiation-selected
`q_out(ω,k)`; (ii) the **first-shape-order** piece as the **two-sided composition** built from the two-face
boundary matching — schematically `N_0 ∘ M_{h_s} ∘ N_0` together with the `Div(h_s ∇)` and `κ² h_s` terms
(`κ ≡ ω/c_s0`), where `M_{h_s}` is multiplication by the face height and `N_0` the flat normal-derivative/DtN
factor — retaining `W_bg(x)` and its in-plane jets as functions of `x`; and (iii) the **same** object as the
two-momentum kernel `Z_s^{(1)}(ω;k,k′)` carrying `Ŵ_bg(k−k′)` and **both** branch legs `q_out(ω,k)`,
`q_out(ω,k′)` **explicit**. ⛔ **Do not** emit the first-shape-order piece as a single-`k` multiplier
`Z(k;∇W_bg(x))` (the WKB/local freeze) **nor** as a left-quantized product `a(x,k)=W_bg(x)·σ(q_out(k))` carrying
only the **one** leg `q_out(k)`: both diagonalize the operator (delete the `q_out(k′)` output leg and the
profile's propagating↔evanescent mode mixing) and are exactly the freeze the non-uniform background forbids
(rule 17; `N5`: the object is the operator/kernel, its spectrum is S11c-d's). A rigid vertical translation of a
face (`k=k′`, `Ŵ_bg(0)`) must **cancel** from the operator — emit that cancellation as a **named residual**
`S11CC1_DTN_RIGID_SHIFT_OPERAND`/`S11CC1_DTN_RIGID_SHIFT_RESIDUAL`, ⛔ not a stated check.

⭐ **Emit, per anchoring and per face, keyed as a CAS map:** the flat symbol operand `Z_s^{(0)}`; the
first-shape-order operator/kernel operand `Z_s^{(1)}`; the assembled curved DtN operator; its structure across
the **regimes on each branch leg** — the ordered `3×3` regime pairs `(q_out(k)², q_out(k′)²)∈{+,−,0}²`, including
the grazing cases where the input leg, the output leg, or both vanish (with the behaviour **at** each grazing
zero). On real `ω`, split the operator into its **radiative/dissipative** part — the **Hermitian part** under the
true-area boundary pairing, `H_a = (Z + Z^{†_a})/2` (§3b), ⛔ not the entrywise `Re/Im` of a symbol — and its
**reactive** part, the anti-Hermitian `K_a = (Z − Z^{†_a})/(2i)`; report a **per-face added mass** `m_add`
(`p=m_add·∂_t²ζ`) **only** on a named block where `Z` is purely reactive and that local relation actually exists,
reported **per face** with the §1c sign convention against outward acceleration, ⛔ not a two-face sum. ⛔ Do
**not** report the sign of an individual kernel entry's real or imaginary part — off-diagonal entries change
phase under translation of the profile/Fourier basis and are not physical; branch/sheet signs are carried as
separate `q_out` data (§1b), not as kernel-entry signs. Emit the **face-parity** structure of the first-shape-order
correction (its parity under `s→−s`, and whether it couples the `δW` and `ζ_c` combinations of the two per-face
impedances — a **computed** block). ⛔ Report the operator as what it is — a permeability/memory kernel wherever
the algebra gives one — and ⛔ do not force it into a single local velocity-impedance if the algebra does not
give one.

```text
⇒ S11CC1_DTN_FLAT_SYMBOL , S11CC1_DTN_OPERATOR (composition form) , S11CC1_DTN_KERNEL (two-momentum k,k′) ,
  S11CC1_DTN_RIGID_SHIFT_OPERAND , S11CC1_DTN_RIGID_SHIFT_RESIDUAL ,
  S11CC1_DTN_BY_REGIME_PAIR , S11CC1_DTN_BY_PARITY , S11CC1_DTN_HERMITIAN_PART , S11CC1_DTN_REACTIVE_PART ,
  S11CC1_DTN_INERTIAL_LOADING (purely-reactive block only) , S11CC1_DTN_GRAZING_BEHAVIOUR ,
  S11CC1_DTN_TERM_ORIGINS .
```

### 3b · The permeable curved face response, its Fredholm loci, and the dissipation audit

Close the face by solving the supplied relations together, per anchoring and face, keeping the three `τ_I`
independent and `Λ_I(ω)`/`μ_θ` symbolic: the bulk relation `δp_s = Z · v_bulk,s` (§3a, the operator), the
interfacial mass balance `v_bulk,s = V_s + J_s/ρ_m` (§1d, T-c′), and the §1d flux closure and traction. Solve for
the face pressure `δp_s`, the relative flux `J_s`, and the traction `t_s` (carrying `Λ_X`) as functions of the
face velocity `V_s` and `μ_θ`, to first shape order, retaining the curved-face measure and normal. Because `Z` is
an **operator**, the elimination is an **operator inverse** `[I + (Λ_A/ρ_m²) Z]^{-1}`, ⛔ not a scalar division;
report the response as what it is — per face and per parity combination — and ⛔ do not force it into a single
scalar velocity-impedance if the algebra does not give one.

⭐ **Degenerate loci — the Fredholm condition, then the finite-dimensional CAS loci.** The response's loss of
solvability is failure of the operator `[I + (Λ_A/ρ_m²) Z]` to be invertible — a **Fredholm/resolvent
condition** on the function space and profile, ⛔ not the vanishing of a single scalar coefficient. Emit the
**formal noninvertibility condition** as a CAS object (the operator whose invertibility is in question, and its
symbol where it is diagonal). Then apply the §6 algebraic **locus protocol only** to (i) the flat
Fourier-diagonal symbol's degenerate loci and (ii) any genuinely finite-dimensional algebraic degeneracy of the
coefficients; ⛔ do not force the operator's profile-conditioned singular set through the algebraic protocol — the
profile-conditioned resolvent/spectrum is **S11c-d's** object (`N5`), named here as reserved, not solved.

⭐ **Dissipation audit — three distinct objects, on real `ω`, with the background measure `a_s^{0,α}`.** The
bare bulk `Z` carries **none** of the `Λ_I`, so its Hermitian part audits bulk **radiation only** and is a
different object from the permeable-response dissipation (which is where the `Λ_I`/`ωτ_I` dependence lives). Emit,
per regime pair and parity combination:

1. **Bulk radiation** — the Hermitian part of the bare DtN operator under the true-area pairing
   `⟨f,g⟩_∂Ω = Σ_s ∫ a_s^{0,α}(x) f_s^*(x) g_s(x) d³x`, `H_a[Z] = (Z + Z^{†_a})/2` — this is the §3a
   `S11CC1_DTN_HERMITIAN_PART` (the radiative part of the DtN split), reused here, ⛔ not a new object.
2. **Permeable-response dissipation** — the **two-port power-conjugate form** (S11b `:705-717`), the block
   Hermitian form of the closed port map `(V_s, μ_s) → (δp_s + Λ_X𝒜_s, J_s)` under the true-area pairing:
   whether it is sign-definite, on which `Λ_I` it depends, and how it varies with each `ωτ_I` (both limits). ⛔ Do
   **not** use `H_a[Z]` (bulk-only) for this. ⇒ `S11CC1_PERMEABLE_PORT_HERMITIAN`,
   `S11CC1_PERMEABLE_DISSIPATION_VS_OMEGA_TAU`.
3. **The independent energy-balance route** (same sub-case: real `ω`, propagating, impermeable, `Λ_X⁰=0`, so it
   is a genuine second construction of the sign). ⛔ The face operand is **not** `½Re(δp_s V_s^*)` — that equals
   the bulk Poynting flux by the acoustic identity and never sees `t_s`. Instead: the **face operand** is the
   true-area **traction** pairing `½ Re Σ_s ∫ a_s^{0,α} t_s·v_{face,s}^* d³x`, built from the §3b `t_s` object
   (in this slice it reduces to `−½ Re Σ_s ∫ a_s^{0,α} δp_s V_s^* d³x`); the **bulk operand** is the outgoing
   **Poynting flux** computed from `φ` on a control surface **in the half-space / at `|w|→∞`**, ⛔ not `δp` at the
   face. Emit both operands and their residual; a one-sided reversal of the `t_s` sign must move **only** the face
   operand (making the residual nonzero) — that is the shared-sign test. ⛔ This route uses **no** slab EOM row
   and does **not** import `S11CB_SLAB_OPERATOR` (that pairing is S11c-c2's, after the fold); a `t_s`-sign or
   outward-acceleration error is caught here, in c1, by the traction-vs-far-field comparison.

⛔ **Sign-definiteness at first shape order (a truncation caveat).** A flat evanescent channel is a nullspace of
the zeroth-order Hermitian form; the curvature-induced radiated power there is `O(η²)`, exactly the term the
first-shape-order truncation omits, so a **passive** operator can look indefinite purely from the truncation, and
an ordered propagating↔evanescent regime block is not itself Hermitian. ⇒ assemble **adjoint-related regime
pairs** before forming the Hermitian part; test the power identity through `O(η)`; restrict any sign-definiteness
claim to the subspace where the zeroth-order Hermitian form is **nondegenerate**, and on its nullspace emit the
typed token `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` (the `O(η²)` leakage there belongs to S11c-e). ⛔ Attribute any
loss to bulk radiation vs an interfacial channel by computation, and ⛔ do not read a propagating-regime `Re Z` as
interfacial transfer (§1c). ⛔ State no sign in prose.

```text
⇒ S11CC1_FACE_RESPONSE , S11CC1_FACE_RESPONSE_COEFFS , S11CC1_NONINVERTIBILITY_CONDITION ,
  S11CC1_DEGENERATE_LOCI_* (locus protocol, flat/finite-dim only) ,
  S11CC1_PERMEABLE_PORT_HERMITIAN , S11CC1_PERMEABLE_DISSIPATION_VS_OMEGA_TAU ,
  S11CC1_ENERGY_FACE_TRACTION_OPERAND , S11CC1_ENERGY_BULK_FARFIELD_FLUX_OPERAND , S11CC1_ENERGY_RESIDUAL .
```

The response objects `(δp_s, J_s, t_s)(V_s, μ_θ)`, the DtN operator, and the flat symbol are the c1 exports
S11c-c2 consumes.

---

## 4 · Objects to compute and emit

Every item is computed for both anchoring branches and both density representatives wherever density enters,
carries the object's computed `(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component,
term, order, coefficient, parity, survival, regime, branch, or cancellation in prose.

- **The curved-face DtN operator** (§3a) ⇒ `S11CC1_DTN_FLAT_SYMBOL`, `S11CC1_DTN_OPERATOR` (composition),
  `S11CC1_DTN_KERNEL` (two-momentum `k,k′`), `S11CC1_DTN_RIGID_SHIFT_OPERAND`, `S11CC1_DTN_RIGID_SHIFT_RESIDUAL`,
  `S11CC1_DTN_BY_REGIME_PAIR`, `S11CC1_DTN_BY_PARITY`, `S11CC1_DTN_HERMITIAN_PART`, `S11CC1_DTN_REACTIVE_PART`,
  `S11CC1_DTN_INERTIAL_LOADING`, `S11CC1_DTN_GRAZING_BEHAVIOUR`, `S11CC1_DTN_TERM_ORIGINS`.
- **The permeable curved face response, loci, and dissipation** (§3b) ⇒ `S11CC1_FACE_RESPONSE`,
  `S11CC1_FACE_RESPONSE_COEFFS`, `S11CC1_NONINVERTIBILITY_CONDITION`, `S11CC1_DEGENERATE_LOCI_*`,
  `S11CC1_PERMEABLE_PORT_HERMITIAN`, `S11CC1_PERMEABLE_DISSIPATION_VS_OMEGA_TAU`,
  `S11CC1_ENERGY_FACE_TRACTION_OPERAND`, `S11CC1_ENERGY_BULK_FARFIELD_FLUX_OPERAND`, `S11CC1_ENERGY_RESIDUAL`.
- **The representation-invariance package** (§5a) ⇒ `S11CC1_REP_INVARIANCE_EULERIAN_OPERAND`,
  `S11CC1_REP_INVARIANCE_HANZAWA_OPERAND`, `S11CC1_REP_INVARIANCE_RESIDUAL`.
- **The one-sided independence control** (§5a — mandatory `W_bg` tilt; `Σ_E`/`μ_R,bg` reserved for c2) ⇒
  `S11CC1_CONTROL_INDEPENDENCE_BASE_OPERAND`, `S11CC1_CONTROL_INDEPENDENCE_CORRUPTED_OPERAND`,
  `S11CC1_CONTROL_INDEPENDENCE_RESIDUAL`.
- **The source-level form ablations** (§5b) ⇒ `S11CC1_CONTROL_FORM_BASE_OPERAND`,
  `S11CC1_CONTROL_FORM_ABLATED_OPERAND`, `S11CC1_CONTROL_FORM_RESIDUAL`.
- **The uniform-limit regression** (§5c) ⇒ `S11CC1_UNIFORM_LIMIT_S11CC1_OPERAND`,
  `S11CC1_UNIFORM_LIMIT_S11B_OPERAND`, `S11CC1_UNIFORM_LIMIT_RESIDUAL`.
- **The zero-jet regression** (§5d) ⇒ `S11CC1_ZERO_JET_S11CC1_OPERAND`, `S11CC1_ZERO_JET_S11B_OPERAND`,
  `S11CC1_ZERO_JET_RESIDUAL`.
- **The branch/momentum liveness control** (§5e — both legs) ⇒ `S11CC1_CONTROL_BRANCH_BASE_OPERAND`,
  `S11CC1_CONTROL_BRANCH_SIGNFLIP_INPUT_OPERAND`, `S11CC1_CONTROL_BRANCH_SIGNFLIP_OUTPUT_OPERAND`,
  `S11CC1_CONTROL_BRANCH_MOMENTUMFREEZE_OUTPUT_OPERAND`, `S11CC1_CONTROL_BRANCH_MOMENTUMFREEZE_INPUT_OPERAND`,
  `S11CC1_CONTROL_BRANCH_RESIDUAL`.
- **Dimensions and homogeneity** (§6) ⇒ `S11CC1_DIMENSIONS`, `S11CC1_HOMOGENEITY_BASE_OPERAND`,
  `S11CC1_HOMOGENEITY_CONTROL_OPERAND`, `S11CC1_HOMOGENEITY_RESIDUAL`.

No task may be replaced by a hand-typed expression carrying its anticipated content. Every comparison emits
operand A, operand B, and `A−B`; ⛔ no residual value is supplied.

---

## 5 · Independent routes and controls

Every comparison emits operand A, operand B, and the computed residual. No residual value is supplied. An honest
zero from a symmetric cancellation is a valid finding (rule 6); ⛔ a control whose baseline and corrupted operands
are identical calls (`A−A≡0`) is **not** a control (rule 2 corollary 3) and must never be emitted as a structural
zero.

### 5a · Representation-invariance routes (`N4`/`N6`) — the genuine control

For each physical anchoring, compute the curved-face DtN operator (§3a) and the permeable face response (§3b) by
two independent routes:

1. **direct Eulerian** boundary-perturbation of the half-space bulk problem about the flat face, from the §§1b–2
   setup, evaluating at the level-set face by the T-a/T-a′/T-e substrate;
2. **radiation-preserving boundary flattening (Hanzawa / layer-potential)**: a **cutoff/Hanzawa** change of
   coordinates that equals the face map near the boundary and the **identity at infinity**, pushing the boundary
   perturbation into the bulk operator's coefficients only within a bounded layer — with its transformed
   radiation condition stated explicitly. ⛔ Do **not** use the global scaling `w′=[w−ζ_c]/[W_bg+δW]` over the
   whole half-space: it is secular at infinity (its transformed operator has coefficients growing with the normal
   coordinate, and an outgoing wave acquires a secular first variation `∝ w′`), so it does not preserve the
   radiation branch and the residual would not compare two well-defined constructions. A **boundary-integral /
   layer-potential** construction of the DtN is an acceptable alternative second route provided its outgoing
   kernel is the §1b radiation-selected one. ⚠ This is a **bulk** problem: `ρ_m` is constant and `Σ_E` is absent
   from `Z`, so ⛔ do **not** invoke S11c-a's slab-density map `N4: Δρ=δρ_E+u·∇ρ⁰` here (it belongs to the slab
   fields, not the bulk field `φ`), and ⛔ do not carry an `S11c-a`-style `MATERIAL_ADVECTED` **density** operand
   into this route — the two routes are Eulerian-direct vs Hanzawa/layer-potential of the **same** bulk problem.

The **anchoring** (LAB_HELD vs MATERIAL_ADVECTED face geometry, S11c-a §2c) is held fixed across the two routes,
and the branch selection of §1b is applied identically in both. Emit both operators and their difference under
`S11CC1_REP_INVARIANCE_EULERIAN_OPERAND`, `S11CC1_REP_INVARIANCE_HANZAWA_OPERAND`,
`S11CC1_REP_INVARIANCE_RESIDUAL`, keyed by object and anchoring.

For the one-sided independence control (`N6`), the load-bearing same-order channel that reaches c1's objects is
the **tilt** of the face (T-a/T-e). **The tilt mutation is MANDATORY on the DtN**: reverse only the `x¹` first jet
of `W_bg` in the upper-face direct-route level set `F_+^α` / boundary data, leaving the Hanzawa route and every
other source unchanged, and recompute the DtN operator and the response (the DtN's first-shape-order piece carries
this source; omitting it would leave the §3a parity untested). ⚠ **The profile-advection channel (`Σ_E`) and the
modulus channel (`μ_R,bg`) are structurally absent from c1's `(V_s,μ_θ)` objects** and are **reserved for
S11c-c2**: with `μ_θ` symbolic, density advection does not enter `μ_s=μ_θ/ρ_br,bg⁰` at linear wave order, `Σ_E`
is not in the B0c algebra (`ρ_m` constant), and `μ_R,bg` lives inside the supplied `μ_θ` operator — so mutating
either on c1 is a forbidden `A−A` control or forces a premature expansion of `μ_θ` into slab DOFs (c2's work). If
a density-representation test on the **response** is wanted in c1, mutate `ρ_br,bg` in `μ_s=μ_θ/ρ_br,bg`, ⛔ never
`Σ_E`. An object is omitted from a mutation only when its computed and emitted dependence shows the source
**structurally absent**; a residual `A−A≡0` from re-running an unaffected object is not a control. Emit the
uncorrupted-route operand, the corrupted-route operand, and their residual under
`S11CC1_CONTROL_INDEPENDENCE_BASE_OPERAND`, `S11CC1_CONTROL_INDEPENDENCE_CORRUPTED_OPERAND`,
`S11CC1_CONTROL_INDEPENDENCE_RESIDUAL`, keyed by object and anchoring. The control does not authorize editing an
already-emitted DtN, response, or residual.

### 5b · Source-level form ablations, one direction at a time

The form control in c1 is the **`W_bg` tilt** ablation. For each direction `i∈{1,2,3}` **separately**, create a
formal first-jet ablation in which only `∂_{xᶦ}W_bg` is set to zero, retaining the other two directions, both
anchorings, both density representatives, and every constitutive law; the density gradient definitionally induced
by `W_bg` (via S11c-a §2b) co-varies. Recompute the §3a DtN and the §3b response from that source. ⛔ Do not
ablate an already-computed operator, ⛔ do not drop all three directions simultaneously, and ⛔ do not use an `η`
rescaling as this control (a coefficient rescale tests arithmetic; only a form change tests physics). ⚠ **The
`μ_R,bg` jet is structurally absent from c1's `(V_s,μ_θ)` objects** — it lives inside the supplied `μ_θ` operator,
which c1 carries as an opaque symbol — so its ablation on c1 would be a forbidden `A−A`; the `μ_R,bg` form control
is **reserved for S11c-c2**, where `μ_θ` is composed with the slab variables. Emit the baseline operand, the
independently recomputed ablated operand, and their residual under `S11CC1_CONTROL_FORM_BASE_OPERAND`,
`S11CC1_CONTROL_FORM_ABLATED_OPERAND`, `S11CC1_CONTROL_FORM_RESIDUAL`, keyed by
`{object,anchoring,density,direction}`.

### 5c · Uniform-limit regression smoke test

For the curved-face DtN operator and the permeable face response, independently obtain the `(η,σ_W)→(0,0)`
operand and the corresponding S11b flat object (§1c). Emit both and their residual under
`S11CC1_UNIFORM_LIMIT_S11CC1_OPERAND`, `S11CC1_UNIFORM_LIMIT_S11B_OPERAND`, `S11CC1_UNIFORM_LIMIT_RESIDUAL`, keyed
by object, regime pair, and parity. ⛔ This is a **regression smoke test only**: it checks the curved construction
reproduces the flat bulk closure at zeroth shape order; it does **not** validate the first-shape-order curvature
coefficient/sign/parity/branch structure, which vanish at `(η,σ_W)=0`. `(η,σ_W)→0` is **not** an accepted
corruption of §5a/§5b/§5e.

### 5d · Zero-jet regression (`σ_W→0` at finite `η`)

Operand A: c1's DtN in the limit `σ_W→0` with `η` **retained** and `w₁` a **constant symbol** (`∇W_bg=0`, so the
two faces are parallel planes at `±W̄₀(1+η w₁)/2`). Operand B: the **UNMODIFIED** S11b half-space impedance
`S11B_Z_IMPERMEABLE` — ⛔ **with no `W_0→W̄₀(1+η)` substitution and no two-face re-solve**. Emit both and their
residual under `S11CC1_ZERO_JET_S11CC1_OPERAND`, `S11CC1_ZERO_JET_S11B_OPERAND`, `S11CC1_ZERO_JET_RESIDUAL`. ⛔ Do
**not** state the residual is zero. ⚠ **The target must be the bare half-space `Z`, not B0b re-solved at gap
`W̄₀(1+η)`** — the flat half-space impedance `Z_0=ρ_m ω/q_out` is thickness-**independent** (a uniform thickness
change is a rigid face shift; on-shell the rigid-shift kernel vanishes, §3a), so re-solving a finite-gap **cavity**
as operand B would reproduce the very `O(η)` cavity error this control exists to catch and drive the residual to a
false zero. The genuine test is whether c1's `σ_W→0` operand acquires a spurious `η`/thickness dependence — the
disconnected-half-space physics (§1b) forbids it, so a nonzero residual is a finding.

### 5e · Branch / momentum liveness control (rule 17)

The branch object `q_out` and the two-momentum operator structure are varying physics, not values to freeze.
Construct the DtN and response with `q_out(ω,k)` and the regime discriminant kept **live and symbolic**, then
form **same-domain** corrupted operands and recompute the DtN and response downstream at the source:

Fix the kernel convention once: in `Z_s(ω;k,k′)` with `Ŵ_bg(k−k′)`, `k` is the **output** leg and `k′` the
**input** leg (§3a). Then form:

- **sign flip** (the S11b opposite-sheet error): flip the branch sign on **one** leg only — the **input-leg** flip
  `q_out(ω,k′)→−q_out(ω,k′)` and, separately, the **output-leg** flip `q_out(ω,k)→−q_out(ω,k)` — recomputing each;
  ⇒ `S11CC1_CONTROL_BRANCH_SIGNFLIP_INPUT_OPERAND`, `S11CC1_CONTROL_BRANCH_SIGNFLIP_OUTPUT_OPERAND`;
- **momentum freeze** (collapsing the two-momentum kernel to a single-`k` multiplier, the §3a WKB freeze): the
  **output-leg freeze** `q_out(ω,k)→q_out(ω,k′)` and, separately, the **input-leg freeze**
  `q_out(ω,k′)→q_out(ω,k)`, recomputing each; ⇒ `S11CC1_CONTROL_BRANCH_MOMENTUMFREEZE_OUTPUT_OPERAND`,
  `S11CC1_CONTROL_BRANCH_MOMENTUMFREEZE_INPUT_OPERAND`. Running both one-leg freezes is the strongest coverage and
  removes the convention ambiguity that would otherwise let the two engines corrupt different legs.

Emit the live baseline operand and both corrupted operands and their residuals against the baseline under
`S11CC1_CONTROL_BRANCH_BASE_OPERAND` and `S11CC1_CONTROL_BRANCH_RESIDUAL` (keyed by corruption, object,
anchoring). ⛔ The corruption is applied **only** in the corrupted operand, never in an emitted physics object.
Independently check, on the real `ω` axis, that the baseline branch carries energy outward / decays (the §1b
radiation condition), then continue the **already-selected** branch without re-selection.

---

## 6 · Method, dimensions, locus protocol, and script obligations

The derivation method is S11b's: balance laws, the binding virtual-displacement rule, variational derivatives
with held-fixed fields named, the supplied sign conventions, the supplied radiation condition and branch object,
and prescribed external virtual work. It is not an ordinary action-principle derivation of retarded response.

The structural rule is binding: physical symbols may be combined by hand only in the supplied setup, background
ansatz, face maps, bulk acoustics, radiation condition, and supplied laws in §§1–2. Every §3–§5 expression is
reached by computation. Every control re-enters at an ansatz, map, level set, bulk operator, radiation condition,
or supplied law, **never at a result**.

Restore units and compute the `[L,T,M]` dimension of every emitted object. `η`, `σ_W`, `w₁`, `m₁`, `θ`, `e_W`
are dimensionless; `W_0`, `L_W` have length dimension; `ρ_m`, `c_s0`, `ρ_br,bg⁰`, the `Λ_I⁰`, and the `τ_I` carry
their S11b dimensions. Each impedance/response is a **ratio** — ⛔ state what of what before assigning a
dimension. For every homogeneity comparison emit both operands and the residual, and include a source-level
dimension corruption demonstrating that the check can fail ⇒ `S11CC1_DIMENSIONS`, `S11CC1_HOMOGENEITY_BASE_OPERAND`,
`S11CC1_HOMOGENEITY_CONTROL_OPERAND`, `S11CC1_HOMOGENEITY_RESIDUAL`.

**The locus protocol (inherited from S11b §8, `N8`) — for FINITE-DIMENSIONAL algebraic loci only.** Wherever §3b
asks for a genuine equation-system locus — the flat Fourier-diagonal symbol's degenerate loci and any
finite-dimensional algebraic coefficient degeneracy (`S11CC1_DEGENERATE_LOCI_*`) — emit all of: `_EQUATIONS`
(the defining system as CAS relations before any solve), `_SOLUTION` (the solution set as the CAS returns it,
variables named, domain unrestricted), `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` (each with ordered fields
`STATUS_TOKEN∈{PROVED_TRUE·PROVED_FALSE·UNDECIDED}`, `TEST_OBJECT` the live re-parseable CAS boolean,
`OPERANDS`), and `_REAL_ADMISSIBLE` (per `_SOLUTION` branch, in branch order, with
`STATUS_TOKEN∈{ADMISSIBLE·EXCLUDED·UNDECIDED}`, `TEST_OBJECT`, `OPERANDS`). ⛔ A branch is never silently
dropped; `UNDECIDED` is an explicit coverage finding. ⛔ "Solve over the reals" is not a specification. ⛔ Do
**not** force the operator's profile-conditioned singular set through this protocol — that is the Fredholm
condition `S11CC1_NONINVERTIBILITY_CONDITION` (emitted formally) and the deferred S11c-d resolvent (§3b). ⛔ A
regime classification (propagating/evanescent/grazing) is reported in whatever form it takes, not forced through
`_EQUATIONS`/`_SOLUTION`.

The three script clauses are exact obligations:

1. A script prints computed CAS objects and never states conclusions.
2. For every comparison it emits operand A, operand B, and `A−B` before any guard. A physics disagreement emits
   and continues with exit status 0; nonzero exit is reserved for operational failure only.
3. Interpretation belongs to the step record, not either engine.

Emission is never conditional on a payload's value; a boolean is a typed CAS object with its operands, never a
host-language native boolean. ⛔ No script emits a terminal judgement, `VERDICT`, `PASS`, or `FAIL`. ⛔ A growing
or decaying root, a real or imaginary impedance part, is emitted as the object the computation returns — never
suppressed, re-branched, or reclassified to reach a passive or lossless form (§1b).

---

## 7 · Names, F9 reservations, chain output, and parallel tag grammar

Every spatially varying object and every new kernel/operator has a **fresh** standard name. The varying names
reserved by S11c-a/S11c-b are inherited unchanged (`W_bg`, `w1_profile`, `L_W`, `sigma_W`, `mu_R_bg`,
`m1_profile`, `rho_4D_bg_rho4_constant`, `rho_br_bg_rho4_constant`, `rho_4D_bg_rhobr_constant`, `e_W_bg`,
`eta_bg`, the S11c-b operator/kernel/`μ_θ` names). S11c-c1's new objects — the curved DtN operator, its flat
symbol, the permeable face response, the noninvertibility condition, and any new permeability/memory-kernel
constant emitted by §3 — take fresh names and ⛔ must **not** reuse the imported constant keys `W_0`, `mu_R`,
`rho_br`, `e_W`, or `v_0` (`N14`/`G3`). `v_bulk_normal_0` remains reserved for the drain and is never aliased to
`v_0`. The importing SymPy engine applies the inherited F9 collision check against the S11b, S11c-a, and S11c-b
keys; a disagreement emits and continues.

**Chain output (`N1`/`N8`; the topology is the two-leg-gated `directives/export_ledger_bind_closure_design.md`,
§D1–§D3, which supersedes `F10` — pointed at, not restated).** The SymPy engine imports the inherited model via
`load_model` (`scripts/ledger_fold.py`) over the atomic frozen base `scripts/S11c_b_exports.py`, binding only its
declared `IMPORT_KEYS`, and writes `scripts/S11c_c1_exports.py` as its **own-rows delta** (§D2 — ⛔ **not** the
accumulated whole-model file): the flat plain-git LEDGER of the rows c1 **defines or re-derives that some later
step binds** (§D1 bind-closure membership) — the §3b consume-set (the DtN operator `dtn_operator`, its flat
symbol `dtn_flat_symbol`, the two-momentum kernel `dtn_kernel`, and the closed face response `(δp_s, J_s,
t_s)(V_s, μ_θ)`, F9c-written `s11c_c1_face_response`/`s11c_c1_face_response_coeffs`) with the **new-symbol**
closure those rows introduce (the inherited symbols they reference already live in the frozen base and are ⛔ not
re-emitted). The §4 dissipation, energy, loci/noninvertibility, and Hermitian/reactive and
regime/parity/term-origin structural-view diagnostics are **emit-only**, in the annexed `.out`, ⛔ not in the
LEDGER. Its `BUILD_INPUT_DIGESTS` pins
`{this sub-step's SymPy audit, `scripts/S11b_exports.py`, `scripts/S11c_a_exports.py`,
`scripts/S11c_b_exports.py`, this spec, and `scripts/ledger_fold.py`}` (the fold module is a shared executable
input, §D3); before publication the delta is validated by the §D3 guard —
`check_consumer(load_model(base), IMPORT_KEYS)` resolves the manifest's recursive closure, the bidirectional
smoke-test asserts the recorded lookups equal `IMPORT_KEYS`, and the minimum-mode check asserts the delta is its
own bind-closure plus named infrastructure — alongside the `D3` in-run round-trip and the `_RELATIONALS` reviver.
⛔ Never `git add -f` a big `.out`; ⛔ never annex an `*_exports.py`. The S11c-c1 comparator is this sub-step's own
frozen `T7` join (`N8`, inherited verbatim — join
by object name with the axis-typed keys of the S11c-a reconciliation schema
`steps/S11c_a_interface_shape_derivatives.md:233-253`, pair residual operands, reject a native boolean as a
residual operand, three-valued, repoint ablation). It computes and prints, deciding nothing (rule 2). ⛔ No
representational fold is pre-registered into the comparator; any cross-engine representative difference is a
computed residual **adjudicated after the run**. ⚠ A branch or regime the two engines select or key differently
is a **computed residual to adjudicate**, not a fold to pre-register (§1b: branch selection is physics).

Both engines use the exact grammar

```text
<ENGINE>_S11CC1_<QUANTITY>
```

where `<ENGINE>` is `PY` or `WL`. Each base name written as `S11CC1_<QUANTITY>` above is emitted by replacing its
leading `S11CC1` with `PY_S11CC1` or `WL_S11CC1`; do not duplicate the `S11CC1` component. Emit one tag per named
object; a single object's branch/face/DOF/density/regime-pair/parity/direction/mutation cases are a keyed CAS map
in that object's payload, not separately invented tag names. Any unavoidable engine-local tag has `_LOCAL_`
immediately after `S11CC1`, and each engine emits one local-tag inventory. The SymPy engine imports the inherited
model via the fold over the atomic frozen base `scripts/S11c_b_exports.py` (the whole F9-resolved S11b + S11c-a +
S11c-b model), binding only its declared `IMPORT_KEYS`; the Wolfram engine re-derives the supplied §§1–2 inputs
and the S11c-a substrate it consumes and the S11c-b `μ_θ` operand, importing nothing.

---

## 8 · Supplied versus computed; builder report

**Supplied and unfalsifiable here:** all of §1 (the inherited DOF and sector split; the bulk acoustics, radiation
condition, and branch object; the flat B0b/B0c reference used only as the §5c regression operand; the face
closure, affinity, kernels, traction with its `Λ_X` channel, and kinematic balance); the background ansatz,
anchoring, and curved half-space setup of §2; the standing rest-frame limit and its non-uniform validity domain
(§2b); and the S11c-a shape-derivative substrate (T-a..T-i) and the S11c-b `μ_θ` operand — imported by SymPy and
re-derived by Wolfram, per-engine-verified with the S11c-b cross-engine residual deferred.

**Computed here:** the curved-face DtN operator (the composition form and the two-momentum `k,k′` kernel, the
flat symbol, the rigid-shift residual, the regime-pair/parity/branch structure, the Hermitian and reactive parts)
(§3a); the permeable curved face response, its formal noninvertibility condition and finite-dimensional loci, and
its dissipation audit as three distinct objects — bulk-radiation Hermitian part, the two-port permeable-port
Hermitian form, and the independent traction-vs-far-field-flux energy balance (§3b); and both operands and the
residual for the representation-invariance (Eulerian vs Hanzawa second route), the mandatory `W_bg`-tilt one-sided
independence control, the `W_bg` per-direction form-ablation, the uniform-limit, the zero-jet, the branch/
momentum-liveness (both legs), and the homogeneity packages (§§5–6) — each with its computed multigrade and
dimension. (The `Σ_E`-advection and `μ_R,bg` channels are structurally absent from c1's `μ_θ`-symbolic objects
and are reserved for S11c-c2.)

The builder's report is under 25 lines and gives the files written, line and tag counts, tasks run or skipped,
runtime, all emitted tag names, and any ambiguity or non-computable requested object. It states that §§1–2 and
the S11c-a/S11c-b substrate were supplied and unfalsifiable in this build, that the rest-frame limit and its
non-uniform validity domain condition every result (§2b, grazing = strict `v_bulk_normal_0=0`), and reports no
computed value and no conclusion about the physics.
