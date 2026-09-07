# Question-vet — is this the RIGHT N6 reconcile question, at the retained order?

You are Codex (`gpt-5.6-sol`, xhigh), helping me **vet a question before I commit to an instrument**. ⛔ Document /
reasoning only; ⛔ do NOT modify the tree; ⛔ run no CAS. Working dir `/var/projects/toy_physics`; paths under
`research/pde_ledger_v3/` unless absolute. This is the **E1 question-vet**: I own the reconcile question and its
adjudication; the reconcile INSTRUMENT will be astra-authored and G1-reviewed afterward. Your job here is only to
confirm the **question is the right one at the retained order — not a convenient proxy** — and to correct the
justified-channel physics if I have it wrong. ⛔ Do NOT tell me the answer to the reconcile (there is no expected
remainder value); tell me whether the QUESTION and its operationalization are correct.

## The object (grounded — read these, ⛔ do not trust my paraphrase)
- Spec: `directives/S11c_c2_SHARED_PHYSICS.md` §5c (≈303-395), plus the cleared route-2 construction
  `_measurements/S11c_c2_N6_route2_spec_astra.md`.
- Diagnostic: `scripts/S11c_c2_N6_diagnostic_sympy.py`. Load-bearing lines:
  - `run` builds the two routes per `(α,ρ)` case: `E,_=build_increment(...e_coeff,es...)` (:850),
    `M,_=build_increment(...m_coeff,ms...)` (:851), `REP_INVARIANCE_RESIDUAL = residual(E,M)` (:853).
  - `residual(a,b)={k:a[k]−b[k]}` (:576) ⇒ **`R_N6 = E − M`**, element-wise over `(block,grade,signature,face)` keys.
  - `build_increment` (:476-513) sums products `carrier_grade × source_grade × kernel_grade × measure × amplitude`
    — bilinear in (carrier `coeffs`, `sources`) modulo `sp.cancel(coefficient/eps)` (:479) and grade truncation.
  - Route inputs differ in THREE places: **carrier** (`e_coeff` = imported Eulerian `inputs.slab[α,ρ]` (:796-797)
    vs `m_coeff` = native MATERIAL face-factory (:806,:809)); **μ** (`mu_e` (:799) vs `mu_m=material.subs(t,1)`
    (:800)); **face velocity/normal** (Eulerian `ev=DELTA_W` (:841) vs material `m_v` (:806)).
  - `constitutive` (:318-360): `mu_e=[variation(monomial)]` (:348, no pullback); `mu_m=[variation(pulled)]` where
    `pulled=material_pullback(monomial.subs(theta→theta+shift))`, `shift=(t-1)·adv`, `adv=u·∇ρ₄/ρ₄` (:328,:334,:342).
    At `t=1`: `shift=0` ⇒ `mu_m` = variation of the material_pullback with no θ-shift. The N4 advection tag `t`
    toggles the θ-advection (`t=0` ⇒ θ→θ−adv).
  - PIT emit (:769): each object → `columns` (the keys) + `numerator_denominator` (per-sample `[num,den]` mod p) +
    `nonzero_modular_numerator` (per column: certified-nonzero iff any sampled numerator ≠ 0). ALL objects in a case
    share identical sample points (:758-759), so cross-object linear relations can be PIT-tested over shared samples.
- **The finding (⛔ NOT pre-adjudicated):** uncorrupted `R_N6` is **certified-nonzero in ~18 forward-block columns**
  (`THETA`/`E_W`/`TRANSVERSE_TO_THICKNESS`), both densities, `LAB_HELD`; reverse blocks zero; controls bite.

## The physics claim N6 asserts
E (Eulerian route: Eulerian carrier + Eulerian μ) and M (material route: native-material carrier + material μ, its
builders already covector-mapped to the common Eulerian face basis) are the **same self-energy increment in two
coordinate representations**. §5c: "N6 is the physics requirement that these two are the same operator in two
representations"; the two routes deliberately bind DIFFERENT μ (Eulerian vs material) so the residual carries the N4
constitutive channel — "N4 enters through the closed pressure source μ." The sanctioned field redefinition is
`Δρ = δρ_E + u·∇ρ⁰`, relating two descriptions of ONE perturbation of a fixed background, **at FIXED anchoring**
(⛔ never bridging the two anchorings — that is the N4 category error §5c forbids).

## My candidate reconcile question + method (VET THIS — correct it if wrong)
**Question:** Does `R_N6 = E − M`, at the retained order `(η^{≤1}, σ_W^{≤1})`, reduce EXACTLY to the justified
material↔Eulerian representation content, with a **certified-zero remainder** (⇒ representation invariance holds; the
nonzero residual is the expected representation content) — or does an **unexplained remainder** survive a
certified-nonzero column (⇒ a genuine N6 failure to adjudicate)?

**Candidate method:** By `build_increment` bilinearity, `R_N6 = Increment(e_coeff,es) − Increment(m_coeff,ms)`
splits into a **carrier channel** `Increment(e_coeff−m_coeff, es)` and a **source channel**
`Increment(m_coeff, es−ms)`. Build the *justified* channel `J` by an INDEPENDENT route (⛔ never by subtracting from
`R_N6`), form `remainder = R_N6 − J` as a new PIT object over the SAME samples, and test `nonzero_modular_numerator`
per column at the retained order. Certified-zero ⇒ invariance holds; certified-nonzero ⇒ real finding.

## The sub-questions I actually need vetted
1. **Operationalization.** Is "`R_N6` reduces to the justified channel with certified-zero remainder" the correct
   operationalization of N6 representation invariance at the retained order, or a proxy that could read as invariance
   when it is not (or vice versa)? If the correct statement is different, state it.
2. **What is `J`?** Is the WHOLE material_pullback μ difference (`mu_m − mu_e`) justified representation content, or
   is only the `u·∇ρ⁰` advection piece sanctioned by `Δρ=δρ_E+u·∇ρ⁰` justified — so that any material_pullback
   content BEYOND `u·∇ρ⁰` would itself be an unexplained remainder (a real finding), not part of `J`? Which is it,
   and why? This is the crux — getting `J` wrong makes the instrument a proxy.
3. **Carrier channel.** Should `Increment(e_coeff−m_coeff, es)` vanish independently (the material factory's covector
   map reconstructs the imported Eulerian carrier), or does the carrier legitimately carry justified representation
   content that belongs in `J`? (Note the diagnostic emits `CARRIER_RECONSTRUCTION_RESIDUAL` = imported-Eulerian vs
   Eulerian-*factory*, :836 — an Eulerian-vs-Eulerian reconstruction check, NOT Eulerian-vs-material.)
4. **Decomposition legitimacy.** Is the bilinear carrier/source split exact at the retained order given
   `sp.cancel(coefficient/eps)` (:479) and grade truncation, or does cancel/truncation break additivity so the
   reconcile must build `J` a different way (e.g. a single full material→Eulerian transform of route 1)?
5. **Defining relation.** Per the reconcile method (`[[feedback_reconcile_representational_bridge]]`: "collapse to
   the DEFINING relation = strong"), what is the defining relation whose collapse would be the STRONGEST evidence of
   invariance here — and is a certified-zero PIT remainder that relation, or merely necessary-not-sufficient?
6. **Leakage / iterate-to-zero.** Any way this question or method leaks an expected value or would let an instrument
   "iterate to zero"? (There is no supplied remainder target; the disposition is adjudicated on my side.)

## Output
For each of 1-6: is my framing right, and if not, the corrected statement. End with the **corrected reconcile
question** (one paragraph) I should hand to the instrument's decision legs, or **QUESTION SOUND AS STATED** if no
change. Brief, evidence-first, grounded in the cited lines.
