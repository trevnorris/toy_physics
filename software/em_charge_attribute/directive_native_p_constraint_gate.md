# Directive (v6) — Stage-1 constraint-class kill-gate: does the native little-arrows `Pᵃ` sector develop an EMERGENT first-class `U(1)` Gauss constraint (that a conserved charge sources), or not?

**Make-or-break gate for the EM sector's "medium-first generative" pass** (`docs/em_medium_first_generative_plan.md` §4; unanimous 3/3 `SOUND_BUT_REFRAME`). **v6 fixes the ONE ill-posedness the v5 build hit: the "unique exact-quartic basis modulo unrestricted nonlinear field redefinitions" is ill-defined (a nonlinear point map sends quartic→degree-8), so the exhaustive classification returned `GATE_ILL_POSED` even though both engines found all constraints second-class / no `G` on the regular stratum. v6 replaces it with the standard well-posed ORDER-BY-ORDER basis + a decisive-linearization lever.** (v1–v5 review folded; physics confirmed correct since v2.) Requirement + acceptance; Codex builds the scripts, freezes+documents the flagged implementation choices (reviewable), and chooses no physics that affects the verdict. Reasoning effort: **xhigh**.

## Goal (one line)
By a **full nonlinear Dirac–Bergmann analysis**, decide whether the frozen theory develops an **ADDITIONAL first-class local Gauss constraint** (an emergent `U(1)` acting nontrivially on physical fields, sourced by a conserved charge), or does not (no native EM). The pre-existing radial pair is second-class; that is NOT the test.

## ⚠️ Central point (no tautology)
`{P²−1, P·π}=2P²≈2≠0` is second-class **for the BARE regular sigma model only**; in coupled theories re-derive tangency from the full Legendre map. Gauging adds a **NEW first-class primary→Gauss chain** (requiring a genuine kinetic-Hessian null direction — compute all rank strata; do not assert rank preservation). **The gate SEARCHES for that additional chain.**

## FROZEN fields (real 4D objects; diagonally-locked rotations — internal index IS spatial)
- **`Pᵃ`, a∈{x,y,z,w}** — real 4-vector, hard `PᵃPᵃ=1` via multiplier `λ` (retain `λ,π_λ`; full chain, do NOT pre-register `C₂`). **Easy-axis:** potential energy `V_easy=−½κ(P·ŵ)²`, `κ>0` — minimized at `P=±ŵ` (easy-AXIS, matching the vacuum), i.e. Lagrangian term `L⊃+½κ(P·ŵ)²`. Residual symmetry **`SO(3)_brane × Z₂(w→−w)`**; vacuum `P=±ŵ`. Parities: **`P_w` ODD, `P_∥` EVEN**.
- **`u_a`, a∈{x,y,z}** — brane MacCullagh shear, `w→−w` EVEN. **Transverse via a LOCAL multiplier `σ` enforcing `∇·u≈0`** (frozen realization); the longitudinal/compression mode is the separate `c_s` density sector, absent here.
- **`χ_B`** — brane scalar, `w→−w` EVEN. **THEORY-A:** independent canonical field, `½(∂_tχ_B)²−½c_χ²(∂_iχ_B)²−V(χ_B)`, `V=a_Bχ_B²(χ_B−1)²`. **THEORY-C:** `χ_B` is the composite — an **auxiliary field with NO independent kinetic term**, constrained `χ_B−|P_∥|²≈0` by multiplier `ξ` (frozen realization, not elimination); `V(χ_B)` + gradient retained (become `|P_∥|²` functionals on the constraint surface). Run the full multiplier chain.
- **Charge = EXTERNAL conserved current `j^μ` (`∂_μj^μ=0`), `w→−w` ODD** — added ONLY in the source test. (Native `±w` supply is a downstream bill; see verdicts.)

## FROZEN finite operator basis (closure rules → unique enumeration; Codex EXHIBITS the explicit list + completeness)
- **Base kinetic (frozen):** `½(∂_tP)²−½c_P²(∂_iP)²`; `½ρ_u(∂_tu)²−½μ_u(∇×u)²`; `χ_B` per theory above; `V_easy` above.
- **Cross-couplings `{g_a}` — ORDER-BY-ORDER basis (v6 fix).** Expand in a combined power counting: **at most two derivatives**, and **field-amplitude order-by-order in the fluctuations about the vacuum, starting at quadratic**. At each order, enumerate all independent operators invariant under **`SO(3)_brane × Z₂(w→−w)` × brane-parity × time-reversal × `u`-shift** (**incl. Levi-Civita `ε_{ijk}` parity-allowed terms**), `P` via `PᵃPᵃ=1`-respecting contractions, reduced **modulo IBP + PERTURBATIVE (order-by-order) field redefinitions + the holonomic constraints**. This quotient **closes at each order and IS well-defined and finite** (a perturbative redefinition `φ→φ+δφ` with `δφ` higher-order only shifts operators upward, so each order's basis is unambiguous — unlike the exact-finite quotient the v5 build correctly rejected). Multipliers `λ,σ,ξ` occur ONLY in their frozen linear constraint terms; coefficients are polynomial constants in the stability/positivity domain (Hamiltonian bounded below). Report the order at which the verdict is decided.

## What you must COMPUTE (requirement, not route), for EACH theory (A, C)
1. **Full functional (nonlinear) Dirac–Bergmann chain** from the complete Legendre map: all primaries (incl. `π_λ`, the `σ`/`∇·u`, THEORY-C `ξ`), tangency chain → secondaries, each constraint's class via **weak** brackets with **every** constraint. Reject trivial multiplier null-directions as gauge generators.
2. **Search for an ADDITIONAL first-class Gauss generator `G`:** first-class combination (weakly commuting with ALL constraints) generating a **nontrivial local `U(1)`** on physical fields, primary→secondary Gauss chain, genuine Hessian null direction — vs `{g_a}`. **⭐ Decisive-linearization lever (v6):** a first-class Gauss constraint of the FULL nonlinear theory necessarily has a **nontrivial linearization** (its leading, quadratic-order form is itself a first-class constraint of the linearized theory). Therefore **ABSENCE of any first-class Gauss at the leading (quadratic-in-fluctuation) order is DECISIVE for `NATIVE_P_NO_EMERGENT_GAUSS`** — no nonlinear gauge symmetry exists without a linear part. **PRESENCE at leading order is NOT yet decisive** (an accidental linearized gauge symmetry may not survive): confirm nonlinear persistence order-by-order before any `EM_CANDIDATE`/`TUNED`/`SYMMETRY_PROTECTED` verdict; genuine ambiguity in that nonlinear persistence → `GATE_ILL_POSED` (that branch ONLY). A per-`k` principal-symbol screen is the natural leading-order instrument.
3. **Source test (frozen):** search **source-free** (`j=0`) first. Only if `G` exists, add **one** linear coupling of conserved `j^μ` to the emergent connection `A_μ` that `G` generates; "sourcing" = `j⁰` enters the Gauss constraint (`G→∂·E−j⁰`), preservation consistency ⇒ `∂_μj^μ=0`. No arbitrary `j`-couplings.
4. **Shear-duplicate check:** compare the **reduced PHYSICAL MODE SUBSPACE** (propagating DOF after full constraint reduction) to the `pathA_36` MacCullagh transverse modes (NOT "`G`'s content" — a Gauss constraint doesn't propagate). Identical → not new EM.
5. **Generic/tuned/symmetry-protected + rank:** full functional independence + weak brackets; **`k=0`/global zero modes flagged separately** (not local first-class); a claimed gauge chain must exist on an **open regular phase-space stratum**, not only at the vacuum/singular configs. **Generic** = open set of `{g_a}`; **tuned** = measure-zero **unless an exact frozen symmetry enforces it** (`SYMMETRY_PROTECTED`).

## Implementation choices Codex FREEZES & DOCUMENTS in the report (reviewable; must be standard + not verdict-affecting)
Spatial topology + boundary conditions + smearing/test-function class + puncture treatment (the functional phase space); the explicit gauge-fixing pair + boundary prescription for control 5; the explicit global-only action for control 6; the exhibited finite operator list. State each; a choice that changes the verdict is itself a `GATE_ILL_POSED` finding.

## Able-to-fail acceptance (MANDATORY — SIX controls, frozen)
1. **Free (ungauge-fixed) Maxwell**, `A₀` RETAINED as a Lagrange multiplier (NOT Weyl/temporal gauge): MUST recover `π⁰≈0` primary → Gauss **secondary**, first-class, conserved current.
2. **Cartesian gauged hard-unit model** — a `U(1)`-gauged 2-component real unit field `(φ₁,φ₂)`, `φ²=1`, `Dφ=(∂−A)φ`, `A₀` retained: MUST detect the **MIXED** first-class(Gauss)+second-class(`φ²=1` radial) structure.
3. **Bare O(4) hard sigma, no couplings:** radial pair second-class, **no** additional Gauss.
4. **Non-conserved external current on Maxwell:** MUST fail by **inconsistent preservation** `∂ₜG+{G,H_T}∝−∂_μj^μ≠0`.
5. **Fully gauge-fixed Maxwell** (Coulomb `∇·A≈0` + its secondary = a complete second-class pair, with the documented boundary prescription): no local gauge — distinct negative.
6. **Global-`U(1)`-only model** (a complex scalar, global charge, NO `A` field): MUST classify as **no local Gauss**.
- **Per-tooth ablation:** every assert fires at its own point. **No source-grep acceptance;** computed runtime guards on the actual bracket/rank/generator output.

## Rigor / tooling
- **Dual engine (REQUIRED):** SymPy runs the Poisson-bracket Dirac algebra + classification + `G`-search; a Mathematica `.wl` **independently** re-derives the chain, classes, `G` verdict. `ENGINE_AGREE` per theory.
- Scripts under `software/em_charge_attribute/`; runners `timeout 600`; no script > 10 min; timeout = reformulate, never raise the cap.
- **Honest scope:** Stage-1 constraint classification only. A pass earns ONLY the next gates (compactness, integer charge, deconfinement, dynamical `±w`-binding — deferred).

## Output → VERDICT decision table (evaluate in THIS order)
`software/em_charge_attribute/reports/native_p_constraint_gate.md`, per theory (A, C): frozen setup + EXHIBITED basis + completeness + the documented implementation choices; full nonlinear Dirac chain (constraints + classes + matrix); the `G`-search (existence vs `{g_a}`, generic/tuned/symmetry-protected, `j`-sourced?, shear-duplicate via reduced-mode compare); ALL SIX controls; dual-engine logs (`ENGINE_AGREE`). **Branch:**
1. any control fails / engine disagreement / a leading-order first-class *candidate* whose NONLINEAR PERSISTENCE is genuinely ambiguous over the order-by-order basis → **`GATE_ILL_POSED`**. **NOT triggered by the basis quotient (now well-posed order-by-order, v6) and NOT by a decisive quadratic-order ABSENCE of first-class (→ branch 2).**
2. else **no additional first-class Gauss ANYWHERE in the admissible `{g_a}` domain** → **`NATIVE_P_NO_EMERGENT_GAUSS`**
3. else (`G` exists somewhere) — evaluate a→e in order:
   a. `G`'s reduced physical modes = MacCullagh modes → **`FIRST_CLASS_IS_SHEAR_DUPLICATE`**
   b. `j·A` symmetry-forbidden OR `j⁰` cannot enter the Gauss constraint → **`FIRST_CLASS_BUT_CHARGE_DECOUPLED`**
   c. `G` only on measure-zero `{g_a}`, not symmetry-enforced → **`FIRST_CLASS_TUNED_INVERSE_DESIGN`**
   d. `G` enforced by an exact frozen symmetry → **`FIRST_CLASS_SYMMETRY_PROTECTED`**
   e. else (generic, sourced, non-duplicate) → **`FIRST_CLASS_GENERIC_EM_CANDIDATE`** — **with the native-charge-supply bill** (the external `j` sources it, but the native `±w` label supplies no continuous current: a flagged downstream requirement, NOT ill-posed).
