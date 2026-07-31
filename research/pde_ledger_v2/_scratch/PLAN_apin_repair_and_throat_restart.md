# PLAN — retire the `a`-pin, then restart on the throat interior

> ⛔⛔ **HISTORY — NOT INSTRUCTION.** This is the step-① plan, completed 2026-07-31 (`407eed94`). Its "OPEN
> pending a user decision" language about `R2.a_pin` describes a question that **no longer exists** —
> the relation was retired by removal. ⛔ Do not execute anything here.
> ⇒ **Live workstream: `research/pde_ledger_v3/` (start at `NEXT_SESSION.md`).**

**Status: DRAFT FOR REVIEW. Nothing has been edited in the repository.**
Author: orchestrating session, 2026-07-30/31. Branch `ledger-v2-rebuild` @ `978bf0a2`, tree clean.

---

## §0. What I am asking you for

Three independent reviewers are reading this (Codex, Grok, GLM). I want **independent judgement**, not
confirmation. In particular:

1. **Is the physics reasoning in §6 correct?** It is the load-bearing new claim and I derived it quickly.
2. **Is step ① (§3) correctly scoped** — does removing the pin break something I have not listed?
3. **Is step ② (§7) the right restart**, given §5? I have already been wrong once this session about
   what the corpus contains (§9). Assume I am wrong again somewhere.
4. **Where does this plan rebuild apparatus instead of doing physics?** That failure has happened twice
   in this project and is the single thing most likely to kill it again.

⛔ Do not be agreeable. A finding that breaks the plan is the most valuable output. If the model itself
is falsified by something here, say so plainly — that is a first-class result in this project, not a
problem to be softened.

---

## §1. The model, for context

One compressible superfluid — a Gross–Pitaevskii / nonlinear-Schrödinger (GNLS) condensate in a **4D**
bulk. Order parameter `ψ`, number density `ρ = |ψ|²`, phase `θ`, flow `v = (ħ/m)∇θ`, closed by a stiff
polytropic EOS `P = Kρ⁵`.

That one medium exists in two material states, distinguished by an order field `χ_B = r_B e^{iθ_B}` with
`r_B = |χ_B| ∈ [0,1]`: **brane** (`r_B=1`, ordered, shear-supporting — our 3D space, a domain wall at
`w=0`) and **bulk** (`r_B=0`, same medium, de-structured, shear-free).

A **particle is a throat**: a puncture of the brane whose wall structure extends into the bulk as a tube
(the **sleeve**), mouth radius `a`, `w`-extent `L`, on one side (`s = ±1`, which is the charge sign —
a Z₂ orientation label, explicitly **not** a winding). Bundled with it is a trapped nonlinear core (the
**geon**) that constitutes the mass, and a bound flow field.

The four forces are modes of that one medium: **gravity** = the steady drain flow between throats;
**light** = in-plane shear of the brane; **charge** = the static ±w throat; **magnetism** = the moving
±w throat.

**Goal of the project:** a mathematical *bridge* between the EM and gravity formalisms — an analog, not
an ontology claim. Judged by calibrate-then-predict surplus on **dimensionless** held-out ratios.

---

## §2. What this session established

The session began as a repair task and turned into a foundations question. Sequence of findings:

### 2.1 The `a`-pin is a units artifact, not physics

`stage004:124-130` imposes **four** natural-unit pins `{a, c_s0, ħ, m_GNLS}` on **three** base dimensions
`{L,T,M}`, reports "rank 3, nullity 1", and reads off the surviving monomial as

```
a·c_s0·m_GNLS/ħ = 1    i.e.    a = ħ/(m_GNLS·c_s0)
```

**Claim (mine, needs checking — see §8 Q1).** The remaining three pins already form a complete,
non-degenerate unit system:

```
c_s0   = L T⁻¹    → ( 1, −1, 0)
ħ      = M L²T⁻¹  → ( 2, −1, 1)      det = 1  ⇒  rank 3, nullity 0
m_GNLS = M        → ( 0,  0, 1)
```

⇒ The fourth pin was **redundant**, and *any* fourth pin on three dimensions must produce exactly one
dimensionless residue. Because `a` is a pure length, the residue it produces is forced to be
`a·m·c_s0/ħ`. So `a = ħ/(m c_s0)` carries **zero information about throats** — it is the unique monomial
that had to exist once the unit system was over-specified.

### 2.2 The damage this caused, triaged

Triage rule used (from the governing decision doc): a mention is not damage. Damage is only
**(a)** the pin presented as determining a physical radius, **(b)** the pin conflated with the throat
radius, **(c)** their dimensional coincidence recorded as evidence of consistency.

| file | verdict |
|---|---|
| `stage004` note `:134`, `:197` | **DAMAGED (a) ×2** — `:134` calls the pin *"a mouth-radius collective moment"*; `:197` says the *"pin geometry `a`"* **follows** from the primitives alongside `c_s0` and `ξ_h` |
| `stage004` `.py` / `.wl` | ✅ clean — both print *"a is a branch collective moment, **not** a base invariant"* |
| `parameter_register.md:132` (the pin row) | ✅ benign and correct — `CONV`, *"never free"* |
| `parameter_register.md` — **13 loci** | ⭐ **DAMAGED (b)** — see §2.3 |
| `stage016:184` | **DAMAGED (b)** — *"`a` is the same throat radius as stage018's — a physical-radius identity"*. `:181` immediately above is **correct**: *"an identical dimension, so no dimensional check here can support or refute"* |
| `stage023:481`, `:407` | **DAMAGED (b)** and **(c)**. `:407` was not on the known list |
| `STATUS.md` | ✅ **already repaired** — `:40` now reads *"a DIMENSIONAL COINCIDENCE, not agreement"*. Only its own metadata is stale (still calls the repair "step ②") |
| `stage005` ×2, `docs/pathA_preregistration.md` | ✅ clean **and correct** — they carry the anti-tautology caveat. Cite, do not touch |
| `stage043` `.py` / `.wl` | ✅ zero damaged occurrences |

### 2.3 ⭐ The load-bearing form of the damage — it understates the count

The damage is **not** mainly narrative. In 13 places `parameter_register.md` uses `a` as a genuine
physical radius inside a physics formula and tags it `CONV` (*"`a` is the `CONV` pin"*), so that it is
**not counted** as a knob:

| locus | the formula `a` appears in |
|---|---|
| `:314` | `K̄₀a⁵/(27c_s⁵) − 2G/(5c⁵) = 0` — the **2.5PN Burke–Thorne matchback** |
| `:306` | `54Gc_s⁵/(5a⁵c⁵)` — the ℓ=2 quadrupole magnitude |
| `:308` | `z = aω/c_s` — the DtN response at the throat boundary |
| `:310` | `g_base ∝ a^(−7/2)`, `N0_den ∝ a⁻⁷` |
| `:186` | `ϖ_Φ2 = (c_s/a)²·(6+(ma)²)` — bulk-Helmholtz geometry |
| `:200` | `ℓ/a = 1/20` — the embedding/PT length ratio |
| `:629, :669, :710, :731, :745` | `a` cited as `CONV` **provenance** for stage-level claims |

Several sit directly under headline count claims — *"ZERO new counted knobs"* (`:308`), *"Part-II counted
CALIB set UNCHANGED"* (`:710`, `:731`). A **calibrated physical input is being excused as a units
convention**, so the irreducible count is understated.

### 2.4 The acceptance check cannot see the question it is being asked

`stage043`'s medium `AlgebraBlock` (`:606-615`) carries 10 symbols and 5 constraint polynomials, one of
which is `scale_a*mass*cs0 − hbar` — **the pin relation, used as a content constraint**. That block is
identical to the v2 reduction registry's medium block (ambient 10, `dim_after` 5, `delta_r` 5), which is
exactly `acceptance_check.py`'s `EXPECTED_MEDIUM_PAYLOAD["baseline"]`.

⇒ **Both sides of the acceptance comparison admit the pin as a constraint.** They agree because they
share the assumption under question. The `MATCH` therefore carries **no information** about the pin's
classification. It tests the rank machinery, not the physics.

⚠ `stage043` also tags the same pin `REG:conv:a_pin → CAT_CONVENTION` in its register (`:334-337`) while
consuming its equation as content 270 lines later. **Classified as convention, used as content, in one
artifact.**

### 2.5 The pin is a dangling leaf

`Q.medium.a_pin` appears in **no** relation's `input_qids` — only as `R2.a_pin`'s own `designated_output`
and in its own positivity guard. Nothing in the medium block consumes it.

Also: `R2.xi_h` and `R2.a_pin` have **identical inputs** `[ħ, mass, c_s0]` and differ only by `√2`. That
is the source of the `ξ_h² − 2a² ≡ 0` entailment found earlier — not a physics coincidence.

### 2.6 ⭐ Why the `ξ_h` identification is a category error

`ξ_h = √2·ħ/(m_GNLS·c_s0)` is built **only** from medium primitives — it is **one number for the whole
medium**. The throat radius is a property of **a defect**, and there are at least three leptons. A single
medium constant cannot equal a species-indexed family. `a = ξ_h/√2` can hold for at most one lepton, by
accident of units.

⇒ This is stronger than "two lengths agreeing on being lengths". It is a **medium constant equated with a
per-particle quantity**.

⭐ **Why it happened is instructive.** In a *standard* superfluid the topological defect is a **quantized
vortex**, and a vortex core size **is** the healing length — universal across all vortices. That is almost
certainly the intuition behind the identification. But this model's defect is explicitly **not** a
winding (charge is a Z₂ orientation), so the model is not entitled to the standard result.

---

## §3. STEP ① — retire the pin (USER DECISION, already made)

**Decision (user, this session):** remove all instances of `a_pin`; `a` means the **throat radius**,
everywhere, one quantity. Anywhere that used the pin part must be re-thought. The throat radius is
**TBD** — calibrated/undetermined for now, its derivation being the program's goal.

⇒ Note this **dissolves** rather than answers the previously-open question of `R2.a_pin`'s registry class.
The relation is deleted, not reclassified.

### 3.1 The edit list

| # | file | change |
|---|---|---|
| 1 | `stage004` note §4 (`:122-137`) | Replace the four-pin imposition with the three-pin system `{c_s0, ħ, m_GNLS}` (rank 3, nullity 0, complete). **Delete** the relation `a = ħ/(m c_s0)`. `a` is no longer defined here at all. ⚠ **Line-count-preserving** or update loci — see 3.2 |
| 2 | `stage004` note `:197` | `a` no longer "follows" from the primitives; remove it from that list |
| 3 | `stage004` `.py` / `.wl` | pin-relation derivation + `A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT` emission removed; the pin-corruption ablation probe retired with it |
| 4 | `parameter_register.md:132` | the `a` (pin) row becomes `a` = **throat mouth radius**, class `CALIB`/TBD, **no defining equation** |
| 5 | `parameter_register.md` ×13 | `CONV` → the physical radius, `CALIB`, **counted**. Affected headline claims re-opened (see 3.3) |
| 6 | `stage016:184`; `stage023:407`, `:481` | the conflation statements removed; `stage016:181`'s correct caveat kept |
| 7 | `reduction/relations.yaml` | **delete** `R2.a_pin` |
| 8 | `reduction/quantities.yaml` | `Q.medium.a_pin` → the throat radius: no defining equation, **must-supply**, `scope` no longer `branch-pin` |
| 9 | `stage043` `.py` `:611` / `.wl` `:543` | remove `scale_a*mass*cs0 − hbar` from the medium `AlgebraBlock`; retag `REG:conv:a_pin` |
| 10 | `acceptance_check.py` | all four numbers **re-established by independent derivation**, ⛔ never preserved |
| 11 | `STATUS.md` | stale repair metadata only (it calls the repair "step ②"; the archive step is done) |

### 3.2 ⛔ A hazard the repair must respect

The registry pins loci **by line number** into both prose files being edited — **16 ranges** into
`stage004` (`65 | 76-83 | 85-92 | 102-120 | 122-137`) and **10** into `parameter_register.md`. A repair
that shifts lines silently invalidates them, in a corpus where *"zero cross-artifact citations resolve to
a locus"* is already the measured defect. ⇒ Edit **in place, line-count-preserving**, or update every
locus in the same commit.

### 3.3 What the numbers do

- One fewer admitted constraint ⇒ baseline moves from `dim_before 10, dim_after 5, Δ 5` to
  `dim_before 10, dim_after 6, Δ 4`. All four medium cases change.
- The **sim-input set grows by one**: the throat radius becomes a quantity you must supply. This is the
  **anti-flattering** direction — the model looks *less* determined — which is the direction to trust.
- `ξ_h² = 2a²` stops being an identity. It becomes a *claim*, and per §2.6 a false one.

### 3.4 ⛔ Scope escalation the reviewers must check

Removing the pin is **not contained to the v2 ledger**. Computational consumers of `scale_a`/`a_pin`:

- `research/pde_ledger_v2/scripts/ledger_stage028_2_5pn_matchback_sympy_audit.py` (+ `.wl`)
- `software/stage1_solver/tools/pathA_2_5pn_matchback_sympy.py` (+ `.wl`)
- `research/4d_1pn_bridge/scripts/throat_kappaPV_closure.py`
- `research/pde_ledger_v2/scripts/midway_knob_audit_codimension_sympy.py`
- `research/pde_ledger_v2/reduction/{acceptance_check.py, able_to_fail.py}`
- `software/em_charge_attribute/u1_*`, `u2_*` (~8 files) — **deferred**, that workstream's charge model
  changed recently

⭐ **The first three are the gravity ladder.** `INV2: K̄₀a⁵/(27c_s⁵) − 2G/(5c⁵) = 0` consumes `a`, and
`parameter_register:314` currently tags that `a` as `CONV`. **Whether that `a` is the pin (≡1,
contentless) or the physical throat radius changes what the GR match means.** This is a physics question
inside the repair, not bookkeeping. → §8 Q2.

---

## §4. Why this mattered more than a symbol clash

`a` is the anchor for the sleeve's other lengths — `L/a = 37/20`, `ℓ/a = 1/20` — and `a` itself was
anchored to nothing. Every geometric fact about the particle is a **frozen ratio to an unanchored
length**. The pin was someone terminating that chain by reaching for the one length that had a formula.

Inventory of lengths in the throat region and what anchors each:

| length | anchored to | real? |
|---|---|---|
| `ξ_h` healing length | core balance from `ħ, m, c_s0` — **derived** | ✅ the distance the condensate heals over |
| `δ`, `σ_wall` (wall) | `δ = √(κ_B/2a_B)`, `σ_wall = √(2a_Bκ_B)/6` — **derived**, but from two **postulated** constants `{a_B, κ_B}` | ⚠ real, but rests on new postulates |
| `W_slab` | ⛔ **nothing** — `FREE-UNREDUCED`, *"double-well selects NO width"* | the width the model actually needs |
| `a` mouth radius | ⛔ **nothing** | real object, no equation |
| `L`, `ℓ` | frozen ratios to `a` | only as real as `a` |
| `r_e` | external electrostatic estimate | ⛔ contradicted — see §6 |

---

## §5. What the corpus actually contains on the throat interior

Three independent surveys were run this session. Summary of what exists, because the plan turns on it.

### 5.1 The brane / wall — attempted, and it hit a no-go

- ✅ **The wall potential IS written down with coefficients.**
  `software/em_charge_attribute/g0_closure_card_v0.md:115-124`:
  `S_χ = ∫ ½Z_χ(∂_t r_B)² − ½κ_χ|∇₄r_B|² − (λ_χ/4)r_B²(1−r_B)²`, with `λ_χ = 2κ_χ/ℓ²` giving a planar
  wall of logistic width `ℓ`. Star values at `:131` (`V_χ = (800/4)r_B²(1−r_B)²`, `[POSTULATE]`).
- ✅ **The kink is solved.** `ledger_stage006_two_phase_chiB_ontology.md:136-141` — Part-I twin
  `f_B = a_B χ_B²(1−χ_B)²`, EL kink solved, *"width SOLVED from the EL equation"*; `δ`, `σ_wall` derived
  (`parameter_register.md:287`, R20 `DERIVED`).
- ⛔ **But it does not select the width the model needs.** `parameter_register.md:163,286` — `W_slab` is
  `FREE-UNREDUCED`, *"double-well selects NO width"*, R19 **PENDING**. The kink gives one interface; the
  brane is a finite slab.
- ⛔ **The medium program returned a NO-GO.** The prior-art survey chose the **GNLS polar-smectic**
  (`medium_requirements_and_prior_art.md:128,143`) and named Gate L *"Highest risk; the most likely
  no-go"* (`:172-183`). Gate L ran: **`FAIL_COUPLE_STRESS_NOGO`**
  (`stage030_pathA35_gateL_source_map.md:16`). `conceptual_history.md:340` records the program as *"now
  SUPERSEDED at the brane-existence level."*
- ⚠ **Another unearned identification.** `ledger_stage044_parent_action_reconciliation.md:104`:
  `ℓ = √(2κ_χ/λ_χ) = √(κ_χ/2a_B) = δ` — the closure-card length and the ledger's kink width are
  **identified, not independently derived**. Same pattern as the `a`-pin.
- ⚠ `model_map.md:26` justifies the postulate as *"the parent potential `U(ρ)` is single-well"*. That
  reasoning is **wrong** — the wall comes from `V_χ(r_B)`, a potential on a *different field*. The
  conclusion (postulated) survives for a different reason: `V_χ` is itself postulated with two free
  constants and selects no slab width.

### 5.2 The geon / mass mechanism — a proxy exists, the geon does not

- ⛔ **"Geon" is always an open declaration.** `two_throat_simulation_handoff_spec.md:219` — *"its profile
  is a declared OPEN input [POSTULATE: mass mechanism; profile OPEN]"*. The string "geon" appears nowhere
  in `notes/`.
- ✅ **But its proxy — the "trapped support wave" — has explicit energy functionals.**
  `notes/lepton_mass_notes.md:50-58` gives a reduced rest-energy ansatz:
  ```
  F(a,ρ) = A(ρ)/a + B(ρ)/a² + C(ρ)a³
      A/a  — trapped support-wave ledger  E_w
      B/a² — self-flow / feed ledger      E_f
      C a³ — pressure-volume ledger       E_PV
  ```
  and `:110` — `F_* = 18u = (18/11)E_w`, so the isolated particle mass is `(18/11)` times the support-wave
  ledger on the same branch.
- ✅ **A self-selection balance exists and was SOLVED.**
  `software/stage1_solver/reports/pathA_26_derrick.md:19` —
  `E*(a,L) = Ω_fluid,excess + I·ω(a,L) + P_vac·V + σ·A`, with the fixed-action wave branch
  `ω = √(c_w²(π²/a² + (π/2)²/L²) + μ_in²)`.
  `:43` — sharp-wall stationary point **`a* = 1.8924`, `L* = 1.8121`, `L/a = 0.9575`**, Hessian
  eigenvalues `[0.311, 1.616]` (positive), virial residuals exactly `0`.
  Also `notes/inner_throat/inner_throat_4d.md:878,1204` — `∂_a(E_fluid+U_wave+E_geom)=0`,
  `∂_a H_tot = 0, ∂_L H_tot = 0`.

### 5.3 ⛔⛔ Three negative results that are not in the front-door docs

1. **The drain destabilizes the throat.** `pathA_26_derrick.md:58-72` — Phase C, with the drain slaved in:
   *"Robust stable region exists: `False`"* (`THROAT_DRAIN_DESTABILIZED`). The failure is **static
   divergence, not flutter**: at the required gain `g_open = 88.83` (vs threshold
   `gcrit_plus = 0.00131`) the stiffness has negative determinant and a negative symmetric-part
   eigenvalue, and *"passive damping `C≥0` cannot remove the negative-stiffness direction"*. Stable fixed
   points exist *"only in a tiny near-zero-drain corner."*
2. **The lepton mass tower is falsified.** `conceptual_foundation.md:487-488` — *"the trapped-wave **mass
   TOWER is falsified** — it predicts mass ratios 1:9:25, reality is 1:207:3477 — and the absolute mass
   scale (the wave amplitude) is undetermined. So this gives **one** soliton's mass, **not** a predictive
   lepton spectrum."*
3. **The solved `L/a` contradicts the frozen one.** `pathA_26_derrick.md:43` solves `L*/a* = 0.9575`.
   Everything downstream uses the frozen ansatz `L/a = 37/20 = 1.85` (`STAGE1_VERDICT.md:13`,
   `two_throat_simulation_handoff_spec.md:73`). **A factor of ~1.9 disagreement between the solved value
   and the value in use.**

Also: `STAGE1_VERDICT.md:13,18` — GATE A ran with `a=1, L=37/20` **frozen**, and returned
`R_norm = −10.7999993`, ~7 orders of magnitude below the `54/5` target — a robust MISS.
`conceptual_history.md:138,142` — *"a frozen input — self-selection explicitly forbidden"*, a fold at
`τ≈0.029`, and *"nothing in the built models self-selects the throat size/shape"*; the candidate fix (a
4D Derrick/virial check) is marked **UNTESTED / never run**.

### 5.4 Gravity — and the good news

- **Part II is CLOSED** (ledger stages 001–002, 008–029), every stage dual-engine verified (SymPy +
  independent Mathematica, both exit 0) + per-tooth ablation + review legs
  (`STAGE_PROVENANCE_INDEX.md:36-38`).
- ⭐ **The PN ladder is explicitly independent of the throat interior.** `research/4d_2_5pn/paper/`
  `4d_2_5pn.tex:613-614` — the dynamics are *"controlled by the worldline/worldtube multipoles rather
  than by arbitrary internal details of the defect"*; `4d_1pn_full.tex:110,145` — *"we do not attempt to
  solve the fully dynamical moving-throat problem in the bulk"*. Earned target-blind: `1/r²` + attractive
  sign, the ℓ=2 DtN fingerprint `{1/9, 4/81, 1/27}`, `χ_Q=+1`, cross-ℓ `{1, 1/2, 1/27}`, SO(3) `λ_m=6`.
  A live falsifier exists (stage009 `RETURN_RESIDUAL_PREDICTION`, ℓ=0/1 return, which GR forbids).
  ⇒ **This survives the restart untouched.**
- **Blocked on the throat interior:** `{μ_R, ρ_br}` (R10, hence `c_γ` and every cone lock), the
  frozen-wall packet `{μ_η, T_w, β}` (R30), the breathing drive `{Vp0/ℓ_c}` (R33), and the density-port
  magnitudes (SIM-deferred).
- ⛔ **`m_defect` is a GAP.** `stage004:181` — `INFLOW_MASS_SOURCE_MISSING` / `BLOCKS_MASS_EMERGENCE`. The
  only relation is a **dimensional bridge** `m_defect = α_J ħJ/c_γ²` with `α_J` a `CANDIDATE`, *"labeled
  candidate; not load-bearing"*. pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED` — *"No action-level,
  boundary-source, Noether, or Hamiltonian equation"*.
- ⭐ `pathA_19_dimensional_foundation.md:32` already names **`J` the better invariant label** over the
  mouth radius `a` — i.e. the corpus had already concluded `a` was the wrong handle.

---

## §6. ⭐⭐ The central physics question — does a heavier lepton have a bigger or smaller throat?

This is the question that broke the session open, and the corpus contains **two opposite answers**.

### 6.1 The external estimate says smaller

`two_throat_simulation_handoff_spec.md:222` — *"the throat's nonzero size may be identified through the
self-energy balance `k e²/a = m_e c²`, giving `a ∼ k e²/(m_e c²) = r_e`"*, tagged `[CALIBRATED size
comparison]` and explicitly de-rated: *"never an interior-PDE or force-sign input"*.

All three leptons carry the **same** charge `e`, so this gives `a ∝ 1/m` — heavier ⇒ **smaller**, spanning
`a_e : a_μ : a_τ = 1 : 1/207 : 1/3477`.

⛔ **I believe this estimate does not apply to this model**, for two reasons:
1. It is the classical electron radius, which assumes the entire rest mass is electrostatic self-energy —
   known to be false even in standard physics (`r_e` is not the electron's size).
2. In *this* model charge is the **Z₂ ±w puncture orientation**, identical across leptons, while mass is
   drain / standing-wave energy. Same charge, different mass ⇒ electrostatics cannot be what sets the
   size.

### 6.2 ⭐ The model's own reduced energy ledger says bigger

Take the corpus's own ansatz (`lepton_mass_notes.md:50-58`), `F(a) = A/a + B/a² + Ca³`, and impose
stationarity:

```
dF/da = −A/a² − 2B/a³ + 3Ca² = 0
```

`F → ∞` as `a → 0` (the `A/a` and `B/a²` terms) and as `a → ∞` (the `Ca³` term), so a **finite minimum
exists** — the size is genuinely self-selected. Dropping the sub-leading `B` term for scaling:

```
3Ca² ≈ A/a²   ⇒   a⁴ ≈ A/(3C)   ⇒   a ∝ A^{1/4}
m = (18/11)E_w = (18/11)(A/a)  ⇒  m ∝ A/A^{1/4} = A^{3/4}
```

⇒ `A ∝ m^{4/3}` and therefore

```
                    a ∝ m^{1/3}
```

**Heavier lepton ⇒ BIGGER throat.** And `m^{1/3}` is exactly the Q-ball / liquid-drop scaling
`R ∝ Q^{1/3}`, which is what `conceptual_foundation.md:484` says the object is (*"a self-bound
Q-ball/soliton"*).

### 6.3 What this means

- The two routes disagree not in magnitude but **in the sign of the slope**. They cannot both be right.
- The model's own equations back the "bigger" answer. The `r_e` estimate is the outlier, and it is
  already de-rated in the very line that states it.
- ⚠ **Caution.** `pathA_26_derrick.md:21` uses a *different* wave branch —
  `ω = √(c_w²(π²/a² + (π/2)²/L²) + μ_in²)`, a cavity mode, which for small `a` gives `ω ∝ 1/a` and hence
  a Compton-like `E ∝ 1/a`. **The two reduced models in the corpus are not obviously the same model.**
  Whether they agree is an open question and is exactly the sort of thing that produced the `a`-pin.
- ⚠ The `1:9:25` of the falsified tower is a **harmonic** pattern (`ω ∝ n`), not this balance. So the
  falsification of the tower does **not** automatically falsify §6.2 — different mechanism. But it does
  show this ground has been walked and produced a negative result.

⇒ **This is the single most checkable physics claim in the plan. → §8 Q3.**

---

## §7. STEP ② — the proposed restart

⛔ **Not** "go write the wall potential" — that exists (§5.1). The honest statement of where the model is:

> The medium is standard and solid. The **particle** is a postulate whose mass mechanism is undeclared,
> whose size is set by a balance that has been solved once and found **unstable under its own drain**,
> whose mass tower is **falsified**, and whose geometry is a chain of frozen ratios to an unanchored
> length. The brane it punctures rests on a postulated double-well with two free constants that selects
> no slab width, and the medium program that would have grounded it returned a **no-go**.

That is the real front. Proposed order, each step stopping for the user:

**R-0. Bank step ①.** Retire the pin. Independent of everything below.

**R-1. Reconcile the two reduced throat models.** `lepton_mass_notes.md`'s `F = A/a + B/a² + Ca³` versus
`pathA_26_derrick.md`'s `E*(a,L) = Ω_fluid + Iω(a,L) + P_vac V + σA`. Are they the same balance in
different variables, or two different physical pictures? Settle §6.2 vs §6.3's caution. **This is cheap
and it is a physics step, not apparatus.**

**R-2. Resolve the `L/a` contradiction.** Solved `0.9575` vs frozen `37/20 = 1.85` (§5.3.3). Everything
geometric downstream uses the frozen value. Either the solve is wrong, the ansatz is wrong, or they are
different quantities — and the third possibility is the `a`-pin pattern again.

**R-3. Face `THROAT_DRAIN_DESTABILIZED` directly.** `pathA_26` says the object is unstable at the drain
strength the model requires, by static divergence, and damping cannot fix it. This is either (i) a real
falsification of the throat-as-modeled, (ii) an artifact of the sharp-wall / slaved-drain reduction, or
(iii) evidence the stabilizing physics is in the sector that was frozen out. ⭐ **A model whose particles
are unstable is a bigger problem than an unanchored length**, and this result is not in `model_map.md`
or `STATUS.md`.

**R-4. Then, and only then, the mass–radius relation `R*(E)`** and the held-out lepton test.

### 7.1 The held-out test this sets up

If `R*(E)` yields `a ∝ m^{1/3}`, calibrating `a` on the electron **predicts** `a_μ` and `a_τ`. It is a
**dimensionless** prediction (`a_e : a_μ : a_τ`), which is the only kind that can test this model, and it
fails loudly — the alternatives span 3.5 orders of magnitude in the opposite direction.

### 7.2 What is explicitly NOT being redone

The bulk GPE/Madelung (textbook); the PN gravity ladder (audited, GR-matched, and **worldtube-reduced so
it does not depend on the interior** — §5.4); the target-blind sector structures; the dimensional
foundation minus the pin. The EM/charge cluster is **deferred** — that workstream's charge model changed
recently and is the wrong place to spend attention now.

---

## §8. Questions I want answered

**Q1 — the unit-system claim (§2.1).** Is `{c_s0, ħ, m_GNLS}` on `{L,T,M}` rank 3 / nullity 0? Verify the
determinant independently. If it is rank 3, does it follow that the fourth pin's residue is contentless,
and is my stronger claim correct — that *any* fourth pin of dimension `L` yields `a·m·c_s0/ħ` and nothing
else? ⛔ Do not take my arithmetic on trust.

**Q2 — the gravity consumers (§3.4).** In `ledger_stage028_2_5pn_matchback` and
`software/stage1_solver/tools/pathA_2_5pn_matchback`, is the `a` in `INV2: K̄₀a⁵/(27c_s⁵) − 2G/(5c⁵) = 0`
being used as the **pin** (≡1, contentless) or as the **physical throat radius**? What happens to
`MATCHBACK_CONSISTENT` under each reading? This decides whether step ① touches an audited GR match.

**Q3 — the scaling (§6.2).** Is `a ∝ m^{1/3}` right, given `F = A/a + B/a² + Ca³` and `m = (18/11)A/a`?
Redo it keeping `B`. Does the conclusion survive? Then: is it consistent with `pathA_26`'s cavity-mode
branch, or do the two reduced models genuinely disagree?

**Q4 — the instability (§5.3.1).** Read `pathA_26_derrick.md` Phase C. Is `THROAT_DRAIN_DESTABILIZED` a
real physical result or an artifact of the sharp-wall / slaved-drain reduction? How much of the model
does it threaten if real?

**Q5 — the restart order (§7).** Is R-1 → R-2 → R-3 → R-4 right? Should R-3 come first, on the grounds
that an unstable particle invalidates the rest?

**Q6 — apparatus.** Where does this plan rebuild census-grade machinery instead of doing physics?

**Q7 — what have I missed?** Assume there is at least one more `a`-pin-shaped identification in this
corpus — two are already known (`ℓ = δ` at `stage044:104`, and the pin itself). Name candidates.

---

## §9. Errors I made this session, recorded so you can calibrate

Stated so reviewers can see the failure pattern rather than trust the author.

1. **I claimed the sonic horizon gives the throat radius** — `r ~ (J/ρc_s)^{1/3}` from continuity. That
   was **the same move as the `a`-pin**: the sonic horizon is a bulk-flow surface, the mouth radius is a
   wall-sector property, and nothing identifies them. Retracted.
2. **I claimed `V_χ` was an unfilled slot and proposed "write the wall potential" as the restart.**
   It is written down with coefficients and the kink is solved (§5.1). Wrong target.
3. **I doubted a correct sub-agent finding** because my own `grep` window truncated the evidence, and had
   to reverse. The reviewer was right about the class *and* the specific loci.
4. **I over-warned about circularity** in calibrating `a` against a lepton mass. Using lepton masses as
   *input* to test whether the resulting radii are sane is legitimate; the circularity only bites if the
   ladder is later claimed as a *prediction*.

⇒ Pattern: **premature identification of two quantities that share a dimension**, and **asserting absence
from a partial search**. Both are exactly what this plan is trying to stop.
