# Defect Interaction & Regime Map

**Purpose.** A living map of the 4D superfluid-analog particle (a finite brane→bulk *throat defect*): its
physical attributes, which attributes form **feedback loops** (bidirectional coupling) vs **one-way** vs
**decoupled**, and how the attributes behave across **regimes** (rest/ground, moving, excited, breakup, shedding,
spin). It records where existing work lives (file:line / file:section + derivation-atlas node IDs), what is solid
vs open, and the net-new gaps. It is the roadmap for the calibrate-and-predict *held-out surplus* and the place to
consult whenever a prediction **misses** (a miss is a missing-physics signal, not a cue to add a fit parameter).

**How to read.** Atlas node IDs (e.g. `PHYS_BRANE_BULK_THROAT_DEFECT`) resolve via
`python3 graph/query_graph.py show <ID>`. Confidence is flagged where a claim is inferred vs explicit in a source.
Scope note: framing here is the private engineering picture; the public/paper framing stays strict toy-analog.

---

## 0. One-line orientation

The particle is a **finite-radius open conduit (throat) through the brane into a 4+1 bulk**, not a dimple or a
capped pocket (`PHYS_BRANE_BULK_THROAT_DEFECT`; `notes/pde_audit_full.md:V2-28`, `notes/atom_work.md`). Its
field content is a gauged superfluid (GNLS) medium `ψ`, a localized Maxwell field `A_M`, and the throat shape
`R(Ω,w,t)`. At rest+ground state the throat is a "hyper-trumpet": a narrow brane **mouth** of radius `a` belling
out to a wider open **exit** `R(L)>a` (`PHYS_OPEN_FINITE_EXIT`: `R(0)=a, R(L)>0`).

---

## 1. Nodes — the physical attributes (grounded in the atlas `physical_ontology` layer)

| node (attribute) | what it is | math | source |
|---|---|---|---|
| **Throat defect** `PHYS_BRANE_BULK_THROAT_DEFECT` | the particle = finite open conduit | `Σ=0`, `r=R(Ω,w,t)` | `notes/pde_audit_full.md:V2-28`, `notes/atom_work.md` |
| **Superfluid medium** `PHYS_SUPERFLUID_MEDIUM` | density / current / EOS substrate (the "substance") | `ρ=|ψ|²`, `P=Kρ⁵` | `research/4d/paper/4d.tex`, `notes/pde_audit_full.md` |
| **Intake/output (inflow)** `PHYS_SUPERFLUID_INTAKE_OUTPUT` | flux-like throat intake → detected mass | leakage/work/export kernels | `notes/pde_audit_full.md`, `research/4d_plasma/paper/4d_plasma.tex` |
| **Mouth** `PHYS_MOUTH_CROSS_SECTION` | brane entrance cross-section (the aperture) | `a(t)=⟨R(Ω,0,t)⟩_{S²}` | `notes/pde_audit_full.md:V2-28`, `notes/moving_throat_notes_full.md` |
| **Interior support** `PHYS_INTERIOR_SUPPORT` | axial cavity carrying wall/BdG/Maxwell/port modes | `0≤w≤L`, D/N ladder | `notes/pde_audit_full.md:V2-28`, `notes/moving_throat_notes_full.md` |
| **Open conduit / exit** `PHYS_OPEN_CONDUIT`, `PHYS_OPEN_FINITE_EXIT` | open finite-radius branch (no hard cap) | `R(L)>0` | `notes/pde_audit_full.md:V2-28` |
| **Magnetic/vortical circulation (swirl)** `PHYS_MAGNETIC_VORTICAL_CIRCULATION` | quantized azimuthal circulation → magnetism | `∮(∂_iθ − q_*A_i/ℏ)dl^i = 2πn` | `research/4d/paper/4d.tex`, `notes/pde_audit_full.md:V2-30`, `notes/circulation/step_01_fluxoid_firewall.md` |
| **Localized EM sector (gauge)** `PHYS_LOCALIZED_EM_SECTOR`, `PHYS_MIXED_EM_CORE` | 4+1 Maxwell field + mixed (a–w) channels | `A_M, F_MN, Z(w)`; `A_w, J^w, F_{μw}` | `research/4d_em_fields/paper/4d_em_fields.tex`, `notes/pde_audit_full.md:V2-30` |
| **Charge** `PHYS_CHARGE_BRANCH` | oriented puncture label (electric sign) | `η_Q=±1`, `q_*=η_Q e_*`, `q_eff=q_*/√Z_int` | `research/4d/paper/4d.tex`, `research/4d_em_fields/paper/4d_em_fields.tex` |
| **Quadrupole / shape response** `PHYS_GROUPED_P2_SHAPE_RESPONSE`, `PHYS_OUTGOING_QUADRUPOLE_PORT` | l=2 STF shape mode + outgoing radiative port | grouped `{20,21,22}`; `iω⁵`, `Γ₂^port=a⁵/27c_s⁵` | `research/4d_2_5pn/paper/4d_2_5pn.tex`, `research/4d_4pn/paper/4d_4pn.tex` |
| **Wall constitutive law** `S_Σ` (Path-A) | tension/inertia/potential of the throat wall (a fixed law) | `μ_Σ,T_{w,Σ},T_{Ω,Σ},U_Σ` | `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md`, `OPEN_NONLINEAR_S_SIGMA` |
| **Response readouts (outputs)** `PHYS_RESPONSE_READOUTS` | low-order projected branch outputs (NOT the throat) | `D0,C,P0,N2,N4` | `notes/pde_audit_full.md` |

Output observables read off the readouts: **mass** (`R_norm` = GR quadrupole normalization), **magnetic moment**
(`μ=q_*S/2m_ψ`), **charge** (`χ_Q`).

---

## 2. Coupling structure — feedback vs one-way vs decoupled

Legend: ⇄ bidirectional feedback · → one-way · ⊥ decoupled (selection rule).

| edge | type | mechanism | source / confidence |
|---|---|---|---|
| density `ρ` ⇄ inflow `v_r` | ⇄ | GNLS amplitude↔phase (continuity/Euler) | `research/4d/paper/4d.tex`; high |
| inflow `(ρ,v_r)` ⇄ throat `R0` | ⇄ **(the Path-A loop)** | inflow pressure vs confinement `V_conf(R0)`; the throat shape sets how fast it can flow in, the inflow loads the wall | decision-08; `software/stage1_solver/`; high |
| wall law `S_Σ` → throat `R0` | → | constitutive law fixed; `R0(w)` solved from the static balance | decision-10/11; `OPEN_NONLINEAR_S_SIGMA`; high |
| chem. potential `μ` → `ρ` (+ weak back) | → (constraint) | sets global density scale; solved jointly for mass-normalization | decision-10; medium |
| gauge `A` ⇄ matter `ψ` | ⇄ | gauged GNLS covariant derivative + Maxwell source | `notes/moving_throat_pde_program_compact.md:556-630`; high |
| swirl `A` → magnetic moment | → (output) | circulation/twist *is* the moment | `notes/4d_lepton_notes.md:406-440`; high |
| **circulation `n_Γ` ⊥ throughput `Φ` (mass)** | ⊥ **(theorem)** | topology fixes the 1-cycle winding (magnetism) but NOT the 3-sphere throughput amplitude (mass); distinct invariants on distinct cycles | `notes/lepton_mass_notes.md:594-638`; **high** |
| gauge moment `Z0` → mass denom. `D0` | → (weak, residual) | the *even/conservative* Maxwell moment enters `D0=K−B0−Z0` | `notes/moving_throat_pde_program_compact.md:540-552`; high |
| **motion ⊥ quadrupole (l=2)** *linearly* | ⊥ **(theorem)** | the wall operator is O(3) block-diagonal in `l`, so different-`l` sectors don't linearly couple — no linear scalar (l=0) feed-down into l=2, and a uniform translation (l=1) cannot *linearly* source l=2 | `notes/moving_throat_pde_program_compact.md:2489` (block-diagonal; l=0 doesn't linearly enter grouped l=2), `notes/moving_throat_notes_full.md:29520` (no linear scalar feed-down); high |
| branch bundle `(K,B,Z,N)` → observables | → (extraction) | `D0,N0→R_norm`; `χ_Q` | `TARGET_PACKET_A`, `OPEN_QUAD_NORMALIZATION`; high |

**The two load-bearing decouplings** (both are *derived theorems*, not assumptions):
1. **"Magnetism doesn't change mass."** The circulation/throughput orthogonality theorem makes the azimuthal
   swirl (magnetism, winding `n_Γ`) orthogonal to the radial throughput (mass, amplitude `Φ`). Quantum-number
   firewall: `η_Q` (charge) / `n_Γ` (circulation·magnetism) / `Φ` (throughput·mass) / `τ` (spin) are all distinct
   (`notes/4d_lepton_notes.md:43-67`). **Caveat:** this is exact for the *topological* facet of the swirl; the
   *conservative* Maxwell moment `Z0` still feeds `D0` (a small residual gauge→mass coupling), so the swirl→mass
   edge is "⊥ for the winding, weak-→ for the even Maxwell moment."
2. **Motion doesn't linearly source radiation.** A uniformly moving (translating) defect cannot linearly excite
   the l=2 quadrupole — whether it does so *non*linearly is the unwritten translation↔throat bridge (§4, §6).

**A subtlety that ties it together (the near-pole amplification).** `B0+Z0 ≈ 0.0047 ≪ K = 4.06`, so at generic
parameters the back-reaction terms (including the swirl's `Z0`) are negligible and the decouplings hold cleanly.
But the GR target needs `D0 = K − B0 − Z0 → 0`, i.e. `K → B0+Z0` — and in that near-pole regime the tiny
back-reaction terms become co-determinative. *A coupling that is negligible at rest becomes load-bearing exactly
where we operate.* (decision-08; `PHYS_RESPONSE_READOUTS`.)

---

## 3. The "hyper-trumpet" and why mass is mouth-controlled

The rest/ground throat is a flared trumpet: narrow mouth `a`, belling exit `R(L)>a` (`PHYS_OPEN_FINITE_EXIT`).
The user's physical claim — *inflow (hence mass) is set by the **mouth**, so injecting energy expands the throat
interior but the leading mass stays invariant* — is **encoded in the model at leading order** (not exact total
mass invariance): with the harmonic wall
(`g=0`) the mass-relevant stiffness `K = τ·κ̂` is **independent of the throat shape `R0(w)`** (it is fixed by the
mouth aperture `a` and the material scale `τ`); the eigenmode `χ` is shape-independent too. Throat expansion
changes the sub-dominant overlaps `(B0,Z0,N0)` only weakly. So "mass invariant under throat expansion" ≈ "K is
mouth-set, throat-shape-independent" — the same statement two ways
(`software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md` §5a; `PHYS_MOUTH_INTERIOR_FIREWALL`).

---

## 4. Regime axis

| regime | technically | status | where / refs |
|---|---|---|---|
| **rest + ground (mass)** | stationary self-consistent background `[ψ,A,R0,μ]`; Schur `D0=K−B0−Z0`, `R_norm` | **building (Path-A)** | `software/stage1_solver/` (chunks 1a–1c done); decisions 08/10/11; `OPEN_QUAD_NORMALIZATION`, `TARGET_PACKET_A` |
| **shedding (radiation)** | radiation-reaction = odd-ω dissipative response; only the STF quadrupole survives universally | **normalization derived; branch-realized rate net-new** | `notes/5pn/`, `notes/25pn_notes.md`; `PHYS_OUTGOING_QUADRUPOLE_PORT`; `R_norm=0` *enforces* the 2.5/4/5PN Burke–Thorne quadrupole normalization |
| **excited / breathing** | higher-energy stationary branches + dynamical (ω≠0); throat breathes; `l=2` internal mode | **engine spec'd, no solver** | `notes/moving_throat*`, `notes/inner_throat/`; `PHYS_GROUPED_P2_SHAPE_RESPONSE` |
| **breakup / stability** | linearized spectrum; instability = eigenvalue→0 (`α_crit` buckling) | **partial; two conflicting pictures (§5)** | `notes/moving_throat_notes_full.md` (`α_crit`); `notes/session_barrier/` (reduced 1D); `inner_throat_hard_mode` |
| **spin / magnetism** | circulation winding `n_Γ` + mixed (a–w) Maxwell moment; magnetic moment | **kinematic yes; g=2 & spin-½ open** | `notes/circulation/`, `notes/g2/`, lepton notes; `OPEN_LEPTON_SPIN_DISCRETIZER`, `OPEN_G2_COMMON_LAYER` |
| **moving (translating)** | bulk translation through the medium; far-field wake | **exists but DISJOINT from the throat corpus** | 1PN/EIH wake (`notes/summaries/`): wake-mixing `α²=3/4`, SR `E=γE₀`; bridge to throat unwritten |

**Terminology landmine:** the "moving-throat" corpus models a **breathing/deforming wall** with the defect
spatially fixed — *not* translation. Genuine translating-defect physics is the separate 1PN wake corpus.

---

## 5. Breakup — experiment vs the two model pictures (the open question for the user)

The model currently has **two incompatible breakup pictures**, and **experiment helps adjudicate**:

- **User's mental model:** inject energy → throat expands → if it expands past the mouth, the particle **breaks
  apart** (predicts an *upper* energy bound / a closed stable island).
- **Reduced-1D `proton_proxy`:** breakup = *dressing-leg collapse* in a timescale race `S=t_cross/t_collapse<1`;
  found **one-sided** (higher energy → *more* stable, no upper edge) — the opposite topology
  (`notes/session_barrier/proton_proxy_stability_report.txt`). But this is a crude `V_eff(r)` closure, not the PDE.
- **Throat-PDE `α_crit`:** a genuine quadrupole **buckling** instability of the throat shape
  (`notes/moving_throat_notes_full.md`) — closer to "throat deforms until it can't be supported."

**What experiment says** (energy injection into a particle):
- **Elementary particles do NOT fragment.** An electron is point-like; you cannot shatter it. Dumping energy in →
  it **radiates** (bremsstrahlung), and **above threshold the energy converts into NEW particles (pair
  production)** — which requires a third body (a nucleus or external field) for momentum conservation; a single
  photon cannot pair-produce in vacuum. The original particle persists.
- **Composite particles DO break up.** A proton/nucleus fragments via deep-inelastic scattering / spallation once
  the momentum transfer `Q²` (or deposited energy) exceeds the binding scale → "current + target fragments" that
  hadronize.
- **"Too fast" maps to sudden vs adiabatic.** What matters is whether energy is deposited *faster than the system
  can shed/relax it*; a sudden transfer exceeding a threshold pushes the system past breakup before it radiates
  the energy away.

**Implication (CONCEPTUAL — for the user to weigh, see §7).** If the defect models a **lepton (elementary)**,
experiment says the reduced-1D fragmentation picture is *wrong*: the correct high-energy behavior is **shed
(radiate — the very `R_norm`/quadrupole channel) + pair-produce a defect–antidefect** above threshold, NOT a
shattering throat. The user's "sheds so fast that breakup needs very rapid bombardment" then fits the elementary
picture exactly: the defect radiates energy faster than you can pile it in, so it *resists* breakup, and the real
escape valve at extreme energy is **pair production**, not fragmentation. "Breakup" for an elementary-like defect
is best modeled as (a) the **`α_crit` stability bifurcation** (the excited static branch ceases to exist → decay)
and (b) **pair production**, not the reduced-1D collapse. Whether the defect is elementary-like or composite-like
is a model-level question (§7).

---

## 6. The big unification — three programs bottleneck on the Path-A solve

5PN radiation, the moving-throat radiative engine, and g−2 all symbolically reduce to the **same** open
computational gate = solve the self-consistent throat branch and read off `χ_Q = N_Q = 1` (≡ `R_norm = 0`):

- **5PN / shedding:** `R_norm=0` *enforces* the 2.5/4/5PN radiation-reaction Burke–Thorne quadrupole normalization
  (`R_norm = m̂0²P0 − 54Gc_s⁵/(5a⁵c⁵)`; the physical Burke–Thorne coefficient is the invariant product `m̂0²P0`,
  not `R_norm` itself) — same `D0=K−B0−Z0`, same `N0`; the 4PN tail is fixed by the same `γ₅`
  (`notes/5pn/5pn_notes_full.md`; `OPEN_QUAD_NORMALIZATION`, `OPEN_5PN_EVEN_GATES`).
- **moving-throat:** lane operator `D_{A,0}=K_A−B_{A,0}−Z_{A,0}` is the per-`l=2`, ω→0 instance of Path-A's `D0`;
  Path-A = static-monopole core, moving-throat = its ω-expanded l=2 generalization. Handoffs are symbolic specs,
  no implemented solver; the V2-23 demo missed ~7 orders (`P0/P0_target≈0.00183`) (`OPEN_EXECUTABLE_BRANCH_SOLVER`,
  `TARGET_PACKET_A/B`).
- **g−2:** reproduces `g_loc`, misses the electron by `2.27e-12`, pending the same canonical-branch datum
  (`notes/g2/g2_full_output.md`; `OPEN_G2_COMMON_LAYER`).

**Consequence:** finishing the Path-A `R_norm` calibration is the shared computational gate for all three — full
completion of the 5PN / moving-throat normalization, and the common branch datum for g−2 (a necessary
prerequisite, not the whole anomaly theorem: `OPEN_G2_COMMON_LAYER` is `conditional_open` and the exact electron
sliver is not naturally forced).
Open gates feeding it: `OPEN_PARENT_PROMOTION_S_SIGMA`, `OPEN_NONLINEAR_S_SIGMA` (the `S_Σ` promotion, =
decision-11), `OPEN_SOURCE_PORT_NORMALIZATION`, `OPEN_MATERIAL_CLOSURE`, `OPEN_WEAK_AXISYM_ORBIT_LOCK`.

---

## 7. Net-new gaps + open conceptual questions

**Net-new (not yet done anywhere):**
- **Shedding *rate* — the branch-realized luminosity / inspiral decay.** The instantaneous-power + Schott
  decomposition identities exist (`research/4d_2_5pn/paper/4d_2_5pn.tex:4688-4702`); what's missing is the realized
  energy-loss rate computed from the actual normalized branch.
- **The translation↔throat bridge** — does a moving defect's wake *non*linearly source the throat quadrupole? (The
  block-diagonal `l`-decoupling only kills the *linear* route; a nonlinear wake→l=2 feed is not excluded.)
- **Throat-level breakup** — `R` as an evolved field with a mouth-overrun / `α_crit` bifurcation condition, to test
  the user's picture vs reduced-1D.
- **Pair production** — defect–antidefect creation as the genuine high-energy channel (per §5).
- **g=2 Dirac factor and spin-½ quantization** (`OPEN_LEPTON_SPIN_DISCRETIZER`); **repeated-injection/bombardment**.

**Open conceptual questions for the user (would change the model's character — do NOT resolve unilaterally):**
1. **Is the defect elementary-like or composite-like?** This sets what "breakup" means (pair production +
   bifurcation vs fragmentation). §5 argues the lepton case is elementary-like.
2. **Should the swirl→mass coupling be exactly zero?** The topological facet is (theorem); the conservative `Z0`
   facet is not. Is the residual `Z0→D0` physical (EM self-energy) or an artifact to remove?

---

## 8. Path-A path forward (assessment)

The map sharpens, rather than changes, the plan:
1. **Finish the Path-A `R_norm` solve first** — it is the shared gate for 5PN / moving-throat / g−2 (§6). Concrete
   next step: USER `frozen: YES` on `decisions/11` → build the §5 extraction module → the unique `R_norm(τ)=0`
   root-find (task #62).
2. **Treat a miss as missing physics (κ_PV discipline), not a rescue DOF** — and this map is where to look: the
   richest candidates are the **shedding rate**, the **translation↔throat bridge**, and **pair production**.
3. **The held-out predictive surplus** (the thing that makes the analog credible) lives in the **moving**,
   **breakup**, and **shedding-rate** regimes — sequence them after the rest test, cheapest first
   (breakup/stability spectrum ≈ already partly computed → moving/bridge → shedding rate).

---

## 9. Next steps (resume-here after /compact)

**Immediate — the computational-physics block (Path-A, the shared gate):**
1. USER `frozen: YES` on `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md` + compute
   `candidate_freeze_hash` over its §8 (before any calibration touches the anchor).
2. Build + validate the §5 field→D0 extraction module (Codex codes / Claude reviews + transliteration-fidelity
   audit); re-solve the closed system with the frozen `homogeneous_isotropic_hooke_v1` (`g=0`) forms.
3. The unique `R_norm(τ)=0` root-find (stable-side `D0>0`) — the decisive Path-A test; report `τ*` + naturalness +
   held-out `R_pole/P2/P4` + §J calibration-covariance / Schur-margin error bars. (task #62)

**If it misses:** κ_PV discipline — diagnose + **derive** the missing physics, never add a rescue DOF. Richest
candidates from this map: the shedding *rate*, the translation↔throat bridge, pair production.

**Parked for the user (conceptual, §7):** elementary-vs-composite defect; exact-zero vs residual swirl→mass.

**Held-out surplus roadmap (after the rest test):** breakup/stability spectrum (cheapest, partly computed) →
moving / translation↔throat bridge → shedding rate.

---

## References

Internal: derivation atlas `graph/` (`query_graph.py show <NODE_ID>`); `notes/5pn/`, `notes/moving_throat*`,
`notes/session_barrier/`, `notes/circulation/`, `notes/g2/`, `notes/lepton_mass_notes.md`,
`notes/4d_lepton_notes.md`, `notes/pde_audit_full.md`; `software/stage1_solver/decisions/08,10,11`.

External (experimental grounding, §5): pair production / elementary particles don't fragment —
[Britannica: Pair production](https://www.britannica.com/science/radiation/Pair-production),
[Physics LibreTexts 4.3 Pair Production](https://phys.libretexts.org/Bookshelves/Modern_Physics/Spiral_Modern_Physics_(D'Alessandris)/4:_The_Photon/4.3:_Pair_Production);
composite-particle breakup —
[Wikipedia: Deep inelastic scattering](https://en.wikipedia.org/wiki/Deep_inelastic_scattering).
