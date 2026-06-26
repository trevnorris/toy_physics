# Directive pathA_35 — The light-confining shear-surface brane: a consistency-gate test (Gate L first)

**Status:** DRAFT v3 (Codex design-review xhigh + GLM tertiary both folded — see §14; GLM foresaw a four-way couple-stress no-go now
encoded as §2.6 + `FAIL_COUPLE_STRESS_NOGO`; awaiting final Codex confirm → user gate before any execution).
**Supersedes the *route* of:** `directives/pathA_25_gnls_polar_smectic_consistency_gates.md` (the GNLS polar-**smectic/density** route —
now CLOSED, see §0).
**Keeps the *content* of:** `reports/pathA_24_T0_freeze.md` (the frozen GNLS + polar-OP action, hash `8fa41ac51e88`).
**Conceptual source (read first):** `docs/conceptual_foundation.md` (§3 light deep-dive, §4 throat, §6.1 prior internal work) and
`docs/medium_requirements_and_prior_art.md` (requirements A/B/C, prior-art survey, the Gate L spec, the inherited walls).

---

## 0. Provenance & what changed (why this is a NEW directive, not a pathA_25 rung)

The **density-smectic route is CLOSED**, earned and tri-reviewed (ledger `reports/pathA_25_STATUS.md`):
- **`B4 = FAIL_NOT_CODIM1`** — the baseline layering driver makes a finite-`k` softening, but the GNLS cubic `U'''=15Kρ0²>0` makes a
  rank-2 equilateral **triad** beat the codim-1 lamella. No codim-1 density layer forms.
- **`R/C = RC_DENSITY_SMECTIC_LIGHT_NOGO`** — the only admitted driver that *does* open a genuine codim-1 density window (`Cpin`,
  `χ_Cpin<0`) **pins `P` along the layer normal** (`P_∥=0` at the static minimizer) → starves the in-plane shear that carries light.
  So **a density modulation can form a brane only by killing light** — a genuine no-go between *brane-exists* and *light-exists*.

The conclusion both results point to (and the user-chosen next front, 2026-06-26): **the brane is not a density layering; it is a
light-confining SHEAR surface.** The candidate is the **4D-throat-light hypothesis** (`docs/conceptual_foundation.md` §4; memory
`project-light-is-4d-throat-hypothesis`): light is intrinsically a (3+1)D in-plane shear field on the brane worldvolume (→ 2 transverse
polarizations automatically), confined because the **bulk is shear-free and light *is* shear** (no medium for it off the brane), and the
brane curves into `w` at throats (extrinsic curvature, not a bulk excursion).

**The methodology choice (approach A, user-decided 2026-06-26).** Rather than re-fight the *emergence* of the codim-1 surface (the
"why-3D" leg the prior-art survey twice flagged as the hardest, and which T1 already killed for the static domain wall), we **postulate
the codim-1 shear surface directly** — axis imposed by hand, labeled as a conceded inherited wall — and spend the budget on the **real
crux: is the postulated MacCullagh shear sector *dynamically consistent* on it (Gate L)?** (Phrasing it as "does light live on it" overstates
it — the MacCullagh form is postulated; the test is its dynamical consistency.) This is squarely the analog methodology
(`feedback-analog-find-consistent-structure`):
postulate the structure freely; the only test is internal consistency; falsification = a no-go between requirements.

---

## 1. The methodology (BINDING)

**Specify the FULL candidate structure (postulated freely — it is an analog), then test internal CONSISTENCY by HUNTING FOR A NO-GO.**
A no-go = **two or more requirements are mutually incompatible in this structure** (here the prime candidate: "any in-plane stiffness
clean enough to give curl-only 2-transverse light is *not* bounded-below / not angular-momentum-closed" — the Kelvin-gyrostat wall).
This is the adversarial target. "The minimal postulate wasn't enough" is a retired failure mode — we postulate the full package.

Binding consequences for every gate:
- **The FULL structure — brane localization, MacCullagh package, couple-stress sector, `P–u` coupling — is frozen UP FRONT in G0
  (§7), before any gate is computed.** Postulating an *ingredient* is allowed; postulating an *outcome* is not. **Adding an ingredient
  after seeing a gate result is `AD_HOC_RESCUE`** and requires a fresh G0. We may freely choose the action's terms (the MacCullagh
  stiffness, the couple-stress family, the `P–u` operator family, the confinement profile). We may **not** assert a mode is
  bounded-below, assert the modes are traction-carrying, assume the longitudinal sector is gauge, or assume the leak is small — those
  are exactly the **consistency questions the gates decide**.
- **Conditional-verdict rule (inherited from pathA_23/24/25).** This whole directive *relies on* two postulates: (i) the imposed axis /
  postulated codim-1 surface, and (ii) the postulated MacCullagh light-sector package (pathA_23 Stage 2 established the medium does not
  *derive* `μ_R` — `FAIL_UNSPECIFIED_SUBSTRUCTURE`). Therefore **the best attainable Gate-L verdict is `FREE_LIGHT_OK_CONDITIONAL`**, with
  both postulates named in the verdict. No unconditional "the medium derives light" claim is reachable here, and **no `decisions/14` /
  paper update follows a CONDITIONAL pass without explicit user sign-off.**
- **Every postulated knob is counted (§7 G0.4).** Because we postulate freely, the only guard against drifting into a second medium is a
  complete parameter ledger: every independent new constant AND function (confinement width, MacCullagh modulus/moduli, couple-stress
  inertia/stiffness, `P–u` coupling, `u_w` gap, brane inertia) is counted. ≥2 independent new inputs ⇒ `SECOND_MEDIUM_DRIFT` pressure,
  reported plainly (pathA_25 G0 cost 5; this package will cost some — state it, don't bury it).
- **A no-go is a first-class success** (`feedback-falsification-is-the-goal`). Never rescue or soften it. A clean "every hurdle passes"
  is the *suspicious* outcome and triggers extra adversarial scrutiny (the standing tri-review backstop, §11).

---

## 2. The candidate structure under test (the freeze content — concrete, so the gates are testable not generic)

The kept compressible superfluid `ψ=√ρe^{iθ}` on `X^i=(x,y,z,w)` provides the bulk, gravity-drain, magnetism, and `c_s`. A **postulated
codim-1 surface** at `w=0` (axis imposed) carries an **in-plane material displacement** `u` with a **rotational-elastic (MacCullagh)
stiffness**, anchored by a **couple-stress sector built from the T0 arrows `Pⁱ`**. Working name: **the shear-surface brane.**

### 2.1 KEPT, byte-for-byte (the existing program rides on this)
- **GNLS medium:** `ψ=√ρe^{iθ}`; quantum pressure `(ℏ²/8mρ)(∇ρ)²`; **single-well** EOS `U(ρ)=(K/4)ρ⁵` → `P(ρ)=Kρ⁵`,
  `c_s²(ρ)=5Kρ⁴/m` (B3); flow `v_i=(ℏ/m)∂_iθ−(q_*/m)A_i`; circulation/vortices → gravity-drain (B1) + magnetism (B2). The **bulk is
  shear-free** — load-bearing (B2 magnetism), and the native reason light is confined (it *is* shear, so it has no bulk medium).
- **The T0 polar OP** `L_pol` (`reports/pathA_24_T0_freeze.md` §2.2, frozen SHA-256 `8fa41ac51e88…`), **unchanged**:
  `L_pol = ½ mρ a²(D_t^v P)² − ½ mρ c_s²(ρ) a²(∂_j P^i)² − ¼ mρ c_s²(ρ)(|P|²−1)²`, `D_t^v P^i=∂_t P^i+v^j∂_j P^i`, `P^i∈ℝ⁴`, polar
  (`P≠−P`), carried by the medium (ρ-weighted, advected), one-constant Frank, O(4)-isotropic, no easy axis. **0 independent OP params.**

### 2.2 The POSTULATED codim-1 shear surface (replaces the §2.2/2.4 smectic driver of pathA_25 — conceded, labeled)
- The brane is a 3-surface at `w=0`. **The axis `w` and the surface's existence are POSTULATED, not derived** — the conceded inherited
  wall ("emergent axis / why-3D"; `docs/medium_requirements_and_prior_art.md` §inherited-walls). The verdict carries this conditionality.
- It carries an **in-plane displacement** `u^a(x,y,z,t)`, `a∈{x,y,z}` (tangent to the brane) — a **surface collective / material DOF of
  the SAME single medium** (the structured surface region of the one substance), **carried by the medium** so that tilting an arrow *is*
  a real material rotation (the traction condition of `conceptual_foundation` §3: only then does the stiffness become
  rotational-*elastic*, Frank→MacCullagh). It is **NOT** a second decoupled field — **but it is also NOT identified with the bulk
  advection velocity** `v^a`: `u^a` is **tangentially free-slip** from bulk mass transport (`u̇^a ≠ v^a`), because the shear-free bulk
  offers no tangential resistance (`decisions/15` §16). These two facets are not in tension; they are **jointly tested** — **L(a-i)**
  needs `u` *carried by the medium* (traction), **L(c)** needs `u` *free-slip from bulk transport* (no-leak). The out-of-plane / bending
  component `u_w` is the drumhead scalar (its gap is an **earned** hurdle, L(d)).
- **Localization & the shear-free bulk:** `u` is localized to `w≈0` by a frozen confinement profile `g(w)` (peaked at 0, decaying into
  the bulk); for `|w|≫` brane thickness the medium is the kept shear-free GNLS. The brane↔bulk interface coupling is the kept **full
  projected interface traction** `T_na = T_wa + (T_ww δ_ab − T_ab)∂_b u_w`, with the convective normal→tangential part `T_wa = m ρ v_w v_a`
  (**mass density `m ρ`** — T0's `ρ` is number density, so the `m` rides along; the pressure part is purely normal ⇒ free slip;
  `decisions/15` §16). These are the **two DIRECT leak channels** — convective (`v_w v_a`) and slope-mixing (`∂_b u_w`); Gate L(c) tests
  these **plus an indirect `u→Pⁱ→v` channel** (`Pⁱ` advected by the bulk flow, §8 L(c)) — all **not** new postulates.
- **The `u_w` bending scalar must be gapped** (a massless out-of-plane mode = an excluded fifth force; `decisions/15`/`conceptual_foundation`
  §3). The gap term + its scale are frozen in G0. (Honest hypothesis to keep visible, NOT to assume: per `project-light-is-4d-throat-hypothesis`,
  `u_w`/the longitudinal sector may be the *curvature-activated throat-binding channel* — non-propagating on the FLAT brane, active only
  at throats. Gate L tests only the flat brane; any throat activation is Gate T's.)

### 2.3 The POSTULATED light-sector constitutive package (frozen in G0)
pathA_23 Stage 2 established the medium does **not** derive a clean MacCullagh law (`FAIL_UNSPECIFIED_SUBSTRUCTURE`). Under the analog
reframe we **postulate** the package and freeze it up front:
- an **in-plane rotational-elastic (MacCullagh) stiffness** — energy in the curl of the in-plane material displacement, `∝ ½μ_R(∇×u)²`,
  modulus `μ_R>0` postulated — **and** a stated decision on the longitudinal sector: **either** a frozen scalar-potential analog `φ` /
  explicit constraint that removes it, **or none** (in which case the **C5 obstruction** below is an expected able-to-fail outcome of
  Gate L(a-iii)). **Provenance rule (anti-impose-the-answer):** any `φ`-analog or constraint must be a **local variational structure with
  independent mechanical provenance** (a real term in the frozen action), counted in G0.4 — a raw kinematic `∇·u=0` projector with no
  provenance **FAILS** (`FAIL_C5_LONGITUDINAL_ZERO_MODE`) or, if retained, must be **labeled an additional postulate in the verdict**
  (degrading the conditionality);
- a **spin / couple-stress sector** with **specified rotational inertia and stiffness**, frozen as a *physical sector* — **not** as "a
  reservoir that achieves closure" (whether it *achieves* angular-momentum closure is Gate L(b)'s to decide). The **preferred
  single-medium form reuses the T0 `Pⁱ` rotational inertia** as the Cosserat micro-rotation (**zero new DOF** — the arrows ARE the
  gyrostat elements, `conceptual_foundation` §3); any independent micro-rotation variable is a **new DOF counted in G0.4**;
- a **specified `P–u` coupling operator** tying `Pⁱ` rigidly to the in-plane displacement `u` (exact admissible operators frozen in G0,
  a finite target-blind family). Frozen as an *operator* — **not** as "a tie that makes the modes traction-carrying" (whether they ARE
  traction-carrying vs Frank torque-only is Gate L(a-i)'s to decide). **Honesty fork (decided by L(a-i)):** if the frozen `P–u` operator
  genuinely **sources** the rotational stiffness `μ_R` from the `Pⁱ` sector, the verdict may read *arrows supply MacCullagh traction*
  (`ARROWS_SUPPLY_TRACTION`); if `μ_R` is instead an **independent postulated surface modulus** (traction follows from the standalone `u`
  energy regardless of `Pⁱ`), the verdict must honestly read **`POSTULATED_SURFACE_ELASTICITY`** (`μ_R` an independent input counted in
  G0.4) — *not* "the arrows supply it."

### 2.4 The C5 MacCullagh gauge obstruction (named — this program's own verified failure mode, `decisions/15` §11)
A curl-only potential `½μ_R(∇×u)²` is invariant under `u→u+∇χ`, but the kinetic term `½ρ_br(∂_t u)²` is **not** invariant for
time-dependent `χ`; the EOM then forces `∂_t²(∇·u)=0` — the longitudinal mode is a **constrained physical zero mode, not a removable
gauge artifact** (Maxwell escapes via the scalar potential `φ→φ−∂_tχ`; MacCullagh has no `φ`). G0 must state whether the frozen package
carries a `φ`-analog/constraint; if not, **C5 is an expected able-to-fail outcome of Gate L(a-iii)**, not a surprise.

### 2.5 The Kelvin-gyrostat bounded-below hazard (named — the most-likely no-go)
**Positive curl energy is necessary but NOT sufficient** (`decisions/15` §11; `reports/pathA_25_STATUS.md`). The antisymmetric MacCullagh
stress carries angular momentum; without a couple-stress sector that *absorbs* it, the historical rotational ether suffered a
**negative-energy / angular-momentum-non-conservation instability** (the Kelvin gyrostat) — the system unwinds even with `μ_R>0`. Gate
L(b) must verify the **full Hamiltonian is bounded-below AND angular momentum closes through the frozen couple-stress sector.** This is
the gate's teeth and this program's honest prior for where it dies.

### 2.6 The couple-stress sector — parity, mode-count, and the foreseen FOUR-WAY no-go (GLM tertiary, folded v3)
The GLM tertiary pass foresaw a concrete **four-way incompatibility** the gate ladder must be able to detect and **attribute correctly**
(a no-go here is a first-class success, not a directive flaw — but misattributing it to L(b) alone would be a reporting failure). The four
hurdles L(a-ii) / L(b) / L(a-iii) / L(d) are **linked through the `Pⁱ` couple-stress / C5 / `u_w` structure:**
- **Closure needs LIVE `Pⁱ` modes; clean mode-count needs them GONE.** L(b-ii) angular-momentum closure needs a *dynamical* couple-stress:
  the closure divergence `∂_b m_{ab} ∝ γk²P` **vanishes as `k→0` if `Pⁱ` is gapped** (frozen micro-rotation → micropolar reduces to
  Cauchy → no reservoir for the antisymmetric MacCullagh stress `2μ_R(∇×u)`). But L(a-ii) (no hidden propagating modes) needs the `Pⁱ`
  spin waves **gapped or gauge-redundant.** Massless `Pⁱ` → 3 extra propagating modes (`FAIL_HIDDEN_PROPAGATING_MODE`); gapped `Pⁱ` → no
  low-`k` closure (`FAIL_GYROSTAT_NO_CLOSURE` at `k→0`).
- **The C5 `φ`-analog may collide with the `u_w` gap.** The only on-brane scalars are `u_w`, `ρ`, `θ` (the latter two used by GNLS). If
  the C5 `φ`-analog **is** `u_w` it must be **massless** (gauge fields are) — contradicting L(d) (`u_w` gapped, no fifth force). If `φ` is
  a new field, "zero new DOF" erodes. If absent, C5 fires (L(a-iii) FAIL).
- **Parity (polar vs axial).** The Cosserat micro-rotation `φ_a` is **axial**; `Pⁱ` is **polar** (T0: `P≠−P`). A linear `P·(∇×u)`
  coupling is therefore a **pseudoscalar (parity-odd)** — it cannot sit in a parity-even energy. The viable fix uses the brane normal:
  `ϖ_a = ŵ×P` (axial micro-rotation; polar×polar = axial — distinct from the C5 scalar potential `φ`), coupling `ŵ·(P×(∇×u))` —
  parity-even, **but it drags the conceded axis `ŵ` into the couple-stress
  sector and re-admits exactly the ε-contracted / chiral class T0 EXPLICITLY EXCLUDED.** So the honest structural cost is **0 new fields +
  1 new parity-constrained, axis-dependent coupling postulate** (count it in G0.4 as a structural constraint; the "arrows ARE the gyrostat,
  zero new DOF" framing must say *the coupling* — not the field — is the load-bearing new ingredient, with weaker provenance than T0).
- **The one possible escape — rigid coupling.** Slaving `P = ŵ×(∇×u)` (no independent `Pⁱ` modes → mode-count clean, closure automatic)
  is the candidate way through, **but it introduces `ω² = c_γ²k² + γk⁴/ρ_br` dispersion** — so L(a-ii)'s "non-dispersive" must explicitly
  decide whether a `k⁴` correction is tolerated, and at what scale.

**G0 must make these decisions UP FRONT** (§7 G0.3): the `Pⁱ` **mass/gap status** (massless | gapped | slaved-rigid, the last a named
alternate branch), the **`φ`-analog identity + provenance** (and whether it equals `u_w`), and the **parity structure** of the `P–u`
coupling. **Verdict:** if the four hurdles cannot be jointly satisfied for any frozen choice, the result is **`FAIL_COUPLE_STRESS_NOGO`**
(the four-way) — reported as a genuine no-go with the full chain named, *not* misattributed to a single hurdle.

---

## 3. Falsification stance & able-to-fail discipline (BINDING)

Every hurdle must be able to return BOTH outcomes from a real computation, with a negative control that flips it:
- **L(a-i) traction-vs-torque:** via a **virtual-work / momentum-balance test**, establish whether the frozen `P–u` operator converts
  `Pⁱ` orientation stiffness into a genuine **surface traction** (force/area across an in-plane cut). A **Frank-only reference** (`Pⁱ`
  decoupled from `u`) must compute to torque-only / no traction → flagged `FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION`. **And** if the
  traction is found to come only from a standalone `u` elastic energy (independent of `Pⁱ`), the verdict must downgrade to
  `POSTULATED_SURFACE_ELASTICITY` (not "arrows supply it"). If the test cannot distinguish these, it is tautological.
- **L(a-ii) hidden-mode audit:** curl-only 2-transverse is largely by construction, so the teeth are (1) a **Cauchy reference**
  (generic `½λ(∇·u)²+μ(∂u)²`) must compute to a stray propagating longitudinal third mode → `FAIL_CAUCHY_STRAY_LONGITUDINAL`; and (2) the
  **full coupled principal symbol** (incl. `P–u`, couple-stress, any `φ`/constraint, the confinement profile, `u_w`) must show **no extra
  propagating scalar/longitudinal mode** beyond the 2 transverse → else `FAIL_HIDDEN_PROPAGATING_MODE`.
- **L(a-iii) C5:** a package with **no `φ`-analog/constraint** must yield the constrained longitudinal zero mode → flagged
  `FAIL_C5_LONGITUDINAL_ZERO_MODE`; a package *with* the device must remove it (both branches computed); a no-provenance kinematic
  `∇·u=0` projector must itself FAIL or degrade the verdict (the §2.3 provenance rule), not silently pass.
- **L(b) bounded-below + closure (two-pronged — curl energy `≥0` is NOT enough):** the test must fire on **EITHER** (i) the full
  *coupled* system unbounded-below — checked by a **Hamiltonian / energy-matrix eigenvalue analysis, NOT the dispersion relation alone**
  (a gyroscopic `P–u` system can have all-real frequencies — "gyroscopic stabilization" — yet be energetically unbounded; dispersion-only
  would pass the Kelvin-gyrostat death) — **OR** (ii) a nonzero **angular-momentum-balance residual** (antisymmetric MacCullagh stress not
  balanced by the couple-stress divergence), tested **in the `k→0` limit** (where gapped `Pⁱ` freezes out and the reservoir vanishes).
  Two distinct controls: **(c1)** *omit* the couple-stress reservoir (DOF change) → ≥1 prong fires; **(c2)** *retain* `Pⁱ` but with a
  **large gap** (no DOF change) → the low-`k` closure residual must be nonzero. → `FAIL_NOT_BOUNDED_BELOW` / `FAIL_GYROSTAT_NO_CLOSURE`.
- **L(c) leak (direct + indirect channels):** (direct) compute the **full projected traction** `T_na = T_wa + (T_ww δ_ab − T_ab)∂_b u_w`
  with `T_wa = m ρ v_w v_a`; a **bent control** (`v_w v_a ≠ 0` and/or slope `∂_b u_w ≠ 0`) must give nonzero `T_na`, the **flat** wave must
  kill BOTH terms (§16 free-slip). (indirect) **`Pⁱ` is advected by the bulk flow `v` (T0 `D_t^v P`) and coupled to `u`**, so test the
  **`u→Pⁱ→v` channel**: does the `P–u` coupling, via `Pⁱ` advection, source bulk shear/vorticity even when the direct traction vanishes?
  Control: Frank-only (decouple `Pⁱ` from `u`) → the indirect leak must vanish. A nonzero indirect leak is `FAIL_LEAK_BREAKS_MAGNUS`
  (or relocate to Gate S), not to be left untested.
- **L(d) `u_w` gap (full coupled spectrum, not the frozen input):** compute the **full coupled** `u_w` spectrum (incl. `P–u`, the
  `φ`-analog, the confinement profile) and confirm **no massless scalar mode** — neither bare `u_w` nor a coupled descendant inheriting
  its character → else `FAIL_BENDING_MASSLESS_FIFTH_FORCE`. (Reading off the frozen gap is pass-by-construction; the real test is whether
  any coupling **un-gaps** it.) Control: a deliberately ungapped `u_w` must trip it.

A negative verdict (any no-go) gets HARDER verification than a positive (`feedback-negative-verdict-short-circuit`): truth in the
computed OUTPUT, not prose; post-exec adversarial-with-ablation is the non-optional backstop.

---

## 4. Honest priors (calibrate expectations before execution — do NOT bank these)
- **Gate L is the crux and the most-likely no-go.** Prior weight on the killer: **L(b) bounded-below / angular-momentum closure** > **L(a-iii)
  C5** ≈ **L(d) `u_w` gap** > L(a-i) > L(a-ii) (the curl-only form is postulated, so a/ii is largely by construction — its teeth are the
  hidden-mode audit + the Cauchy control).
- **The single most-likely shape of the failure is the §2.6 four-way `FAIL_COUPLE_STRESS_NOGO`** (GLM tertiary prior): closure needs live
  `Pⁱ` modes, clean mode-count needs them gone, the C5 `φ` may collide with the `u_w` gap, and the parity-correct coupling drags in `ŵ`.
  Expect the death to be *linked*, not isolated to L(b); the rigid-coupling branch (`P=ŵ×(∇×u)`, `k⁴` dispersion) is the one escape to test.
- **Best realistic harvest:** `FREE_LIGHT_OK_CONDITIONAL` (clean curl-only 2-transverse light, bounded-below via the `Pⁱ` couple-stress
  closure, free-slip no-leak) — CONDITIONAL on (imposed axis + postulated MacCullagh package). That is a genuine *win* for the analog: it
  would be the first time the light sector is internally consistent on a concrete brane.
- **A no-go (most likely L(b)) is equally valuable** — it would say the shear-surface brane cannot carry stable light either, which would
  be a deep statement about the whole light program (and would redirect to the lattice/Wen route we deliberately declined, or kill it).

---

## 5. Verdict grammar (every gate result carries one provenance tag + one outcome label)
**Provenance:** `DERIVED_FROM_KEPT_ACTION` (follows from GNLS+T0 alone — rare here) | `CONDITIONAL_ON(imposed_axis)` |
`CONDITIONAL_ON(MacCullagh_package)` | `CONDITIONAL_ON(both)` | `CALIBRATION_GAP`.
**Gate-L outcome labels:** `FREE_LIGHT_OK_CONDITIONAL` (best case — **free-field surface light only**, NOT a full Maxwell-dictionary /
charge / cone-lock / packet claim) **[sub-tag:** `ARROWS_SUPPLY_TRACTION` | `POSTULATED_SURFACE_ELASTICITY`**]** |
`FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION` | `FAIL_CAUCHY_STRAY_LONGITUDINAL` | `FAIL_HIDDEN_PROPAGATING_MODE` |
`FAIL_C5_LONGITUDINAL_ZERO_MODE` | `FAIL_BENDING_MASSLESS_FIFTH_FORCE` | `FAIL_NOT_BOUNDED_BELOW` | `FAIL_GYROSTAT_NO_CLOSURE` |
`FAIL_LEAK_BREAKS_MAGNUS` | **`FAIL_COUPLE_STRESS_NOGO`** (the §2.6 four-way: closure ⊥ mode-count ⊥ C5/`φ` ⊥ `u_w`-gap not jointly
satisfiable — the foreseen no-go; report the full chain, do not collapse it onto one sub-hurdle).
**Cross-cutting guards (any gate):** `FAIL_DIMENSIONAL` (§10) | `FAIL_TAUTOLOGICAL` (a control that cannot flip) | `ENGINE_DISAGREE`
(dual-engine mismatch beyond tolerance).

---

## 6. The gate ladder & ordering
1. **G0 — structure freeze** (§7): freeze the shear-surface brane action + the full parameter ledger + anti-circularity preregistration.
   Its *calculated* claims (dimensions; the flat-brane bulk/brane mode content) get the full dual-engine + arbiter + fidelity treatment
   (T0 lesson: a "freeze" that makes arithmetic claims still needs the engines).
2. **Gate L — light on the surface (THE CRUX, §8)** — the sub-hurdles a(i–iii), b, c, d on the postulated flat brane. **This is the
   make-or-break; everything downstream is gated on it.** Per §2.6, a(ii)/b/a(iii)/d are a **linked chain** through the `Pⁱ`/C5/`u_w`
   structure — a failure is to be reported as the linked no-go (`FAIL_COUPLE_STRESS_NOGO`), not collapsed onto one hurdle.
3. **Gates S / B / Q / T (§9)** — sketched here, fully specified gate-by-gate later (each its own execution-prompt design-review).
   S (magnetism preserved), B (brane↔gravity compat), Q (two `±w` charges), T (throat/mass + the relocated curvature leak).

**Routing:** a Gate-L no-go (esp. L(b)) **stops** the ladder and is the deliverable (report it; do not proceed to S/B/Q/T on a dead
brane). `FREE_LIGHT_OK_CONDITIONAL` proceeds to S.

---

## 7. G0 — structure freeze (anti-circularity preregistration)
G0 freezes, as a single hashed action, **before any gate computes**:
- **G0.1** the kept GNLS + T0 `L_pol` (§2.1), referenced by the T0 hash `8fa41ac51e88` (byte-for-byte; no edits).
- **G0.2** the postulated brane: the confinement profile `g(w)` family (finite-parameter, e.g. a fixed-shape localized profile with ≤1
  width), the in-plane displacement `u^a` as the medium's restricted displacement, the imposed axis `w` (labeled conceded), the brane
  inertia `ρ_br`. **Target-blind:** `g(w)` admitted on locality/minimality grounds only, never because it helps a gate.
- **G0.3** the light-sector package (§2.3) — including the four **§2.6 decisions** made explicit and target-blind:
  - the MacCullagh `½μ_R(∇×u)²`;
  - **`Pⁱ` mass/gap status:** `massless` | `gapped` | `slaved-rigid` (the last = `P=ŵ×(∇×u)`, a named alternate branch carrying a `k⁴`
    dispersion to be checked at L(a-ii)) — with the per-gate consequence of the choice stated;
  - **`P–u` coupling operator + parity structure:** a frozen finite family; because `Pⁱ` is polar and the micro-rotation axial, a
    parity-even linear coupling needs the brane normal (`ϖ_a=ŵ×P`, axial; distinct from the C5 scalar `φ`) — state this, note it re-admits the T0-excluded ε-contracted class, and
    that the couple-stress sector then depends on the conceded axis `ŵ`;
  - **`φ`-analog / longitudinal-sector decision:** present or absent; if present, its **identity** (is it `u_w`? a new field?) and its
    **independent variational provenance** (the §2.3 anti-impose rule) — and an explicit **`φ`-vs-`u_w`-gap consistency check** (a massless
    `φ`=`u_w` contradicts L(d));
  - the **`u_w` gap term + scale.**
- **G0.4 the complete parameter ledger — split into four explicit sub-counts** (each item tagged kept / postulated-ingredient /
  conceded-wall): **(i) fields / DOF** (reuse-`Pⁱ` = 0 new DOF; any independent micro-rotation = new DOF); **(ii) independent constants**
  (`μ_R`, `ρ_br`, couple-stress inertia/stiffness, the `u_w` gap scale, coupling normalizations); **(iii) independent functions** (the
  `g(w)` profile family, any kernel); **(iv) structural constraints / postulates** (the imposed axis, the Cosserat-reuse interpretation,
  the `P–u` operator choice, any `φ`-analog/constraint). **≥2 independent new inputs (constants + functions + structural) ⇒
  `SECOND_MEDIUM_DRIFT_AT_FREEZE(n)`, reported plainly** (judged later vs held-out surplus, not disqualifying — pathA_25 G0 cost 5 and
  proceeded by user gate).
- **G0.5 the flat-brane mode content** (a *calculated* claim → dual-engine): linearize the frozen action on the flat brane; report the
  in-plane spectrum's DOF count and the bulk's shear-free status. (This is the input Gate L dissects; G0 only *states* it, Gate L
  *interrogates* it.)

**G0 verdict labels:** `T0_SHEAR_FROZEN(<hash>)` | `SECOND_MEDIUM_DRIFT_AT_FREEZE(<n inputs>)` (record-and-proceed) | a fresh-G0 trigger
if any later gate needs an un-frozen ingredient.

---

## 8. Gate L — light on the shear surface (THE CRUX)
Run on the **postulated flat brane** (G0's frozen action, linearized). Each sub-hurdle is able-to-fail with the §3 negative control.

- **L(a) — the right wave (curl-only, 2-transverse, longitudinal not a stray physical mode):**
  - **L(a-i) traction-not-torque (+ provenance of the stiffness).** Via a **virtual-work / momentum-balance test**, establish whether the
    frozen `P–u` operator converts `Pⁱ` orientation stiffness into a genuine **traction** (force/area across an in-plane cut) vs only a
    **torque** on `Pⁱ`'s orientation (Frank → not light). Negative control: Frank-only reference → torque-only
    (`FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION`). **Provenance fork:** if the traction comes only from a standalone `u` elastic energy
    (independent of `Pⁱ`), the verdict sub-tag is `POSTULATED_SURFACE_ELASTICITY` (`μ_R` an independent input), **not**
    `ARROWS_SUPPLY_TRACTION`.
  - **L(a-ii) hidden-mode audit.** Solve the in-plane dispersion: exactly **2 transverse modes**, **non-dispersive `ω=c_γ k`** (state
    whether this means exactly linear or linear-at-low-`k`; for the slaved-rigid branch report the `k⁴`-correction scale and whether it is
    tolerated). Two teeth (the bare count is largely by construction): (1) a **Cauchy reference** (`½λ(∇·u)²+μ(∂u)²`) must yield a stray
    longitudinal third → `FAIL_CAUCHY_STRAY_LONGITUDINAL`; (2) the **full coupled principal symbol** (`u`, `Pⁱ`, couple-stress, any
    `φ`/constraint, the confinement profile, `u_w`) must show **no extra propagating scalar/longitudinal mode** — in particular the `Pⁱ`
    spin waves must be gapped/gauge-redundant, NOT 3 extra massless modes → else `FAIL_HIDDEN_PROPAGATING_MODE`. This is the mode-count leg
    of the §2.6 chain (it pulls against L(b) closure). `c_γ²(μ_R,ρ_br)` recorded for Gate K `λγ` later (flag, don't bank).
  - **L(a-iii) C5 longitudinal sector.** Show the longitudinal sector is **non-propagating on the flat brane** — removed by the frozen
    `φ`-analog/constraint (which must satisfy the §2.3 provenance rule), OR shown to be a pure (non-physical) gauge mode. A **constrained
    physical zero mode is a FAIL** (`FAIL_C5_LONGITUDINAL_ZERO_MODE`); a no-provenance `∇·u=0` projector FAILS or degrades the verdict. **If
    the `φ`-analog is `u_w`**, the gauge field must be massless — flag the collision with L(d)'s required gap (the §2.6 chain); a new-field
    `φ` is counted in G0.4. (The throat-activation hypothesis is explicitly OUT of scope — flat brane only.)
- **L(b) — bounded-below WITH angular-momentum closure (the teeth, two-pronged).** Curl energy `(∇×u)²≥0` is **necessary but NOT
  sufficient**. Establish **both:** (i) the full *coupled* system is bounded-below — by a **Hamiltonian / energy-matrix eigenvalue
  analysis**, NOT a dispersion-only check (a gyroscopic `P–u` system can be all-real-frequency yet energetically unbounded — gyroscopic
  stabilization); and (ii) the antisymmetric MacCullagh stress's angular momentum is **balanced by the couple-stress divergence**, tested
  **in the `k→0` limit** (where a gapped `Pⁱ` reservoir freezes out). Two controls: omit the reservoir (DOF change) → ≥1 prong fires; OR
  retain `Pⁱ` with a large gap (no DOF change) → nonzero low-`k` closure residual. (`FAIL_NOT_BOUNDED_BELOW` / `FAIL_GYROSTAT_NO_CLOSURE`.)
  **This is the most-likely no-go (see §2.6).**
- **L(c) — no leak into the shear-free bulk (direct + indirect).** *Direct:* compute the **full projected interface traction** `T_na =
  T_wa + (T_ww δ_ab − T_ab)∂_b u_w` with `T_wa = m ρ v_w v_a`; the flat wave (`v_w=0`, `∂_b u_w=0`) must kill **both** terms (free-slip,
  §16), a bent control must give nonzero `T_na`. *Indirect:* test the **`u→Pⁱ→v` channel** (`Pⁱ` advected by bulk `v`, coupled to `u`) —
  does the coupled dynamics source bulk shear/vorticity even when the direct traction vanishes? Control: Frank-only → indirect leak must
  vanish. Either channel nonzero → `FAIL_LEAK_BREAKS_MAGNUS` (or relocate to Gate S). The curvature² residual is **relocated to Gate T**.
- **L(d) — `u_w` bending scalar gapped, in the FULL coupled spectrum (no fifth force).** Compute the **coupled** `u_w` spectrum (incl.
  `P–u`, the `φ`-analog, the confinement profile) and confirm **no massless scalar mode** — neither bare `u_w` nor a coupled descendant
  inheriting its character → `FAIL_BENDING_MASSLESS_FIFTH_FORCE`. (Reading the frozen gap is pass-by-construction; the real test is whether
  any coupling un-gaps it.) Control: ungapped `u_w` must trip it.

**Pass = `FREE_LIGHT_OK_CONDITIONAL`** (with the L(a-i) provenance sub-tag) only if a/i, a/ii, a/iii, b, c, d **all** hold on the frozen
package. Any single isolated failure → the corresponding no-go label. **If the failure is the §2.6 linked chain** (closure ⊥ mode-count ⊥
C5/`φ` ⊥ `u_w`-gap), report **`FAIL_COUPLE_STRESS_NOGO`** with the full chain named — do NOT collapse it onto one hurdle. Either way the
ladder stops and the no-go is the deliverable (report it).

---

## 9. Gates S / B / Q / T (summary specs — each fully specified + design-reviewed when reached)
- **Gate S — magnetism preserved.** In-plane stiffness confined to the brane; inter-layer/bulk shear-free (Magnus/vortices intact).
  Tightly tied to L(c). Outcome: `MAGNETISM_OK` | `FAIL_BULK_SHEAR_BREAKS_MAGNUS`.
- **Gate B — brane↔gravity compatibility.** The brane sector must NOT disturb long-wavelength `c_s`/flow/drain or the existing
  GR-quadrupole bundle (`χ_Q`, `P0`; the gravity arm Gates 1–5). Outcome: `GRAVITY_COMPAT_OK` | `FAIL_DISTURBS_CS_OR_QUADRUPOLE`.
- **Gate Q — two charge signs.** Throat puncture direction `±w` (which side of the postulated surface) → two mirror, mass-independent
  charge signs. Outcome: `TWO_SIGNS_OK` | `FAIL_NOT_TWO_SIGNS`.
- **Gate T — throat / mass + the relocated leak.** A defect puncturing the surface → finite trapped-wave throat (consistency with
  pathA_26 `THROAT_DRAIN_DESTABILIZED`-not-a-kill); the curvature-localized leak from L(c) and the throat-activation hypothesis for the
  longitudinal/`u_w` channel are tested HERE. Outcome: `THROAT_OK_CONDITIONAL` | `FAIL_NO_FINITE_THROAT` | `FAIL_LEAK_AT_THROAT`.
- **Gate K — cone-lock `c_γ≈c_s`** is a **deferred downstream gate** (calibration-gap expected, per `medium_requirements`; flag
  `c_γ²(μ_R,ρ_br)` from L(a-ii), do not bank). It is **NOT** an inherited wall — the true inherited walls (conceded, not gates) are
  **dynamics / `G` / `α`**, the **emergent axis / why-3D**, and **Lorentz / preferred-frame**.

---

## 10. Dimensional firewall (user-mandated — "no LLM establishes a dimensional fact")
Every dimensional claim in G0 and Gate L is established by **dual-engine `dim_of`/`dimOf` on the REAL expressions, units restored,
able-to-fail** (`feedback-dimensional-consistency-check`): the check must FIRE on a deliberately corrupted term and must be **free of
any back-solved free carrier** (the Gate-4 trap: do not let a free constant absorb the homogeneity). No prose dimensional argument is
accepted as evidence. SymPy + Mathematica must agree (mismatch = `ENGINE_DISAGREE`). **Comprehensive + inline (binding):** the scripts
must dim-check **EVERY expression they construct, as they build it** — not a single end-pass on the action. That means each action term,
every derived coupling/stiffness/inertia, the full traction `T_na`, and **every quantity in each gate's computation** (linearized
operators, dispersion relations, `c_γ²`, energy/Hamiltonian matrices, closure residuals, etc.) — and they must all be **mutually
consistent** (the same constant carries the same dimensions everywhere it appears; combinations resolve to the right dimensions at every
use). A script that computes any quantity it has not dim-checked is non-compliant. **Note (correction propagated from GLM):** `decisions/15`
§16 writes the traction `T_wa = ρ v_w v_a` with `ρ` = T0 **number** density (dimensionally missing `M`); this directive uses `T_wa =
m ρ v_w v_a` (mass density) — the dim check must carry the `m` and FIRE if it is dropped. **Note (correction propagated from GLM):** `decisions/15`
§16 writes the traction `T_wa = ρ v_w v_a` with `ρ` = T0 **number** density (dimensionally missing `M`); this directive uses `T_wa =
m ρ v_w v_a` (mass density) — the dim check must carry the `m` and FIRE if it is dropped.

---

## 11. Discipline (BINDING)
- **Codex codes / designs the scripts; Claude reviews only** (`feedback-claude-reviews-codex-codes`). This directive states
  requirement + acceptance + verdict ladder; it does **not** pre-design the route.
- **Dual-engine** wherever Mathematica can independently verify (test = "is it possible," `feedback-dual-engine-required`): Mathematica
  leads / SymPy cross-checks (or vice-versa); Codex needs `--sandbox danger-full-access` to RUN `.wl`; the orchestrator can run `math`
  directly as arbiter; ≤2 concurrent `math -script`.
- **Reports-only**; **one gate at a time** with an explicit user gate between gates (`feedback-sequential-audit-chunks`); the **per-gate
  EXECUTION PROMPT gets its own design-review** before the expensive run (pathA_25 T1 lesson).
- **Tri-review every compute gate:** orchestrator arbiter re-run (refreshes committed output) + transliteration-fidelity audit (code-vs-
  equations, term-by-term) + **adversarial-with-ablation on a clean agent** (the non-optional backstop that caught pass-by-construction
  in the gravity Gates 2/3/4/5). Review rounds offloaded to a gauntlet-runner agent (`feedback-offload-review-gauntlet`).
- **No fake scripts** (no `python3 -c` commentary; read and reason); **YAML/markdown** for any LLM-read/written file; audit/verify
  prompts rendered **under the project root** (`software/stage1_solver/_scratch/`), not `/tmp`.
- Codex `exec … -c model_reasoning_effort=xhigh`, backgrounded with `< /dev/null`, **never wrapped in shell `timeout`**; the scripts it
  runs get `timeout 600` (a 124 = reformulate the math, never raise the cap).

---

## 12. Deliverables
- `reports/pathA_35_G0_freeze.md` (+ `_results.yaml`) — frozen action, parameter ledger + count, flat-brane mode content, dim-check.
- `reports/pathA_35_gateL_light.md` (+ `_results.yaml`) — the six sub-hurdle results (a-i, a-ii, a-iii, b, c, d) with their negative
  controls, the verdict.
- `tools/pathA_35_*_sympy.py` + `tools/pathA_35_*.wl` — the dual-engine scripts (Codex-authored).
- Resume/ledger update in `reports/pathA_25_STATUS.md` (or a fresh `pathA_35_STATUS.md`) + `STATUS.md` + `decisions/13` §0 at the gate.

---

## 13. Review plan (the gauntlet for THIS directive)
1. **Codex design-review (xhigh), via an agent** — iterate to GREEN (`feedback-review-ordering-codex-then-glm`, `feedback-offload-review-gauntlet`).
2. **One GLM tertiary pass** — fold its concerns.
3. **Codex confirm-pass (xhigh)** to GREEN again.
4. **User gate** → then the G0 execution-prompt design-review → execute G0 → Gate L.

---

## 14. Changelog
- **v1 (2026-06-26):** initial draft. Pivots pathA_25's closed density route to the postulated light-confining shear surface (approach
  A); axis imposed (conceded wall); Gate L is now the first gate; smectic/B4 dropped; MacCullagh + couple-stress (`Pⁱ`) + `P–u` package
  kept and made the crux.
- **v2 (2026-06-26):** folded Codex design-review (xhigh, SOUND-WITH-CONCERNS, 9 FIX + 2 NIT). Key changes: (1) `u^a` redefined as a
  surface collective DOF of the one medium — carried (traction) yet free-slip from bulk transport (`u̇^a≠v^a`, no-leak) — resolving the
  §16 leak-premise tension; (2) L(c) now tests the **full** projected traction `T_na = T_wa + (T_ww δ−T_ab)∂_b u_w` (both channels), with
  `T_wa = m ρ v_w v_a` (mass-density, carry the `m`); (3) C5 `φ`/constraint must have independent variational provenance (no
  impose-the-answer projector); (4) L(a-i) adds a virtual-work test + the `ARROWS_SUPPLY_TRACTION` vs `POSTULATED_SURFACE_ELASTICITY`
  honesty fork; (5) L(a-ii) reframed as a hidden-mode audit of the full coupled principal symbol (`FAIL_HIDDEN_PROPAGATING_MODE`); (6)
  L(b) made two-pronged (unbounded/complex-frequency OR angular-momentum-balance residual); (7) new hurdle **L(d)** = `u_w` gapped, else
  `FAIL_BENDING_MASSLESS_FIFTH_FORCE`; (8) Gate K relabeled a deferred downstream gate (not an inherited wall); (9) G0.4 ledger split into
  fields/constants/functions/structural sub-counts; (10) best-case pass renamed `FREE_LIGHT_OK_CONDITIONAL` (free-field surface light
  only, not a full Maxwell/charge claim).
- **v3 (2026-06-26):** folded GLM tertiary (SOUND-WITH-CONCERNS, 7 FIX + 4 NIT). Headline: GLM **foresaw a four-way no-go** now encoded as
  **§2.6** + verdict `FAIL_COUPLE_STRESS_NOGO` — closure (needs live `Pⁱ` modes) ⊥ mode-count (needs them gone) ⊥ C5/`φ` (may collide with
  the `u_w` gap) ⊥ parity (polar `Pⁱ` vs axial micro-rotation forces a `ŵ`-dependent, T0-excluded coupling); plus the rigid-coupling escape
  (`P=ŵ×(∇×u)`, `k⁴` dispersion). Also: L(b-i) now a **Hamiltonian-eigenvalue** check (gyroscopic stabilization ≠ bounded-below), L(b-ii)
  tested at **`k→0`** with a large-gap control; L(c) adds the **indirect `u→Pⁱ→v`** leak (Frank-only control); L(d) reframed to the **full
  coupled** `u_w` spectrum (frozen-gap readout was pass-by-construction); G0.3 now fixes the `Pⁱ` mass/gap status + `φ`-identity/provenance +
  `P–u` parity structure up front; the §16 `T_wa` number-density form corrected to `m ρ` (dim firewall). The "does light live on it" framing
  sharpened to "is the postulated MacCullagh sector dynamically consistent." Post-confirm cleanup (flagged by both Codex confirms):
  renamed the 2 straggler best-case labels to `FREE_LIGHT_OK_CONDITIONAL`, deliverable count five→six, and §2.2 "two leak channels" → "two
  direct + indirect `u→Pⁱ→v`" (consistency with the v3 L(c) edit).
- **v3.1 (2026-06-26, post-sign-off touches):** (a) §10 dim firewall strengthened to **comprehensive + inline** — dim-check every
  constructed expression as it is built (not an end-pass) and require mutual cross-consistency (user directive). (b) Symbol rename: the
  axial Cosserat micro-rotation `ŵ×P` is now **`ϖ_a`** (was `φ`, which collided with the C5 scalar-potential `φ`) — §2.6 + §7 G0.3. Both
  are clarity/rigor touches, no change to any gate's substance or verdict ladder.
