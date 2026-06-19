# Decision 13 — B2c verdict is UNDETERMINED; the real next chunk is the EMERGENT-CONSTANTS derivation

**Date:** 2026-06-19
**Status:** DECIDED direction (user-driven, 2026-06-19). This is the **resume-here record after /compact** for the
Path-A build. Supersedes the "rigorous MISS" reading of B2c (decision-12 B2c STATUS block — now flagged superseded).
**Mechanism:** B2c 3-round build/review → two audits (`pathA_17` validity, `pathA_18` dimensional) → two independent
verification agents → user methodology call (derive the emergent constants before `m̂0²·S_port`).

---

## 0. STATUS / NEXT ACTION (resume here after /compact — 2026-06-19)
**Where we are:** `pathA_19` (foundation), `pathA_20` (`c_s`+`c`), and `pathA_20b` (`c_γ/c_s`) are ALL EXECUTED +
REVIEWED (fidelity-clean dual-engine; all verdicts adversarially verified HONEST) + COMMITTED — ledgers in §8/§9/§10.
**The emergent-constants program has now extracted everything derivable from the symbolic action + dimensions** (the
convergence, §10): base `{L,T,M}` retained; `c_s∝ρ²`, `[c_s]=[c_γ]=L/T`, `[J]=T⁻¹` derived; **`c_γ/c_s` is a CALIBRATION
KNOB** (photon rides the bulk cone, sound the emergent acoustic cone; the non-relativistic matter sector doesn't lock
them); the flux-law value, brane photon cone, `G`, and `α_J` all need the **SOLVED THROAT PROFILE**; `ħ`-emergence needs
new substrate physics. `pathA_21` (the SYMBOLIC spec-completion step: derive `G` + mass-bridge FORMS + re-test M-collapse
+ produce the PROFILE-SOLVE SPEC) is **WRITTEN + design-reviewed**: Codex SOUND-WITH-FIXES → all fixes applied →
confirm-pass **SOUND-AS-IS** (reviews in `_scratch/pathA_21_directive_review*.md`). Execute prompt STAGED at
`_scratch/pathA_21_execute_prompt.md`. Tree clean (everything committed).
**NEXT ACTION (resume here after /compact — user chose option A, 2026-06-19): FIRE the `pathA_21` EXECUTION run.**
Command (backgrounded, NEVER shell-timeout-wrapped):
`codex exec --sandbox workspace-write -m gpt-5.5 -c model_reasoning_effort=xhigh - < software/stage1_solver/_scratch/pathA_21_execute_prompt.md > software/stage1_solver/_scratch/pathA_21_execute.log 2>&1`
Then REVIEW per the directive's Review section (transliteration-fidelity on the new module + the `.wl`; an adversarial
pass with **distrust-restated-target** — for each of P1/P2/P4, is there a real source-equation chain or a restatement?
is `G` extracted only after a `G`-free force law? is EP derived from two SEPARATE masses?; plus a check that the P5
profile-solve spec is concrete + complete); Claude reads only residuals → gate to `pathA_22`.
**Discipline reminder:** Codex derives/codes, Claude reviews. `pathA_21` is SYMBOLIC — the win is the rigorously-DERIVED
FORMS (`C_F`, `m_defect↔J`, `G`) + the P5 profile-solve spec, NOT numbers. Expected honest outcomes that are VALID
PASSES: `{L,T,M}` retained (`ħ` undetermined), symbolic/profile-dependent `α_J`/`G`, `MASS_BRIDGE_FORM_NOT_DERIVED`,
`EP_NOT_DERIVED`, `NEWTON_G_FORM_NOT_DERIVED`/`FORCE_NOT_NEWTONIAN`. The DERIVED-FORM GATE forbids restating any target
(`α_J:=m_defect c²/(ħJ)` or `G` by rearranging `F=Gm₁m₂/r²` are FAILs). After `pathA_21`: option **C** = the throat-profile
SOLVE (closes the numbers, driven by the P5 spec); `pathA_22` = scale-map → `m̂0²·S_port` → re-run B2c.

---

## 1. The headline: B2c does NOT give a validated hit-or-miss
B2c computed, on the self-consistent Path-A background, that the model's `P0` is **6.7–9.6 orders of magnitude below**
the GR quadrupole target `54/5 = 10.8` (measured at every converged τ: `P0=2.795e-9` at τ=1 … `1.22e-6` at
τ=0.029). The numbers are real. **But the comparison is `R_norm = m̂0²·S_port·P0 − 54Gc_s⁵/(5a⁵c⁵)`, and B2c pinned
`m̂0²·S_port = 1`.** Two audits + two independent verification agents established:

- **NOT a dimensional/units bug.** With units restored, the dimensional gaps the audit found are **order-neutral**
  (they equal `×1` under the natural-unit pins `a=c_s=ħ=m=1`); restoring them shifts NOTHING numerically. (This
  *corrected* `pathA_18`, which had overstated "8 inconsistencies" → really 2, and wrongly assumed a 4-spatial `G`.)
- **NOT a validated "missing physics" miss either.** `m̂0²·S_port` is **dimensionful** — it is the conversion factor
  `T_target/[P0]`, NOT a free dimensionless number (both verification agents, independently). Pinning it to `1` is a
  **scale choice** that *forces* the `D0→0` knife-edge. The pre-reg bakes this in: `T := 54Gc_s⁵/(5a⁵c⁵)/(m̂0²S_port)`,
  `D0 = N0/T` (`docs/pathA_preregistration.md` §253). So `m̂0²·S_port` directly sets whether calibration lands at a
  knife-edge (miss) or naturally (match).
- **The ~9-order "miss" is exactly the magnitude `m̂0²·S_port` would supply** (`10.8/P0 = 3.86e9` at τ=1). Whether the
  *derived* `m̂0²·S_port` is ≈1 (→ real miss) or carries large factors (→ match) is **plausible either way** (the
  healing length `a` is microscopic while the GR comparison is at the orbital/radiation scale, so a large
  `(scale-ratio)^n` is on the table) and **cannot be settled by dimensional analysis** — only by deriving it.

**Verdict: UNDETERMINED, hinging on the un-derived dimensionful normalization `m̂0²·S_port`.** B2c is NOT committed as
a miss. The B2c machinery (assembly/dual-engine/verdict-logic/warm-start/negative-controls) was fidelity-clean and is
retained; only the *interpretation* changed.

## 2. The foundational reframe (user-driven): the constants are EMERGENT, derive them first
`c`, `c_s`, `G` are NOT fundamental inputs in this model — they are **outputs of the 4D superfluid PDE**, and their
familiar "3+1" dimensions are the **brane (observed) dimensions**, while the **bulk PDE works with 4D-dimensioned
quantities**. The bulk→brane reduction (integrating out the transverse `w`-direction, width ~`a`) is where the
scale/dimension factors live — and `m̂0²·S_port` is almost certainly a piece of *that same bulk→brane normalization*.
So **deriving the emergent constants is logically PRIOR to `m̂0²·S_port`**; once the constants + reduction are in
closed form, `m̂0²·S_port` should largely fall out.

What the verification already established about the constants:
- **`c_s` (sound speed):** emergent from the EOS, `c_s² = 5Kρ⁴/m`; `[c_s]=L/T`; **value drifts with the background
  `ρ`** (not constant).
- **`G`:** emergent from defect back-reaction. Bulk 4-spatial Laplacian → `1/r²` (`[G_bulk]=L⁴T⁻²M⁻¹`); but the
  observed sector is **brane-projected to 3 spatial dims** → `1/r`, standard Burke-Thorne/Peters, `[G]=L³T⁻²M⁻¹`
  (confirmed from the parent action: defects live on the 3D brane, `w` transverse). **User hypothesis (likely right,
  must DERIVE): `G` reduces to the standard 3D form after the brane reduction, with the transverse-scale (`a`)
  factors being the content of the emergent `G`.**
- The GR target `54Gc_s⁵/(5a⁵c⁵)` is NOT dimensionless on its own (it carries `M⁻¹L⁻²T⁻²` with 3D `G`); it is used as
  the Peters **pure number `54/5`** with constants pinned to 1, and `m̂0²·S_port` silently absorbs the dimension.

## 3. The velocity structure (refined with the user, 2026-06-19): THREE scales + the standing-wave bridge
The model has **three distinct velocities, all `L·T⁻¹`** (so dimensional analysis alone cannot separate them — only the
dynamics can). The earlier "c vs c_s = wave-vs-terminal-velocity" framing is superseded by this richer picture:
- **`v_b` — background condensate flow velocity** `v_b=(ħ/m)∇θ` (drift of the medium). This is the **gravitational-sector
  variable**: gradients of `v_b` and `ρ` ARE the analog field; `v_b=c_s` is the acoustic horizon. Varies in space, not
  constant. **`ρ` is ALSO not constant and is coupled to `v_b`** via stationary continuity `∇·(ρv_b)=0` + quantum-
  Bernoulli `½mv_b²+μ(ρ)+V+Q=const` (flow speeds up → `ρ` drops → `c_s` drops). So `ρ, v_b, c_s` are ONE coupled
  profile, not independent dials; genuine constants are `K, m, ħ, a`, the Bernoulli head, the mass flux, and the
  asymptotic `ρ₀/c_{s,0}` (the `c_s=1` pin = `c_{s,0}`). (Its full `G` role → `pathA_20`.)
- **`c_s` — bulk sound speed** (density/phonon waves through the medium). EOS-set, `∝ρ²`.
- **`c_γ` — photon/gauge-wave speed** (massless gauge excitation on the brane); the brane light cone.

**The standing-wave bridge (the key insight).** A massive particle (throat) is a **standing wave of the photon/gauge
field** — two counter-propagating `c_γ`-waves. Consequences, to be derived: (a) rest internal oscillation = Compton
frequency `ω₀=m*c_γ²/ħ`, so `E_rest=m*c_γ²` is trapped-wave energy; (b) driving the envelope at `v` Doppler-shifts the
components and slows the internal clock as `ω₀/γ`, **freezing at `v→c_γ`** (= relativistic time dilation, the user's
"oscillations stop at light speed"); (c) the envelope cannot outrun its constituent waves, so the **terminal-velocity
ceiling is `c=c_γ`** — `c` is the terminal velocity AND the photon speed, unified because matter is standing light.
This makes the old `c≠c_s` hypothesis and the existing `c=c_s` result COMPATIBLE: `c=c_γ` is forced;

**The one surviving open number: `c_γ/c_s`.** Particle ceiling = photon speed is forced; whether the photon (brane
gauge) speed equals the bulk sound speed depends on whether gauge + density share the acoustic metric, or the
localization profile `Z(w)`/width `a` rescales the brane gauge cone (the EM sector already shows `μ₀^eff=μ₀/Z_int`).
**`c_γ/c_s` (closed form) IS the `(c/c_s)³` tail factor `R_tail`** (=1 iff they coincide). Let the math decide.

**Is DEFECT mass emergent? (user idea, 2026-06-19; sharpened by Codex directive-review into a falsifiable test.)** The
idea: a defect is a SINK pulling superfluid inward, and that inflow strength is the defect's mass; the same inflow
shared between two defects is their attraction, so **`m_defect` and `G` are two faces of one quantity** (→ `pathA_21`).
**Two masses must NOT be conflated:** `m_GNLS` (the constituent mass in the parent action's kinetic/EOS/current terms)
is part of the EXACT action and stays `[M]` unless an action-level rewrite keeps every term homogeneous; only the
DEFECT mass `m_defect` (a throat branch property) might be emergent. The ontology says defect mass ↔ drainage/VOLUME
deficit (volume flux), charge ↔ vorticity flux — so do not assume mass = number flux; frame-tag number flux
`J=T⁻¹`, volumetric `Q_vol=ρ⁻¹J`, mass flux. Honest relation: `m_defect=α_J ħ J/c_γ²` (`α_J` dimensionless/branch
data); `m_defect=J` only after `ħ=c=1` (→ using the not-yet-derived `c` to drop `M` is CIRCULAR; de Broglie route is a
note to `pathA_20`). **Falsifiable:** `{L,T}` (mass emergent) is allowed ONLY if the whole action rewrites homogeneously
AND a boundary/source/Noether/Hamiltonian derivation ties `m_defect` to the inflow; else RETAIN `{L,T,M}` and record
`J` as a conserved rate — a valid negative result. **`a`-pin:** `a` is a mouth-radius collective MOMENT (not a
fundamental coordinate), particle-dependent + deformation-fragile; the conserved invariant is the flux `J` (Gauss, in a
no-leakage region — but the throat bottom may be open, so no-net-accretion must be a derived BC or a logged gap). Assess
re-pinning `a→J`, but do NOT claim it changes `m̂0²·S_port` (may be dimensionally neutral) — that value is `pathA_22`.

## 4. NEXT CHUNK — the emergent-constants derivation (the agreed next step)
**Execution split (user call, 2026-06-19 — foundation first):** the derivation is run as FOUR gated sub-directives, not
one: **`pathA_19` (FOUNDATION** — base set / mass fork M-vs-`J` / flux-invariant re-pin of `a` / pin over-determination
/ dictionary / paper-prose reconciliation; this is split out first because it can change the base dimensional system
everything else sits on; ✅ DONE 2026-06-19, see §8) → **`pathA_20`** (`c_s` + `c`, the velocity structure; FINALIZED
against the pathA_19 base — now also carries **S2b** the throat as a transonic/choked drain → the flux law
`J_crit(ρ₀, geometry)` set by the sonic point = acoustic horizon, replacing the "constant `J`" assumption, and **S3**
the `ħ`-provenance fork [`ħ=m c_s0 a/√2` ⇒ `ħ` may be EMERGENT, not fundamental — decides whether `M` ever collapses])
→ **`pathA_21`** (`directives/pathA_21_emergent_G_mass_bridge.md`: `G` from defect back-reaction + 4D→3D reduction;
the **mass-bridge** `m_defect=α_J ħ J/c²` [= `E=mc²` with the inflow as the energy ≡ the standing-wave `ħω₀`] DERIVED +
the **M-collapse re-test** of the pathA_19 negative + the equivalence-principle check + the `m↔G` unification) →
**`pathA_22`** (natural-unit→physical scale map → derive `m̂0²·S_port` → re-run B2c → real verdict). Each gated by the
user. Conceptual basis for S2b/S3/the bridge co-developed with the user 2026-06-19 (see §3 + the directives).

A full **reasoned + dimensional derivation** of the model's constants, from the parent action outward (NOT assumed —
derived; much already exists scattered across `research/pde_ledger/paper/parts/part01_parent_geometry.tex`, the
sound-speed ledger, and `research/4d_2_5pn/paper/4d_2_5pn.tex` — so partly consolidate + make rigorous + dimension-
check, partly fill gaps):
1. **`c_s`** from the EOS (`c_s²=5Kρ⁴/m`); its dimension + state-dependence.
2. **`c`** — what it physically IS (test the `c≠c_s` / terminal-velocity hypothesis), its definition, and its
   relation to `c_s` (the `c/c_s` ratio).
3. **`G`** — DERIVE from the defect back-reaction + the 4D→3D brane reduction: closed-form in superfluid params
   (`K,a,ρ,ħ,m`), bulk vs brane dimensions, the transverse-scale factors; TEST the "reduces to 3D `L³T⁻²M⁻¹`"
   hypothesis explicitly.
4. **The natural-unit → physical-scale map:** what the pins `a=c_s=ħ=m=1` actually stand in for — which determines
   whether `m̂0²·S_port` is genuinely ≈1 or carries orders.
Run everything through the committed harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
(machine-verify); produce a consolidated "constants & dimensions" reference doc (4D-bulk vs 3D-brane frame tags +
reduction maps). Discipline: Codex derives/codes, Claude reviews; dimensional-check + fidelity + adversarial.

**THEN:** derive `m̂0²·S_port` (should largely follow from #3/#4) → re-do the B2c `R_norm(τ)=0` calibration with the
derived normalization → real hit-or-miss verdict.

## 5. GATE-A freeze implication (methodology decision, do NOT do unilaterally)
`m̂0²·S_port=1` is currently part of the **GATE-A freeze** (decision-11 §2/§3, hash `ed3585…`). Deriving it means
**un-pinning a frozen convention → a documented freeze amendment (new hash)**. This is a methodology call to settle
with the user when the constants derivation is in hand (it may show the pin was fine, or that it must be replaced by
the derived value).

## 6. Standing step locked in (user request)
**Dimensional-consistency check (units restored, SymPy harness) BEFORE trusting any number** is now a STANDING gate:
`docs/pathA_preregistration.md` "STANDING VERIFICATION STEPS"; memory `[[feedback-dimensional-consistency-check]]`;
harness `dimensional_check.py`. Pairs with `[[feedback-transliteration-fidelity-audit]]`.

## 7. Audit trail (distilled; raw outputs regenerable, in gitignored `_scratch/`)
- `pathA_17` (validity): the "miss" = the pinned `m̂0²·S_port`; `54/5`←Burke-Thorne is correctly derived; no dropped
  4π/Jacobian. Verdict: ARTIFACT for the full claim, REAL only inside the pinned convention.
- `pathA_18` (dimensional, units restored): built the reusable harness; flagged dimensional gaps — BUT overstated
  (8→2 distinct; 4-spatial-`G` wrong). Honest caveat: dimensional analysis fixes no decimal order (factors =1 in
  natural units).
- Verification agent A (adversarial): 8→2 real issues, both order-neutral natural-units conventions, not code bugs;
  "if anything STRENGTHENS the real-missing-physics reading" (i.e. dimensions don't rescue the miss).
- Verification agent B (independent re-derivation): `[G]=L³T⁻²M⁻¹` effective-3D (brane projection, from the parent
  action); `m̂0²·S_port` is the dimensionful conversion factor `T_target/[P0]`, not dimensionless; `D0=K−B0−Z0` is a
  legal subtraction (`K`=modal stiffness eigenvalue) within the frozen normalization conventions.

## 8. pathA_19 (FOUNDATION) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh (`_scratch/pathA_19_execute.log`, 162k tokens). Deliverables: extended harness group in
`src/stage1_solver/dimensional_check.py` (+509, side-by-side), `reports/pathA_19_dimensional_foundation.md` (F5 ref
doc), `tools/pathA_19_foundation_dimensional_crosscheck.wl` (dual-engine). Python 17/17 algebraic checks pass;
Mathematica 20 checks PASS; full suite 92 passed; `git diff --check` clean. Acceptance = **PASS_WITH_NAMED_RESIDUALS**
(exit-0 NOT treated as acceptance).

**F1 mass-fork verdict — RETAIN `{L,T,M}` (honest negative result):** `m_GNLS` stays an explicit parent-action mass
`[M]` (appears in kinetic operator, current, Madelung velocity, Euler eq, sound-speed law). `m_defect` emergence is
**NOT derived** — the ontology gives drainage/volume-deficit *scaling* (`brane_bulk_ontology.tex:1267-1297`) but no
boundary source / Noether charge / Hamiltonian energy theorem tying `m_defect` to the inflow rate. `ħJ/c_γ²=M` is a
*dimensional conversion only*, explicitly NOT a derivation. So mass-as-inflow is rejected **for this gate** and carried
to pathA_21 as `INFLOW_MASS_SOURCE_MISSING` (BLOCKS_MASS_EMERGENCE) with a concrete derivation target.

**F2 flux + a-pin:** `[J_number]=T⁻¹` in BOTH 4D-bulk (closed 3-surface) and brane (2-surface) frames; `Q_vol`=`L⁴T⁻¹`
(bulk)/`L³T⁻¹` (brane); `m_GNLS J`=`MT⁻¹`. Gauss shape-independence holds ONLY with no enclosed source/leakage — but
projection makes `S_leak` and the throat bottom is open/closed/connected-undecided → `NO_NET_ACCRETION_BC_UNDERIVED`
carried. `a` confirmed a mouth-radius collective MOMENT (`a0=R0(0)`, `a(t)`=mouth average), not fundamental →
`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`; use `J` as the invariant scale-map label.

**F3 pins / healing / dictionary:** 4 pins (`a=c_s=ħ=m=1`) on 3 base dims ⟹ one null relation
**`a=ħ/(m_GNLS·c_s0)`**. GNLS core balance ⟹ **`ξ_h=√2·ħ/(m_GNLS·c_s0)`** with `h0=(5K/4)ρ0⁴=m_GNLS c_s0²/4`; so if
`a`≡healing core then `a=ξ_h/√2` (convention/branch factor). Independent set = {`ħ`, `m_GNLS`, `K`, chosen `ρ0`};
derived = {`c_s0`, `ξ_h`, `a`(if core-identified), `m_defect`(blocked)}. 4D dictionary homogeneous; the 3D GR target is
a downstream conversion problem, not a base-system change.

**F4 paper-prose reconciliation:** `part01`/`pde.tex` 4D action/EOS/current/projection AGREE with the harness
dictionary; `em_fields.tex:1717-1786` flagged **WRONG-3D-CONVENTION** (ρ₀ as kg m⁻³, pressure/enthalpy-per-mass,
`V=πa²L` throat volume) — legacy 3D/SI prose, not the 4D number-density dictionary.

**Three R_norm dimensional flags (carried, NOT repaired here):** formal 4D R_norm target is NOT dimensionless (needs
`L T² M`); observed 3D GR target needs `L² T² M`; a TRUE `{L,T}` base FAILS the R_norm gate (needs `L² T²`). These
corroborate the B2c-undetermined verdict (§1) and belong to pathA_21/pathA_22 — pathA_18 behavior preserved.

**Review:** two clean transliteration-fidelity agents (Python module; Mathematica `.wl`) → both **FIDELITY-CLEAN**;
`.wl` confirmed INDEPENDENT (no `Import`/`Get` of Python results; own representation); no tautology/can't-fail gate (the
Maxwell `c²` factor and the `{L,T}`-gate rejection both demonstrably *can* fail); side-by-side scoping CLEAN. Two
prose-only (non-machine-checked) claims — the `√2` in `ξ_h` and `h0=(5K/4)ρ0⁴=m c_s0²/4` — hand-verified correct. One
harmless dead helper (`_dim_dict`).

## 9. pathA_20 (c_s + c velocity structure) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh. Deliverables: harness group in `dimensional_check.py` (+478, side-by-side,
`--patha20-velocity`), `reports/pathA_20_velocity_constants.md`, `tools/pathA_20_velocity_constants_crosscheck.wl`.
Python 21 dim + 5 alg checks pass; Mathematica PASS; full suite 92 passed. Acceptance = `PASS_WITH_NAMED_RESIDUALS`.

**Three `UNDETERMINED` verdicts — ALL VERIFIED HONEST WALLS (adversarial premature-punt audit), not punts:**
- `C_GAMMA_RATIO_UNDERDETERMINED`: the parent gauge action gives the photon the **bulk Minkowski cone** (`Z(w)` is an
  overall prefactor on `F_{MN}F^{MN}` → cancels from the principal symbol; it renormalizes the coupling `μ0_eff`/`q_eff`,
  NOT the cone), while `c_s` is the **emergent acoustic cone** — and the sources NEVER relate them (`em_fields.tex` only
  ASSERTS `c=c_s` weak-field by fiat). Carried: `λ_γ=c_γ/c_s`, tail `(c/c_s)³=λ_γ³`. [Caveat: the report's stated reason
  over-credits "`Z(w)` unsolved"; the real reason is structural (gauge-bulk-cone vs emergent-acoustic-cone) — sharpen
  this when pathA_21 carries it forward.]
- `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`: the model-faithful throat profile (with `Q`, `V_conf`, throat
  geometry `R0`, bottom-topology BC) is explicitly unsolved (`pde.tex:2845-2849`). Codex DID derive the conditional
  ideal-Euler nozzle (`c_*/c_s0=1/√3`, `ρ_*/ρ0=3^(-1/4)`, `J_crit/(ρ0 c_s0 A_*)=3^(-3/4)`, independently verified
  correct for the `Kρ⁵` EOS) and correctly REFUSED to promote it (dropping `Q` fails at the throat). [Caveat: report
  prose overstates "ODE cannot be closed" — the classical one WAS closed; the wall is the `Q`/`V_conf`/geometry/BC profile.]
- `HBAR_PROVENANCE_UNDETERMINED`: NO `ħ`-free substrate relation anywhere in the 4 cited papers (grep-confirmed);
  circulation `Γ∈ℤ` is an IMPOSED classical input label, not a derived substrate quantum. Refusing emergence = correct
  anti-tautology; `UNDETERMINED` over `FUNDAMENTAL` is the honest conservative call.

**Clean derivations:** `[c_s]=[c_γ]=L/T`, `c_s∝ρ²` (`d ln c_s/d ln ρ=2`), `[J]=T⁻¹` (bulk+brane). **Standing-wave
`c=c_γ` GENUINELY NON-CIRCULAR** (`ω₀=c_γ k_⊥` from the trapped-mode BC; `ω₀/γ` clock from a boosted wave operator;
NO `E=mc²`/Compton premise) — a sketch (asserts boost-covariance) but non-circular. Mass-bridge recorded candidate-only;
`M` not collapsed; scope CLEAN.

**Review:** Python + `.wl` transliteration both FIDELITY-CLEAN; `.wl` confirmed INDEPENDENT. Non-blocking: the
"group-velocity" check is dimensionally vacuous (doesn't differentiate the dispersion); the sonic `1/√3` is asserted as
a literal not derived-in-code (Bernoulli is prose-only → dual-engine agreement on it isn't meaningful) — both labeled
conditional, neither changes a verdict.

**STRATEGIC FINDING:** the emergent constants CANNOT be closed by dimensional analysis + the existing source equations —
the NUMBERS (`c_γ/c_s`, the flux law, and downstream `G`/`α_J`) all need (a) the SOLVED stationary defect profile and
(b) for `c_γ/c_s`, the coupled GNLS+Maxwell linearization (does the photon ride `g_μν(c_s)` or `η_MN`?); `ħ` needs NEW
substrate microphysics. This is the SAME bottleneck as the core Path-A solver → the emergent-constants chunk and the
main solver work converge here. FORK presented to the user (§0).

## 10. pathA_20b (c_γ vs c_s coupled linearization) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh. Deliverables: harness group in `dimensional_check.py` (+486, side-by-side,
`--patha20b-cgamma-cs`), `reports/pathA_20b_cgamma_cs_linearization.md`, `tools/pathA_20b_cgamma_cs_crosscheck.wl`.
Python 11 dim + 7 alg; Mathematica PASS; full suite green. Acceptance = `PASS_WITH_NAMED_RESIDUALS`.

**Verdicts (both UNDETERMINED, VERIFIED HONEST + sharper than pathA_20 — adversarial audit):**
- `bulk_verdict = C_GAMMA_BULK_UNDERDETERMINED` (`BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`). The coupled principal
  symbol FACTORIZES — off-diagonal GNLS↔gauge couplings are lower-order London/plasma terms (`δJ⁰=q★δρ`,
  `−(q★/m)ρ0 δA^i` algebraic in `δA`) below the `∂²` cone — giving `det P = ħ(ω²−c_s²k²)·(C_E(ω²−c_bulk²k²))²`. Photon
  cone `c_bulk²=C_B/C_E` (`Z/μ0` cancels). DERIVED, verified (not asserted).
- **KEY PHYSICS: `c_γ=c_s` is NOT forced.** The parent gauge metric is bulk Minkowski `η_MN` (no `K/ρ0/m`); `c_s` is an
  emergent Bogoliubov speed from a DIFFERENT sector; the matter sector is NON-RELATIVISTIC (parabolic, 1st-order in `t`)
  so it imprints NO coordinate light cone on the shared `t` → sharing `t` does NOT lock the cones. **`c_bulk/c_s` is a
  CALIBRATABLE normalization freedom (the `η` time-vs-space ratio), not a derived number and not a permanent gap.** If
  `c_bulk` fixed, `λ_γ=c_γ/c_s ∝ ρ0⁻²`.
- `brane_verdict = C_GAMMA_RATIO_STILL_UNDERDETERMINED` (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED` +
  `BRANE_PHOTON_CONE_REQUIRES_PROFILE`) — the observed photon cone needs the zero-mode reduction / solved profile (parent
  paper declares the Maxwell/mixed spectrum OPEN: part01:1502,944; pde.tex:541-565). pathA_21 consumes this, `λ_γ` symbolic.

**Background legality:** legal only with a neutralizing external source `J_ext0⁰=−q★ρ0` (jellium; sourced
pde.tex:370-374) — handled, not punted.

**Review:** Python + `.wl` both FIDELITY-CLEAN; `.wl` INDEPENDENT (own CAS, computes the determinant natively); the
factorization + `c_bulk²=C_B/C_E` + the "unsourced `c_γ=c_s`" claim VERIFIED DERIVED; negative control BINDING (a forced
`c_γ=c_s` genuinely fails); coupled-symbol GENUINELY-COUPLED; scope CLEAN. **Non-blocking caveats (recorded):** (1) the
`lambdaLogSlope=−2` harness check is TAUTOLOGICAL in code (hardcodes `ρ0⁻²` rather than deriving it from `c_s∝ρ0²`) — the
CLAIM is hand-verified correct but the check doesn't back it; (2) the `forced-equals` check is cosmetic (binding teeth =
the residual-independence check) and the `coupled-det` check restates the block factorization → "7/7 algebraic" is NOT 7
independent confirmations; (3) the bulk residual is better framed as a CALIBRATABLE FREEDOM than a gap. [If we revisit
this harness, make `lambdaLogSlope` a real derivation.]

**STRATEGIC STATE:** pathA_19/20/20b now CONVERGE — the emergent-constants program has extracted everything derivable
from the symbolic action + dimensions. `c_γ/c_s` is a CALIBRATION KNOB; the flux law, brane cone, `G`, and `α_J` all need
the SOLVED THROAT PROFILE; `ħ`-emergence needs new substrate physics. Remaining moves: (A) `pathA_21` SYMBOLIC (complete
the unknowns/knobs map) or (C) the throat-profile solve (closes the numbers). User decision (§0).
