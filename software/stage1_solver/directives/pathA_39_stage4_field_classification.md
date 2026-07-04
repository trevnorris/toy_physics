# pathA_39 Stage 4 — field-coupling classification: assemble the linearized EM multiplet + classify from a computed DOF/constraint/residue analysis (sector close)

**Status:** ⭐ v2 — folds the Codex design-review (`NOT_SOUND` → v0 relabel; 4 blockers + 6 concerns) AND the GLM tertiary
(`SOUND_WITH_CONCERNS`; scalar-sector stability sub-gate + wording sharpeners). The genuine content is the MAIN computation:
assemble ONE combined linearized quadratic action for the multiplet from imported BLOCKS, and — using the longitudinal
constraint CLASS imported from pathA_36 (earned there, NOT re-derived here) — COMPUTE from that one system the physical DOF
consequences, the scalar-sector stability (§3.2b), `h`'s propagating status, and per-eigenmode charge residues; classify from
the computed features. **Next:** dual-engine execution (SymPy+Mathematica) → tri-review. **Parent:**
`pathA_39_magnetism_close_maxwell.md` §3/§4. **Author:** orchestrator (scaffolding: requirement + acceptance; imports action
BLOCKS + source vectors, NOT the final feature vector; does not re-derive Stages 0–3).

> **⚠️ HONEST FRAME (§6).** Exact Maxwell is EXCLUDED up front (`M_h>0`, `B_eff>0`). This gate CLASSIFIES the assembled
> linearized EM field content and is able-to-fail via ACTION-LEVEL counterfactuals whose class is DERIVED from the same
> constraint/rank/residue extractor (not assigned by the counterfactual label): decouple `h` + impose the first-class chain
> → the extractor must RETURN 2 physical DOF + a first-class generator + transverse-only charge = Maxwell structure;
> zero scalar source (charge kept transverse) → clean coexistence; the real branch → the departure. Expected =
> `FIELD_SCALAR_VECTOR_DEPARTURE` (2 transverse vector modes + a propagating charge-coupled `h`-branon; "photon" is reserved
> for the Maxwell counterfactual branch, where the first-class gauge chain is imposed). Per
> [[feedback-negative-verdict-short-circuit]], the expected departure matches our prior, so the Maxwell-recovered-when-`h`-
> decoupled and clean-recovered-at-zero-residue counterfactuals are the discrimination proofs. **This is a SPEC close
> (a classified field-content skeleton with conditional flags), NOT a solved-parameter close, and NOT the `λγ` knit.** Do not
> rescue, do not oversell.

---

## §0. The crux in one paragraph
Stages 0–3 each examined one facet. Stage 4 assembles the COMPLETE linearized EM sector as ONE quadratic system and asks the
field-theory question — **what field content is this?** Using the longitudinal constraint CLASS imported from pathA_36 (that
first-/second-class classification is earned there, NOT re-derived here), it COMPUTES from the one assembled system: the
physical propagating DOF (the Dirac count given that class + the kinetic rank), the scalar-sector STABILITY (§3.2b), whether
`h` is a physical propagating scalar, and each eigenmode's charge residue from the one inverse operator — then classifies.
The provenance-fixed inputs force a scalar-vector departure, but the extractor must be able to RETURN Maxwell / clean
coexistence under action-level counterfactuals, or it is a rubber stamp.

## §1. Imported BLOCKS + source vectors (structured values, not prose; not the feature vector) — C5
Import the quadratic-action BLOCKS and the charge/current SOURCE vector, bind each into the assembled operator (below).
Do NOT import the Stage 0–3 verdicts/DOF-counts/residues as the answer — compute those here from the assembled system.
- **Transverse (pathA_36):** `L_T = ½ρ_br u̇_T² − ½μ_R(∇×u_T)²`; `c_γ²=μ_R/ρ_br`. (2 in-brane transverse polarizations.)
- **Longitudinal/density (pathA_36 + Stage 0+1):** `L_L = ½ρ_br u̇_L² − ½B_eff k² u_L²`, `B_eff=ρ_B0²/χ_c>0`; AND the
  pathA_36 constraint/first-class data that distinguishes the finite-compressibility (second-class physical) branch from the
  tuned first-class-Maxwell locus (import the constraint CLASS, not a DOF number).
- **Branon `h` (pathA_38):** `L_h = ½M_h ḣ² − ½K_h k² h²`, `K_h=M_h c_E²`, `M_h>0`; `c_E` = the dynamic-Green speed.
- **`h`–`u_L` mixing (Stage 0+1):** off-diagonal `−C_hu k² u_L h` (the `D_x` block). `C_hu` is a **sim-deferred mixing
  coefficient** (Stage 0+1 computed residues only in the `C_hu→0` decoupled limit) — so the "healthy scalar spectrum" is NOT
  automatic; it is a COMPUTED condition on `C_hu` (see §3.2b, the stability sub-gate).
- **Source/current vector `J_q`** (charge couplings, the objects the residue computation uses): `q_A^T` (→ `u_T`),
  `q_L=Nu·a_L·s` (→ `u_L`, `a_L` sim-deferred), `q_h=2Q_E tanh(b/ℓ)/b` (→ `h`), and the mass source `q_M` (→ `u_L`, for the
  clean/mass channels). **SIGN source of truth = `reports/pathA_39_magnetic_force.md`** (both transverse and longitudinal
  like-current channels ATTRACT; the stale `_scratch/…_resolution_v5.md §4` "longitudinal repels/`q_L²/q_A^T²` crossover"
  line is SUPERSEDED — do not import it).
- **Stage 3 (operator parity):** `FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING`. Enters as a **sim-deferred contamination FLAG**
  (`operator_parity_contamination = true, magnitude sim-deferred`), NOT as an added propagating DOF (C3). Do not let it
  change the DOF count.

## §2. Setup — assemble ONE combined quadratic operator (the genuine-gate core, B1)
Build a single inverse-propagator matrix `Q(ω,k)` over the named canonical coordinates
`Φ = (u_T1, u_T2, u_L, h)` (name them explicitly so `h`, `u_L`, `u_T` are not accidentally overlapping descriptions of the
same mode — physics-fidelity note), by placing the §1 blocks: transverse 2×2 (decoupled), the `(u_L,h)` 2×2 scalar block
with the `−C_hu k²` mixing, and the charge/current source vector `J_q(k)` coupling `Φ` to the external charge. `Q` is the
one object; all features below come FROM it (and its constraint structure), never from an imported answer.

## §3. The executable gate (dual-engine; features COMPUTED from the one system; able-to-fail via action-level counterfactuals)
1. **Assemble** `Q(ω,k)` + `J_q` from the §1 blocks (bind every imported value into `Q`/`J_q`).
2. **Physical DOF (B1/B4):** compute the Hessian/kinetic rank + the Dirac constraint count → the number of physical
   propagating DOF, per-sector: transverse (expect 2), longitudinal `u_L`, `h` (physical propagating vs absent/decoupled).
   **The constraint CLASS (first- vs second-class) for the longitudinal sector is IMPORTED from pathA_36 (earned there), NOT
   re-derived here** (GLM concern 2/4): Stage 4 computes the DOF CONSEQUENCES of that imported class + the kinetic rank, and
   the residues (§3.4) — the able-to-fail content is at the extractor/counterfactual level, not a re-classification of
   constraints. On the real branch `u_L` is second-class-physical (finite compressibility); under the Maxwell counterfactual
   the first-class chain is imposed (§5).
2b. **Scalar-sector STABILITY sub-gate (GLM concern 1 — the one place the computed DOF count can fail the builder).** The
   mixed `(u_L, h)` potential is `V = [[B_eff k², C_hu k²],[C_hu k², K_h k²]]` (`K_h=M_h c_E²`). The two scalar eigenmodes are
   BOTH healthy propagating DOF **iff** `det V/k⁴ = B_eff K_h − C_hu² > 0` AND both eigenvalues `> 0`. COMPUTE this
   (`scalar_sector_stable`) on the finite-`k` stiffness matrix — the kinetic block `diag(ρ_br, M_h)` stays positive, so a
   violation is a GRADIENT/tachyonic instability, NOT a kinetic ghost; `det=0` is the marginal/degenerate edge. If violated
   (`C_hu² ≥ B_eff K_h`) one scalar eigenmode is gradient-unstable → the "4 healthy DOF" claim FAILS → emit
   `FIELD_SCALAR_SECTOR_UNSTABLE` (§4). Since `C_hu` is sim-deferred, the real-branch healthy-departure verdict is CONDITIONAL
   on the stability bound `C_hu² < B_eff K_h`, reported as a computed condition (not assumed).
3. **Poles/speeds:** `det Q(ω,k)=0` → the propagating-mode speeds (transverse `c_γ`, scalar eigenpoles of the `(u_L,h)`
   block, incl. `c_E` and the density `c_s`-branch).
4. **Charge residues from the SAME inverse (B1):** compute `R = J_q† Q⁻¹ J_q` and its residue at each physical pole — the
   charge coupling to each eigenmode (NOT copied from Stage 0+1). Record which poles carry nonzero CHARGE residue beyond the
   2 transverse vector modes.
5. **Classify** from the computed features {physical DOF; longitudinal class; `h` status; per-pole charge residues} per §4
   (deterministic schema).
6. **Counterfactual branches (B2/B4):** run each control (§5) through the SAME assembly + extractor, each as an ACTION-LEVEL
   modification (e.g. Maxwell = REMOVE the `h` block + `q_h` from the action AND impose the pathA_36 first-class chain — NOT
   `M_h=0` substituted into a residue, which is singular). The class must be DERIVED from the branch's own
   constraint/rank/residue output.
7. **State deferrals:** `c_E=c_γ` (Lorentz/`E–B` form) and `c_γ=c_s` (`λγ`) are NOT tested — the separate knit.

## §4. Verdict taxonomy — deterministic schema (primary + flags; B3)
Emit a `primary` class + explicit `flags`; the same feature vector maps to exactly one report (no double-reporting).
- **`primary`** (from computed DOF/constraint/`h`-status/transverse-only-charge):
  - `FIELD_EXACT_MAXWELL_STRUCTURE` — 2 physical transverse DOF, a first-class (gauge) generator present, no propagating
    physical `h`, no charged longitudinal/density pole, charge residue ONLY in the transverse sector. **Named "…STRUCTURE"
    (DOF/gauge level), NOT full Maxwell** — the `c_E=c_γ`/`λγ` Lorentz-cone is deferred (C6). KNOWN UNREACHABLE on the real
    branch (`M_h>0`); reachable only under the action-level Maxwell counterfactual = the primary able-to-fail target.
  - `FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY` — 2 transverse + a longitudinal `B_eff>0` density/gravity DOF that
    carries ZERO charge residue, and (if `h` propagates) `h` also carries zero charge residue. Reachable only when all
    scalar/density charge residues vanish (`q_h=0 ∧ q_L=0`, `q_A^T≠0`).
  - `FIELD_SCALAR_VECTOR_DEPARTURE` (EXPECTED, conditional on `scalar_sector_stable`) — 2 transverse VECTOR modes
    (physical polarizations, NOT gauge "photons" on this branch — no first-class chain here, GLM concern 3) + a PROPAGATING
    physical `h` carrying a nonzero charge residue (robust `h`-floor: `M_h>0`, `q_h≠0`, independent of `a_L`).
  - `FIELD_SCALAR_SECTOR_UNSTABLE` (GLM concern 1) — the mixed `(u_L,h)` stiffness is NOT positive-definite
    (`C_hu² ≥ B_eff K_h`): one scalar eigenmode is gradient-unstable (tachyonic; kinetic block stays positive) → the healthy
    4-DOF classification FAILS. A qualitatively worse departure (or a sign the linearized description / imported `C_hu` regime
    breaks down). Reachable via the large-`C_hu` control (§5) and, on the real branch, if the sim-deferred `C_hu` violates the
    bound.
- **`flags`:**
  - `scalar_sector_stable` — the COMPUTED positivity `B_eff K_h − C_hu² > 0 ∧ both eigenvalues > 0`; the real-branch healthy
    departure is CONDITIONAL on it (`C_hu` sim-deferred). Expected true, but computed — not assumed.
  - `density_charge_coupled` — true iff the density (`c_s`) pole carries a nonzero charge residue (the `q_L∝a_L` admixture;
    rests on the sim-deferred `a_L≠0`). (This is v0's `FIELD_CHARGE_COUPLED_DENSITY_DEPARTURE`, now a flag on the departure
    primary, not a competing primary — resolves the B3 ambiguity.)
  - `operator_parity_contamination` — true (sim-deferred magnitude), from Stage 3.
- **`FIELD_CLASSIFICATION_UNDERDETERMINED`** (integrity override): the assembled system is inconsistent with a source block
  (import-fidelity fails) → no verdict until reconciled.

## §5. Controls — each RUN through the same assembly+extractor; each records {assembled Q, constraint/rank, poles, J_q, residues, class} (B4)
| control | ACTION-LEVEL modification | derived expectation |
|---|---|---|
| **real provenance-fixed branch** | `M_h>0`, `B_eff>0`, `J_q` as earned | **IF `scalar_sector_stable` (`C_hu²<B_eff K_h`):** 4 physical DOF (2 `u_T` + `h` + density); `primary=FIELD_SCALAR_VECTOR_DEPARTURE`; `density_charge_coupled=(a_L≠0)`; `operator_parity_contamination=true`. **ELSE** (`C_hu` violates the bound): `primary=FIELD_SCALAR_SECTOR_UNSTABLE` overrides |
| **Maxwell counterfactual (B2)** | REMOVE the `h` block + `q_h` from the action; impose the pathA_36 first-class Maxwell chain on `u_L` | extractor must DERIVE: 2 physical DOF, a first-class generator, no charged scalar/longitudinal pole, charge residue only transverse → `FIELD_EXACT_MAXWELL_STRUCTURE` (the primary discrimination proof) |
| **clean-coexistence (C2)** | `q_h=0 ∧ q_L=0` (zero scalar/longitudinal charge source) but `q_A^T≠0` (transverse charge kept), `M_h>0`, `B_eff>0` | diagonalize the `(u_L,h)` block with zero scalar source, VERIFY zero scalar/density charge residues (do not bypass the residue computation) → `FIELD_TRANSVERSE_EM_PLUS_CLEAN_GRAVITY_DENSITY` |
| **`a_L→0` (h-floor robustness)** | `q_L=0`, keep `M_h>0`, `q_h≠0`, `q_A^T≠0` | `density_charge_coupled` drops to false; `primary=FIELD_SCALAR_VECTOR_DEPARTURE` stands — the ROBUSTNESS proof that the `h`-floor survives losing the density admixture (distinct from the clean row's coexistence proof) |
| **large-`C_hu` (stability discriminator, GLM 1)** | set `C_hu² ≥ B_eff K_h` (strong scalar mixing) | `scalar_sector_stable=false` → must emit `FIELD_SCALAR_SECTOR_UNSTABLE` (the healthy-DOF count fails); restoring `C_hu²<B_eff K_h` recovers the departure — the able-to-fail-the-builder proof |
| **import-fidelity (C5)** | corrupt one bound imported value (`c_E≠` pathA_38 Green speed; `B_eff≠ρ_B0²/χ_c`; `K_h≠M_h c_E²`) | `FIELD_CLASSIFICATION_UNDERDETERMINED` (not a silent classify) |
| **DOF-count discriminator** | (main gate) count physical DOF from the constraint/rank on each branch | 4 on the real branch, 2 under the Maxwell counterfactual — DERIVED, the discriminator |

## §6. Honest scope + expected landing
- **In-scope:** assemble the one quadratic system, COMPUTE DOF/constraint/`h`-status/residues, classify, run action-level
  counterfactuals. NO nonlinear throat solve; NO new force/mode derivation (import blocks only).
- **Sim-deferred / conditional:** `a_L` (hence `density_charge_coupled`); the Stage-3 contamination MAGNITUDE
  (`operator_parity_contamination` is a flag); the `c_E=c_γ` (Lorentz form) and `c_γ=c_s` (`λγ`) — NOT tested (the knit).
- **Expected landing:** `primary=FIELD_SCALAR_VECTOR_DEPARTURE`, `density_charge_coupled=(a_L≠0)`,
  `operator_parity_contamination=true`. The EM sector is a transverse-vector + propagating charge-coupled-scalar multiplet —
  NOT exact Maxwell — the coherent field-content statement of the Stage 0–3 departures. **First-class for a unified toy
  analog; a clean coexistence would be stronger-than-expected. Do not rescue, do not oversell.**
- **This is a SPEC close, not a solved-parameter close (C4):** a classified field-content SKELETON with conditional
  density/parity flags; the `h`-floor is robust, the density departure and the contamination magnitude are sim-deferred.
- **Does NOT:** un-earn Stages 0–3; do the `λγ`/`c_E=c_γ` knit; claim full Lorentz-Maxwell (only the DOF/gauge STRUCTURE
  class is defined here).
- **Sector close:** the 4th/last magnetism stage. After it the sector is complete as a spec (a real fourth force + a
  characterized scalar departure); the program moves to the consistency knit ([[project-brane-existence-defect-structure]])
  and the post-Stage-4 "extra scalar" follow-up (task #110, `notes/pathA_39_charge_coupled_scalar_followup.md`).

## §7. Integrity (mandatory)
Dual-engine SymPy + Mathematica `ENGINE_AGREE` on the intermediate COMPUTED FEATURE VECTOR (assembled `Q`, constraint/rank,
DOF per sector, poles, `J_q† Q⁻¹ J_q` residues) — NOT only the final enum. Import-fidelity: bind every §1 value into `Q`/`J_q`
and check against the source reports (corrupt→`UNDERDETERMINED`). Every class empirically reachable via its action-level
branch — the classifier MUST EMIT it from the branch's own computed features, not raise, not read the label
([[feedback-negative-verdict-short-circuit]]): the Maxwell-counterfactual DOF-drop-to-2 + first-class generator, and the
zero-scalar-residue clean branch, are the discrimination proofs. Dimensional firewall on the assembled action. No target
readback (the class is a pure function of the computed features).

## §8. References (read first)
- Codex review this folds: `_scratch/pathA_39_s4_dirreview.md`. Parent: `directives/pathA_39_magnetism_close_maxwell.md`
  §3/§4; resolution `_scratch/pathA_39_s2_resolution_v5.md` §4 (⚠️ sign line superseded — see §1).
- Imports: `reports/pathA_39_scalar_admixture_screen.md` (Stage 0+1: `(u_L,h)` block, `C_hu`, residues),
  `reports/pathA_39_magnetic_force.md` (Stage 2: signs = source of truth), `reports/pathA_39_stage3_operator_parity.md`
  (Stage 3 flag), `reports/pathA_36_c5_phase_potential*` (transverse + longitudinal constraint class, `B_eff`, `c_γ`),
  `reports/pathA_38_results.yaml` (`M_h`, `c_E`, `q_h`).
- `notes/pathA_39_charge_coupled_scalar_followup.md` (task #110).
- Memories: `[[project-brane-existence-defect-structure]]`, `[[feedback-native-em-mechanisms]]`,
  `[[feedback-negative-verdict-short-circuit]]`, `[[feedback-dual-engine-required]]`.

## §9. Changelog
- **v2 (2026-07-04)** — folds the GLM tertiary review (`SOUND_WITH_CONCERNS`; legitimate sector-close, no kill-level
  defects). MATERIAL fix: added the scalar-sector STABILITY sub-gate (§3.2b) — the mixed `(u_L,h)` potential positivity
  `B_eff K_h − C_hu² > 0` is now COMPUTED (`scalar_sector_stable`), with a new `FIELD_SCALAR_SECTOR_UNSTABLE` primary + a
  large-`C_hu` control; the real-branch healthy departure is conditional on this bound (`C_hu` sim-deferred) — the one place
  the computed DOF count can fail the builder (GLM 1). Wording: the longitudinal constraint CLASS is IMPORTED from pathA_36,
  Stage 4 computes only its DOF consequences + residues (GLM 2/4); "photon" reserved for the Maxwell-counterfactual branch,
  real-branch = "transverse vector modes" (GLM 3); `a_L→0` annotated as the `h`-floor robustness proof vs clean's coexistence
  proof (GLM 5).
- **v1 (2026-07-04)** — folds the Codex design-review (`NOT_SOUND`): promotes the combined-quadratic-action assembly +
  Dirac DOF/constraint computation + `J_q† Q⁻¹ J_q` residues from a control row into the MAIN gate (B1 — v0 was a relabel);
  Maxwell counterfactual is now an ACTION-LEVEL branch (decouple `h` + impose the first-class chain, DERIVE 2 DOF), not a
  singular `M_h=0` substitution (B2); deterministic `primary`+`flags` schema, `density_charge_coupled` as a flag on the
  departure primary (B3); controls run through the same feature extractor recording all artifacts (B4); `FIELD_EXACT_MAXWELL`
  → `FIELD_EXACT_MAXWELL_STRUCTURE` (DOF/gauge only, `λγ` deferred, C6); clean control keeps `q_A^T≠0` + verifies zero scalar
  residues (C2); Stage-3 mixing = sim-deferred flag not a DOF (C3); "spec close not solved-parameter close" (C4); import
  structured values not prose (C5); reconciled the stale resolution-v5 sign line (C1). Dual-engine on the computed feature
  vector, not just the enum.
- **v0 (2026-07-04)** — DRAFT (relabel-risk; superseded).
