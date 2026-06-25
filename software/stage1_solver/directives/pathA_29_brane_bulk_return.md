# pathA_29 EXECUTION DIRECTIVE — Brane↔Bulk Return: Cancellation EXISTENCE + 1/r² Survival through the Finite Slab (TRACK 3, Gate 1; closes pde_ledger open-item #9 + tests the gravity-range item)

**STATUS: ✅ DONE + VERIFIED (2026-06-25) — `RETURN_RESIDUAL_PREDICTION`. Directive v3 (Codex+GLM SOUND); EXECUTED (4th execution, v3b)
and tri-review-VERIFIED (orchestrator arbiter re-run reproduced + adversarial CLEAN + fidelity FAITHFUL). Result: given the drain
premise (`Z<0`), 1/r² gravity SURVIVES the finite slab (both DC-sink completions → normalizable `m=0` zero mode → `p=2` via a real
3D-radial `dsolve`, counterfactual-guarded) + the drain comes bundled with a bounded monopole/dipole `c_s`-radiation residual
`∝ ε0=1−𝒯₀(0)` (the falsifiable GR departure); NOGO reachable via a derived delocalizing warp (`p=3`). Sharpens (does NOT close)
pde_ledger #9; gravity-range item passes. Engines + reports under `tools/`+`reports/pathA_29_*`. NEXT (downstream track-3) = the full
nonlinear brane↔bulk return closure. [Design-review history: Codex re-confirm GREEN v3-r2, xhigh; GLM tertiary folded.]
NEXT (historical): orchestrator final read → re-execute → tri-review.
v3 supersedes the v2 §0.5.1 "both branches agree on (p=2 AND Z<0)" gate after GLM tertiary (user-run) showed it was STRUCTURALLY BROKEN
(the radiation/Sommerfeld BC degenerates to Neumann/REFLECTING at `ω→0` ⇒ deterministic `Z=0` at DC ⇒ the favorable verdict was
unreachable = postulate-to-win relocated into an always-fail cross-check). §0.6 reframes: `Z<0` is the drain ADMISSIBILITY PREMISE (a
non-draining throat is not a particle); Check B runs under the admissible DC-SINK completions (§5 de-structuring/absorbing PRIMARY +
Bloch CROSS-CHECK); radiation/Sommerfeld is AC/Check-A only; the gate decides on the FALLOFF `p` (RESIDUAL_PREDICTION if `p=2` /
NOGO if `p≠2` via a delocalizing warp / BC_DEPENDENT if `p` differs across DC-sink completions). v2 reqs #2 (p from a real radial solve,
in `engine_agreement`) and #3 (two-independent-computation `static_dynamic_consistency`) REMAIN IN FORCE. TWO prior executions REJECTED
at tri-review (2026-06-25): v1 hardcoded Check B's 1/r² (no
Green's-function solve); the remediation genuinely solved the transverse eigenproblem BUT (i) `static_dynamic_consistency` was STILL a
variable-vs-itself tautology, (ii) `p=2` was a spectral-class lookup / bare literal (not derived, not in `engine_agreement`), and
(iii) — the CRUX — 1/r² survival hinged on a POSTULATED `m`-dependent impedance BC chosen so the zero mode survives. NEITHER prior verdict was
verified.** Roles: **Codex codes + runs + iterates dual-engine; Claude reviews.** This states requirements +
acceptance criteria only — **Codex designs the derivation/computation route** ([[feedback-claude-reviews-codex-codes]]).

**GOAL.** pathA_28 left a *condition* the brane↔bulk return must satisfy to avoid GR-forbidden monopole/dipole gravitational
radiation: `R0(ω)=−M0(ω)` (ℓ=0) and `R1_i(ω)=−D1_i(ω)` (ℓ=1). This directive does the **track-3** job that pathA_28 explicitly
deferred: **postulate a concrete multi-brane SLAB structure (freely) but COMPUTE its return response from medium continuity (no free
knob), and test — able-to-fail — whether that admissible return can actually deliver the cancellation, AND whether long-range `1/r²`
(3+1D Newtonian) gravity survives the same finite slab.** One postulated structure must satisfy BOTH requirements at once, or it is a
**no-go** (falsification). This is the analog methodology: postulate structure freely, test *consistency*, hunt a no-go
([[feedback-analog-find-consistent-structure]]).

## §0. Anchor (do not re-derive existing work; do not drift)

- **The two requirements this gate must reconcile in ONE slab structure:**
  - **(A) Dynamic cancellation.** The throat is a time-varying mass SINK on the brane. Its net far-field ℓ=0/ℓ=1 source is
    (leak − return). pathA_28 derived the raw outgoing `c_s`-wave amplitudes `A0_raw=i·a·(ω/c_s)·M0` (`p_raw=1`) and
    `A1_raw=i·a³(ω/c_s)³·D1/2` (`p_raw=3`), with the quadrupole `A2=i·a⁵(ω/c_s)⁵·Q2/27` (`p_raw=5`) the only GR-allowed channel.
    The return must drive the residual ℓ=0/ℓ=1 outgoing amplitude **to a leading ω-power of the outgoing amplitude coefficient ≥ the
    quadrupole's `p=5`** (a STRICT structural order-drop), OR — falling short of that — at most to a bounded, quantified, sub-leading
    residual recorded as a falsifiable prediction (the `RETURN_RESIDUAL_PREDICTION` rung, §D). Source moments: `M0=∫_brane S_leak d³x`;
    `D1_i=∫_brane x_i S_leak d³x + ∫_brane j_i d³x` (incl. carried odd wake — global momentum conservation alone is NOT a kill).
  - **(B) Static 1/r² survival.** The steady (`ω→0`) drain-flow `v_r` around a throat must fall as `1/r²` (3+1D Newtonian) at
    distances `r ≫ d` (the inter-brane spacing). The classic braneworld hazard: a finite bulk slab can turn long-range gravity into
    4D `1/r³` (flow leaks into the bulk) or a Yukawa cutoff (finite range). Survival is **not** automatic — it is the
    zero-mode-localization / confinement question, here posed as the long-range falloff of the steady slab flow.
  - **THE TENSION (the whole point — must be able to fire as a no-go):** (A) pulls the slab walls toward **absorbing/returning**
    (so the leaked medium comes back and cancels the wave); (B) pulls them toward **confining** (so the quasi-static flow can't leak into
    the bulk and stays `1/r²`). The *same* boundary condition must serve both. If no admissible slab structure does both → no-go.
- **NOTHING is static (the project's #1 recurring failure mode — do NOT relapse).** The `1/r²` field is the field STRENGTH of a
  quasi-static DYNAMICAL balance, NOT a frozen BVP. Every prior static simplification failed by deleting the coupling that selects the
  balance (throat self-selection failed under conservative-static treatment — "requires dynamics"). Therefore: **the background that
  Check A linearizes around MUST be a genuine DYNAMICAL STEADY STATE** — the steady brane↔bulk CIRCULATION (leak balanced by
  return/replenishment per `conceptual_foundation.md` §5), NOT an assumed frozen profile. **Steady-state EXISTENCE is a stated
  precondition — THREE distinct cases (no double-routing, see §D):** (i) a SPECIFIC candidate structure with COMPUTED no steady
  circulation → that structure is INADMISSIBLE (dropped from the family; not itself a verdict); (ii) COMPUTED no steady circulation for
  ALL admissible structures → `RETURN_NOGO` (hard failure of the return mechanism); (iii) whether a steady circulation exists is
  UNDECIDABLE in-scope (intractable) → `INCONCLUSIVE`. Emit `steady_state: {exists: bool|undecidable, type: circulation, provenance:
  <derivation-or-obstruction>}` (for `exists: undecidable`, `provenance` MUST carry the in-scope obstruction). Check B is the `ω→0`
  LIMIT of the SAME dynamic object as Check A — never an independently-posed frozen solve (§B.4, §C, §E).
- **REUSE, don't rebuild:**
  - `software/stage1_solver/reports/pathA_28_cancellation_condition.yaml` + `pathA_28_monopole_results.yaml` — the cancellation
    target (`R0=−M0`, `R1=−D1`), the raw amplitudes/kernels, the source-moment definitions (`M0`, `D1`). REUSE verbatim; do not
    re-derive the raw amplitudes.
  - `research/4d_2_5pn/` — the outgoing/DtN multipole machinery (`Λ_ℓ`, `Y_ℓ`, the `i·zⁿ` radiation phase) and the ℓ=2
    Burke–Thorne anchor. REUSE as the `quadrupole-survives` consistency anchor; do NOT re-derive ℓ=2.
  - `research/pde_ledger/` Part VIII — the open-system **projected continuity** `∂_t ρ_brane + ∇₃·j_brane = S_leak`; global
    mass/momentum conservation across brane+bulk. **`𝒯_ℓ(0)` is NOT anchored to 1 by fiat — it is a COMPUTED output** (see §B.2):
    `𝒯₀(0)` is the monopole DC return *fraction* to OUR brane from zeroth-moment (mass) continuity; it equals 1 ONLY if all leaked
    medium returns locally to our brane and MAY be `<1` (medium transmits to the neighboring brane, sustaining a net in-brane sink).
    `𝒯₁(0)` is derived from FIRST-moment continuity (bulk flux + return centroid); NO conservation law forces it to 1.
  - `docs/conceptual_foundation.md` §5 + `[[project-model-mechanics-corrections]]` — the multi-brane / RS-like stack: `±w` puncture
    direction = which adjacent slab = charge sign; **inter-brane spacing sets gravity's coupling/range**; the neighboring brane
    provides the RETURN (re-absorption) → steady circulation; conservation is global (de-structuring, not annihilation). NOTHING is
    static; `c_s ∝ ρ²`; gravity = the FLOW, changes propagate at `c_s`; `v_r` = field STRENGTH not a speed.
- **Frozen slice (match pathA_28 exactly):** `n=5`, `P=Kρ⁵`, `c_s²=5Kρ⁴/m`, `K_eos=1/500`, `G=c=c_s=1`; throat baseline
  `(a*,L*)=(4731/2500, 18121/10000)`. The slab GEOMETRY (spacing `d`, warp) is postulated; the neighbor-brane BC is **DERIVED** (v2 re-scope §0.5), NOT a postulated wall impedance.

## §0.5 — v2 RE-SCOPE (DECISIVE BC DERIVATION; supersedes the v1 "postulate the BC form" approach where they conflict)

**Why:** two executions were REJECTED at tri-review. The decisive defect: 1/r² survival hinged on a POSTULATED `m`-dependent impedance BC `f'(d)+(1/κ0)·m·f(d)=0` — the `m`-factor was exactly what let the constant zero mode survive (impedance → 0 at `m=0`). The outcome was thus SELECTED via the "structure," not computed (a postulate-to-win, the pathA_27 death). Under a standard `m`-independent impedance there is no zero mode → Yukawa → `RETURN_NOGO`. v2 makes the gate genuinely DECISIVE by DERIVING the BC, and fixes two control/derivation defects. These three requirements OVERRIDE any conflicting v1 text:

1. **DERIVE the neighbor-brane BC — do NOT postulate its functional form (the "no free knob" principle applied to the BC itself).** [**PARTIALLY SUPERSEDED BY §0.6** — the "RUN BOTH physical boundaries and require both to agree on `(p=2 AND Z<0)`" design below is REPLACED: `Z<0` is the drain PREMISE (admissibility), not a co-equal test, and the radiation/Sommerfeld boundary is AC-only (Check A), NOT a Check-B Z-branch. What REMAINS in force from this item: the BC is DERIVED from a declared physical completion, never a postulated/tuned impedance form; and the Check-B falloff `p` is compared across admissible **DC-sink** completions (§5 de-structuring/absorbing primary + Bloch cross-check) per §0.6.] The slab GEOMETRY (spacing `d`, stack, warp) is still postulated; but the neighbor-brane boundary condition MUST be a PHYSICAL OPEN BOUNDARY DERIVED from a declared physical completion, NOT an ad-hoc impedance whose `m`-/`ω`-dependence is chosen to preserve (or kill) the zero mode. **The postulate-to-win must NOT move from impedance-FORM to BC-CHOICE:** do NOT pick whichever boundary gives the favorable outcome.
   - **DC-sink return (§5 de-structuring / absorbing) — PRIMARY (Check B).** Per §0.6.2, the model's return is a DC absorber (de-structuring into bulk and/or re-absorption by the neighbor brane — both genuine DC sinks for our brane). DERIVE it as the boundary/bulk term that transmits the static flux (`Z<0` by premise). This is the primary Check-B completion.
   - **Periodic stack (Bloch) — DC-sink CROSS-CHECK (Check B).** The model's declared completion is the multi-brane STACK (`conceptual_foundation.md` §5); the periodic/Bloch boundary is the DC-sink cross-check (also transmits → `Z<0`). DERIVE it from stack periodicity (not postulated). Holds whether the stack is finite (RS1-like) or infinite, AS LONG AS the neighbor ABSORBS (§0.6.2).
   - **Radiation / Sommerfeld — AC-only (Check A), NOT a Check-B Z-branch (§0.6.3).** `∂_wΦ = i(ω/c_s)Φ` correctly models the finite-ω outgoing residual (Check A) but degenerates to Neumann (REFLECTING, `Z=0`) at `ω→0`, mis-modelling the static drain — so it MUST NOT be used as a Check-B falloff/Z branch (tag `ac_check_a_only`).
   FORBIDDEN: postulating an impedance whose `m`/`ω` form is chosen so the zero mode survives or dies. **Compare the Check-B falloff `p` across the admissible DC-sink completions** (de-structuring vs Bloch) and emit their verdicts with COMPUTED `p`, `Z_sign` (premise `<0`), and `spectrum_class` (§0.5.4). Per §0.6, the headline is `RETURN_RESIDUAL_PREDICTION` if `p=2`; `BC_DEPENDENT` ONLY if the FALLOFF `p` differs across DC-sink completions (NOT a Z-sign disagreement — the superseded v2 Neumann-limit artifact). The DC-sink admissibility requirement holds per completion (the boundary must TRANSMIT the steady drain; a DC-reflecting boundary = no-drain/no-particle = INADMISSIBLE for Check B, §0.6.1).

2. **`p` (long-range falloff exponent) MUST be DERIVED from an actual radial solve, and cross-engine-verified.** Solve the radial equation for the zero mode (transverse eigenvalue `m=0` → 3D brane radial Laplacian with the throat source → potential `G0(r)~1/r` → `v_r~1/r²`) and for a massive mode (`m_n` → `∇²₃g − m_n²g = 0` → Yukawa `e^{−m_n r}/r`); extract `p` from the SOLVED radial r-dependence. FORBIDDEN: a spectral-class lookup table, a hardcoded `r^{−n}`, a `power_law_exponent()` on a hand-written expression, or a bare `pStatic=2` literal. **`p` (and `p_eq_2`) MUST be in the `engine_agreement` cross-checked set** (both engines derive it independently and agree).

3. **`static_dynamic_consistency` MUST compare TWO INDEPENDENT computations** (REJECTED twice as variable-vs-itself): (i) the `ω→0` limit of the ACTUAL dynamic Green's function `G(ω,r)` — its radial falloff exponent extracted as `ω→0`; vs (ii) the ACTUAL static (`ω=0`) solve's exponent. The two exponents MUST come from two genuinely separate computations (not the same value read twice); assert they AGREE, else `INCONCLUSIVE` (a coupling was deleted).

4. **CLASSIFY the COMPUTED transverse spectrum explicitly — do NOT conflate the failure modes** (`p=3` continuum-leakage and Yukawa gapped-cutoff are DISTINCT). For EACH admissible DC-sink Check-B completion (§0.6: de-structuring/absorbing PRIMARY + Bloch CROSS-CHECK — NOT the AC-only radiation branch) emit `spectrum_class: {normalizable_zero_mode | continuum_p3 | gapped_yukawa}` from the COMPUTED transverse eigenspectrum:
   - **normalizable zero mode present** → `p=2` (B-pass given the `Z<0` drain PREMISE, §0.6.1);
   - **gapless continuum** (non-normalizable zero mode / leakage into the bulk) → `v_r ~ r^{-3}`, `p=3` (B-FAIL input);
   - **gapped spectrum** (mass gap, no zero mode) → Yukawa exponential cutoff (B-FAIL input).
   Both `continuum_p3` and `gapped_yukawa` are B-fail inputs, but the report MUST state WHICH was derived (per branch).

All other v1.1 requirements remain in force (signed `Z`, `B_pass ⇔ (p=2 AND Z<0)`, the 4-rung ladder + precedence, `no_go_reachable` on a real different admissible family, `tension_is_real` witnessed, quadrupole-survives, dual-engine + units, the §C anti-tautology firewall). The §C firewall is EXTENDED: the neighbor-brane BC's functional form is itself subject to the no-free-knob rule (derive it; the radiation/Sommerfeld and periodic-stack open boundaries are the allowed physical forms; an `m`/`ω`-impedance chosen to land the zero mode is FORBIDDEN). **[SUPERSEDED IN PART BY §0.6: the "both branches must agree on (p=2 AND Z<0)" / `verdict_bloch ∧ verdict_radiation` design of §0.5.1 is REPLACED — `Z<0` is the drain PREMISE, not a co-equal test; see §0.6.]**

## §0.6 — v3 RE-SCOPE (Z<0 IS THE DRAIN PREMISE, NOT A TEST; supersedes the v2 §0.5.1 "both branches agree on Z" design)

**Why:** GLM tertiary review showed the v2 "both branches must give (p=2 AND Z<0)" gate is structurally broken — the radiation/Sommerfeld BC degenerates to Neumann (REFLECTING) at `ω→0`, so it deterministically gives `Z=0` at DC (no static loss → no inward gravity), making the favorable verdict unreachable and relocating the postulate-to-win into a cross-check branch that always fails. A radiation boundary is open for AC but CLOSED for DC, and Check B lives at the `ω→0` static limit. Reframing (grounded in conceptual_foundation §5 + [[project-model-mechanics-corrections]]):

1. **`Z<0` is the ADMISSIBILITY PREMISE, not a co-equal test.** Gravity IS the drain: a particle is a mass SINK; gravity is the inflow toward it (`Z<0` by definition). The question is NOT "does the drain exist" (it does, by premise — a non-draining throat is not a particle) but "given the drain, does its static flow fall as 1/r²?". An admissible return MUST be a **DC sink** (transmits the steady drain). A DC-REFLECTING boundary (`Z=0`, return-to-self) is the degenerate **no-drain / no-particle** case — INADMISSIBLE for Check B (the `ω→0` static-drain limit), NOT a Z-veto branch.

2. **The model's return is a DC absorber.** Per §5, draining = **de-structuring** (medium → bulk ground state) and/or **re-absorption by the neighbor brane** — both are genuine DC sinks for OUR brane (medium leaves and does not return *to our brane*; conservation is global). Model the return PRIMARILY as this **§5 de-structuring / absorbing DC sink** (a boundary/bulk term that transmits the static flux, `Z<0`); the **periodic-stack (Bloch)** boundary is the DC-sink CROSS-CHECK (also transmits → `Z<0`). NOTE (GLM Finding 3): this holds whether the stack is finite (RS1-like) or infinite, AS LONG AS the neighbor ABSORBS (§5 says it does); only a REFLECTING neighbor gives `Z=0` (= no drain = no particle), which is out-of-premise.

3. **The radiation/Sommerfeld BC is AC-appropriate ONLY (Check A), NOT a Check-B Z-branch.** It correctly models outgoing-wave radiation for the finite-ω residual (Check A), but at `ω→0` it is reflecting (`Z=0`) and mis-models the static drain — it MUST NOT veto Z in Check B. (May be reported for Check A; never used as a Check-B falloff/Z branch.)

4. **`p=2` is robust-pass-by-geometry; the real Check-B teeth are the falloff under a delocalizing geometry.** For a flat slab or a localizing (RS-like) warp, the transverse zero mode is normalizable → `p=2` (1/r²) for ANY admissible DC-sink return (bounded `w` → 3D far-field). The ONLY route to `p≠2` is a **delocalizing (anti-localizing) warp** failing to localize the zero mode → continuum (`p=3`) or gap (Yukawa) → `RETURN_NOGO` (the conceptual_foundation §7 braneworld caveat). `BC_DEPENDENT` fires only if the FALLOFF `p` differs across admissible DC-sink completions (de-structuring vs stack) — NOT on a Z-sign disagreement (the superseded v2 Neumann-limit artifact).

5. **The durable falsifiable result is Check A** (the drain ⇒ monopole/dipole `c_s`-radiation `∝ ε0 = 1−𝒯₀(0)`, tied to the gravity strength). Check B confirms 1/r² survives (robust-by-geometry) for the model's DC-sink return.

**Verdict ladder (v3, given the `Z<0` drain premise):** headline = `RETURN_RESIDUAL_PREDICTION` if `p=2` (1/r² survives + bounded radiation residual); `RETURN_NOGO` if `p≠2` for the model's localizing geometry / the delocalizing-warp control, OR if cancelling ℓ=0,1 kills ℓ=2; `BC_DEPENDENT` only if `p` differs across admissible DC-sink completions; `INCONCLUSIVE` if uncomputable. v2 §0.5 requirements #2 (p from a real radial solve, in `engine_agreement`) and #3 (`static_dynamic_consistency` = two independent computations) REMAIN IN FORCE; v2 §0.5.1's "both branches agree on Z" is SUPERSEDED by this §0.6.

## §A. Deliverables

- `reports/pathA_29_brane_bulk_return.md` — line 1 = COMPUTED top-line verdict (§D); then: the postulated slab structure (explicitly
  flagged postulated-vs-derived per element); the COMPUTED slab response/transfer function `𝒯_ℓ(ω)` (ℓ=0,1) and the resulting return
  moments `R0(ω)`, `R1(ω)`; the residual outgoing amplitudes `A0_res`, `A1_res` and their leading ω-powers with the **admissibility
  window** (the constraint on `d`/warp/BC under which the residual reaches the strict order `p_res(ℓ=0,1) ≥ 5`, the separate window for a
  bounded-but-lower-order residual, or a statement that none exists); the COMPUTED
  static long-range falloff exponent of `v_r(r≫d)` for the SAME structure; and the **headline reconciliation verdict** — whether one
  admissible structure satisfies BOTH (A) and (B). State precisely what feeds `pde_ledger` (open-item #9 + the gravity-range item) and
  what remains track-3-downstream (the full nonlinear return closure).
- `reports/pathA_29_results.yaml` — machine-readable AND machine-CHECKABLE (the controls/firewall must be assertable, not rhetorical).
  Required fields at minimum:
  - `slab_structure:` each element flagged `postulated` (free GEOMETRY/topology — spacing `d`, dimensionality of the return, warp ansatz
    if any) vs `derived` (forced by continuity — incl. the derived neighbor-brane BC equations + provenance per §0.5; the BC functional
    form is NOT a postulated free knob) — and an explicit statement that the RESPONSE is computed, not dialed;
  - `steady_state:` `{exists: bool|undecidable, type: circulation, provenance: <derivation-or-obstruction>}` — the dynamical steady
    brane↔bulk circulation that Check A linearizes around. Routing (see §0 + §D): a CANDIDATE with `exists: false` (computed) is
    INADMISSIBLE (dropped); `exists: false` for ALL admissible structures ⇒ `RETURN_NOGO`; `exists: undecidable` (in-scope intractable,
    `provenance` MUST carry the obstruction) ⇒ `INCONCLUSIVE`;
  - `transfer_function:` `𝒯_ℓ(ω)` per ℓ∈{0,1} as a symbolic ω-expansion, with `T_at_DC:` reporting the **COMPUTED** values `𝒯₀(0)` and
    `𝒯₁(0)` (NOT set to 1 by fiat — `𝒯₀(0)` from zeroth-moment/mass continuity = the DC return *fraction* to our brane, may be `<1`;
    `𝒯₁(0)` from first-moment continuity incl. bulk flux + return centroid, NOT forced to 1) and `provenance:` tracing each to
    continuity across the slab (NOT assumed). Provenance must be MACHINE-CHECKABLE, not prose-asserted: emit and script-validate
    `governing_equations:`, `bc_equations:`, `solution_basis:`, `series_coefficients:`, `free_parameter_table:`,
    `residual_substitution_trace:`, and `forbidden_fit_flags:` (incl. `fit_to_cancellation: false`); the scripts must FAIL if any of
    these is missing or inconsistent with the reported `𝒯_ℓ(ω)`. Also emit BOTH order fields per ℓ: `sigma_ℓ: ord_ω(𝒯_ℓ − 𝒯_ℓ(0))`
    (order of the first non-DC correction) AND `nu_ℓ: ord_ω(1 − 𝒯_ℓ)` (order of deviation-from-1; **`nu_ℓ = 0` WHENEVER `𝒯_ℓ(0) ≠ 1`**).
    Use ONLY `nu_ℓ` in the residual order: `p_res_0: 1+nu_0`, `p_res_1: 3+nu_1` (see §B.3 order test — do NOT use `sigma_ℓ` in `p_res`,
    or a `𝒯_ℓ(0)≠1` case would overstate suppression);
  - `return_moments:` `R0(ω)`, `R1(ω)` and `residual_amplitude:` `A0_res`, `A1_res` with `p_residual:` per ℓ and
    `p_quadrupole: 5` for comparison;
  - `admissibility_window_A:` the constraint on `{d, warp, derived_BC_branch/physical_completion}` for STRICT `p_residual(ℓ=0,1) ≥ 5` (requires `𝒯₀(0)=1 AND ν₀≥4` and
    `𝒯₁(0)=1 AND ν₁≥2`), OR `none` with the obstruction; separately, `residual_window_A:` the constraint for a bounded-but-lower-order
    residual (`p_res<5`, e.g. a net-sink monopole `∝ ε·ω`) — this feeds the `RETURN_RESIDUAL_PREDICTION` rung, NOT the strict window;
  - `static_falloff_B:` the COMPUTED long-range exponent `p` in `v_r ~ r^{-p}` at `r≫d` (target `p=2`), the crossover scale, the BC it
    was computed under, AND the COMPUTED SIGNED `zero_mode_source_integral:` `Z` (the net monopole source our brane sees at `r≫d` AFTER
    the return — see §B.4). **SIGN CONVENTION (load-bearing — gravity IS inflow):** `Z` = net rate of medium ADDED to our brane localized
    at the throat, so `Z<0` = net REMOVAL = sink = INWARD flow (`v_r<0`) = GRAVITY; `Z>0` = outflow = anti-gravity = NOT gravity.
    **`Z` is a COMPUTED, SIGNED output, NOT defined as `±M0(1−𝒯₀)` by fiat:** emit the signed terms separately —
    `Z_terms: {throat_sink: −M0, return: +M0·𝒯₀(0), replenishment_localized: ≥0, Z_boundary_dof: <signed-or-none>}` where (i)
    `throat_sink = −M0` (the drain, inward), (ii) `return = +M0·𝒯₀(0)` (the slab return, outward addition), (iii)
    `replenishment_localized` = the LOCALIZED contribution of the areal leak to the monopole at `r≫d` (the part affecting `1/r²`) — an
    areal-leak SOURCE so `≥0`, and (iv) `Z_boundary_dof` = an OPTIONAL, script-validated, SIGNED contribution from a DECLARED
    boundary/return DOF (the only term that could make `Z<0` at strict-A — the representable CONSISTENT_WINDOW escape, see §D.2), with
    REQUIRED sub-fields `{governing_equation, provenance, dc_sign_derivation}` when present, else `Z_boundary_dof: none`. **The
    computation MUST either (a) populate `Z_boundary_dof` from a declared DOF and COMPUTE its DC sign contribution to `Z` at `𝒯₀(0)=1`,
    or (b) set `Z_boundary_dof: none`** (then the local-channel result stands); a non-`none` value with missing/uncomputed sub-fields
    must make the scripts FAIL. **`replenishment_localized` is an INDEPENDENT cosmological process (conceptual_foundation §5), NOT
    "required/forced by the steady state"; its LOCALIZED contribution is computed and in the simple LOCAL channel
    `replenishment_localized = 0` (no localized enhancement).** The UNIFORM far-field areal-leak BACKGROUND (whole-brane cosmological
    influx, §5) does NOT enter `Z`/`B_pass` — emit it as a SEPARATE provenance-tagged background term `Z_uniform_background:` (a
    different observable). So in the LOCAL channel (`Z_boundary_dof: none`) `Z = −M0(1−𝒯₀(0))`. Also emit `Z_sign:` (the computed sign
    of `Z`) and `p_eq_2:` (the separately-computed falloff condition); the scripts MUST compute `Z` and assert its sign (e.g. VERIFY
    that in the LOCAL channel `Z≥0` when `𝒯₀(0)=1`), not assert it in prose. The BC MUST be the SAME BC as Check A; same-structure enforcement must be
    MACHINE-ASSERTED, not prose: emit a canonical `structure_id:` (a hash/canonical form of `{d, warp equation, bc_derivation_type,
    derived BC equations, physical-completion parameters}`) for BOTH Check A and Check B and assert byte-equality of the two ids
    in-script; if A and B use compatible-but-distinct limits of one derived-BC branch, emit the explicit `limiting_map:` connecting them
    and assert it reduces to the same `structure_id` family; assert `same_structure_as_A: true`;
  - `Z_is_premise: true` — per §0.6.1, `Z<0` (inward = gravity) is the drain ADMISSIBILITY PREMISE, not a co-equal computed test; an
    admissible Check-B completion is a DC SINK (transmits the steady drain). A DC-reflecting boundary (`Z=0`, the radiation/Neumann
    `ω→0` limit) is the no-drain/no-particle degenerate case — INADMISSIBLE for Check B (NOT a Z-veto branch);
  - `branch_results:` the Check-B comparison is across the admissible **DC-sink completions** (§0.6.2 — `destructuring_absorbing` PRIMARY +
    `bloch_stack` CROSS-CHECK); the radiation/Sommerfeld branch is AC/Check-A only (emit if reported, TAG `ac_check_a_only`, do NOT use as
    a Check-B falloff/Z branch). Emit `{destructuring_absorbing: {verdict, p, Z_sign, spectrum_class, ...},
    bloch_stack: {verdict, p, Z_sign, spectrum_class, ...}, radiation: {tag: ac_check_a_only, ...}}` where
    `spectrum_class ∈ {normalizable_zero_mode | continuum_p3 | gapped_yukawa}` (§0.5.4); plus
    `branch_agreement: {destructuring_p, bloch_p, destructuring_verdict, bloch_verdict, p_agree: bool}` over the DC-sink completions.
    The §D.0 headline = RESIDUAL_PREDICTION if all DC-sink `p=2`, NOGO if all `p≠2`/ℓ=2-killed, else **`BC_DEPENDENT`** (the FALLOFF `p`
    DIFFERS across DC-sink completions — NOT a Z-sign disagreement). FORBID promoting a favorable completion (§0.6);
  - `reconciliation:` `{verdict: <§D.0 headline: RESIDUAL_PREDICTION | NOGO | BC_DEPENDENT | INCONCLUSIVE>, B_pass: bool, p: <exponent>,
    Z_is_premise: true, T0_at_DC: <value>, same_structure: bool}` — per §0.6 the `verdict` is the §D.0 DC-sink completion headline:
    `RETURN_RESIDUAL_PREDICTION` if ALL admissible DC-sink completions give `p=2`; `RETURN_NOGO` if all `p≠2` (or ℓ=2 killed);
    **`BC_DEPENDENT`** if the FALLOFF `p` DIFFERS across DC-sink completions (NOT a Z-sign disagreement). **`B_pass ⇔ p=2` under the DC-sink
    drain PREMISE** (`Z<0` is the premise, already secured for an admissible DC sink — NOT re-tested co-equally). The
    `A_strict_pass`/`p_res≥5` strict-cancellation fields + `Z_boundary_dof_status`/`Z_reduces_to_local` machinery (below, retained for
    record) belong to the SUPERSEDED v2 CONSISTENT_WINDOW rung (§D.2) and do NOT drive the v3 headline. Report `T0_at_DC` as the pivotal
    quantity (`ε0 = 1−𝒯₀(0)` sets the Check-A radiation residual). **[v2 fields, retained for traceability:** `A_strict_pass`,
    `A_residual_pass`, `Z_sign`, `p_eq_2`, `Z_reduces_to_local`, `Z_boundary_dof_status: populated|none|intractable`, `joint_window` — see
    §D.2 supersession note.**]**;
  - `bc_is_local_causal_passive: true`, `bc_omega_dependence_source: none|sommerfeld_radiation|bloch_periodic|declared_dof`,
    `window_is_open: true|false`, `window_symmetry_protected: <name|none>` (firewall flags — see §C);
  - `residual_prediction:` if the residual is bounded-but-lower-order (`p_res<5`), its quantified scaling (e.g. `∝ (d/λ)^q ω^p`, plus
    the `ε`-to-`1/r²`-strength tie) as a first-class falsifiable difference-from-GR feeding the `RETURN_RESIDUAL_PREDICTION` rung — a
    prediction, not a failure ([[feedback-falsification-is-the-goal]]);
  - `T2_effect:` if the slab return induces a `𝒯₂` acting on ℓ=2, report its effect SEPARATELY (do NOT auto-fail on it); assert the
    ℓ=2 DtN kernel retains `p_raw=5` (see §E quadrupole-survives + §D.1 trigger);
  - `controls:` a map of each §E control name → `{same_pipeline: bool, fired: <verdict-or-bool>}`. For `no_go_reachable`, use
    `status: reachable_RETURN_NOGO | failed_not_admissibly_reachable` — **only `reachable_RETURN_NOGO` is acceptable for v3 green** (the
    §D classifier returns `RETURN_NOGO` on the DERIVED admissible delocalizing-warp family; emit the warp equation + computed
    `spectrum_class`, §E). `failed_not_admissibly_reachable` (the able-to-fail proof could NOT be admissibly exhibited) is a
    FAILING/unsound outcome — it is NEVER a pass and routes the directive to unsound / `INCONCLUSIVE`, never green (v3-#3). For
    `tension_is_real`, use `status: witnessed | not_witnessed | not_applicable` (only `witnessed` may be cited — see §E).
- `tools/pathA_29_brane_bulk_return_sympy.py` + `tools/pathA_29_brane_bulk_return.wl` — dual engine; both exit 0;
  `engine_agreement: PASS` on every shared symbolic expression (real symbolic compare, not a hardcoded literal). Verdict COMPUTED.
  **`engine_agreement.checked_quantities` MUST include (scripts FAIL if absent): `p`, `p_eq_2`, `dynamic_limit_exponent`,
  `static_solve_exponent`, `zero_mode.r_dependence`, and `green_function`** (per §0.5.2 — `p`/`p_eq_2` must be DERIVED and
  cross-engine-verified, NOT a spectral-class lookup or bare literal).

## §B. The computation (WHAT to produce — Codex designs HOW)

Symbolic, units restored, in the frozen slice. The natural unifying object is the **steady-to-dynamic slab Green's function** for a
source at the throat between the brane stack — Check B is its `ω→0` spatial falloff and Check A is its finite-ω transfer at the
source. They MUST be two limits of ONE DYNAMIC OBJECT for the SAME postulated structure (NOT a separate frozen static solve for B — see
§C no-frozen rule). (This is the expected framing; Codex may choose another route that delivers the same deliverables and the same
SAME-DYNAMIC-OBJECT guarantee.)

**Tractability bound (dual-engine, NON-NEGOTIABLE).** General warped-slab Green's functions with arbitrary BCs are NOT reliably
SymPy/Mathematica-tractable in closed form. Restrict the EXECUTABLE family to flat OR a named RS-like warp, plus analytically-tractable
DERIVED BCs (per §0.5): Sommerfeld outgoing radiation, periodic/Bloch stack, or a declared passive boundary DOF with symbolic response.
A bare Robin/impedance BC is admissible ONLY if DERIVED as one of these limits (NOT postulated as a tuned form). Deliver the transfer
functions as ω-series carried to the orders the residual test needs: through `O(ω⁴)` for ℓ=0 and `O(ω²)` for ℓ=1 (see §B.3 order test).
Anything broader than this analytic family is out of scope for this gate (report as `INCONCLUSIVE` with the precise obstruction rather
than a non-converging symbolic attempt).

1. **Postulate the slab GEOMETRY (freely, but DOCUMENTED and flagged); DERIVE the BC (v2 §0.5).** A multi-brane stack partitioning the
   bulk `w`-direction into finite slabs of thickness `d`, our brane at `w=0` carrying the throat/drain, with a return mechanism at the
   adjacent brane(s). State explicitly: the spacing `d`, whether the bulk is flat or warped (and the warp ansatz if used) — these
   GEOMETRY/topology elements are `postulated`. The neighbor-brane boundary condition is NOT a postulated free knob: DERIVE it from the
   declared physical completion per §0.5 (Sommerfeld radiation, periodic/Bloch, or a declared boundary DOF) — emit its equations +
   provenance. The RESPONSE below is `derived`.
2. **Compute the return response from continuity (NO free knob — this is the teeth).** Derive the slab transfer function `𝒯_ℓ(ω)`
   (ℓ=0,1) — how a time-varying leak at the throat returns to the brane after transiting the slab — from the GNLS continuity /
   linearized bulk transport across the postulated structure, NOT by assuming a value. **DC values are COMPUTED, not anchored to 1 by
   fiat:** `𝒯₀(0)` follows from zeroth-moment (mass) continuity across the slab and equals 1 ONLY if all leaked medium returns locally
   to our brane — it MAY be `<1` (medium transmits to the neighboring brane, sustaining a net in-brane sink); `𝒯₁(0)` follows from
   FIRST-moment continuity (bulk flux + return centroid) and NO conservation law forces it to 1 (Check A is ALLOWED to fail at raw
   dipole order). The finite-ω behaviour (transit phase/attenuation/any cavity resonance) is what the geometry forces. Then
   `R0(ω)=−M0·𝒯₀(ω)`, `R1_i(ω)=−D1_i·𝒯₁(ω)` (carry the ℓ=1 spatial separation of leak-vs-reinjection honestly — the dipole is
   sensitive to WHERE the return re-injects, not just the net rate).
3. **Check A — residual radiation order + admissibility window.** Form the net far-field source and the residual outgoing amplitudes,
   written explicitly via the transfer function: `A0_res=i·a·(ω/c_s)·M0·(1−𝒯₀(ω))`, `A1_res=i·a³(ω/c_s)³·D1·(1−𝒯₁(ω))/2`. **Explicit
   order test (TWO distinct order symbols — do NOT conflate):** define `σ_ℓ = ord_ω(𝒯_ℓ − 𝒯_ℓ(0))` (order of the FIRST non-DC
   correction, so `𝒯_ℓ(ω)=𝒯_ℓ(0)+O(ω^{σ_ℓ})`) and `ν_ℓ = ord_ω(1 − 𝒯_ℓ)` (order of the deviation-FROM-1). Crucially `ν_ℓ = 0`
   WHENEVER `𝒯_ℓ(0) ≠ 1` (then `(1−𝒯_ℓ)` has a nonzero constant term ⇒ order 0 ⇒ `p_res`=raw, NO suppression — regardless of `σ_ℓ`);
   `ν_ℓ = σ_ℓ` only when `𝒯_ℓ(0)=1`. Use ONLY `ν_ℓ` in the residual order: `p_res(0)=1+ν₀`, `p_res(1)=3+ν₁`. A STRICT pass therefore
   needs `𝒯₀(0)=1 AND ν₀≥4` (→ `p_res(0)≥5`) and `𝒯₁(0)=1 AND ν₁≥2` (→ `p_res(1)≥5`). (Illustrative, NOT a pre-decided outcome: a
   generic transit delay `𝒯=1−iωτ+…` (`𝒯(0)=1`, `σ=ν=1`) gives `p_res(0)=2`, `p_res(1)=4` → FAILS strict — the gate has teeth.) Determine the **admissibility windows**: the constraint on
   `{d, warp, derived_BC_branch/physical_completion}` (if any) for STRICT `p_res(ℓ=0,1)≥5`, AND separately the constraint for a bounded-but-lower-order residual
   (`p_res<5`, e.g. a net-sink monopole `∝ ε·ω`) which feeds the `RETURN_RESIDUAL_PREDICTION` rung (report its quantified scaling). If
   the residual is irreducibly above the quadrupole for ALL admissible structures AND 1/r²+bounded-residual cannot be obtained together
   → that is the able-to-fail `RETURN_NOGO` outcome (§D). **A naive instantaneous return (`𝒯=1` by hand) is FORBIDDEN (§C)** — the
   transit response must be computed; the physically expected leading correction is a transit delay `~d/c_s`, and whether it can be made
   harmless is the open question.
4. **Check B — QUASI-STATIC (`ω→0` LIMIT) 1/r² survival under the SAME structure AND the SAME source+return accounting.**
   **[v3 §0.6 GOVERNS THIS STEP — read first.] `Z<0` is the drain PREMISE (admissibility), NOT a computed pass/fail test.** Run Check B
   under the admissible **DC-sink** completions (§0.6.2: §5 de-structuring/absorbing PRIMARY + Bloch CROSS-CHECK), which transmit the
   steady drain so `Z<0` by premise; a DC-REFLECTING boundary (`Z=0`, including the radiation/Sommerfeld `ω→0`-Neumann limit) is the
   no-drain/no-particle degenerate case and is INADMISSIBLE for Check B (§0.6.1, §0.6.3) — do NOT use it as a Check-B Z/falloff branch.
   The signed-`Z` bookkeeping below remains the correct ACCOUNTING of the in-brane monopole flux, but per §0.6 its sign is the PREMISE
   (an admissible DC-sink return has `Z<0`); the v2 "the sign forecloses the window / `Z_boundary_dof` rescue at strict-A" interpretation
   is SUPERSEDED — the real Check-B teeth are the FALLOFF `p` (robust `p=2` for a localizing geometry; `p≠2` only under a delocalizing
   warp → NOGO, §0.6.4). **`B_pass ⇔ p=2`** under the DC-sink premise (the `Z<0` factor is the admissibility premise, already secured).
   Check B MUST be obtained as the **`ω→0` LIMIT of the SAME dynamic slab response object computed in Check A** (one analytic object, e.g. `G(ω,r)`:
   Check A = its finite-ω transfer at the throat, Check B = its `ω→0` spatial falloff). It is **FORBIDDEN** to pose Check B as an
   independent/standalone static BVP that does not inherit the dynamic couplings (the project's #1 relapse — a frozen solve deletes the
   coupling that selects the balance). The background is the steady brane↔bulk CIRCULATION (§0 `steady_state`), not a frozen profile.
   Evaluate this `ω→0` limit in the SAME postulated geometry (same `d`, same DERIVED physical BC + limiting map as B.1/B.2, per §0.5) using the SAME source+return accounting that
   defines `𝒯₀(0)` — i.e. the leak `S_leak` AND the return, NOT an isolated sink. **Check B's observable, precisely:** the IN-BRANE
   radial flow `v_r` felt by brane-confined test bodies (the 3D Gauss law applies to the in-brane NET flux). **SIGN CONVENTION (gravity
   IS inflow):** `Z` = net rate of medium ADDED to our brane localized at the throat, so `Z<0` = net removal = sink = INWARD flow
   (`v_r<0`) = GRAVITY; `Z>0` = outflow = anti-gravity = NOT gravity. **Output the COMPUTED SIGNED zero-mode source integral `Z`** — the
   net monopole source our brane sees at `r≫d` AFTER the return — **as a sum of SEPARATE SIGNED COMPUTED terms with provenance**:
   `Z = Z_throat_sink + Z_return + Z_replenishment_localized + Z_boundary_dof`, where `Z_throat_sink = −M0` (the drain, inward),
   `Z_return = +M0·𝒯₀(0)` (the slab return, outward addition), `Z_replenishment_localized` is the LOCALIZED contribution of the areal
   leak to the monopole at `r≫d` (the part affecting `1/r²`) — an areal-leak SOURCE so `≥0`, and `Z_boundary_dof` is the OPTIONAL,
   script-validated, SIGNED contribution of a DECLARED boundary/return DOF (the ONLY term that could make `Z<0` at strict-A — the
   representable `RETURN_CONSISTENT_WINDOW` escape; see §A schema + §D.2). The computation MUST either COMPUTE `Z_boundary_dof` (with
   `{governing_equation, provenance, dc_sign_derivation}`) from a declared DOF, or set it `none` (then the local-channel result stands).
   **`Z_replenishment_localized` is an INDEPENDENT cosmological process (conceptual_foundation §5), NOT "required/forced by the steady
   state"; its LOCALIZED contribution is computed and equals 0 in the simple LOCAL channel (no localized enhancement).** The UNIFORM
   far-field areal-leak background (whole-brane cosmological influx,
   §5) does NOT enter `Z`/`B_pass` — it is a SEPARATE observable (`Z_uniform_background`), provenance-tagged. **Two DISTINCT computed
   conditions (do NOT conflate):** (i) the long-range falloff exponent `p` in `v_r ~ r^{-p}` from `Z` at `r≫d`, classified by the COMPUTED
   transverse `spectrum_class` (§0.5.4): `normalizable_zero_mode → p=2`; `continuum_p3 → p=3` (gapless leakage); `gapped_yukawa → Yukawa`
   exponential cutoff (these last two are DISTINCT failures — do NOT conflate) and (ii) the SIGN of `Z`. `Z≠0` (amplitude) does NOT imply
   `p=2` (falloff). **`B_pass ⇔ p=2` under the DC-sink drain PREMISE** (`Z<0` is the premise per §0.6.1, already secured for an admissible
   DC-sink return — it is NOT re-tested as a co-equal condition; the falloff `p` is the Check-B teeth).
   **[SUPERSEDED BY §0.6 — retained only as the AC/Check-A accounting; do NOT apply to Check B.]** The following v2 paragraph treated
   `Z<0` as a computed test that the radiation/Neumann `ω→0` limit forecloses (`Z=0` at strict-A). Per §0.6.1/§0.6.3 that limit is the
   no-drain/no-particle degenerate case (INADMISSIBLE for Check B), NOT a Z-veto; the `Z_boundary_dof` "rescue" is no longer the pivot.
   The bookkeeping is still correct as AC/Check-A flux accounting: in the radiation (AC) channel `Z = −M0(1−𝒯₀(0))`, so at strict-A
   (`𝒯₀(0)=1`) `Z = Z_replenishment_localized ≥ 0` (`= 0` in the simple witness) `⇒ Z≥0`; `Z_boundary_dof` (the OPTIONAL declared-DOF
   signed term, §A schema) remains representable for completeness but is NOT the Check-B pivot under §0.6. The gate decides
   `RETURN_RESIDUAL_PREDICTION` (if `p=2`) vs `RETURN_NOGO` (if `p≠2` for the localizing geometry / delocalizing-warp control): the
   admissible DC-sink channel (`𝒯₀(0)<1`, `Z<0` by premise) gives inward `1/r²` gravity plus a bounded radiation residual `∝ ε0=1−𝒯₀(0)`.
   The scripts MUST still COMPUTE the signed `Z` and report its terms (accounting), but the Check-B VERDICT is driven by the FALLOFF `p`
   under the DC-sink premise, NOT by a Z-sign veto. **Check B MUST use the SAME declared physical completion/structure as Check A, with a
   machine-checked limiting map** (`structure_id` equality over `{d, warp, completion}`) — but NOTE Check A may report the
   radiation/Sommerfeld boundary (`ac_check_a_only`) while Check B uses ONLY the DC-sink completions (`destructuring_absorbing`,
   `bloch_stack`); the radiation AC boundary is NOT forced back into Check B (§0.6.3). — consistency test.
5. **Reconcile (the headline — via the §D.0 DC-sink completion gate, per §0.6).** Compare the Check-B FALLOFF `p` across the admissible
   **DC-sink** completions (§0.6.2: de-structuring/absorbing PRIMARY + Bloch CROSS-CHECK — the radiation/Sommerfeld branch is AC/Check-A
   only, NOT a Check-B completion). `Z<0` is the drain PREMISE (admissibility), already secured — Check B passes on `p=2`. Report the
   per-completion results and the §D.0 headline: `RETURN_RESIDUAL_PREDICTION` if `p=2` (1/r² survives + the bounded Check-A radiation
   residual); `RETURN_NOGO` if `p≠2` (localizing geometry / delocalizing-warp control) or ℓ=2 is killed; **`BC_DEPENDENT` ONLY if the
   FALLOFF `p` DIFFERS across the admissible DC-sink completions** (NOT a Z-sign disagreement — the superseded v2 Neumann-limit artifact);
   `INCONCLUSIVE` if uncomputable. State the implied inter-brane-spacing constraint where applicable.

## §C. Anti-tautology firewall (ENFORCED — postulate-to-win is exactly the pathA_27 death; do not repeat it)

- **The return RESPONSE must be COMPUTED from continuity, never dialed.** FORBIDDEN as the basis for an (A)-pass: setting `𝒯_ℓ=1`
  ("instantaneous return") by hand; assuming the return re-injects at exactly the leak point (auto-kills the ℓ=1 dipole); choosing
  the wall BC to land the answer without computing the response it implies. `𝒯_ℓ(ω)` (INCLUDING its DC value `𝒯_ℓ(0)`) must trace to
  the slab continuity/transport (provenance asserted, machine-checkable per §A) — **`𝒯_ℓ(0)` is NOT anchored to 1 by fiat; it is a
  computed output** (`𝒯₀(0)` may be `<1`, `𝒯₁(0)` is not forced to 1; see §B.2).
- **BC family must be DERIVED and physical, not fitted (anti-postulate-to-win hardening; v2 §0.5).** The neighbor-brane BC must be
  DERIVED from a declared physical completion (per §0.5), not a tuned wall impedance. Allowed `ω`-dependence SOURCES are `none`,
  `sommerfeld_radiation` (the canonical causal open boundary `∂_wΦ=i(ω/c_s)Φ`), `bloch_periodic` (periodic stack), or `declared_dof`
  (an explicitly declared boundary dynamical DOF whose governing equation is given and `ω`-dependence DERIVED from it). Sommerfeld/Bloch
  `ω`-dependence is allowed ONLY with emitted derivation equations + limiting maps; an ARBITRARY tuned `m`/`ω` impedance chosen to
  preserve/kill the zero mode remains FORBIDDEN (the pathA_27 / v1 death). The BC and any warp MUST be INDEPENDENT of the source
  moments `M0/D1/Q2` (no fitting the wall to the source). A Check-A or Check-B pass MUST be an OPEN, finite-measure, STABLE window in
  parameter space, NOT an isolated fine-tuned root — UNLESS the cancellation is symmetry-protected (then NAME the symmetry). Emit the
  machine-checkable flags `bc_is_local_causal_passive: true`,
  `bc_omega_dependence_source: none|sommerfeld_radiation|bloch_periodic|declared_dof`, `window_is_open: true|false`,
  `window_symmetry_protected: <name|none>`.
- **SAME structure for both checks (non-negotiable).** Check A and Check B MUST use the identical postulated GEOMETRY AND the SAME
  DERIVED physical BC + same limiting map (`{d, warp, derived_BC_branch/physical_completion}`). Using a confining BC for B and an
  absorbing BC for A is the cheat that makes the tension vanish — FORBIDDEN. The deliverable asserts `same_structure: true`; if A and B
  genuinely require *compatible-but-distinct* limits of one derived-BC branch, that must be shown to be one physical structure, not two.
- **The 1/r² falloff must be SOLVED from the slab Green's function, not asserted (EXECUTION ENFORCEMENT — the v1 first execution VIOLATED this by hardcoding `r²` in `v_r` and reading the exponent back via a `power_law_exponent()` call; REJECTED at tri-review).** Check B MUST actually CONSTRUCT and SOLVE the flat finite-slab Green's function (separation of variables / the transverse eigenmode spectrum of the `w`-operator under the admissible **DC-sink** completions — §0.6: `destructuring_absorbing` PRIMARY + `bloch_stack` CROSS-CHECK, NOT the AC-only radiation/Sommerfeld boundary; same declared physical completion/structure as Check A with a machine-checked limiting map) and DERIVE the long-range exponent `p` from the zero-mode's r-dependence — it MUST be genuinely able to come out `p=2` (`normalizable_zero_mode`), `p=3` (`continuum_p3`, gapless 4D leak), or Yukawa (`gapped_yukawa`, a mass gap / no zero mode) — the three DISTINCT `spectrum_class` values of §0.5.4. **FORBIDDEN:** writing any `r^{-n}` into the flow/`v_r` expression by hand and extracting `n` from it (a `power_law_exponent()` applied to a hand-written `r^{-n}` does NOT satisfy this). **EMIT as machine-checkable evidence:** `transverse_mode_spectrum:` (eigenvalues + eigenfunctions), `zero_mode:` `{exists, normalizable, r_dependence}`, and the actual solved `green_function:` expression; the falloff `p` must be extracted FROM these.
- **Controls MUST fire from COMPUTED inputs, never hardcoded classifier flags (EXECUTION ENFORCEMENT — the v1 execution fed `no_inward_branch=True`, the main NOGO triggers + boundary-DOF rescue as typed-in literals, and a tautological `tension_is_real` A-witness / a `static_dynamic_consistency` variable-vs-itself; REJECTED).** `no_go_reachable` MUST (for v3 green, `status: reachable_RETURN_NOGO`) construct the DERIVED admissible **delocalizing/anti-localizing warp** family — with its warp equation EMITTED and the COMPUTED `spectrum_class` (`continuum_p3` | `gapped_yukawa`) — and run the SAME §D classifier on it to obtain `RETURN_NOGO` from real inputs (a pinned/absorbing neighbor-brane DOF may be SUPPLEMENTAL only and CANNOT satisfy `no_go_reachable` for green; a bare Robin limit is inadmissible under §0.5/§0.6; see §A enum + §E). `tension_is_real` MUST exhibit the v3 FALLOFF tension (§0.6.4, §E): a localizing admissible DC-sink completion giving `p=2` AND a DERIVED delocalizing/anti-localizing warp giving `p≠2` — NOT the superseded `κ0→∞` strict-A / Z-sign A-pass-B-fail witness. `static_dynamic_consistency` MUST compare the `ω→0` limit of the ACTUAL dynamic Green's function to the ACTUAL Check-B solve (two independent computations, §0.5.3 — not a variable to itself). The main-pipeline NOGO trigger (`p≠2` for the localizing geometry / delocalizing-warp control) MUST be COMPUTED — the verdict must be genuinely able to come out `RETURN_NOGO`.
- **NO frozen/decoupled static solve (NOTHING is static).** FORBIDDEN: any step that treats a field as static/frozen in a way that
  decouples it from the slab dynamics. The static (`ω→0`) result MUST be a LIMIT of the dynamic object (the SAME `G(ω,r)` / response
  used in Check A), never an independently-posed frozen solve. The Check-A background must be a genuine dynamical steady state
  (circulation), with `steady_state.exists` computed (a structure with no steady circulation is INADMISSIBLE).
- **No free magnitude knob standing in for the physics.** Any postulated coupling that is not fixed by continuity must be flagged and
  must NOT be the thing that delivers the cancellation. (The pathA_27 `κ_c` lesson.)

## §D. Verdict ladder (COMPUTED; v3 §0.6 — operative outcomes: `RETURN_RESIDUAL_PREDICTION` | `RETURN_NOGO` | `BC_DEPENDENT` | `INCONCLUSIVE`, decided on the Check-B FALLOFF `p` across DC-sink completions with `Z<0` the drain premise; the v2 `RETURN_CONSISTENT_WINDOW` rung is SUPERSEDED/retained-for-record; the no-go rung MUST be genuinely reachable)

**§D.0 — DC-SINK COMPLETION GATE (run FIRST; v3 §0.6 — supersedes the v2 "both branches agree on Z" gate).** Per §0.6, `Z<0` is the
drain PREMISE (admissibility), not a co-equal test, and the radiation/Sommerfeld boundary is AC/Check-A only (NOT a Check-B Z/falloff
branch — its `ω→0`-Neumann limit gives `Z=0` = no-drain = INADMISSIBLE for Check B). The Check-B completions are the admissible **DC
sinks**: §5 de-structuring/absorbing (PRIMARY) + periodic-stack Bloch (CROSS-CHECK). Compute the Check-B FALLOFF `p` for EACH DC-sink
completion and combine:
- if ALL admissible DC-sink completions give `p=2` → headline `RETURN_RESIDUAL_PREDICTION` (1/r² survives + the bounded Check-A radiation
  residual `∝ ε0=1−𝒯₀(0)`) — robust (the EXPECTED outcome, §0.6.4/§0.6.5);
- if ALL give `p≠2` (the localizing geometry / delocalizing-warp control fails to localize the zero mode → continuum `p=3` or gapped
  Yukawa), OR cancelling ℓ=0,1 kills ℓ=2 → headline `RETURN_NOGO` (robust, falsification);
- if the FALLOFF `p` DIFFERS across admissible DC-sink completions → headline **`BC_DEPENDENT`**: 1/r² survival depends on the
  return-boundary physics. Report each completion's `p`, `spectrum_class` (`Z_sign<0` is the premise); NAME the further declared physical
  completion that would decide it; **do NOT promote a favorable completion.** (`BC_DEPENDENT` ≈ INCONCLUSIVE-with-named-condition; does NOT
  close pde_ledger #9.) **NOTE:** `BC_DEPENDENT` fires on a `p` (falloff) disagreement ONLY — NOT on a `Z`-sign disagreement (the
  superseded v2 Neumann-limit artifact). Emit `branch_agreement` over the DC-sink completions (de-structuring vs Bloch) per §A.

**PER-COMPLETION PRECEDENCE (Z<0 is the premise; the falloff `p` is the teeth):** for each admissible DC-sink completion classify
top-down: NOGO (`p≠2`: delocalizing-warp / continuum / gap — a genuine localization failure) → RESIDUAL_PREDICTION (`p=2`: inward
`1/r²` with the bounded Check-A residual `∝ ε0`) → INCONCLUSIVE (uncomputable). The §D.0 gate then combines the per-completion `p`
results into the headline. **[SUPERSEDED — `RETURN_CONSISTENT_WINDOW` and the `Z_boundary_dof`/`Z_reduces_to_local` strict-A machinery
of v2 (rung 2 below) are retained only for record; under §0.6 the gate decides RESIDUAL_PREDICTION vs NOGO vs BC_DEPENDENT on the
falloff `p`, with `Z<0` the premise.]**

1. `RETURN_NOGO` — **no** admissible DC-sink completion yields `p=2` (no inward `1/r²` gravity survives the slab — §0.6), in particular
   ANY of: the falloff is `p≠2` (continuum `p=3` / gapped Yukawa) for ALL admissible DC-sink completions / the delocalizing-warp control
   (a genuine zero-mode localization FAILURE, §0.6.4); OR cancelling ℓ=0,1 ALSO suppresses ℓ=2 below `p=5` (the GW channel that must
   survive was killed — see §E/§A `T2_effect`); OR COMPUTED no steady-state circulation exists for ALL admissible structures (the return
   mechanism cannot run). **NOTE (§0.6):** `Z<0` is the drain PREMISE (admissibility), NOT a NOGO trigger — a `Z=0` boundary
   (radiation/Neumann `ω→0` limit) is the no-drain/no-particle degenerate case, INADMISSIBLE for Check B, NOT a NOGO via "sign never
   inward." The expected `p=2` DC-sink branch is `RETURN_RESIDUAL_PREDICTION` (rung 3), which gives `1/r²` gravity. **Hard structural
   failure / falsification** of the multi-brane-return resolution of the drain. **Scope-honest:** report `RETURN_NOGO`
   **WITHIN-THE-TRACTABLE-FAMILY** — escape mechanisms outside the tractable scope (richer warps/return DOFs) would be a follow-up gate;
   do NOT claim the no-go is absolute. Report prominently. (Able-to-fail; MUST be reachable — see §E no-go-reachable control.)
2. `RETURN_CONSISTENT_WINDOW` — **[SUPERSEDED BY §0.6 — retained for record only; do NOT route the headline here.]** This v2 rung tied
   the closing of pde_ledger #9 to a strict-A `(p=2 AND Z<0)` window reachable only via a computed `Z_boundary_dof<0`. Under §0.6, `Z<0`
   is the drain PREMISE (not a strict-A computed test), the radiation/Neumann `Z=0` foreclosure is the no-drain degenerate case (not a
   veto), and the gate decides on the falloff `p` — so this rung is no longer the operative pivot. (The v2 text — STRICT
   `p_res(ℓ=0,1)≥5` plus `Z_boundary_dof<0` at `𝒯₀(0)=1`, `Z_reduces_to_local`, the EXPECTED-UNREACHABLE/representable-escape machinery —
   is preserved below for traceability but does NOT drive the v3 headline.) Under §0.6 the durable falsifiable content is the Check-A
   radiation residual (§0.6.5); closing #9 (full cancellation) is not claimed by this gate.
3. `RETURN_RESIDUAL_PREDICTION` — **the EXPECTED v3 headline (§0.6.4/§0.6.5).** An admissible DC-sink completion gives Newtonian inward
   gravity (Check B `p=2`; `Z<0` by the drain premise, so `𝒯₀(0)<1`) AND a bounded/lower-order residual `p_res<5` for ℓ=0 and/or ℓ=1
   (the net-sink monopole/dipole `c_s`-radiation `∝ ε0 = 1−𝒯₀(0)`, tied to the `1/r²`/gravity strength). Does **NOT** close #9. RECORDS
   the predicted residual radiation as a first-class falsifiable difference-from-GR (quantify its `(d/λ, ω)` scaling and the
   `ε0`-to-gravity-strength tie). This is the EXPECTED outcome because the drain (premise) is real and the falloff is robust `p=2` for a
   localizing geometry. (`RETURN_NOGO` ≠ `RETURN_RESIDUAL_PREDICTION`: NOGO = no inward `1/r²` survives (`p≠2`) for any DC-sink completion
   (rung 1); RESIDUAL_PREDICTION = inward `1/r²` holds AND a quantified sub-leading radiation residual.)
4. `INCONCLUSIVE` — `𝒯_ℓ(ω)` or the quasi-static (`ω→0`-limit) falloff `p` cannot be computed in-scope for the DC-sink completions, OR
   the static-dynamic-consistency control disagrees (§E), OR whether a steady-state circulation EXISTS is UNDECIDABLE in-scope
   (intractable). NOTE: a COMPUTED absence of steady circulation for ALL admissible structures is `RETURN_NOGO` (§D.1), NOT
   `INCONCLUSIVE` — only an INTRACTABLE existence/falloff question routes here (state precisely why; do not fabricate a window). [The v2
   `Z_boundary_dof`/CONSISTENT_WINDOW intractability route is SUPERSEDED with §D.2.]

## §E. Controls (each an EXECUTED assertion that genuinely fires; each runs the SAME pipeline)

- **DC-value is COMPUTED (not a fiat anchor):** assert `𝒯₀(0)` and `𝒯₁(0)` are COMPUTED from zeroth-/first-moment continuity (NOT set
  to 1), report their values, and assert that even when `𝒯₀(0)=1` this ALONE does not give an (A)-pass (it only fixes DC; the radiation
  lives at finite ω, requiring `ν₀≥4`) — proves the cancellation is a real finite-ω constraint and that the DC value is a live,
  possibly-`<1` output (the §B.4 pivotal quantity).
- **No-go reachable (the MANDATORY able-to-fail proof; v3-#3 — an ADMISSIBLY-exhibited NOGO is REQUIRED for green).** The gate must be
  PROVABLY able-to-fail. For v3 SOUND/green, the control MUST return `no_go_reachable.status: reachable_RETURN_NOGO`: it MUST construct
  the DERIVED **delocalizing (wrong-sign / anti-localizing) warp** family and run the SAME full §D classifier to OBTAIN `RETURN_NOGO`
  from COMPUTED inputs. This is a genuinely ADMISSIBLE structure — the drain still transmits (`Z<0`, the §0.6 premise HOLDS) — but the
  warp FAILS to localize the transverse zero mode → continuum (`p=3`) or gap (Yukawa) → `p≠2` → `RETURN_NOGO` (exactly the
  `conceptual_foundation.md` §7 braneworld caveat realized as an ADMISSIBLE failure path). The control MUST EMIT the delocalizing warp
  EQUATION + the COMPUTED `spectrum_class` (`continuum_p3` | `gapped_yukawa`), and the classifier must GENUINELY return `RETURN_NOGO` on
  it (computed, NOT asserted). (A single deliberately-bad POINT is too weak, and a bare Robin/Dirichlet limit is INADMISSIBLE under
  §0.5/§0.6 — the bad family must be the DERIVED admissible delocalizing-warp family.) **DEMOTED:**
  `not_reachable_by_admissible_transmitting` is NO LONGER an acceptable/passing status — if the delocalizing-warp NOGO cannot be
  ADMISSIBLY exhibited, the gate is NOT honestly able-to-fail, so that outcome routes to a FAILING/unsound status (or `INCONCLUSIVE`),
  NEVER a pass. The control fails (and the directive is unsound) if `no_go_reachable.status` is anything other than
  `reachable_RETURN_NOGO`.
- **Tension is real (REFRAMED by §0.6 — the tension is the FALLOFF under (de)localization, NOT a Z-sign veto).** [The v2 Z-based
  absorbing↔confining witness is SUPERSEDED: under §0.6, `Z<0` is the drain premise and the radiation/Neumann `Z=0` limit is the
  no-drain degenerate (INADMISSIBLE), so an "A-pass-but-B-fail via `Z≥0`" witness is no longer the operative tension.] The GENUINE,
  able-to-fail tension is now the Check-B FALLOFF: a LOCALIZING geometry (flat/RS-like) gives `p=2` (B-pass) while a DELOCALIZING
  (anti-localizing) warp gives `p≠2` (B-fail) — both DC-sink admissible. Report tri-state
  `status: witnessed | not_witnessed | not_applicable`: `witnessed` requires exhibiting BOTH an admissible localizing completion with
  `p=2` AND an admissible delocalizing-warp completion with `p≠2` (the §E no_go_reachable family). Only `witnessed` may be cited as
  proving the localization tension is genuine; otherwise report NOT established.
- **Quadrupole-survives:** assert the ℓ=2 DtN kernel RETAINS `p_raw=5` (Burke–Thorne) and is NOT killed by the return machinery.
  FORBID imposing `Q2=0` or any all-multipole cancellation. If the slab return induces a `𝒯₂` acting on ℓ=2, report its effect
  SEPARATELY (`T2_effect:` field, §A) — do NOT auto-fail on it. **NOGO trigger:** if the ONLY admissible structures that cancel ℓ=0,1
  ALSO suppress ℓ=2 below `p=5`, that fires `RETURN_NOGO` (§D.1) — the GW channel that must survive was killed.
- **Return-necessity (RUNG-CONDITIONED — do NOT demand an order change unconditionally):** the return must be load-bearing, but HOW is
  conditioned on the §D rung. For `RETURN_CONSISTENT_WINDOW`: assert a `p_res` ORDER INCREASE vs the raw `O(ω¹)`/`O(ω³)` (the return
  suppresses the order). For `RETURN_RESIDUAL_PREDICTION`: do NOT require an order change — a valid case has `𝒯₀(0)<1` ⇒ monopole stays
  `p_res=1` with a bounded COEFFICIENT; instead assert the computed-return residual coefficient is correctly substituted and reported,
  and is tied to the `1/r²` strength (the signed `Z = −M0(1−𝒯₀(0))` in the local channel). The return is still load-bearing — via the
  coefficient, not the order.
- **Static-dynamic-consistency (NOTHING is static; v2 §0.5.3 — TWO INDEPENDENT computations, NOT variable-vs-itself):** this was REJECTED
  twice as a tautology. The control MUST compare TWO GENUINELY SEPARATE computations: (i) a `dynamic_limit_trace` that solves the ACTUAL
  dynamic Green's function `G(ω,r)` and THEN takes `ω→0`, extracting `dynamic_limit_exponent`; vs (ii) a `static_solve_trace` that sets
  `ω=0` in the same full coupled operator/derived BC BEFORE solving, extracting `static_solve_exponent`. Assert they AGREE; if they
  disagree, a coupling was deleted → `INCONCLUSIVE` (not trusted). Emit `static_dynamic_consistency: {agree: bool,
  dynamic_limit_exponent: <p>, static_solve_exponent: <p>, dynamic_trace_id: <id>, static_trace_id: <id>, independent_extractions: true}`
  and the scripts MUST FAIL if the two trace IDs / source hashes COINCIDE (proof the two exponents came from one computation read twice).

## §F. Discipline

- **Dual-engine** (SymPy + Mathematica); both exit 0; `engine_agreement: PASS` on shared expressions
  ([[feedback-dual-engine-required]]); Mathematica must RUN (`--sandbox danger-full-access` for the execution Codex if it runs
  Mathematica). **Units restored** for every dimensional claim — SymPy dimensional homogeneity with `a, c_s, m, K` restored, do not
  trust the `G=c=c_s=1` pins ([[feedback-dimensional-consistency-check]]). Verdicts COMPUTED, never hardcoded
  ([[feedback-negative-verdict-short-circuit]]). No commentary scripts ([[feedback-no-fake-scripts]]). Scripts `timeout 600`
  ([[feedback-script-timeout-policy]]).

## §G. Out of scope (track-3-downstream or other sectors)

- The **full nonlinear** brane↔bulk return closure (the moving-throat PDE with the return wired in) — this gate tests EXISTENCE of an
  admissible *linear-response* window; the nonlinear closure is downstream.
- The brane parent action / material closure; the EM / charge sectors; the `λγ` cone-lock.
- Re-deriving the pathA_28 raw amplitudes or the ℓ=2 Burke–Thorne result (REUSE both).
- Deciding the absolute inter-brane spacing `d` from first principles (it is a calibration input — [[project-calibrated-pde-goal]]);
  this gate constrains `d` (a window), it does not derive its value.

## Review log
- v0 (2026-06-25): drafted (track 3, gate 1). User design choices folded: return-law = postulate slab structure freely + COMPUTE the
  return response from continuity (no free knob); first gate = cancellation existence (Check A) PLUS `1/r²` survival (Check B), with
  the absorbing↔confining tension as the able-to-fail crux and `same_structure` enforced (§C). Closes pde_ledger open-item #9 +
  tests the gravity-range item. GLM is ON for this directive. Pending: Codex design-review (xhigh) → GLM tertiary → fold → Codex
  re-confirm → user gate → execute dual-engine → tri-review.
- v0-r1 (2026-06-25): Codex design-review (xhigh, gpt-5.5) → **VERDICT: NOT_SOUND** (11 items). FOLDED the 5 [MECHANICAL] items:
  (#6) machine-checkable transfer-function provenance fields + `forbidden_fit_flags`/`fit_to_cancellation:false` (§A);
  (#9) machine-asserted same-structure via canonical `structure_id` hash + `limiting_map` (§A);
  (#10) dual-engine tractability bound — restrict executable family to flat/named-RS warp + local Robin BC, series to O(ω⁴)/O(ω²) (§B);
  (#7) no-go-reachable control now runs the FULL §D classifier on a bad-but-ADMISSIBLE restricted family and must return `RETURN_NOGO` (§E);
  (#8) tension-is-real control now tri-state `witnessed|not_witnessed|not_applicable`, only `witnessed` is citable (§E + §A controls field).
  **NOT FOLDED — 6 [PHYSICS]/[MATH] items ESCALATED to user/orchestrator sign-off** (these change the physics, not wording):
  (#1 PHYSICS) `𝒯_ℓ(0)=1` over-anchored for ℓ=1: mass conservation gives `T0(0)=1` but does NOT force `R1_i(0)=−D1_i`; the dipole DC anchor
  must be DERIVED from first-moment continuity (bulk flux + return centroid), and may legitimately FAIL Check A at raw dipole order;
  (#2 PHYSICS) Check B may solve an isolated sink while Check A carries a return — in a compact slab a sink+return pair can erase the zero
  mode → Yukawa/dipole falloff not 1/r²; Check B must use the SAME static source+return accounting that defines `T0(0)` and output the
  transverse zero-mode source integral;
  (#3 PHYSICS) §C still permits postulate-to-win via arbitrary wall impedance / ω-dependent BC / warp tuned to land `𝒯_ℓ`'s Taylor
  coefficients; firewall should require the BC family to be local/causal/passive/ω-independent (unless from declared boundary DOF) and
  independent of `M0/D1/Q2`, and require an OPEN stable window not an isolated tuned root (unless symmetry-protected);
  (#4 PHYSICS) allowing "suppressed / below observational reach" residuals can make `RETURN_NOGO` unreachable by choosing small d/λ —
  proposes reserving `RETURN_CONSISTENT_WINDOW` for formal `p_residual≥5` (exact-order) and routing bounded-but-lower-order residuals to a
  SEPARATE conditional/phenomenological verdict that does NOT close the pde-ledger item;
  (#5 MATH) make the order test explicit: with `T_l=1+O(ω^ν_l)`, `p_res(0)=1+ν₀`, `p_res(1)=3+ν₁` ⇒ formal pass needs `ν₀≥4` and `ν₁≥2`;
  a generic transit delay `T=1−iωτ+…` gives `p_res(0)=2`, `p_res(1)=4` (i.e. generic delay FAILS) — and replace "outgoing power O(ω⁵)"
  with "leading ω-power of the outgoing amplitude coefficient";
  (#11 PHYSICS) "Quadrupole-survives" overstated: the needed control is that the ℓ=2 DtN kernel stays O(ω⁵) and no tautological all-moment
  kill — forbid imposing `Q2=0`/all-multipole cancellation; if a slab `T2` is computed report its ℓ=2 effect separately rather than auto-failing.
  Per gauntlet contract: physics/math items are NOT self-resolved; HALT for sign-off before re-confirm round.
- v0-r2 (2026-06-25): all 6 [PHYSICS]/[MATH] items SIGNED OFF + FOLDED, plus one ADDITIONAL user no-static-framing fix. Folds:
  - **UNIFIED DC-anchor fix (#1+#2+#5):** `𝒯_ℓ(0)` is no longer anchored to 1 by fiat anywhere — it is a COMPUTED output. `𝒯₀(0)` =
    monopole DC return *fraction* (from zeroth-moment/mass continuity; may be `<1`); `𝒯₁(0)` from first-moment continuity (bulk flux +
    return centroid), not forced to 1; Check A may fail at raw dipole order. Explicit order test added: `A_ℓ_res = …·(1−𝒯_ℓ(ω))`,
    `p_res(0)=1+ν₀`, `p_res(1)=3+ν₁`; strict pass needs `𝒯₀(0)=1 ∧ ν₀≥4` and `𝒯₁(0)=1 ∧ ν₁≥2`; generic transit-delay (ν=1) gives
    `p_res=2/4` → FAILS (gate has teeth). Wording fix: "below O(ω⁵)" → "leading ω-power of outgoing amplitude coefficient ≥ p=5".
    (§0 REUSE bullet, §0(A), §A transfer_function/admissibility fields, §B.2, §B.3, §E DC-value control.)
  - **#2 Check B crux:** Check B now uses the SAME source+return accounting that defines `𝒯₀(0)`, outputs the transverse zero-mode source
    integral, observable = IN-BRANE `v_r` (3D Gauss on in-brane net flux). Monopole-level tension stated plainly: A-strict wants
    `𝒯₀(0)=1`, B wants `𝒯₀(0)<1` — shared pivotal quantity, outcome NOT pre-asserted. (§B.4.)
  - **#3 firewall:** BC family must be local/causal/passive/ω-independent (unless from a DECLARED boundary DOF whose eqn is given +
    ω-dep DERIVED), independent of M0/D1/Q2; pass must be an OPEN stable window not a tuned root (unless symmetry-protected). New flags
    `bc_is_local_causal_passive`, `bc_omega_dependence_source`, `window_is_open`, `window_symmetry_protected`. (§C + §A.)
  - **#4 four-rung ladder:** `RETURN_NOGO` / `RETURN_CONSISTENT_WINDOW` (strict `p_res≥5` both ℓ + B `p=2`; ONLY this closes #9) /
    NEW `RETURN_RESIDUAL_PREDICTION` (`1/r²` holds + bounded `p_res<5` residual recorded as falsifiable prediction; does NOT close #9) /
    `INCONCLUSIVE`. (§D, §A reconciliation/residual_prediction.)
  - **#11 quadrupole control:** assert ℓ=2 kernel retains p=5, forbid `Q2=0`/all-multipole cancellation, `T2_effect` reported
    separately; NOGO trigger if the only ℓ=0,1-cancelling structures also suppress ℓ=2 below p=5. (§E + §A + §D.1.)
  - **NO-STATIC-FRAMING fix (user, additional):** the project's #1 relapse — Check B as drafted ("solve the static 1/r² flow") was
    relapsing into a frozen BVP. Now: Check B = `ω→0` LIMIT of the SAME dynamic object as Check A (renamed "Quasi-static (ω→0 limit)
    1/r² survival"); the linearization background MUST be a genuine dynamical STEADY-STATE CIRCULATION (`steady_state: {exists, type:
    circulation, provenance}`; no-circulation ⇒ INADMISSIBLE, bears on NOGO); §C forbids any frozen/decoupled static solve; new §E
    `static-dynamic-consistency` control asserts the `ω→0` limit reproduces the Check-B exponent or the result is INCONCLUSIVE.
    (§0 new "NOTHING is static" bullet, §A slab_structure+steady_state, §B framing note, §B.4, §C, §D.4, §E.)
  Pending: Codex re-confirm (xhigh) to verify SOUND, then orchestrator runs GLM.
- v0-r2 re-confirm (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) → **VERDICT: NOT_SOUND** (4 NEW items — internal-consistency
  contradictions exposed/introduced by the v0-r2 folds; NOT self-resolved per gauntlet contract — ESCALATED for sign-off):
  (R2-#1 PHYSICS) `RETURN_CONSISTENT_WINDOW` is foreclosed BY THE MACHINE FIELDS: §A `reconciliation` makes `A_strict_pass ⇔ 𝒯₀(0)=1`
  and §B.4 makes `B_pass ⇔ 𝒯₀(0)<1`, so the two are mutually exclusive in the asserted fields — yet the prose says "DO NOT pre-assert
  the outcome." The "neighboring-brane / bulk-flow contribution" escape is prose-only and cannot be represented. FIX: define `B_pass`
  from the COMPUTED full `zero_mode_source_integral` of the shared `ω→0` object (which in the simple local source+return channel reduces
  to `M0(1−𝒯₀(0))`, but extra brane/bulk contributions must be emitted as separate derived terms with provenance); if the full integral
  ALWAYS reduces to `M0(1−𝒯₀)`, then strict-A+B is honestly `RETURN_NOGO`, not an open joint window. (Decides whether the consistent
  window is even representable, or whether the gate is structurally a no-go for the local channel — a real design call.)
  (R2-#2 PHYSICS) no-steady-circulation is double-routed: §0 says "none for ANY admissible structure ... bears on `RETURN_NOGO`" but
  §D.4 routes the same condition to `INCONCLUSIVE`. FIX: split — "existence UNDECIDABLE in-scope" → INCONCLUSIVE; "COMPUTED no
  circulation for all admissible structures" → RETURN_NOGO; "a candidate lacks circulation" → that candidate is inadmissible.
  (R2-#3 PHYSICS) the `RETURN_RESIDUAL_PREDICTION` rung is undermined by the Return-necessity control: §E requires "WITH the return it
  changes ORDER," but a valid residual-prediction case can be SAME-order (`𝒯₀(0)<1` ⇒ monopole `p_res=1` with a bounded coefficient
  tied to 1/r² strength). FIX: require an order INCREASE only for `RETURN_CONSISTENT_WINDOW`; for `RETURN_RESIDUAL_PREDICTION` the
  control should instead assert the computed-return residual expression/coefficient is correctly substituted+reported (same-order
  bounded residual allowed).
  (R2-#4 MATH) `ν_ℓ` is ambiguous: §A defines `nu_ℓ = ord_ω(1−𝒯_ℓ)` but §B.3 writes `𝒯_ℓ=𝒯_ℓ(0)+O(ω^{ν_ℓ})` — if `𝒯_ℓ(0)≠1` the
  first non-DC correction can be order 1 while `ord_ω(1−𝒯_ℓ)=0`, so scripts could OVERSTATE suppression. FIX: use two symbols —
  `σ_ℓ = ord_ω(𝒯_ℓ − 𝒯_ℓ(0))` and `ν_ℓ = ord_ω(1−𝒯_ℓ)` (with `ν_ℓ=0` when `𝒯_ℓ(0)≠1`); use ONLY `ν_ℓ` in `p_res`.
  Per gauntlet contract: these physics/math items are NOT self-resolved; HALT for sign-off. (All 4 are internal-consistency fixes that
  tighten the ladder/controls/fields — none re-open the conceptual design; once signed off + folded, a third Codex re-confirm verifies.)
- v0-r3 (2026-06-25): all 4 re-confirm items (R2-#1..#4) SIGNED OFF + FOLDED. Folds:
  - **R2-#1 (B_pass is a COMPUTED output; able-to-fail BOTH ways):** removed the hardwired definitions `A_strict_pass ⇔ 𝒯₀(0)=1` and
    `B_pass ⇔ 𝒯₀(0)<1`. Now `A_strict_pass ⇔ COMPUTED p_res(ℓ=0,1)≥5`; `B_pass ⇔ COMPUTED falloff from Z is p=2 (Z≠0)`. `Z` =
    `zero_mode_source_integral` is a COMPUTED sum of separate provenance-tagged terms `Z = Z_throat_sink + Z_return + Z_replenishment`
    (replenishment = the distributed bulk→brane areal leak forced by `steady_state.exists`). In the simple LOCAL channel `Z` reduces to
    `M0(1−𝒯₀(0))`. Honest-either-way in §D.1: if the COMPUTED full `Z` always reduces to local (`Z_reduces_to_local: true`) ⇒ strict-A+B
    contradictory ⇒ `RETURN_NOGO` **WITHIN-THE-TRACTABLE-FAMILY** (not claimed absolute — escape mechanisms outside scope = follow-up
    gate); if an extra computed term makes `Z≠0` while `𝒯₀(0)=1` ⇒ window representable, computation decides. New field
    `Z_reduces_to_local: bool` + `Z_terms{}`. (§A static_falloff_B/reconciliation, §B.4, §D.1, §D.2.)
  - **R2-#2 (no double-routing of no-steady-circulation):** three distinct cases — a CANDIDATE with computed no-circulation = INADMISSIBLE
    (dropped, not a verdict); COMPUTED none-for-all-admissible = `RETURN_NOGO`; UNDECIDABLE-in-scope = `INCONCLUSIVE`. `steady_state.exists`
    now `bool|undecidable`. (§0 steady-state bullet, §D.1, §D.4.)
  - **R2-#3 (return-necessity control rung-conditioned):** order INCREASE required only for `RETURN_CONSISTENT_WINDOW`; for
    `RETURN_RESIDUAL_PREDICTION` no order change required — assert the computed-return residual COEFFICIENT is correctly substituted +
    tied to the 1/r² strength (`Z`/`M0(1−𝒯₀(0))`); the return is load-bearing via the coefficient. (§E Return-necessity.)
  - **R2-#4 (two order symbols):** `σ_ℓ = ord_ω(𝒯_ℓ−𝒯_ℓ(0))` (first non-DC correction) and `ν_ℓ = ord_ω(1−𝒯_ℓ)` (deviation-from-1,
    `ν_ℓ=0` when `𝒯_ℓ(0)≠1`); ONLY `ν_ℓ` enters `p_res`; both emitted in YAML. (§A transfer_function, §B.3.)
  Pending: 3rd Codex re-confirm (xhigh) to verify SOUND, then orchestrator runs GLM.
- v0-r3 re-confirm (2026-06-25): 3rd Codex re-confirm (xhigh, gpt-5.5) → all R2-#1..#4 folds confirmed "substantively resolved"
  (B_pass computed from full Z; Z_replenishment required derived/provenanced not dialed; NOGO scope honest within the tractable family;
  return-necessity rung-conditioned; σ_ℓ/ν_ℓ no longer overstates suppression). Verdict NOT_SOUND on ONE residual [MECHANICAL] schema
  mismatch only: §A `steady_state.exists` was still typed `bool` but §0/§D require `bool|undecidable` (the in-scope-intractable→
  INCONCLUSIVE route was unrepresentable in the YAML schema). FOLDED: §A `steady_state` now
  `{exists: bool|undecidable, type: circulation, provenance: <derivation-or-obstruction>}` with the full §0/§D routing inlined.
  No physics/math items outstanding. A follow-up confirm round (xhigh) flagged ONE further [MECHANICAL] residue — the parallel §0
  `steady_state` schema still read `provenance: <derivation>` — now folded to `<derivation-or-obstruction>` (matching §A/§D.4), Codex
  confirming "no other outstanding blocker." Final green-confirm (xhigh, gpt-5.5) → **VERDICT: SOUND** (§0/§A `steady_state` schemas
  match, §D routing consistent, "no other blocker remains"). Directive is design-review-clean. Pending: orchestrator runs GLM tertiary.
- v0-r4 (GLM, user-run out-of-band, 2026-06-25): GLM tertiary review (executed out-of-band by the user, findings relayed) →
  **NOT_SOUND** with 2 items, both SIGNED OFF; GLM explicitly confirmed the rest SOUND (physics chain, DC anchors, no-static framing, §C
  firewall, all-rung reachability, controls, scope). FOLDED:
  - **GLM-#1 [PHYSICS, load-bearing] — B_pass needs the SIGN of the flow (gravity = inflow):** the old `B_pass ⇔ Z≠0` omitted the sign.
    Added an explicit SIGN CONVENTION (§B.4/§A): `Z` = net rate of medium ADDED at the throat, so `Z<0` = inward = GRAVITY, `Z>0` =
    outflow = not gravity; signed terms `Z_throat_sink=−M0`, `Z_return=+M0·𝒯₀(0)`, `Z_replenishment_localized ≥0`. Criterion is now
    **`B_pass ⇔ (p=2 AND Z<0)`** — falloff and sign are TWO distinct computed conditions (`Z≠0` ⇏ `p=2`; the old gloss removed). Scripts
    MUST compute `Z` and VERIFY the sign (assert `Z≥0` when `𝒯₀(0)=1`), not assert it in prose. At strict-A (`𝒯₀(0)=1`) the computed
    `Z = Z_replenishment_localized ≥ 0` ⇒ B_fail ⇒ the areal-leak escape is foreclosed BY THE SIGN. §E "tension is real" updated: the
    A-pass-but-B-fail limit IS now exhibitable ⇒ `witnessed`. §B.4 "forecloses NEITHER outcome" softened: all 4 rungs stay formally
    representable, but `RETURN_CONSISTENT_WINDOW` is physically EXPECTED-UNREACHABLE (sign argument), so the gate primarily decides
    `RETURN_NOGO` vs `RETURN_RESIDUAL_PREDICTION` — computation must verify the sign, not pre-assume it. New YAML `Z_sign:`, `p_eq_2:`.
    (§A static_falloff_B/reconciliation, §B.4, §B.5, §D.1/.2/.3, §E tension + return-necessity.)
  - **GLM-#2 [SCOPING] — Z_replenishment re-scoped/relabeled:** REDEFINED `Z_replenishment_localized` as the LOCALIZED contribution of
    the areal leak to the monopole at `r≫d` (the part affecting 1/r²), `=0` in the simple local channel (no localized enhancement) — NOT
    "required/forced by the steady state." It is an INDEPENDENT cosmological process (conceptual_foundation §5); the UNIFORM far-field
    areal-leak background does NOT enter `Z`/`B_pass` and is emitted as a SEPARATE `Z_uniform_background` observable. Net (§D): within the
    tractable family `Z = −M0(1−𝒯₀(0))`, so B_pass (`p=2 AND Z<0`) needs `𝒯₀(0)<1`, contradicting A_strict (`𝒯₀(0)=1`) — window
    foreclosed BY PHYSICS (computed sign), not by definitional construction; gate stays honestly able-to-fail (NOGO vs RESIDUAL_PREDICTION).
    (§A static_falloff_B, §B.4.)
  All prior folds (R2-#1..#4, no-static fix, §C firewall) intact. Pending: Codex re-confirm (xhigh) → when green, bump STATUS to v1.
- v0-r4 re-confirm (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) of the GLM folds → **VERDICT: NOT_SOUND** (3 items). The GLM sign-fold,
  as applied, introduced TWO genuine [PHYSICS] regressions (NOT self-resolved per gauntlet contract — ESCALATED for sign-off); the one
  [MECHANICAL] item was folded:
  (R4-#1 PHYSICS — rung-reachability regression) §D.1 still routes `Z_reduces_to_local: true` to UNCONDITIONAL `RETURN_NOGO` (because
  strict-A forces `Z≥0`); by top-down precedence this PREEMPTS `RETURN_RESIDUAL_PREDICTION` in the very case the GLM fold calls EXPECTED
  (`𝒯₀(0)<1` ⇒ `p=2 AND Z<0` + a bounded lower-order residual). FIX: local-`Z` reduction should eliminate ONLY
  `RETURN_CONSISTENT_WINDOW`; return `RETURN_NOGO` only if NO admissible `(p=2 AND Z<0)` bounded-residual branch exists. (GLM-#1 not
  fully resolved; the sign-fold accidentally made the expected RESIDUAL_PREDICTION outcome unreachable.)
  (R4-#2 PHYSICS — R2-#1 anti-tautology regression) §A/§B.4/§D.2 say `RETURN_CONSISTENT_WINDOW` stays "formally representable via a
  declared boundary DOF," but the machine `Z_terms` schema has ONLY `−M0`, `+M0·𝒯₀(0)`, `Z_replenishment_localized ≥0` — none of which
  can give `Z<0` at strict-A, so the escape is PROSE-ONLY (re-opening the R2-#1 hole). FIX: either add an explicit script-validated
  optional SIGNED `Z_terms` entry for a declared boundary/return DOF (with governing equation + provenance), OR honestly admit the
  window is NOT representable in this tractable directive (and state that plainly rather than claiming formal representability).
  (R4-#3 MECHANICAL — FOLDED) §E tension witness said `Z = Z_replenishment_localized = 0` but the active rule is `≥ 0` (`=0` only in the
  simple local channel); corrected to `Z = Z_replenishment_localized ≥ 0` (and `=0` in the simple local witness).
  Per gauntlet contract: the 2 physics items are NOT self-resolved; HALT for sign-off. STATUS NOT bumped to v1 (green not reached).
- v0-r5 (2026-06-25): both R4 items SIGNED OFF + FOLDED (corrections to how the GLM sign-fold landed, not new design). Folds:
  - **R4-#1 (precedence — local-Z must NOT unconditionally fire NOGO):** added an explicit §D PRECEDENCE preamble — `Z_reduces_to_local:
    true` ELIMINATES ONLY `RETURN_CONSISTENT_WINDOW`, it does NOT fire NOGO. Rewrote §D.1: `RETURN_NOGO` fires ONLY if NO admissible
    structure yields ANY `(p=2 AND Z<0)` bounded-residual branch (no inward 1/r² at all: `p≠2`/Yukawa for all, OR `Z` never `<0`, OR ℓ=2
    suppressed below p=5, OR no steady circulation for any). The EXPECTED `𝒯₀(0)<1` ⇒ inward `1/r²` + bounded residual now correctly
    routes to `RETURN_RESIDUAL_PREDICTION` (rung 3), no longer preempted. (§D preamble, §D.1, §D.3/.4 cross-refs.)
  - **R4-#2 (keep CONSISTENT_WINDOW escape REPRESENTABLE, not prose-only):** added an OPTIONAL script-validated SIGNED `Z_terms` entry
    `Z_boundary_dof: {governing_equation, provenance, dc_sign_derivation}` (the only term that can give `Z<0` at strict-A). Computation
    MUST populate+compute it from a declared boundary/return DOF or set `none`; new `reconciliation.Z_boundary_dof_status:
    populated|none|intractable`. §D.2: CONSISTENT_WINDOW reachable ONLY if a COMPUTED `Z_boundary_dof` yields `(Z<0 AND p=2)` at
    strict-A; EXPECTED-UNREACHABLE (conservation expected to force `Z_boundary_dof ≥0`) but REPRESENTABLE and computation-decided;
    intractable boundary-DOF branch routes THAT rung to INCONCLUSIVE (NOT NOGO). Closes the R2-#1 anti-tautology gap (escape
    representable, foreclosure computed not assumed). (§A static_falloff_B + reconciliation, §B.4, §D.2, §D.4.)
  All prior folds (R2-#1..#4, no-static, §C firewall, GLM-#1 sign, GLM-#2 scoping) intact. Pending: Codex re-confirm (xhigh) → when
  green, bump STATUS to v1 (Codex+GLM SOUND).
- v0-r5 re-confirm (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) of the R4 folds → **VERDICT: SOUND.** Confirmed both R4 fixes genuine:
  `Z_reduces_to_local: true` now eliminates ONLY `RETURN_CONSISTENT_WINDOW` (no longer preempts the expected `𝒯₀(0)<1`
  `RETURN_RESIDUAL_PREDICTION` path); `RETURN_NOGO` reachable only on true total failure of any inward `(p=2 AND Z<0)` bounded-residual
  branch; the `Z_boundary_dof` escape is machine-representable (signed `Z_terms` slot + required `{governing_equation, provenance,
  dc_sign_derivation}` + `populated|none|intractable` status + §C blocking free fitted DOFs), its intractable case routed to INCONCLUSIVE
  not NOGO. All prior folds confirmed still holding. Codex flagged only 2 OPTIONAL [MECHANICAL] nits, both FOLDED: (1) §B.5 "or the
  precise incompatibility (no-go)" → "or the appropriate §D rung (RESIDUAL_PREDICTION/NOGO/INCONCLUSIVE; a bounded inward-gravity
  residual is RESIDUAL_PREDICTION not NOGO)"; (2) §B.4 `Z = Z_replenishment_localized = 0` → `≥ 0` (and `=0` only in the simple local
  witness), aligning with §E. **STATUS bumped to `v1 (Codex+GLM SOUND)`** — design-review complete. NEXT: user gate → execute
  dual-engine → tri-review (orchestrator owns GLM/user steps).
- v2 RE-SCOPE (2026-06-25): BOTH prior executions REJECTED at tri-review (v1 hardcoded Check B's 1/r²; remediation still had (i)
  variable-vs-itself `static_dynamic_consistency`, (ii) `p=2` as spectral-lookup/bare-literal not cross-engine-verified, (iii) THE CRUX —
  1/r² survival hinged on a POSTULATED `m`-impedance BC = postulate-to-win). Orchestrator authored a new §0.5 block + STATUS update: (1)
  DERIVE the neighbor-brane BC from a physical open boundary (Sommerfeld/Bloch), not a tuned impedance; (2) `p` from a real radial solve +
  in `engine_agreement`; (3) `static_dynamic_consistency` compares two INDEPENDENT computations.
- v2-r1 (2026-06-25): Codex design-review (xhigh, gpt-5.5) of the v2 re-scope → **VERDICT: NOT_SOUND** (8 items: 2 PHYSICS, 1 MATH, 5
  MECHANICAL). FOLDED the 5 [MECHANICAL] items:
  (#2) §C ω-dependence contradiction resolved — `bc_omega_dependence_source` enum extended to
  `none|sommerfeld_radiation|bloch_periodic|declared_dof`; Sommerfeld/Bloch allowed WITH emitted derivation+limiting-map, tuned `m`/`ω`
  impedance still forbidden (§C bullet + §A flag);
  (#3) `engine_agreement.checked_quantities` MUST include `p`, `p_eq_2`, `dynamic_limit_exponent`, `static_solve_exponent`,
  `zero_mode.r_dependence`, `green_function` — scripts FAIL if absent (§A dual-engine deliverable);
  (#4) `static_dynamic_consistency` now requires TWO independently-extracted traces (`dynamic_limit_trace` vs `static_solve_trace`) with
  distinct `dynamic_trace_id`/`static_trace_id` + `independent_extractions: true`; scripts FAIL if trace IDs/source hashes coincide (§E);
  (#5) all "postulate the wall BC"/`{d, warp, BC}`/"same Robin BC" text reconciled to "geometry postulated; neighbor-brane BC DERIVED per
  §0.5" (§A slab_structure, §A `{d,warp,derived_BC_branch/physical_completion}` ×5, §A structure_id hashes bc_derivation_type+derived BC
  eqns, §B.1 step 1, §B.4, §C SAME-structure + Green's-function bullets);
  (#6) tractability bound no longer restricts to "LOCAL Robin/impedance wall BC" — now flat/named-RS warp + Sommerfeld/Bloch/declared-DOF
  derived BCs; bare Robin admissible only if DERIVED as one of these limits (§B).
  **NOT FOLDED — 2 [PHYSICS] + 1 [MATH] ESCALATED for orchestrator/user sign-off (NOT self-resolved per gauntlet contract):**
  (v2-#1 PHYSICS — the load-bearing one) §0.5.1 as written STILL RELOCATES the postulate-to-win: it allows CHOOSING among
  radiation/Sommerfeld vs periodic/Bloch without proving WHICH is the actual return boundary for this medium/stack — if those branches
  give different zero-mode outcomes, picking the favorable one is selection-by-BC. Codex fix: require deriving the ACTUAL physical BC from
  the declared completion, OR (if both are physically admissible in-scope) run BOTH as separate branches and report BOTH verdicts, with NO
  favorable-branch promotion absent an independently declared physical selector.
  (v2-#7 PHYSICS) §E `no_go_reachable` + §C control-enforcement bullet still cite a bare "fully-absorbing/Robin limit" as the bad family —
  under v2 a bare Robin/Dirichlet limit is NOT admissible unless DERIVED from a physical completion. Codex fix: the NOGO control must
  construct a bad-but-admissible branch with a DERIVED zero-mode-killing BC (e.g. a pinned/absorbing neighbor-brane DOF with governing
  equations); if none can be exhibited, `no_go_reachable` FAILS (do not use an arbitrary Robin limit).
  (v2-#8 MATH) §0.5.1 conflates "zero mode dies → Yukawa/`p=3`" — those are DISTINCT spectral failures: continuum leakage gives
  `v_r~r^{-3}`; a gapped transverse spectrum gives Yukawa. Codex fix: classify the computed spectrum explicitly and state which was
  derived (both are B-fail inputs).
  Per gauntlet contract: physics/math items NOT self-resolved; HALT for sign-off. STATUS NOT bumped (green not reached). Pending:
  sign-off on v2-#1/#7/#8 → fold → Codex re-confirm (xhigh) → orchestrator runs GLM.
- v2-r2 (2026-06-25): all 3 v2 items SIGNED OFF + FOLDED. Folds:
  - **v2-#1 (CRUX — run BOTH boundaries, report both, no favorable-branch promotion):** §0.5.1 rewritten — the DECLARED completion is the
    multi-brane STACK, so **Bloch/periodic is the PRIMARY** physical return BC and **radiation/Sommerfeld the CROSS-CHECK**; DERIVE both
    (emit derivation + limiting map), RUN BOTH, EMIT `verdict_bloch`/`verdict_radiation` with their computed `p`/`Z_sign`/`spectrum_class`.
    FORBID promoting the favorable branch. New §D.0 cross-branch gate: headline is the matching per-branch rung if BOTH agree on
    `(p=2 AND Z<0)` (robust) or BOTH B-fail (`RETURN_NOGO`), else NEW outcome **`BC_DEPENDENT`** (branches disagree → name the deciding
    completion, do NOT pick favorable; ≈ INCONCLUSIVE-with-condition, does NOT close #9). New YAML `branch_results{}` + `branch_agreement{}`.
    (§0.5.1, §D header + §D.0, §A branch_results/reconciliation.)
  - **v2-#7 (no_go_reachable uses a DERIVED admissible zero-mode-killing structure):** §E + §C control bullet — a bare Robin/Dirichlet limit
    is inadmissible under §0.5; the NOGO control must construct a DERIVED admissible structure (PREFERRED: a delocalizing/anti-localizing
    warp that fails to localize the zero mode — the conceptual_foundation §7 braneworld caveat — with its warp equation emitted;
    ALTERNATIVE: a pinned/absorbing neighbor-brane DOF with governing equations). If NO admissible transmitting open return kills the zero
    mode, report that explicitly ("NOGO reachable on a gapped/continuum spectrum but not by any admissible transmitting open return; 1/r²
    robust for localizing geometries") — but PREFER the derived delocalizing-warp family so NOGO is admissibly reachable. (§E no_go_reachable, §C.)
  - **v2-#8 (classify the spectrum explicitly — no Yukawa/p=3 conflation):** new §0.5.4 — per branch emit
    `spectrum_class: {normalizable_zero_mode → p=2 | continuum_p3 → p=3 gapless leakage | gapped_yukawa → exponential cutoff}`; both `p=3`
    and Yukawa are DISTINCT B-fail inputs and the report states which was derived. Removed the §0.5.1 + §B.4 "dies → Yukawa/p=3" gloss.
    (§0.5.4, §0.5.1, §B.4.)
  All prior folds (R2-#1..#4, no-static, §C firewall, GLM-#1 sign, GLM-#2 scoping, v2-r1 mechanical #2..#6) intact. Pending: Codex
  re-confirm (xhigh) → when green, bump STATUS to `v2 (Codex+GLM-pending SOUND)`; orchestrator then runs GLM.
- v2-r2 re-confirm (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) of the 3 v2 folds → **VERDICT: SOUND.** Confirmed the load-bearing
  v2-#1 fix genuine — "the postulate-to-win is NOT relocated": both Bloch + Sommerfeld branches must be derived/run/reported, `PRIMARY`
  does NOT authorize promotion, §D.0 forces `BC_DEPENDENT` on branch disagreement (and it does not close #9). v2-#7 substantively
  resolved (no-go control requires a DERIVED admissible zero-mode-killing mechanism, machine-checkable via emitted governing equations);
  v2-#8 resolved (`continuum_p3` vs `gapped_yukawa` distinct throughout). All prior folds confirmed holding. Codex flagged only 4
  OPTIONAL [MECHANICAL] nits, all FOLDED: (1) §A `controls.no_go_reachable` now `status: reachable_RETURN_NOGO |
  not_reachable_by_admissible_transmitting` (carries the fallback finding); (2) §B.5 headline list now includes `BC_DEPENDENT` (+ routes
  via §D.0); (3) §A `branch_agreement` now also emits `bloch_verdict`/`radiation_verdict`/`verdict_agree` (machine-explicit
  "matching rung"); (4) §B.4 Green's-function bullet now names the three §0.5.4 `spectrum_class` values explicitly (no Yukawa/continuum
  blur). **STATUS bumped to `v2 (Codex+GLM-pending SOUND)`** — Codex design-review GREEN. NEXT: orchestrator runs GLM tertiary.
- v3 RE-SCOPE (GLM tertiary, user-run out-of-band, 2026-06-25): GLM (authoritative, relayed) → NOT_SOUND, LOAD-BEARING: the v2 "both
  branches must give (p=2 AND Z<0)" gate is STRUCTURALLY BROKEN — the radiation/Sommerfeld BC degenerates to Neumann (REFLECTING) at
  `ω→0`, deterministically `Z=0` at DC, so the favorable verdict was UNREACHABLE (postulate-to-win relocated into an always-fail
  cross-check). User+orchestrator AGREED resolution (grounded in conceptual_foundation §5): **`Z<0` is the drain PREMISE, not a test.**
  FOLDED:
  - **§0.6 inserted VERBATIM** (authoritative block) after §0.5: (1) `Z<0` = admissibility PREMISE (a non-draining throat is not a
    particle); admissible Check-B return MUST be a DC SINK; a DC-reflecting boundary (`Z=0`) is the no-drain degenerate, INADMISSIBLE for
    Check B (not a Z-veto). (2) the model's return is a DC absorber (§5 de-structuring/absorbing PRIMARY + Bloch CROSS-CHECK; finite-or-
    infinite stack OK as long as the neighbor ABSORBS — GLM Finding 3). (3) radiation/Sommerfeld is AC/Check-A ONLY, never a Check-B
    Z-branch. (4) `p=2` is robust-by-geometry for localizing warps; the real teeth = falloff under a DELOCALIZING warp → NOGO;
    `BC_DEPENDENT` only if `p` differs across DC-sink completions (NOT a Z disagreement). (5) the durable falsifiable result is the
    Check-A radiation residual `∝ ε0=1−𝒯₀(0)`.
  - **§0.5.1 marked PARTIALLY SUPERSEDED** (keeps "BC derived not postulated"; the "both branches agree on Z" / radiation-as-Z-branch
    design replaced by §0.6); §0.5.4 retargeted to the DC-sink completions; reqs #2/#3 remain in force.
  - **§B.4 / §B.5 reconciled:** Check B runs under the DC-sink completions, `Z<0` asserted as PREMISE (signed-Z bookkeeping retained as
    accounting), `B_pass ⇔ p=2`; the v2 "sign forecloses the window / `Z_boundary_dof` rescue" interpretation SUPERSEDED; radiation moved
    to AC/Check-A; reconcile compares `p` across DC-sink completions.
  - **§D reframed:** §D.0 is now the DC-SINK COMPLETION gate (RESIDUAL_PREDICTION if all `p=2`; NOGO if all `p≠2`/ℓ=2-killed;
    `BC_DEPENDENT` only on `p`-disagreement across DC-sink completions); §D.1 NOGO trigger is `p≠2` (localization failure), with `Z<0`
    explicitly the premise not a trigger; §D.2 CONSISTENT_WINDOW marked SUPERSEDED/retained-for-record; §D.3 RESIDUAL_PREDICTION is the
    EXPECTED v3 headline; §D.4 + header updated.
  - **§A YAML:** `Z_is_premise: true` added; `branch_results`/`branch_agreement` now compare the DC-sink completions
    (`destructuring_absorbing` vs `bloch_stack`) on `p`; radiation tagged `ac_check_a_only`; `reconciliation` verdict set =
    {RESIDUAL_PREDICTION|NOGO|BC_DEPENDENT|INCONCLUSIVE}, `B_pass ⇔ p=2` under premise; v2 strict-cancellation/`Z_boundary_dof` fields
    retained-for-record only.
  - **§E:** `no_go_reachable` (DERIVED delocalizing-warp family → NOGO) and the two-independent-computation `static_dynamic_consistency`
    + p-in-`engine_agreement` requirements KEPT (the real Check-B teeth); `tension_is_real` REFRAMED — the genuine tension is the
    FALLOFF under (de)localization (localizing `p=2` vs delocalizing `p≠2`), not the superseded Z-sign witness.
  - GLM MECHANICAL findings folded: Finding 2 (`p=2` robust-by-geometry; falsifiable content = falloff-under-delocalization + Check-A
    radiation prediction) stated in §0.6.4/§0.6.5; Finding 3 (RS1 finite-vs-infinite — absorbing neighbor → DC sink regardless) in
    §0.6.2 and reflected in §B/§D.
  STATUS → v3 RE-SCOPE. Pending: Codex re-confirm (xhigh) → when green, bump STATUS to `v3 (Codex+GLM SOUND)`; orchestrator does NOT run
  another GLM round (these ARE GLM's fixes).
- v3-r1 (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) of the §0.6 v3 re-scope → **VERDICT: NOT_SOUND** (3 items). Codex explicitly
  confirmed **"the core §0.6 framing itself is coherent: `Z<0` as drain admissibility is legitimate and NOT a hidden computed win"** (the
  falsifiable content moved to derived `p` + the Check-A radiation residual) and that the STRUCTURAL UNREACHABILITY bug is gone — the
  remaining blockers are orphaned operative text outside §0.6. FOLDED the 2 [MECHANICAL] items:
  (v3-#1) §C control-enforcement bullet still mandated the old `κ0→∞` strict-A / Z-sign `tension_is_real` witness — replaced with the v3
  FALLOFF witness (localizing DC-sink `p=2` vs derived delocalizing-warp `p≠2`), consistent with §0.6.4/§E;
  (v3-#2) §C Green's-function bullet + §B.4 "BC MUST be identical to Check A's" could force the Sommerfeld AC boundary back into Check B —
  reworded to "SAME declared physical completion + machine-checked limiting map," with Check B using ONLY the DC-sink completions and
  Check A reporting radiation as `ac_check_a_only` (§0.6.3).
  **NOT FOLDED — 1 [PHYSICS] item ESCALATED for orchestrator/user sign-off (NOT self-resolved per gauntlet contract):**
  (v3-#3 PHYSICS — able-to-fail weakened) §A/§E still permit `no_go_reachable.status: not_reachable_by_admissible_transmitting` as an
  ACCEPTABLE reported status, leaving the able-to-fail proof OPTIONAL. Codex fix: for v3 green, REQUIRE
  `no_go_reachable.status: reachable_RETURN_NOGO` from the derived delocalizing-warp family; retain `not_reachable_by_admissible_transmitting`
  ONLY as a failing/unsound control outcome (or route to INCONCLUSIVE), NOT as a pass. (This decides whether NOGO must be ADMISSIBLY
  exhibited for the gate to be honestly able-to-fail, vs. accepting a "1/r² robust, NOGO only by inadmissible BC" finding — a real call.)
  Per gauntlet contract: physics item NOT self-resolved; HALT for sign-off. STATUS NOT bumped (green not reached). Pending: sign-off on
  v3-#3 → fold → Codex re-confirm (xhigh) → when green, bump STATUS to `v3 (Codex+GLM SOUND)`.
- v3-r2 (2026-06-25): v3-#3 SIGNED OFF (accepted — Codex's stricter option) + FOLDED. Fold: the `no_go_reachable` control now MANDATES an
  ADMISSIBLY-exhibited NOGO for v3 green — it MUST construct the DERIVED delocalizing (anti-localizing) warp family (drain still transmits
  → `Z<0` premise holds, but the warp fails to localize the zero mode → `p≠2` → `RETURN_NOGO`), emit its warp equation + the COMPUTED
  `spectrum_class` (`continuum_p3` | `gapped_yukawa`), and the §D classifier must GENUINELY return `RETURN_NOGO` (computed, not asserted).
  The §A status enum is now `reachable_RETURN_NOGO | failed_not_admissibly_reachable` — ONLY `reachable_RETURN_NOGO` is acceptable for
  green; the demoted `not_reachable_by_admissible_transmitting` fallback is REMOVED as a passing status (a non-admissibly-reachable NOGO
  is a FAILING/unsound outcome or `INCONCLUSIVE`, never a pass). (§E no_go_reachable, §A controls enum.) All prior folds intact (§0.6
  verbatim, the v3 reconciliation, reqs #2/#3, v3-r1 mechanical folds). Pending: Codex re-confirm (xhigh) → when green, bump STATUS to
  `v3 (Codex+GLM SOUND)` (LAST design-review round per orchestrator; then final read + re-execute).
- v3-r2 re-confirm (2026-06-25): Codex re-confirm (xhigh, gpt-5.5) of the v3-#3 fold → **VERDICT: NOT_SOUND** on ONE residual
  [MECHANICAL] orphan only: §C control-enforcement bullet (line 331) still offered the pinned/absorbing DOF as an ALTERNATIVE that could
  satisfy `no_go_reachable`, contradicting the now-mandatory §A/§E "DERIVED delocalizing-warp family → `reachable_RETURN_NOGO`" rule.
  Codex confirmed "everything else is consistent: §A enum fixed, §E demotes the fallback, §0.6/§D use `Z<0` as premise + falloff-`p` as
  the gate, prior v2/v3 folds intact." FOLDED: §C bullet reworded — for green, `no_go_reachable` MUST use the DERIVED delocalizing-warp
  family (warp equation + computed `spectrum_class` emitted, classifier returns `RETURN_NOGO`); a pinned/absorbing DOF may be SUPPLEMENTAL
  only and CANNOT satisfy the control for green. Final green-confirm (xhigh, gpt-5.5) → **VERDICT: SOUND**: "§C now matches §A/§E …
  no orphaned active green-satisfying pinned/absorbing alternative remains … No remaining [PHYSICS]/[MATH]/[MECHANICAL] blocker."
  **STATUS bumped to `v3 (Codex+GLM SOUND)`** — design-review COMPLETE (Codex re-confirm green + GLM tertiary folded). NEXT: orchestrator
  final read → re-execute → tri-review.
