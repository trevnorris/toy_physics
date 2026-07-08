# II-G1b (ledger_stage012) source map — pathA_30 DtN pole ladder + Robin falsifier

> Running-start prep captured 2026-07-08 (post stage011) so a fresh session can author the reshape directive without
> re-discovery. Verify line refs against the sources before finalizing.
> Companion: `part2_gravity_atomic_split.md` (rows 011/012 + the Cluster-A reshape-cost bullet + the pathA_30 trip-ups)
> and `stage011_pathA30_frozen_reduction_source_map.md` (§1/§2 already delineate this 012 slice; §4 the 011↔012 flow).
> Build-order id **012**, Part II. **pathA_30 SPLIT: 011 (frozen-reduction certificate, DONE `ad469278`) + 012 (this =
> the DtN pole ladder + Robin falsifier).** Same two-stage pattern as pathA_29 → 009/010.
> **Source top-line: `DN_UNITTEST_BC_DEPENDENT`** — the joint pathA_30 verdict; **012 COMPLETES it** (011 carried the
> `REDUCTION_CERTIFIED` component; 012 carries the D/N ladder + the `BC_DEPENDENT` landing).

## ⭐ The headline difference from 011 (READ FIRST)
**011 was a bridge-strip PLUS a de-rig. 012 is a (nearly) PURE bridge-strip** — like 009/010. The physics here is
already genuinely computed in the source (a real D/N boundary-value problem: `dsolve` → coefficient matrix → `LUsolve`
→ DtN → pole ladder → Robin counterfactual), and the `.wl` is a genuine transfer-matrix route that **already carries
the DtN/pole/Robin/static** (unlike 011, where that content had to be EXCLUDED as 012-territory). So for 012 the `.wl`
reshape = strip ONLY the `Get[sympyExprFile]` bridge (L18–19) + the `sympy*`-comparison equalities (L50–57) + the YAML
write (L138–229), and **KEEP** the native transfer-matrix DtN/Robin/pole/static (L25–48) as the independent route.
- **No X≡X de-rig here**, but two honest-labeling watch-items: (1) `dtn_matches_target` (L454) compares the DERIVED
  `dtn` (from `LUsolve`) against the hand-typed `dtn_target=−k·tan(k·L0)` (L451) — a GENUINE comparison (dtn is
  computed), keep it; (2) `bc_derivation_emitted = False` (L542, literal) is an **honest SCOPE flag** (the mouth/cap
  `V_wall` gradient derivation is NOT emitted → the verdict stays `DN_UNITTEST_BC_DEPENDENT`, not `..._PASS`) — label
  it as a banked CALIBRATION input (like stage008's `cancellation_possible`), do NOT dress it up.
- **The able-to-fail burden is the Robin falsifier** (split-note trip-up): "keep the Robin able-to-fail guard (the
  D/D-vs-D/N swap is hardcodable without it) + the K-dim mutation probe." The Robin guard is what makes the DtN
  determination genuinely constrain — it must stay computed, not a stamp.

## File inventory
- **Report:** `software/stage1_solver/reports/pathA_30_dn_unit_test.md` — 012 content: `## DtN Derivation` :26–33,
  `## Pole Ladder` :35–39, `## Static Limit` :41–45, `## Round Trip` :47–49, `## Robin Counterfactual` :51–58,
  `## BC Provenance` :60–64; the `## Dimensional Check` `tan_argument`/`Z00_prefactor` legs :76,:79 are 012 (the
  `cs_squared` leg was 011). `## Engine Agreement` :66–68 = dropped bridge.
- **`.py`:** `software/stage1_solver/tools/pathA_30_dn_unit_test_sympy.py`. **`.wl`:** `…/pathA_30_dn_unit_test.wl`.
  Directive `…/directives/pathA_30_dn_unit_test.md`; `pathA_30_results.yaml`.

## §1 The 012 slice (`.py` line ranges) — the clean cut from 011 is at the `dsolve` (L426)
- **General solution (opens 012):** `dsolve(ode).rhs = C1·sin(ωs/cS) + C2·cos(ωs/cS)` **L426** (`ode` itself =
  stage011's `L_s`, CITED); constants `C1,C2` L427–433.
- **D/N boundary problem:** `mouth_trace` (Dirichlet at `s=0`) L435, `cap_neumann` (Neumann at `s=L0`) L436;
  `dn_matrix` L437–442, `dn_rhs=[ψM,0]` L443, `dn_det = −ω·cos(L0ω/cS)/cS` L444; `LUsolve` coeff L445–446;
  `bc_applied_solution = ψM·(sin(ωs/cS)tan(L0ω/cS)+cos(ωs/cS))` L447.
- **DtN:** `dtn_raw = −∂ₛ(bc)|₀/bc|₀` L448; `dtn_sincos` L449; `dtn` L450; `dtn_target = −k·tan(k·L0)` **L451**;
  `tan_argument = k·L0` L452; `dtn_prefactor = −k` L453; `dtn_matches_target` **L454**.
- **Pole ladder:** `denominator_full`/`factor_list` → `pole_denominator` (the `cos(kL0)` factor) L456–464;
  `pole_equation = Eq(cos(kL0),0)` L465; `pole_ladder = π·cS·(n+½)/L0` **L466**; `pole_residual` L467;
  `halfshift = (pole_residual==0)` **L468**; `denominator_zeros = solveset(cos(x))` L470–471.
- **Static limit:** `static_series = series(dtn,ω,0,6)` L472; `static_series_poly` L473;
  `static_limit = limit(dtn,ω,0⁺)` L474 (report `−L0ω²/cS² − L0³ω⁴/(3cS⁴) + O(ω⁶)`, :43).
- **Robin counterfactual (the falsifier):** `cap_robin = ∂ₛ(sol)|_{L0} − α·sol|_{L0}` L476; `robin_matrix`/`det`/
  `LUsolve` L477–487; `robin_dtn` L488; `robin_denominator_core = k·cos(kL0) − α·sin(kL0)` L489; `dn_from_robin`
  (α→0) L490; `dd_from_robin` (α→∞) L491; `dd_target = k·cot(kL0)` L492; `dd_denominator = sin(kL0)` L493;
  `recovers_dn`/`recovers_dd` L494–495; the **half-shift-destroyed-for-DD** samples L496–506 (D/D poles land on the
  INTEGER ladder `πcS·j/L0`, not the half-shifted `(j+½)`); `dd_zero_mode_removable` L507; numeric `α=2/L0` distinctness
  L509–520; `dtn_mismatch` L521; `counterfactual_guard` dict L522–529.
- **Round trip:** `r_D=−1, r_N=+1`, `round_trip = −exp(2i·k·L0)` L531–533; `round_trip_on_ladder` (sub pole ladder) →
  `R_rt = 1` L534; `round_trip_closes` L535 (report `φ0=0 mod 2π`, :49).
- **Dimensional legs (012 share):** `build_dimensional_check` `tan_argument` L233 + `Z00_prefactor`/`dtn` walk L234,
  L244–245, L259–260,L313 → `[tan_argument]=1`, `[Z00_prefactor]=L⁻¹` (report :76,:79). (The `cs_squared` leg was 011.)
- **BC provenance:** `bc_derivation_emitted=False` L542, `bc_provenance="imposed"` L543; `bc_derivation` dict L571–578.
- **012 verdict rungs** in `compute_verdict` L719–738: `dtn_matches_target`∧`halfshift`→else `FAIL_POLE_LADDER`
  (L732–733); `all(counterfactual_guard.values())`→else `FAIL_COUNTERFACTUAL` (L734–735); `bc_derivation_emitted`
  False → **`DN_UNITTEST_BC_DEPENDENT`** (L736–737); (True path → `DN_UNITTEST_PASS`, L738 — the deferred upgrade).

## §2 The 012 claim-set (derive + assert; report quotes)
- **DtN (EARNED):** solve the D/N BVP on stage011's `L_s` (Dirichlet mouth `s=0`, Neumann cap `s=L0`); the outward-mouth
  DtN is `Z00 = −(ω/cS)·tan(L0ω/cS)` (report :33) — DERIVED via the coefficient-matrix `LUsolve`, then
  `dtn_matches_target` compares it to the typed target (genuine).
- **Pole ladder (EARNED):** the DtN denominator's `cos(L0ω/cS)` factor → pole equation `cos(L0ω/cS)=0` → the
  **half-shifted** ladder `ω_n = πcS(n+½)/L0` (report :37–39); `halfshift` COMPUTED as `pole_residual==0`.
- **Static limit (EARNED):** `static_limit = −L0ω²/cS² − L0³ω⁴/(3cS⁴) + O(ω⁶)` (report :43) = the small-ω expansion of
  the dynamic DtN (NO separate static solve — report :45).
- **Round trip (EARNED):** `r_D=−1, r_N=+1`, `R_rt = −exp(2i L0ω/cS)`; on the D/N ladder `R_rt = 1`, `φ0=0 mod 2π`
  (report :49) — COMPUTED via substitution.
- **Robin counterfactual (EARNED; the FALSIFIER — do NOT drop):** the Robin cap `∂ₛψ − α·ψ = 0` gives
  `Z00^Robin = ω(αcS cos + ω sin)/(cS(αcS sin − ω cos))` (report :56); the guard asserts it **recovers D/N at α→0**,
  **recovers D/D `k cot(kL0)` at α→∞**, the D/D poles are **half-shift-DESTROYED** (land on the integer ladder), a
  numeric `α=2/L0` is **distinct** from both D/N and D/D, and `dtn_mismatch` holds (report :58). This is what proves the
  D/N determination is not hardcodable (the D/D-vs-D/N swap is reachable without it).
- **BC provenance (honest scope):** `bc_provenance=imposed`, `bc_derivation_emitted=False` → the explicit mouth/cap
  `V_wall` gradient derivation is NOT emitted, so the verdict stays **`DN_UNITTEST_BC_DEPENDENT`** rather than
  `DN_UNITTEST_PASS` (report :62–64). Label this a banked CALIBRATION input; earning it → PASS is a deferred upgrade.
- **The COMPLETED joint verdict:** print `DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, cited from
  ledger_stage011) ∧ (012: DtN ladder EARNED + bc_derivation_emitted=False → BC_DEPENDENT landing, computed here)`.

## §3 Reshape cost (the bridge to sever)
Scratch-yaml payload-mirror (Cluster-A variant) — SAME severing as 011. Strip: `write_sympy_exports` L661–670
(`_scratch/pathA_30_sympy_results.yaml` + `_sympy_exprs.wl`); `load_engine_agreement` L673–716 (MMA scratch read +
`digest_matches`/`engine_agreement`); `digest_mapping`/`yaml_read`/`yaml_write` L63–82; `mma_exports`/`expression_digest`
L580–588 (all 6 exported exprs are 012 — sever the digest, keep the exprs as native derivations); report/feed/yaml
writers L768–931. **Zero file I/O.**
- **⭐ The `.wl` is ALREADY the genuine 012 route** (native transfer matrix `transferMatrix` L25–28 → `dtnTransfer`
  L30–36, `poleDenominatorTransfer` L37, `robin*` L39–46, `ddTransfer` L47, `staticSeriesTransfer` L48). **Reshape =
  delete L18–19 `Get[sympyExprFile]` + the L50–57 `sympy*`-comparison lines + the L138–229 `Export`; KEEP the native
  transfer-matrix DtN/pole/Robin/static and print them.** Contrast with 011, where L32–48 had to be EXCLUDED — for 012
  they are the payload. The `.wl` already carries the `tan_argument`/`Z00` dim legs (L121–130) — keep them (012's).
- **Arity discipline (standing, stage007/008/009 lesson):** def/call arity scan + unevaluated-leakage transcript scan.

## §4 Consumed / exported
- **Consumes (cite, dual-site integrity, don't re-derive):** stage011's frozen **`L_s`** (the const-coeff Helmholtz
  operator), the domain **`[0,L0]`** (cap `R0(L0)=0`), and **`c_S`** (`k=ω/cS`). The `dsolve` of `L_s` opens 012 — cite
  `L_s` and `c_S` from stage011 (two independently-typed sites, `consumed − pipeline ≡ 0`, a corruption tooth firing in
  both engines). Also cites stage004's `{L,T,M}` + the `[K]=[P]−5[ρ]` chain for the `tan_argument`/`Z00` dim legs.
  (`c_S² = 5Kρ*⁴/m` is Part-I edge R1 via stage005/011 — do NOT re-derive.)
- **Exports (split-note flow):** `011/012 export the frozen throat packet + D/N provenance + Helmholtz operator → 013
  (β) + 017 (calibration input)`. 012's share = the DtN `Z00 = −(ω/cS)tan(L0ω/cS)`, the resonance ladder
  `ω_n=πcS(n+½)/L0`, and the `BC_DEPENDENT` provenance (the imposed-BC calibration slot).

## §5 Teeth candidates
Keep/assign to 012: **Robin able-to-fail guard** (the falsifier — verify all 6 `counterfactual_guard` booleans are
COMPUTED: recovers-D/N-at-α→0, recovers-D/D-at-α→∞, half-shift-destroyed-for-D/D on the integer vs half-shifted ladder,
numeric-α-distinct, dtn_mismatch; a mutation that hardcodes the D/D swap must be caught); **pole-ladder tooth** (corrupt
the ladder to the integer `πcS·j/L0` → `halfshift` / `pole_residual==0` assert fires); **dtn-matches-target tooth**
(corrupt the derived `dtn` → the `dtn_matches_target` assert fires — proves the DtN is derived-not-typed); **round-trip
tooth** (corrupt `r_D`/`r_N` → `R_rt=1` assert fires); **K-dim tan-arg probe** (the `tan_argument`/`Z00` dim legs +
corrupt-`[K]` → `[Z00_prefactor]≠L⁻¹` → `FAIL_DIMENSIONAL`; the 012 half of the shared dim block — the c_S² leg was
011); **consumed-`L_s`/`c_S` corruption** (dual-site, fires both engines). PLUS the `.wl` def/call arity scan +
unevaluated-leakage transcript scan. ⚠ The `bc_derivation_emitted=False` flag is an honest scope flag, NOT a tooth —
do not fabricate a BC-derivation to force `..._PASS`.

## §6 Register expectation
**Zero new counted knobs.** `α` (Robin cap admittance) is a **control-construction symbol** (tracked, not counted —
like `k_warp` at stage010; it builds the falsifiable counterfactual, not the physics; `[α]=L⁻¹` so `α·cS` matches
`ω`... verify the dim). Cites `L0` (011, ACTION-geometry), `c_S` (R1), `ħ,m,K,ρ*` (Part I). **New edge candidate:**
**R28** — the D/N mouth/cap boundary condition is **IMPOSED** (banked calibration; `bc_provenance=imposed`,
`bc_derivation_emitted=False`) → the `DN_UNITTEST_BC_DEPENDENT` landing; the explicit `V_wall` mouth/cap gradient
derivation is the deferred obligation whose discharge would earn `DN_UNITTEST_PASS` (a `PENDING`/`IMPOSED` edge, NOT a
reduction — analogous to R23's constraint-spec obligation). The half-shifted resonance ladder is a DERIVED consequence
of the imposed D/N pair (record, but it collapses no knob). No *reduction* edge (the operator + speed are banked
011/R1).

## Verdict tokens + honest scope
012 COMPLETES **`DN_UNITTEST_BC_DEPENDENT`**: the DtN `−(ω/cS)tan(L0ω/cS)` + the half-shifted pole ladder + `R_rt=1` are
genuinely derived (EARNED); the Robin falsifier proves the D/N determination has teeth; but the boundary condition
itself is IMPOSED (`bc_derivation_emitted=False`), so the joint verdict lands at `BC_DEPENDENT`, NOT `PASS`. Earning
`PASS` (an explicit mouth/cap `V_wall` derivation) is a deferred upgrade (banked calibration, sim-adjacent). CITED: R1
(speed, stage005/011), stage011's `L_s`/domain. Nothing here selects `L0`. No new counted knobs.

## Process (unchanged, calibrated)
Author reshape directive (§1 bridge-strip + KEEP the `.wl` transfer-matrix 012 route + §2 faithful cover + the
honest-label watch-items for `dtn_matches_target`/`bc_derivation_emitted` + §5 Robin-falsifier teeth) → Codex xhigh
design-review → fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI; 012 is repackaging) → **pre-exec USER GATE** → Codex
builds the two scripts (`--sandbox danger-full-access`, background, `< /dev/null`-equivalent stdin, absolute paths,
xhigh) → dual-engine exit 0 (repo root + foreign CWD) → arbiter re-run → tri-review (fidelity + adversarial-with-ablation,
incl. the Robin-falsifier genuineness + arity scan + tally spot-check) → remediate → fresh-agent re-verify →
registration 11→12 + parameter register (edge R28) + Codex-verify → note/card/`\input{stages/stage_012}` → PDF →
commit + docs/memory sync. Target stem: `ledger_stage012_dtn_pole_ladder_robin` (confirm slug at authoring).
