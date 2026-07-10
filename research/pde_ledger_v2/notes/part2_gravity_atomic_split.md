# Part II — Gravity: FINALIZED atomic-stage split (user-approved 2026-07-08)

> Per-Part user gate (blueprint §7) satisfied 2026-07-08. This is the frozen split for Part II. Assembled from a
> three-cluster source-mapping fan-out (Gates 1–5 / force–return–radiation / port–PN–Gate-6), 2026-07-08.
> Per-stage pre-exec user gate stays in force (unchanged calibrated process).
> Governing: `notes/ledger_v2_blueprint.md` §2 (granularity), §3 (Part II table), §5 (reshape spec), §6 (verification).
> Build order continues from Part I: Part-II new stages = build-order ids **008–029** (22 new stages; stages 001/002
> already built in the pilot belong to Part II, so Part II totals **24** — inside the blueprint's 20–28 estimate).

## Already built (pilot) — no double-count

- **001** solid-angle & second-moment primitives + **002** matter-stress force assembly (`FORCE_ATTRACTIVE_DERIVED`).
  ⭐ Fan-out verdict: **pathA_21c is FULLY folded** — both the reduced-3D `r⁻²` and bulk-4D `R⁻³` lanes, the attractive
  sign chain, and all named residuals live in 001/002. Do NOT re-add a 4D or residual stage. (Only stated-not-derived
  fragment: the Lagrangian→Noether-flux reduction 002 §2 cites — a possible in-place refinement, not a stage.)

## The proposed 22 new stages (dependency order)

### Cluster B — return & radiation constraints (3 stages; sources pathA_28, pathA_29)

| id | Stage | Source | Headline token | Content (atomic step) | Notes |
|---|---|---|---|---|---|
| 008 | II-B1 monopole/dipole constraint spec | pathA_28 | `MONOPOLE_DIPOLE_RETURN_CONDITIONAL` | raw DtN outgoing amplitudes ℓ=0/1/2 (orders p=1/3/5) + the cancellation targets `R0=−M0`, `R1_i=−D1_i` | self-contained constraint-SPEC (not a falsifiable suppression test — honest scope carried); exports the targets 34 consumes + Q₂ as FREE anchor |
| 009 | II-B2 flat-slab return residual | pathA_29 | `RETURN_RESIDUAL_PREDICTION` (part 1) | transport phase → `T0(0)=1/(ε0+1)` (so `1−T0(0)=ε0/(1+ε0)`), residual orders p_res=1,3, signed `Z=−M0·ε0/(ε0+1)` | the falsifiable residual-radiation prediction; `Z<0` labeled PREMISE (v3 relabel history); (fan-out slip `ε0=1−T0(0)` corrected 2026-07-08, fidelity leg catch) |
| 010 | II-B3 localization p=2 + NOGO control | pathA_29 | `RETURN_RESIDUAL_PREDICTION` (part 2) | both DC-sink completions → normalizable zero mode → real 3D-radial dsolve → p=2 (1/r² survives the slab); the anti-localizing warp → p=3 `RETURN_NOGO` folded in as the able-to-fail companion | 1/r² survival; NOGO control stays IN-stage (not split) |

### Cluster A — the reduced-closure gate ladder (13 stages; sources pathA_30–34)

| id | Stage | Source | Headline token | Content (atomic step) |
|---|---|---|---|---|
| 011 | II-G1a frozen-reduction certificate | pathA_30 | `DN_UNITTEST_BC_DEPENDENT` (1/2) | frozen wall → Helmholtz `ψ''+(ω/c_S)²ψ=0`, `c_S²=5Kρ*⁴/m`, projection measure + certificate |
| 012 | II-G1b DtN pole ladder + Robin falsifier | pathA_30 | `DN_UNITTEST_BC_DEPENDENT` (2/2) | dsolve → D/N determinant → DtN `−(ω/c_S)tan(L0ω/c_S)`, half-shifted ladder, R_rt=1, Robin counterfactual, K-dim probe; BC = banked CALIBRATION input |
| 013 | II-G2a harmonic profiles + M/K projection | pathA_31 | `BREATHING_CALIBRATED` (1/3) | α_a/α_L harmonic lifts → `M_AB`,`K_AB` by real ∫dw operator projection → the (a,L) collective closure |
| 014 | II-G2b truncation consistency | pathA_31 | `BREATHING_CALIBRATED` (2/3) | combined-basis generalized eig, β_L0 sweep + N-convergence overlaps → validity window (order-unity wall stiffness caveat carried) |
| 015 | II-G2c legacy-structure recovery + HF force | pathA_31 | `BREATHING_CALIBRATED` (3/3) | legacy-Hessian pattern match + Hellmann–Feynman force, BOTH routes independently emitted (the v1 x−x rig's fix preserved) |
| 016 | II-G3a ℓ=2 SO(3) covariance theorem | pathA_32 | `ISOTROPY_CALIBRATED` (1/2) — EARNED slice | real ℓ=2 harmonics, angular Gram=I₅, computed `−Δ_S²` eigenvalues λ_m=6, K₂ angular stiffness |
| 017 | II-G3b grouped-P2 lane isotropy | pathA_32 | `ISOTROPY_CALIBRATED` (2/2) | grouped {20,21,22} lanes, raw-D defects=0 (PRIMARY), normalized u-defects cross-check, calibration partition; exports the ℓ=2 port kernel |
| 018 | II-G4a DtN Hankel fingerprint | pathA_33 | `QUAD_CALIBRATED` (1/4) — EARNED slice | outgoing series u2=a²/9c_s², u4=4a⁴/81c_s⁴, v5=a⁵/27c_s⁵; χ_Q=+1 outgoing vs −1 incoming |
| 019 | II-G4b prefactor algebra | pathA_33 | `QUAD_CALIBRATED` (2/4) | `P(ω)=D₀N(ω)/D^cons(ω)²`, squared-denominator P2=(D0N2−2D2N0)/D0², P4 |
| 020 | II-G4c 54/5 provenance partition | pathA_33 | `QUAD_CALIBRATED` (3/4) | `54/5=2·27/5` — 27 COMPUTED (`derived_in_gate`), 2/5·G labeled `external_bridge_input`; `G=GENUINE_BLOCKED`; the 4-way provenance partition + its able-to-fail |
| 021 | II-G4d μ̂₀-free dimensional closure | pathA_33 | `QUAD_CALIBRATED` (4/4) | `[P0_phys]=1` μ̂₀-FREE (the natural-units trap catch), drop-normalization + corrupt-`[N₀]` probes; μ̂₀ diagnostic stays non-gate by design |
| 022 | II-G5a cross-ℓ fingerprints + Gate-4 non-regression | pathA_34 | `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (1/2) | outgoing `−(ℓ+1)/Λ_ℓ` for ℓ=0,1,2 (ω¹/ω³/ω⁵) + quadrupole non-regression |
| 023 | II-G5b nullspace underdetermination (departure) | pathA_34 | `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (2/2) | genuine native nullspace (dim 8, return nullity 2), residuals vs pathA_29, `Z0_ret/Z1_ret` untouched → the Gate-6 selector need; selector-equation control → `CROSS_L_RESIDUAL_PREDICTION` (able-to-fail) |

### Cluster C — the density port + consistency caps (6 stages; sources pathA_43, matchback, PN corpus)

| id | Stage | Source | Headline token | Content (atomic step) |
|---|---|---|---|---|
| 024 | II-P1 density-port derivation | pathA_43 | `DENSITY_PORT_HOSTED` (1/4) | N0_den two-port over (q2 wall, Φ2 bulk) — Schur elim / Green-DtN routes; consumes 29's bulk mode + 32's wall mode |
| 025 | II-P2 vector-freedom taint | pathA_43 | `DENSITY_PORT_HOSTED` (2/4) | COMPUTED ancestry/taint over free_symbols + host-set/source_map_complete guard + expression ablation (the anti-rig mechanisms preserved verbatim) |
| 026 | II-P3 continuity lineage | pathA_43 | `DENSITY_PORT_HOSTED` (3/4) | ℓ0(M0)→ℓ1(D1)→ℓ2(`∫Y2* S_leak`) moment ancestry via 29's operator; token-check computed, never a flag |
| 027 | II-P4 port checks + closure slot | pathA_43 | `DENSITY_PORT_HOSTED` (4/4) | the 6 able-to-fail checks (`[N0]=L⁻¹M`, a⁻⁵, `+i z⁵/27`, χ_Q=1, …) + `P0_phys=(c_s/a)²N0/D0` + the K̄ closure slot (the A3 boundary — marked shared, not double-counted) |
| 028 | II-P5 2.5PN match-back | pathA_2_5pn_matchback | (consistency over CALIBRATED moments) | INV1–INV5 (both Γ̄₅ forms → Burke–Thorne, cross-form, K̄₄=4K̄₂²/K̄₀, the INV5 literal anchors) + the 11-mutation matrix; `G=GENUINE_BLOCKED` carried |
| 029 | II-P6 PN corpus DOI-cite | `research/4d_*pn*` + `1pn_orbital_dynamics` | (audited external corpus) | THIN CITE-only stage note: paper/DOI/imported-result/scope per corpus entry (Zenodo DOIs confirmed in the .tex bibliographies); NO scripts, NO re-derivation |

## Settled recommendations (for the gate)

1. **Gate-6 = a Part-VII register item, NOT a Part-II stage** (fan-out recommendation, matches blueprint §8): it
   produces no verdict and no algebra — it is the canonical open-item record (branch realization
   `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵`, the `Z0_ret/Z1_ret` selector need, `L/a` self-selection). Part II carries it only as the
   sim-deferred caveat line on 023/027/028.
2. **PN-cite (029) as a thin Part-II stage** rather than a register item — it is a derivation-chain dependency the
   imports in 027/028 lean on, so it needs visible in-Part provenance.
3. **NOGO warp control folded into 010** (not split) — it is pathA_29's mandatory able-to-fail companion.
4. Totals: **22 new + 2 built = 24 Part-II stages** (blueprint est. 20–28 ✓).

## Reshape-cost map (what §5 work each cluster needs)

- **pathA_30–34 (13 stages):** NO argparse, but a **scratch-yaml payload-mirror variant**: `.py` writes scratch
  yaml/`_sympy_exprs.wl` exports, `.wl` READS the sympy expr export and cross-checks, `.py` re-reads the MMA scratch
  into `engine_agreement`/`digest_matches`. Reshape = strip the whole scratch-file bridge + results-yaml writes;
  KEEP each `.wl`'s already-genuine independent route (transfer-matrix / native DSolve-BVP / Eigensystem) as a truly
  standalone engine (it must derive its own inputs instead of Importing the `.py`'s expressions).
- **pathA_28/29 (3 stages):** JSON digest bridges (`.py` emits json, `.wl` Imports + asserts digest) AND
  **cross-script runtime YAML reads** (29 reads 28's results yamls; both read their own report yamls). Reshape =
  sever all of it; inline the consumed 28 amplitudes/moments as cited symbolic inputs (stage002 consuming pattern —
  008 is the formal home 009/010 cite).
- **pathA_43 + matchback (5 stages):** already contract-clean (no argparse/JSON/reads) — lightest reshapes; the work
  is decomposition + preserving the anti-rig mechanisms.

## Per-gate trip-ups the directives must pin (from the rig histories — do not reintroduce)

- **pathA_31 (015):** v1 was REJECTED for an HF `x−x` tautology + typed counterfactual flags + gamed threshold —
  both HF routes independently emitted; flags/threshold COMPUTED.
- **pathA_32 (016/017):** keep the aggregate probe battery intact (neuter-one → aggregate flips); eigenvalues
  computed (the `tautology_hash` probe); dim-mutation on sourced `T_Ω` must fire `FAIL_DIMENSIONAL`.
- **pathA_33 (020/021):** the 27 stays COMPUTED; `54/5` asserted only as labeled `external_bridge_input`; the v1
  rig (back-solved free-carrier μ̂₀ + constant self_ablation) must not return — μ̂₀-free dim gate + real two-verdict
  self-ablations.
- **pathA_34 (023):** DEFAULT verdict is PASS — `UNDERDETERMINED` must be EARNED from a genuinely-computed native
  nullspace (the v1 rigged-to-UNDERDETERMINED zero-padded-constraint history); back-solving `ε_eff`/`Z` from
  residuals is forbidden (`FAIL_TAUTOLOGICAL` firewall preserved).
- **pathA_43 (025/026):** the computed-taint + host-set guard + lineage token-check are THE earned mechanisms
  (two caught rigs live in their history) — never collapse to name-checks or flags.
- **pathA_2_5pn_matchback (028):** the INV5 literal anchors are the anti-rig teeth vs the coherent-rescale mutation —
  never dedupe them against the closure values.
- **pathA_30 (012):** keep the Robin able-to-fail guard (D/D-vs-D/N swap is hardcodable without it) + the K-dim
  mutation probe.
- **Wolfram arity-mismatch lesson (stage007):** every `.wl` review now includes a def/call arity scan + an
  unevaluated-leakage transcript scan (a mismatched call silently skips its section at exit 0).

## Cross-stage flows (consumed/exported, for the consuming-stage discipline)

001/002 (force) ← stage001 primitives. 008 exports `R0=−M0`/`R1=−D1` + raw amplitudes + Q₂-as-anchor → 009/010/023/026.
011/012 export the frozen throat packet + D/N provenance + Helmholtz operator → 013 (β) + 017 (calibration input).
013 exports the (a,L) closure + M/K → 022/023 (ℓ=0 map). 017 exports the grouped-P2 ℓ=2 port kernel + `K_η+2T_Ω` +
support scalars → 018–021 + 022/023 (ℓ=2 map) + 024 (wall mode). 018–021 export the Λ₂ fingerprint + χ_Q=1 + the
54/5 partition → 022 (non-regression) + 027 (closure). 009/010 export the bulk Helmholtz mode + projected-continuity
operator → 024/026. 024–027 export N0_den + K̄ moments → 028. 029 backs the imports in 027/028.

## Progress

- **II-B1 `ledger_stage008` DONE (2026-07-08)** — pathA_28 constraint-spec reshape. `MONOPOLE_DIPOLE_RETURN_CONDITIONAL`
  with the scope caveat carried verbatim-class (constraint-spec not suppression test; `x−x` bookkeeping labeled;
  `cancellation_possible` = literal flag printed as such). EARNED: the DtN ladder p=1/3/5 SCANNED (not typed) +
  kernels + steady limit + dominance + the 243/244 anchors with computed strict-recovery observations. JSON/YAML
  engine bridge severed both directions; ZERO file I/O. Dual-engine SymPy 54 / Mathematica 57 (3 = arity
  self-checks), CWD-independent. Tri-review CLEAN (`FIDELITY_CLEAN` hand-verified math; `ADVERSARIAL_CLEAN` 14/14
  mutation matrix incl. dual-corruption + planted-arity-mismatch detection) → 5 nits remediated (8 vacuous tally
  stamps de-counted; recovery observations restored computed; tooth 7 pipeline-routed; `real=True` alignment; dead
  param wired) → `REVERIFY_CLEAN` (tallies independently re-counted 54/57 = 60−8+2 / 63−8+2). Register: zero new
  knobs + edge R23 (`PENDING` obligation), `REGISTER_CLEAN` first pass. Registration at count 8.
- **II-B2 `ledger_stage009` DONE (2026-07-08)** — pathA_29 Check-A reshape (`RETURN_RESIDUAL_PREDICTION`, Check-A
  component + explicit stage010 qualifier; no faked standalone headline). EARNED: τ=2d/c_S SOLVED; DC fractions +
  steady balances; ν_ℓ=0 scanned → p_res=1/3 COMPUTED; Z accounting + sign certificate (premise-vs-accounting v3
  labels exact); per-channel strict limits (strict_ν1 computed from the ℓ=1 limit DIRECTLY — tightened over the
  source's ν0 reuse). Heaviest bridge-severing: JSON-digest + pathA_28 YAML reads + trace hashes ALL dead; zero
  file I/O; stage008 kernels consumed via DUAL-SITE citation-integrity. Dual-engine SymPy 48 / Mathematica 52.
  ⚠ Tri-review CAUGHT A RIG-CLASS DEFECT: the v1 kernel citation-integrity was set-then-compare-to-self (single-
  source 27→11 corruption exited 0; adversarial 16-mutant matrix) → remediated (dual-site; all four one-site
  corruptions now fail both engines) + 3 nits → fresh-agent `REVERIFY_CLEAN` (killing mutants re-run
  independently). `Off[Limit::alimv]` unmasked + shown benign. Register: rows {d, ε0, ε1} + edge R24,
  `REGISTER_CLEAN` first pass. Registration at count 9. (Fidelity also caught + fixed a fan-out slip in THIS
  note's row 009: `1−T0(0)=ε0/(1+ε0)`.)
- **II-B3 `ledger_stage010` DONE (2026-07-08)** — pathA_29 Check B reshape (`RETURN_RESIDUAL_PREDICTION` COMPLETED:
  Check B `p=2` computed + Check A cited from stage009). EARNED: two normalizable `m=0` transverse zero modes
  (compact-cell + Bloch) → genuine 3D-radial `dsolve` → **`p=2`** BOTH DC-sink completions; static–dynamic
  consistency (exponent + the stronger Green-equality `1/(4πdr)`, design-review-confirmed); gapped Yukawa contrast;
  counterfactual `r⁻⁴` → residual `5/(π·d·r⁷)` rejected; **NOGO warp control → `p=3` → `RETURN_NOGO`** (able-to-fail);
  computed classifier + falloff-tension witness. CITED dual-site-from-the-start: stage008 `p_raw(ℓ2)=5` →
  `quadrupole_survives` (+ exact-value anchor) & stage009 `A_residual_pass=True`; no ℓ=2 recompute (`T2_applied=false`).
  Same bridge-severing as 009 (zero file I/O; SHA/`structure_id`/`expr_digest` retired). Dual-engine SymPy 71 / WL 81
  (= 71 + 10 arity), CWD-independent. Tri-review CLEAN both legs: `FIDELITY_CLEAN` (all values hand-re-derived; 009/010
  boundary clean) + `ADVERSARIAL_CLEAN` (22-mutant matrix; a hardcoded-wrong-Green mutant proved the counterfactual
  guard genuine; classifier genuinely reaches NOGO). Adversarial's one substantive find — a *coordinated both-sites*
  stage008 citation corruption escaped — folded to stage009 parity via an exact-value anchor (4 nits total: + de-count
  4 vacuous/set-then-compare dim stamps + m=0 seed restatement + a cosmetic `rewrite(sp.exp)` no-op; tallies 76/86 →
  71/81) → fresh-agent `REVERIFY_CLEAN` (anchor proven the sole gate closing the escape; four one-site corruptions
  still fire). Register: **zero new counted knobs** (`k_warp` = control-construction `CANDIDATE`, tracked-not-counted)
  + edge R25 (localization `p=2` EARNED-within-family, structural — discharges NOTHING; R19/`W_slab` & R23 stay
  `PENDING`). Registration at count 10. ⚠ Process note: the remediation Codex was killed mid-final-check but the edits
  had landed + passed at repo root; orchestrator cleaned 2 orphaned `.tmp` files + completed the foreign-CWD
  verification. **⭐ pathA_29 fold COMPLETE (stages 009+010).**
- **II-G1a `ledger_stage011` DONE (2026-07-08)** — pathA_30 frozen-reduction certificate reshape **+ a DE-RIG** (the
  first de-rig since stage007). `FROZEN_REDUCTION_HELMHOLTZ_CERTIFIED` (011 component of the joint
  `DN_UNITTEST_BC_DEPENDENT`; the D/N ladder + `BC_DEPENDENT` landing = stage012). EARNED: `L_s = ψ''+M·ψ'+N·ψ−B·(ħ²/4m²c_S²)ψ''''`
  **assembled** from the reduction (projection measure + every intruding coeff computed with its vanishing/deferral
  condition) → const-coeff Helmholtz `ψ''+(ω/c_S)²ψ=0` on `[0,L0]`; the three source tautologies genuinely
  **de-rigged** (`operator_is_helmholtz` PRODUCED from the assembly, `speed_is_cs` EXTRACTED from the ψ-coeff,
  `domain_is_L0` SOLVED from the pinch-off `R0(L0)=0`; `unsuppressed_operator_intrusion` COMPUTED) — each able-to-fail;
  the BdG `k⁴` term a genuine **4th-order** intrusion (design-review CAUGHT my first draft's coefficient-shift error);
  validity window `{ρ0'/ρ0=0, √γ0 const, δV_conf=0, ∇Q=0, kξ≪1}`; c_S² dim leg `(2,0,-2)` + corrupt-`[K]` probe; ξ≠ℓ_c
  firewall. CITED R1 (`c_S²=5Kρ⁴/m`, stage005) at ρ* dual-site + frozen-export anchor; clean 011/012 cut at the
  `dsolve` (no DtN/pole/Robin here). Zero file I/O (scratch-YAML/`Get`/`sympy*`/digest bridge severed); the `.wl`
  independent (native `D`/`Solve`/`Coefficient`, own trig-basis cross-check, arity self-check). Dual-engine SymPy 61 /
  Mathematica 71 (= 61 + 8 arity + 2 trig-basis), CWD-independent. Codex directive design-review caught 3 BLOCKING
  (BdG operator treatment; dimensionally-safe assembly; R1 dual-site example) → folded → confirm-pass `DIRECTIVE_CLEAN`.
  Tri-review CLEAN both legs (`FIDELITY_CLEAN` all values hand-re-derived + de-rig confirmed produced; `ADVERSARIAL_CLEAN`
  14+ mutant matrix — de-rig proven able-to-fail, BdG decisively a 4th-order term, dual-site R1 anchor catches the
  coordinated-drift escape, no 012 leakage) → 3 teeth/label nits remediated (firewall tooth made genuine — was a vacuous
  `xi≠xi`; witness-path tooth 3b added → the intrusion flag now independently load-bearing, closing the
  hardwire-flag-to-False escape; arity labels corrected 6/8→5/7) → fresh-agent `REVERIFY_CLEAN` (both escapes closed).
  Register: **zero new counted knobs** (`L0`=`ACTION`-geometry like stage009's `d`; `ℓ_c` INERT; `R_mouth` cancels; `ξ`
  R2-family DERIVED) + structural edges **R26** (frozen-reduction validity record) + **R27** (ξ≠ℓ_c firewall) —
  discharge NOTHING; `REGISTER_CLEAN` first pass. Registration at count 11.
- **II-G1b `ledger_stage012` DONE (2026-07-08)** — pathA_30 DtN pole ladder + Robin falsifier reshape (a near-PURE
  bridge-strip, NOT a de-rig — contrast 011). **⭐ COMPLETES the joint `DN_UNITTEST_BC_DEPENDENT`** = (011:
  `REDUCTION_CERTIFIED`, cited) ∧ (012: DtN ladder EARNED + `BC_DEPENDENT` landing, computed here). EARNED: `dsolve` of
  the cited frozen `L_s` → D/N coeff-matrix `LUsolve` → DtN `Z00=−(ω/c_S)tan(L0ω/c_S)` (`dtn_matches_target` = genuine
  derived-vs-typed, NOT X≡X); half-shifted ladder `ω_n=πc_S(n+½)/L0` (`halfshift=pole_residual==0` computed); static
  small-ω series (distinct from `limit=0`); round-trip `R_rt=1`; the **Robin falsifier** — the true 6-member
  `counterfactual_guard` `{robin_determinant_emitted (hardened to the computed det, not bool(hstr)), recovers_DN_at_alpha0,
  recovers_DD_at_alpha_inf, halfshift_destroyed_for_DD, numeric_alpha_distinct, dtn_mismatch}`, each a computed residual
  (`dd_zero_mode_removable`=1/L0 held out as artifact); `tan_argument`/`Z00` dim legs + corrupt-`[K]` (propagates
  through `[c_S²]→[c_S]→[k]`). CONSUMES stage011's `L_s`/`[0,L0]`/`c_S` dual-site — the `L_s` site B = **null-space
  reconstruction** from the `dsolve` `{sin,cos}` pair (posit `y''+ay'+by`, solve `(a,b)=(0,(ω/c_S)²)`; a genuinely
  independent construction, NOT a `k²`-rename — sinh/cosh corruption fires both engines). Zero file I/O; the `.wl` KEEPS
  its transfer-matrix route as the independent engine (native self-checks `robinAlpha0==dtnTransfer`/`robinAlphaInf==ddTransfer`
  preserved; DtN cross-check vs the `.py` `LUsolve` computed-not-hardcoded). Dual-engine SymPy 84 / Mathematica 90
  (+6 = arity block net of 2 SymPy-only DtN-detail checks), CWD-independent. Codex directive design-review caught 3
  BLOCKING (`L_s` dual-site was a rename → null-space route; wrong guard membership + `robin_determinant_emitted`
  hardening; §4.1 tooth was a no-op forced-true → break-the-derivation) + 1 nit (§2d series≠limit) → folded → confirm
  `DIRECTIVE_CLEAN`. Tri-review CLEAN both legs (`FIDELITY_CLEAN` full coverage diff; `ADVERSARIAL_CLEAN` 16-mutant
  matrix — guard un-stampable SM12, sinh/cosh fires SM2/MM3, hardcoded-target caught SM7, all 4 rungs reachable, arity
  clean MM2). One documented non-blocking nit (the `.wl`'s `robin_determinant_emitted` witnesses the typed
  `robinDenominatorCore` — directive-permitted, transfer route has no coeff-matrix det; derived route independently
  guarded) — no remediation (no rig/vacuous-tooth/escape). Register: **zero new counted knobs** (`α`=Robin cap
  admittance, `[α]=L⁻¹`, control-construction `CANDIDATE`, tracked-not-counted like `k_warp`) + edge **R28** (D/N
  boundary IMPOSED → `BC_DEPENDENT`, `IMPOSED`/`PENDING`, discharges NOTHING; deferred discharge = the mouth/cap `V_wall`
  derivation earning `DN_UNITTEST_PASS`), `REGISTER_CLEAN` first pass. Registration at count 12, PDF rebuilt (39pp).
  ⭐ **pathA_30 fold COMPLETE (stages 011+012).**
- **✅ II-G2a `ledger_stage013` DONE (2026-07-08) — the FIRST Part-II CALIBRATED stage.** pathA_31 harmonic-profile +
  M/K-projection reshape (`BREATHING_CALIBRATED` component 1/3): DERIVED harmonic-lift profiles `α_a=sinh(L0β−βw)/sinh(L0β)`,
  `α_L=rAL sinh(βw)/sinh(L0β)` (proven by the `𝓛₀[α]=0` residual); `M_AB`/`K_AB` by real ∫dw operator projection
  (`K_aa=4πT_wβ/tanh(L0β)`, …, `det K=16π²T_w²β²rAL²`) with `forbidden_fit_flags` computed False via **free-symbol-name**
  ancestry (`operator_projection_not_static_Hessian` — NOT the legacy static Hessian); the dynamical-EOM **LHS** (RHS
  `F_A^(HF)` deferred to 015). Clean 3-way cut (L627–679; no 014 truncation, no 015 HF force / legacy Hessian); `c_S` NOT
  consumed (matter-sector deferred, `kξ≪1`). Consumed the frozen wall packet `{L0=37/20,T_w=1,K_η=1,β=1}` via a dual-site
  guard with an **anti-tautology** proviso (packet values held as independent `K_eta_cited` datums, NOT re-expanded via
  the `K_η=T_w β²` alias — else site A collapses to `β−β≡0`). ⭐ **Register: the FIRST Part-II counted knobs — 3 CALIB
  `{μ_η, T_w, β}`** (frozen-wall constitutive packet), with `K_η=T_w β²` a DERIVED manifestation (R29) + R30 the
  named-PENDING nonlinear-throat reduction; ⚠ `β` counted NOT dressed as geometry (source: "geometry alone does not
  derive it"); `r_AL` control-ratio tracked-not-counted. Codex `REGISTER_CLEAN`. **⭐ Directive review used the new
  Codex→Grok→Codex bookend:** Codex `DIRECTIVE_CLEAN` → Grok-4.5 compute-verify pass CAUGHT a kernel-preserving
  residual-tooth defect (my §4.1 example mutants stayed in `ker(𝓛₀)`; verified in SymPy) + 5 hardening nits, all folded →
  Codex confirm-pass (caught 1 dangling `M_posdef` label). Dual-engine SymPy 78 / Mathematica 84, CWD-independent;
  tri-review CLEAN both legs (`FIDELITY_CLEAN` re-derived all M/K to ~1e-124 via a different path; `ADVERSARIAL_CLEAN`
  10-ablation matrix all firing, incl. the AB2b non-kernel-residual proof + the AB3b anti-tautology proof) → 1 non-blocking
  nit remediated (baseline re-integration was value-wise `X−X` → retargeted to the hardcoded closed-forms) → fresh-agent
  `REVERIFY_CLEAN`. Register at count 13, PDF 41pp. **⚠ Also fixed a manifest lag** (the `LINEAR_STAGE_RENUMBERING_MANIFEST.json`
  policy array had stalled at 010 — 011/012 were bumped in coverage but never added to the manifest; added 011/012/013).
- **⭐⭐ TEETH-HARDENING PASS DONE (2026-07-08, post-013, user-directed).** A Grok per-tooth-ablation backstop of the
  committed gravity stages 008–012 found ~11 vacuous/weak able-to-fail constructs the original tri-reviews missed (X≡X
  set-then-compare, dim-tautologies, subsumed/miswired guards) — 009 CLEAN, 013 already hardened. **Falsification-STRENGTH
  defects, NOT earned-physics errors** (no result rested solely on a vacuous tooth; `.wl`s all genuine independent
  engines). Codex remediated (make-genuine or honest de-count); fresh-Grok re-verify: 008/010 CLEAN, 011/012 SymPy CLEAN.
  New tallies 008 53/56 · 010 71/81 · 011 60/70 · 012 82/90. Re-verify also caught 2 WL-only constructs (011 trig, 012
  pole-denom) → Codex WL follow-up fixed; and 1 Grok FALSE-POSITIVE (012 `halfshift_destroyed_for_DD` is genuine — always
  verify Grok negatives by reading code). ⭐ Per-tooth ablation now MANDATORY in the tri-review (blueprint §6 +
  [[feedback-per-tooth-ablation]]). ⚠ OPEN FOLLOW-UP: the `.wl`s were only STATIC-read (2-seat cap) — a dedicated
  WL-ablation pass over all built stages is a tracked Phase-C-adjacent follow-up.
- **✅ II-G2b `ledger_stage014` DONE (2026-07-08) — pathA_31 `BREATHING_CALIBRATED` 2/3 (truncation consistency); the
  FIRST NUMERIC / float-bearing stage.** The combined-basis generalized eigenproblem `K v=ω²M v` over `{α_a,α_L,g_1..g_N}`
  certifies 013's 2-mode closure captures the two lowest generalized modes to the overlap floor `o_k≥0.9` (`o_1=0.99311,
  o_2=0.98776` at `β_L0=37/20`), across a COMPUTED validity window `β_L0∈[0.1,3.0]` (⟺ `K_η/T_w≲2.6`; genuine FAIL rows at
  `β_L0≥5`, sweep spans to 50) and N-converged over `N=4/8/12/16`. ⭐ **The pathA_31 v1 gamed-threshold scar is ABSENT** —
  floor + window genuinely computed + able-to-fail. ⚠ **Honest caveat carried:** the modal overlap does NOT guard
  profile-correctness — `constant_one` (wrong profile) PASSES the overlap (`o_1=1.0,o_2=0.974`); the profile guard is 013's
  residual + 015's HF. Consumed 013's profiles + frozen packet via dual-site integrity (site A `β·L0=37/20` anchor; site B
  residual∧BC — the residual alone misses kernel-preserving corruptions); the M/K seam guarded via the consumed
  profiles/packet, NOT a naïve `numeric M_aa==symbolic M_aa` (normalization differs by `4π L0 μ_η`, cancels in eig/overlaps);
  `c_S` NOT consumed. Bridge severed (numeric-WL scratch export + the `.wl`'s `sympy*`-float diff); the `.wl` KEEPS its
  native `NIntegrate`+`Eigensystem` route, re-targeted to its OWN floor/window; transcript-level dual-engine agreement.
  Dual-engine SymPy/SciPy 93 / Mathematica 100, CWD-independent. **⭐ Directive review = Codex→Grok→Codex bookend:** Codex
  `DIRECTIVE_CLEAN` after 3 BLOCKING + 2 nits folded (site-B residual∧BC, the `min(ω²)>0` tooth, the N-convergence
  non-converging ablation, the identity-sub-Gram projection tooth, the native cited-closed-form `.wl`), Grok-4.5
  compute-verify `DIRECTIVE_CLEAN` (all anchors independently re-derived), Codex final confirm `DIRECTIVE_CLEAN`. Tri-review:
  `FIDELITY_CLEAN` (independent from-scratch Galerkin re-derived every value) + `ADVERSARIAL_ISSUES` — anti-gaming
  floor/window CLEAN (23/23 genuine ablations) but **2 BLOCKING vacuous able-to-fail teeth (4 & 5, `x:=const; expect_fail(x)`
  in both engines)**; remediated to a shared predicate on a mutated copy → fresh-agent `REVERIFY_CLEAN` (deciding
  coupling meta-test: corrupt the predicate → the tooth stops firing → audit fails); a follow-up deleted 3 dead `.wl`
  Module locals. Register: **zero new counted knobs** (numeric controls `{FLOOR,N_FINAL,N_CONVERGENCE,BETA_L0_SWEEP}`
  tracked-not-counted; wall packet consumed from 013, no double-count) + structural edge R31 (truncation validity window,
  discharges nothing), Codex `REGISTER_CLEAN`. Registration at count 14, PDF rebuilt (43pp).
- **✅ II-G2c `ledger_stage015` DONE (2026-07-08) — pathA_31 `BREATHING_CALIBRATED` 3/3 (legacy-Hessian structure recovery +
  the Hellmann–Feynman force); COMPLETES the joint `BREATHING_CALIBRATED` = 013 ∧ 014 ∧ 015; the SECOND calibration-adding
  Part-II stage.** SYMBOLIC / float-free (the inverse of numeric 014). EARNED: (Q1) the legacy-Hessian structure recovery —
  `H_legacy=hessian(E_geom)` own-built; the CITED 013 `M_AB`/`K_AB` match its structural signature (M pos-def by
  **exact-identity Sylvester certificates**, K symmetric, `K_aL<0`, rank + zero-pattern) → the `(a,L)` closure RECOVERED not
  re-postulated (`full_matrix_fit=False`); (Q2) the **Hellmann–Feynman force `F_A^(HF)`** (the RHS 013 deferred) by TWO
  genuinely-different routes — distributed projection vs Hellmann–Feynman parametric derivative — that AGREE
  (`hf_force_reduces=True`) with `unsimplified_routes_identical=False` (⭐ the anti-`x−x`, the pathA_31 v1 rig's OTHER locus,
  ABSENT); (Q3) the static-dynamic limit `Q̈→0 ⇒ K_AB Q = F_A`. ⭐ **The 014↔015 boundary:** the constant profile PASSES 014's
  overlap but FAILS 015's HF (`hf_mismatch=True`) — the HF is the profile guard the overlap could not supply. Consumed 013's
  profiles + M/K + packet (dual-site: site A `β·L0=37/20`, site B residual∧BC, + the M/K det-identities **plus off-diagonal
  sign checks `M_aL>0`/`K_aL<0`** — the det is blind to the off-diagonal sign, a Codex catch) + 014's truncation cert; `c_S`
  NOT consumed. Bridge severed (the `sympy*` export + `.wl` `checks` cross-read); the `.wl` KEEPS its native HF `Integrate` +
  Hessian `D` + structure `Det`/`MatrixRank` route (cites profiles + M/K as literals, NOT `DSolveValue`/re-`Integrate`).
  Dual-engine SymPy 95 / Mathematica 102, CWD-independent. **⭐ Directive review = Codex→Grok→Codex bookend:** Codex 1 BLOCKING
  (the det-identity was blind to the `M_aL→−M_aL` sign flip — det carries `M_aL²`; added explicit sign checks) → confirm
  `DIRECTIVE_CLEAN`; **Grok-4.5 compute-verify caught 1 BLOCKING** (`M_aa>0` is NOT float-free dischargeable — `sp.ask` returns
  `None`, `bool(>0)` raises — so it would silently regress to the vacuous `M_aa−m_aa_positive_form==0` form-equality;
  required exact-identity positivity certificates + a form-equality prohibition) + 1 nit → Codex final confirm
  `DIRECTIVE_CLEAN`. Tri-review: `FIDELITY_CLEAN` (hand-verified all 7 certificate identities + every value) +
  `ADVERSARIAL_ISSUES` — the three central burdens proven genuine by live copy-mutations, but **1 BLOCKING vacuous `x−x`
  tooth (tooth 4, `expect_fail(compact(F_legacy_a−F_legacy_a)!=0)`, a self-compare) in both engines** + 3 lower-severity
  literal/scaffolding flags → remediated (make-genuine ×4: tooth 4 reads the real wrong-profile forces, static-dynamic → a
  genuine `Q̈→0` EOM-residual identity, `full_matrix_fit` → `M_aa≠H_aa`, `structure_from_computed` → firewall-tied; de-count ×2:
  the typed-twice scaffolding + the bypassing-probe meta-asserts) → fresh-agent `REVERIFY_CLEAN` (coupling meta-test: each
  remediated construct fires under a mutation of its named object; tallies 99/106→95/102, net −4 from the honest de-counts).
  **Register: ONE new counted CALIB knob `Vp0/ℓ_c`** (the breathing driving-force scale; `Vp0` + now-live `ℓ_c` its
  manifestations, `ρ*` consumed; `[Vp0/ℓ_c]=M L⁴ T⁻²` Codex-verified) — the SECOND calibration-adding Part-II stage; legacy
  `{κ,χ,σ_a,σ_L}` = pattern basis (NOT counted afresh, edge R32); edges R32 (structure recovery) + R33 (confinement-drive
  reduction debt, sibling of R30). Codex `REGISTER_CLEAN` (1 nit — stale scope line — fixed). Registration at count 15, PDF
  rebuilt (45pp). ⭐ **pathA_31 fold COMPLETE (stages 013+014+015).**
- **✅ II-G3a `ledger_stage016` DONE (2026-07-09) — pathA_32 `ISOTROPY_CALIBRATED` 1/2 (the ℓ=2 SO(3) covariance theorem; the
  EARNED-FIRST leg of the 2-way split).** SYMBOLIC/float-free. EARNED: the five real ℓ=2 harmonics form one orthonormal
  SO(3)-irrep (`Gram=I₅` by genuine S² integrals) with a single COMPUTED `−Δ_S²` eigenvalue `λ_m=ℓ(ℓ+1)=6` for every m
  (Laplace–Beltrami + Rayleigh + eigenfunction residual, NOT typed), and the K₂ ANGULAR stiffness `K₂=K̃+λ_m·T̃_Ω` uses that
  computed eigenvalue (the `forced_eigenvalue_probe` rejects a typed coefficient → `FAIL_NOT_COVARIANT`). Lands the joint
  `ISOTROPY_CALIBRATED` as a PARTIAL component (016 EARNED, 017 PENDING). ⭐ Register: **ZERO new counted knobs** (structural
  covariance edge R34, like 011/012/014); `T_Ω`/`T̃_Ω` + `β₂(w)` first-appear here but their counting is DEFERRED to 017's
  calibration partition; Part-II CALIB set unchanged = 4; Codex `REGISTER_CLEAN`. ⚠ The source's vacuous `k_coeff_equal`
  `λ−λ≡0` self-compare DE-COUNTED (K₂-coefficient computed-ness on a residual-on-the-assembled-K₂-coefficient + live
  `build_K2(lambdas)` + the bare forced probe). ⭐⭐ Directive review = Codex→Grok→Codex bookend: Codex 2 BLOCKING → folded;
  **Grok-4.5 compute-verify caught 1 BLOCKING — the volume-vs-line dimensional-convention mismatch** (pathA_32 = VOLUME
  densities on `a²dwdΩ` / dimensionless `β₂`; stage013 = LINE dims on `4π∫dw` / `β=L⁻¹`; related by `∫a²dΩ`≈L², NOT equal →
  the cross-stage dual-site `K_η=T_wβ²` is NON-TRANSFERABLE → reframed to provenance + self-contained dimensional integrity)
  + 5 nits, all folded → Codex final confirm `DIRECTIVE_CLEAN`. The `.wl` is ALREADY a native independent engine (only the
  scratch-YAML handoff severed). Dual-engine SymPy 82 / Mathematica 91, CWD-independent; tri-review CLEAN (`FIDELITY_CLEAN` +
  `ADVERSARIAL_CLEAN` per-tooth ablation — all 82 SymPy + 4 native `.wl` ablations fire at their own assert) → 1 low-severity
  stamped-literal (`participates_in_verdict`) remediated to a computed verdict-propagation check → fresh-agent
  `REVERIFY_CLEAN`. Registration at count 16, PDF 47pp. ⚠ **Lesson banked:** across a reduction-convention boundary
  (volume↔line densities, dimensionless↔`L⁻¹` profile), a cross-stage relation or dimension-identity may NOT transfer —
  cite as PROVENANCE + rely on self-contained dimensional integrity, not a theatrical cross-stage dual-site.
- **✅ II-G3b `ledger_stage017` DONE (2026-07-09) — pathA_32 `ISOTROPY_CALIBRATED` 2/2 (the grouped-P2 lane isotropy; the
  COMPLETING leg); ⭐ pathA_32 fold COMPLETE (016+017); the THIRD calibration-adding Part-II stage.** MIXED symbolic+numeric. The
  grouped {20,21,22} ℓ=2 lanes — assembled from 016's CONSUMED covariance (`λ_m=6`, `K₂=K̃+λ_m·T̃_Ω`, Gram-diagonal `c_self=1`) via a
  **GENUINE cross-stage dual-site** (same pathA_32 convention → a one-site `λ`/K₂-form corruption is CHECKABLE + fires; the coordinated
  escape closed by a single-`Y20` `(−Δ_S²)Y=λY` echo; ⚠ **the λ-dimensionless trap** — `λ:6→4` is dim-silent, window-positive, AND
  leaves raw-D isotropy incidentally unbroken (all-lanes λ=4 keeps raw-D=0) → an EXPLICIT integrity check, NOT an incidental gate) —
  respond ISOTROPICALLY: **raw-D defects=0 PRIMARY** ∧ normalized-u defects=0 CROSS-CHECK (a pure-prefactor anisotropy MOVES raw-D but
  leaves normalized-u zero → raw-D decisive). The 6-probe aggregate battery {pure_prefactor/sector_selective/m_dependent/
  degenerate_beta/singular_denominator/static_drop_inertia} (neuter-one flips; ⚠ the 3 anisotropy probes' forced verdict is
  vacuous-by-design → load rides the computed move-flags only). Joint `verdict_from_gates` lands `ISOTROPY_CALIBRATED` COMPLETE
  (016∧017). ⭐ **Register: 2 new counted CALIB knobs `{T_Ω, β₂}`** (the ℓ=2 angular-stiffness density + the frozen ℓ=2 radial profile —
  both calibrated NOT derived from the Gate-1 R0 support equation, edge R36 = the `ISOTROPY_PASS` target); `M̃/K̃/T̃_Ω` = DERIVED
  manifestations `∫density·β₂²` (edge R35); the port-kernel support/Maxwell scalars `{B̃,Z̃}` tracked/downstream-pinned (isotropy
  value-independent of them; Z̃ Maxwell couplings pin 018–024); `μ_η/T_w` provenance (013; `K_η=T_wβ²` R29 non-transferable); R34
  backfilled into the edges table; Part-II CALIB set → 6; Codex `REGISTER_CLEAN`. ⭐⭐ **Directive review = Codex→Grok→Codex bookend:**
  Codex `DIRECTIVE_CLEAN` first-pass → **Grok-4.5 compute-verify caught 1 BLOCKING** (§5 EXCLUDE wrongly stripped the harmonics/`intS2`
  infrastructure 017 needs for its anisotropy coefficients + the echo) + 2 nits (export coeff pinned `K̃+6·T̃_Ω`; the anisotropy probes'
  vacuous forced-verdict flagged), all folded → Codex confirm `DIRECTIVE_CLEAN` (Grok compute-confirmed all math incl. all-lanes-λ=4-keeps-raw-D=0).
  Consumes 016 via the genuine dual-site; `c_S` NOT consumed. **The `.wl` is ALREADY a native independent engine** (only the scratch-YAML
  severed). Dual-engine SymPy 118 / Mathematica 127, CWD-independent. Tri-review: `FIDELITY_CLEAN` (independent re-derivation of the
  anisotropy coeffs {2/7,1/7,−2/7}, the Y20 echo, the pure_prefactor discriminator, the six probe verdicts) + `ADVERSARIAL_ISSUES` (42
  per-tooth ablations, 41 firing at their own assert) → **1 finding: the `.wl` Site-B (K₂-form) dual-site tooth was SUBSUMED** (reconstructed
  K₂ from `lambdaByChannel` instead of reading the assembled lane, so an assembly-formula corruption fired only downstream) → **Codex
  remediated** (rewired to read the assembled K₂, mirroring SymPy, + matching assembly-formula ablation teeth both engines) → fresh-agent
  `REVERIFY_CLEAN` (coupling meta-test: the fixed tooth fires at its own named assert on an assembly-formula corruption both engines,
  neutering the fix makes it vacuous; no regression). Registration at count 17, PDF rebuilt (49pp). ⚠ **Lesson banked:** a consuming
  stage's dual-site "second site" must read the ASSEMBLED downstream object (not reconstruct it from the same cited datum) — else it is
  subsumed by the first site and doesn't independently guard the assembly formula (the per-tooth ablation caught the engine asymmetry vs
  the SymPy sibling).
- **✅ II-G4a `ledger_stage018` DONE (2026-07-09) — pathA_33 `QUAD_CALIBRATED` 1/4 (the outgoing ℓ=2 DtN Hankel fingerprint + χ_Q sign;
  the EARNED-FIRST leg of a 4-way split).** SYMBOLIC/float-free (like 016). EARNED: the outgoing ℓ=2 DtN response `Ŷ₂ᵒᵘᵗ=−3/Λ₂ᵒᵘᵗ`,
  `Λ₂ᵒᵘᵗ=z·h₂⁽¹⁾′/h₂⁽¹⁾`, `z=aω/c_s`, `h₂⁽¹⁾=j₂+i·y₂` (explicit rational·sin/cos), series-expands to the DERIVED fingerprint
  `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`, `v₅=a⁵/27c_s⁵` (the `27`=`1/v₅ᶻ` COMPUTED — the `derived_in_gate` `27` that 020's `54/5=2·27/5`
  partition rides); the sign `χ_Q=+1` outgoing / `−1` incoming COMPUTED from `j₂±i·y₂` (a typed χ_Q would be a tautology); the standing
  `j₂` branch → `Λ_stand(0)=+2`, no radiating term (proving `+1/27` is outgoing-BC-selected); passivity (radiating `v₅` from the
  outgoing BC, `3b` DYNAMIC self-ablation — the v1 constant-self_ablation trip-up avoided). Lands the joint `QUAD_CALIBRATED` as a
  **PARTIAL** (018 EARNED; prefactor=019, `54/5` partition + CALIBRATED label=020, μ̂₀-free dim closure=021). ⭐ **Register: ZERO new
  counted knobs** (an EARNED/structural fingerprint slice, like 016/011/012/014). ⚠ **`c_s` (density sound speed) FIRST becomes a LIVE
  object in the Part-II radiative sector here** (a shift from 013–017's `kξ≪1` deferral) — but only as the units-restoring carrier (via
  `z=aω/c_s`): the earned rationals + χ_Q are `c_s`-FREE (Codex re-derived `free=[]`), so `c_s` (= `c_s0`, R1 `c_s²=5Kρ⁴/m`, stage005) is
  cited PROVENANCE, NOT a consumed value (distinct from the frozen-wall `c_S`, 011–017). New structural edge **R37** (fingerprint
  provenance; discharges nothing). ⚠ **Consumption is PROVENANCE ONLY (no dual-site)** — 018 is SELF-CONTAINED exterior spherical-Hankel
  algebra (the fingerprint is built from explicit `j₂`/`y₂`; the literal port-kernel `N_n/D_n` consumption is 019's), so — UNLIKE 017's
  genuine cross-stage dual-site — there is no checkable cross-stage relation to guard (a guard on an unused object would be a vacuous
  tooth). The `χ_Q` reconciliation (`+1` vs pathA_22b `≈0.712`, same name / different computation) flagged as a tracked **Part-VII** item,
  NOT merged. ⭐⭐ **Directive review = Codex→Grok→Codex bookend:** Codex `DIRECTIVE_CLEAN` first-pass (compute-verified the series
  `1+z²/9+4z⁴/81+i·z⁵/27`, incoming `χ_Q=−1`, standing `Λ_stand(0)=+2`, the c_s-free-ness) → Grok-4.5 compute-verify `DIRECTIVE_CLEAN`
  (independently re-confirmed the same) + 3 non-blocking clarity nits folded (the `P₀` name collision between the prefactor `N₀/D₀` and
  the fingerprint static slot; the `passivity_from_source` scope; the 3a/3b self-ablation scoped to 018's local gates) → Codex confirm
  `DIRECTIVE_CLEAN`. The `.wl` is ALREADY a native independent engine (native `j2`/`y2` + `Series`/`Coefficient` + own `dimOf`) — only the
  scratch-YAML `Export` severed. Dual-engine SymPy 59 / Mathematica 65, CWD-independent. Tri-review CLEAN both legs (`FIDELITY_CLEAN`
  hand-re-derived the full series + `χ_Q=±1` + the standing slot; `ADVERSARIAL_CLEAN` — 30/32 per-tooth mutations fired at their own
  assert, the 2 non-fires harness artifacts re-covered; no vacuous tooth) → **NO remediation** (the fingerprint firewall, the
  χ_Q-computed-ness, the DYNAMIC passivity self-ablation, the μ̂₀-free cut, the no-theatrical-dual-site, and the native `.wl` all genuine;
  the two documented caveats are non-defects). Register `REGISTER_CLEAN` (Codex SymPy-re-derived `free=[]`). Registration at count 18, PDF
  rebuilt (51pp). ⚠ Also fixed a coverage-file lag (the by-Part / coverage-class tables were stale at 016; brought to 018). Committed
  `4872e8b7`.
- **✅ II-G4b `ledger_stage019` DONE (2026-07-09) — pathA_33 `QUAD_CALIBRATED` 2/4 (the squared-denominator prefactor algebra; the SECOND leg
  of the 4-way split).** SYMBOLIC/float-free (like 018). EARNED: the squared-denominator prefactor `P(ω)=D₀·N/D^cons²`
  (`N=N₀+N₂ω²+N₄ω⁴`, `D^cons=D₀+D₂ω²+D₄ω⁴`) series-expands (via `(1+x)⁻²=1−2x+3x²`) to `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`,
  `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` — SERIES-EXTRACTED off the actual series (NOT typed), checked vs an independent typed reference; the
  `−2D₂N₀` factor-of-2 is the SIGNATURE of the squared denominator; the N/D self-check (plain `N/D` gives `P₂=(D₀N₂−D₂N₀)/D₀²`, the computed
  gap `D₂N₀/D₀²≠0` → probe `3g` fires `FAIL_PREFACTOR_ALGEBRA`; correct object NO_FAIL) with a DYNAMIC 019-local self-ablation (the v1
  constant-`self_ablation` trip-up avoided). ⭐ **Register: ZERO new counted knobs** (edge **R38**, the squared-denominator prefactor-algebra
  provenance). ⚠ **This is where 018's deferred literal `N_n/D_n` port-kernel consumption lands — at the PROVENANCE level, NOT a
  value-consumption dual-site:** the abstract `D₀,D₂,D₄` are 017's exported D-lanes + the `N₀,N₂,N₄` are `build_port_moments`' concrete port
  N-moments (deferred Gate-6 branch data, **emitted-but-never-checked** — grep-confirmed by both the directive Codex→Grok→Codex pass and the
  tri-review); the algebra is PORT-AGNOSTIC (holds for any nonzero `D0..N4`) → NO checkable cross-stage relation, no dual-site (same landing
  as 018; a guard on the unused/deferred moments would be a vacuous tooth). 019 is **units-FREE** (no `c_s`/`a`/`G`/`μ̂₀`, no dim leg — the
  `[P₀^phys]=1` closure is 021's; enforced by a runtime free-symbol guard `== {ω, D₀..N₄}`). ⭐⭐ **Directive review = Codex→Grok→Codex bookend:**
  Codex `DIRECTIVE_ISSUES` (3 nits, no BLOCKING — a line-ref L512→L401, the D-lane/N-moment provenance split, the sibling-symbol wording; all
  folded) → **Grok-4.5 compute-verify `DIRECTIVE_CLEAN`** (independently re-derived `P₀/P₂/P₄` + the plain-`N/D` gap `D₂N₀/D₀²`, grep-confirmed
  `build_port_moments` emitted-but-never-checked) + 1 genuine nit (D0..N4-only sample subs, dropping the source `.wl`'s `a/cs/c/G` carriers)
  folded → Codex confirm (folds 1/3/4 clean + 2 straggler wording spots fixed). The `.wl` is ALREADY a native independent engine (native
  `serW`/`Coefficient`/`FullSimplify` on its own `prefObj`/`plainObj`, typing its own expected `P₀/P₂/P₄`) — only the scratch-YAML `Export`
  severed. Dual-engine SymPy 18 / Mathematica 24, CWD-independent. Tri-review: `FIDELITY_CLEAN` (independent re-derivation of `P₀/P₂/P₄` +
  the N/D gap) + `ADVERSARIAL_ISSUES` (21 live `.py` mutations + static `.wl`; the central firewall / N/D self-check / dynamic 019-local
  ablation all genuine; **2 subsumed/mirror teeth** — tooth 16 (mirror of the P₂ tooth) → **made-genuine** as a swapped-in correct-object
  positive control, tooth 11 (assert on the constant `rerun_gate_logic:True`) → honestly **de-counted** — + a **token de-obfuscation** (the
  narrative tokens `c_s`/`G`/`54`/`chi_Q`/`mu_hat0` had been string-concatenated/`chr()`-assembled to dodge a source grep — the pathA_41
  anti-pattern; substance was fine, units-freeness is enforced by the runtime symbol guard not a grep) + an `assert_no_float` legibility fix)
  → fresh-agent `REVERIFY_CLEAN` (the coupling meta-test on tooth 16 decisive: mutate its object → fires, neuter the fix → vacuous; no
  live-symbol leak; no regression). Tallies 19/25 → 18/24 (net −1 per engine from the honest tooth-11 de-count). Registration at count 19,
  PDF rebuilt (53pp). ⚠ **Lesson banked:** don't hand Codex a grep-based acceptance it can dodge by string-concatenating a literal (the
  pathA_41 anti-pattern resurfaced — caught by BOTH fresh agents; the genuine units-free enforcement is the runtime free-symbol guard).
- **✅ II-G4c `ledger_stage020` DONE (2026-07-09) — pathA_33 `QUAD_CALIBRATED` 3/4 (the `54/5=2·27/5` provenance partition + the CALIBRATED
  verdict label; the THIRD leg — the one that LANDS the CALIBRATED headline).** UNITS-BEARING but exact symbolic/rational for the earned content
  (a reversal from 019's units-free slice: `{c_s,a,c,G}` live in the ALGEBRA, but 020 does ALGEBRA+PROVENANCE, NOT the dimensional-homogeneity gate
  `[P₀^phys]=1` which is 021's — the 020/021 cut is by OPERATION). EARNED: the assembled ℓ=2 magnitude `54Gc_s⁵/(5a⁵c⁵)` decomposes as `54/5=2·27/5`
  — a **SymPy-VERIFIED rational identity BOUND to `target_rhs`/`v5_slot`** (`mag=target_rhs/(G·c_s⁵/(a⁵·c⁵))→54/5`, `27_from_slot=a⁵/(c_s⁵·v5_slot)→27`,
  `compact(mag−2·27_from_slot/5)==0`), ⭐ NOT the source's typed STRING (`.py` L622–627) NOR a bare-literal `Rational(54,5)−2·27/5` tautology (a
  Grok/Codex genuineness catch); the `27` is `derived_in_gate` (018's `1/v₅ᶻ`, cited NOT re-derived), the `2/5`+`G` are `external_bridge_input` (GR
  Burke–Thorne `2G/5c⁵`, `G=GENUINE_BLOCKED`); the a⁻⁵ scaling DERIVED via a-cancellation from 018's frozen `v5_slot` (`a⁵` typed once,
  `derived_power=a_power(gamma_target)−a_power(v5_slot)=0−5`); the Γ5/χ_Q equivalence `54Gc_s⁵/5a⁵c⁵ ⟺ 2G/5c⁵` closes iff χ=1 ∧ 54/27=2
  (`forward=2G(χ−1)/5c⁵`). The 4-way PROVENANCE partition (`classify_provenance` tag-dominance `deferred>external>derived>convention`) classifies
  the assembled `54/5` as `external_bridge_input` (external DOMINATES the mixed tags), so the source-faithful verdict (`both derived→QUAD_PASS;
  else→QUAD_CALIBRATED`) lands **`QUAD_CALIBRATED` not `QUAD_PASS`**. ⭐ **The verdict is PROVENANCE-driven, NOT `G→λG`-invariance-driven** — a
  SEPARATE g-invariance diagnostic exposes the invariance-only TRAP (`54/5` is G-invariant yet calibrated → an invariance-only test would MISLABEL
  it as earned). ⭐ **Register: ZERO new counted knobs** (`G=GENUINE_BLOCKED` already registered; `c`=the `c_γ` GR-units-bridge cited benchmark
  `P₀∝c_s⁵/c⁵=1/λγ⁵`; `2/5`=GR bridge; `27`=018's derived); new structural edge **R39**; Part-II CALIB set UNCHANGED = 6; Codex `REGISTER_CLEAN`
  (1 wording fold: the `dimensional_ok`-independence ablation was remediated to a structural signature/bytecode cut). ⭐⭐ **Directive review =
  Codex→Grok→Codex bookend:** Codex `DIRECTIVE_ISSUES` (4 BLOCKING genuineness gaps — the strengthened-3c duplicated-literal, the trivially-true
  `partition_ok`, the missing verdict positive control, the unsound no-μ̂₀ mechanism incl. the false "ZERO μ̂₀ refs" claim vs source L525's
  `gamma_quad_eff` string — + 3 nits) → folded → Codex confirm `DIRECTIVE_CLEAN`; **Grok-4.5 compute-verify CAUGHT the CALIBRATED verdict-rule
  INVERSION** (the source default is `else→QUAD_CALIBRATED`; the directive had inverted it to `PASS-unless-both-external` — a shippable rig since both
  endpoint tests still pass) + a bare-literal-identity nit → folded (added a REQUIRED MIXED control one-external→CALIBRATED) → Codex confirm caught 2
  consistency-sweep gaps (stale bare-literal forms + the mixed control unpropagated) → swept → final Codex confirm `DIRECTIVE_CLEAN`. The `.wl` is a
  **GENUINELY AUTHORED** independent route (the source `.wl` had ZERO 020 content — native `Exponent`/`Together`/`Cancel`/`FullSimplify` +
  `Rational`/`Simplify` + a rank-based `MaximalBy` `Association` classifier + `/. G->lambdaG G`, NOT a `.py` mirror, NO `Series`). Consumes 018's
  χ_Q=+1/`27` (enter the self-contained equivalence bridge) + 019's `P0=N0/D0` (`Gamma5`-definitional-only) + 017's D-lanes as PROVENANCE (no
  dual-site). Dual-engine SymPy 74 / Mathematica 82, CWD-independent. Tri-review: `FIDELITY_CLEAN` (independent re-derivation of the Γ5 bridge, the
  bound `54/5=2·27/5`, the classify_provenance dominance, the g-invariance trap; confirmed the `.wl`'s genuinely-different rank-based classifier) +
  `ADVERSARIAL_ISSUES` (~65 `.py` mutations + 7 native `.wl` mutants; the 4 key genuineness risks — bound identity, MIXED-control-catches-inverted-rule,
  proven classifier, genuine `.wl` — all CLEAN; **3 LOW-severity vacuous/subsumed teeth** [a near-tautological `dimensional_ok`-independence `f(x)==f(x)`
  in both engines; a subsumed `P0=N0/D0`-disjoint firewall; a subsumed tag-mutation `≠EXTERNAL` weaker than its `==DERIVED` sibling]) → Codex
  remediated all 3 make-genuine, **no de-counts** (the structural signature/bytecode dim-cut; the `{N0,D0}⊆Gamma5 ∧ ∉residual` before/after
  run-before-form; the independently-classified baseline-external "from" side) → fresh-agent `REVERIFY_CLEAN` (coupling meta-test: each fires at its
  own named assert + goes vacuous when neutered, both engines, no regression). Registration at count 20, PDF rebuilt (56pp). ⚠ Also fixed a
  coverage/manifest YAML lag (the count fields had stalled at 18 while the tables read 19; brought to 20). ⚠ **Lesson banked:** Grok's compute-verify
  catches a *rule inversion* that a reasoning-only review + a self-report can miss — a verdict whose two endpoint tests both pass under an inverted
  default is a shippable rig; pin the default with a MIXED control (the middle case). Committed (see git log).
- **✅ II-G4d `ledger_stage021` DONE (2026-07-09) — pathA_33 `QUAD_CALIBRATED` 4/4 (the μ̂₀-free `[P₀^phys]=1` dimensional closure; the FOURTH,
  COMPLETING leg — LANDS the joint `QUAD_CALIBRATED` COMPLETE); ⭐ pathA_33 fold COMPLETE (018∧019∧020∧021).** EXACT symbolic `(L,M,T)`-triple
  dimensional-vector algebra, float-free. EARNED: the μ̂₀-FREE gate `dimensional_ok=(dim_of(P₀^phys)==ZERO_DIM)` — `P₀^phys=(c_s/a)²·(N₀/D₀)` is
  dimensionless from the SOURCED port dims `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M` (`[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`); the natural-units trap (the
  handoff `P₀=N₀/D₀` drops `(c_s/a)²`) is CAUGHT (`3d` `FAIL_DIMENSIONAL`); the corrupt-dim SCOPE truth-table {`[N₀]`:FAIL, `[D₀]`:FAIL, `[G]`:NO_FAIL,
  `[c_s]`:FAIL, correct:NO_FAIL} — corrupt-`[N₀]`→`[P₀^phys]=(1,−1,0)` (the `(c_s/a)²` factor REMAINS, NOT `−[D₀]`), corrupt-`[G]` NO_FAIL (`G∉free_symbols`,
  a scope diagnostic). ⭐ **The decisive anti-v1 tooth** (the v1 REJECTION locus): the μ̂₀ back-solve is a TAUTOLOGY (re-solving `[μ̂₀]` keeps
  `homogeneity_pass=True` under EVERY corruption → fires on nothing), correctly DEMOTED to a non-verdict diagnostic — so the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]`
  FAIL rows are what reject it (a computed read-set exclusion + a wired back-solve mutant re-run per corruption, all-NO_FAIL). The `Yhat` dimensionless
  check wrapped in a structured `try_dim_of`/`Catch` (a corrupt ω-power fires the NAMED assert, not an uncaught `DimError`); the `3d`/`3d′` self-ablations
  DYNAMIC 021-local re-runs (`rerun_gate_logic` derived from the actual re-run). ⭐ **Register: ZERO new counted knobs** (μ̂₀ = free-carrier NON-verdict
  diagnostic; the SOURCED port dims `[N₀]`=pathA_43 density-port numerator / `[D₀]`=carried reduced static conservative denominator `D₀=K−B₀−Z₀` = dimensional
  PROVENANCE; `c`/`G` already registered) + structural edge **R40**; Part-II CALIB set UNCHANGED = 6; Codex `REGISTER_CLEAN`. ⭐ **021 IS units-bearing AND
  does the `[·]` dim-homogeneity gate** — the OTHER half of the 020/021 operation-level cut (020 did algebra+provenance, 021 does dimensions), the ONLY
  pathA_33 leg with a dimensional gate. Consumption is PROVENANCE (the sourced port dims genuinely ENTER the gate so the corrupt-`[N₀]` tooth is genuine;
  018's `u₂/u₄/v₅` a local frozen fixture for `Yhat`; 019's `P0=N0/D0` enters `P0_raw`; 020's `target_rhs` the μ̂₀ diagnostic's rhs — NO dual-site). ⚠
  **COMPLETE ≠ PASS:** the joint token STAYS `QUAD_CALIBRATED` (calibrated; 020's provenance + `G=GENUINE_BLOCKED`). ⭐⭐ **Directive review = Codex→Grok→Codex
  bookend: Codex CAUGHT 6 BLOCKING** — most notably the ORIGINALLY-PROPOSED corrupt-`[G]` anti-v1 discriminator was BACKWARDS (a back-solved μ̂₀ gate re-solves
  `mu_dim` after every mutation → fires on NOTHING, so it is the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]` FAIL rows that reject it, NOT `[G]`); + the corrupt-`[N₀]`
  dim `(1,−1,0)` not `−[D₀]`; the self-containment local fixture (`build_dimensions`/`YhatPhysical` referenced 018's `u2/u4/v5`); the `[D₀]` provenance; the
  `homogeneityPass`-in-the-`.wl`-guard; the `Yhat` structured-catch — all folded → Codex confirm (1 residual dim-label) → **Grok-4.5 compute-verify
  `DIRECTIVE_CLEAN`** (SymPy-confirmed the μ̂₀-free gate dims + the corrupt-dim truth-table + the back-solve-is-a-tautology crux + the `Yhat` catch;
  reproduced the contrast that the original framing held only for a PINNED μ̂₀) → final Codex confirm `DIRECTIVE_CLEAN`. The `.wl` KEEPS its native dim block
  L101–174 (like 018/019, unlike 020 — a `Which`-based `dimOf` + native `rawDims` Association + `Series`-free dim algebra; severs only YAML, REMOVES the
  018/019 blocks, replaces `Yhat` coeffs with local `u*Sourced`, ADDS the corrupt-`[G]` control + DYNAMIC self-ablation, `yhatOk` split out of the guard).
  Dual-engine SymPy 42 / Mathematica 50, CWD-independent. Tri-review: `FIDELITY_CLEAN` (independent hand re-derivation of the μ̂₀-free gate dims, the corrupt-dim
  truth-table, the back-solve-tautology table; no 018/019/020 leakage; genuine independent `.wl`) + `ADVERSARIAL_ISSUES` (42 `.py` per-tooth ablations + 6 live
  `.wl` mutations; the two KEY teeth — anti-v1 read-set/wired-mutant + corrupt-`[G]` scope control — GENUINE; 5 LOW-severity stamped/subsumed teeth) → Codex
  remediated 2 make-genuine (`rerun_gate_logic` derived from the actual re-run; `participates_in_verdict` derived from the computed read-set) + 3 de-count (the
  QUAD-landing literal==literal tautology + two subsumed aggregate summaries, retained as labeled prints) → fresh-agent `REVERIFY_CLEAN` (coupling meta-test:
  each make-genuine tooth fires under a mutation of its object + goes vacuous when neutered, both engines; the de-counts keep per-row/per-mutant coverage; no
  KEY-tooth/earned-logic regression). Tallies 45/53 → 42/50 (net −3/engine from the honest de-counts). Registration at count 21, PDF rebuilt (59pp). ⚠ Also
  refreshed 3 prior stages' (016/019/020) committed output transcripts to the runner-header convention (they had been committed by direct-run, lacking the
  `# Date`/`EXIT_CODE` stamp; both engine summaries now show 21 pass). ⚠ **Lesson banked:** the Codex design-review caught a directive-level REASONING error
  (the anti-v1 discriminator was backwards — corrupt-`[G]` fires a *pinned*-μ̂₀ gate, but the v1 rig RE-SOLVES μ̂₀ so it fires on nothing); the discriminator
  that a gate is genuinely able-to-fail is which corruptions MAKE IT FIRE, not a single negative control. Committed (see git log).
- **✅ II-G5a `ledger_stage022` DONE (2026-07-10) — pathA_34 `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` 1/2 (the cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the Gate-4
  non-regression; the EARNED-FIRST leg of the pathA_34 2-way split).** SYMBOLIC/float-free/z-space (like 018, generalized cross-ℓ). EARNED: for ℓ=0,1,2 the
  outgoing DtN `Ŷ_ℓᵒᵘᵗ=−(ℓ+1)/Λ_ℓ` series-expands to the DERIVED radiative fingerprint `{1, 1/2, 1/27}` at orders `{ω¹, ω³, ω⁵}` (COMPUTED from the spherical
  Hankel, derive-vs-typed); the static slot `Λ_static=−(ℓ+1)` DERIVED from the Hankel log-derivative (`lam_series.coeff(z,0)`, NOT the source's hand-set-numerator
  X≡X — de-rigged); the radiative order verified by SCANNING the imaginary series for its first nonzero power + all-lower-vanish (NOT the source's dodgeable
  preselected-`2ℓ+1`-nonzero check — a `+i·z` corruption dodges it, the scan catches it); the incoming branch flips ONLY the radiative sign (the genuine sign
  tooth); the ℓ=2 leg (`−3/Λ₂`) reproduces the completed pathA_33 quadrupole `{u₂=1/9, u₄=4/81, v₅=1/27}` (the Gate-4 non-regression, derive-vs-typed vs 018's
  independently-earned literals — NO subsumed X≡X; per-slot u₂/u₄/v₅ mutants since `3e` flips only v₅). ⭐ **Two distinct verdicts printed:**
  `LOCAL_AUDIT_VERDICT=CROSS_L_FINGERPRINT_OK` (the exit-0 gate; read-set = `{cross_l_fingerprints, ell2_non_regression}` provably EXCLUDES any
  nullspace/`base_verdict`, a computed read-tracking guard) + `JOINT_LANDING_LABEL (PARTIAL): FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)` (a printed string). 022 is
  the earned half of a gate that ultimately fails (cf. stage003); the DEFAULT-verdict-is-PASS trip-up avoided (022 owns NO nullspace + does NOT back-solve
  `ε_eff`/`Z` — the `FAIL_TAUTOLOGICAL` firewall preserved for 023). ⭐ **Register: ZERO new counted knobs** + edge **R41** (the cross-ℓ fingerprint + non-regression
  provenance; discharges nothing); `c_s`=R1 units-carrier provenance, `a`=`CONV`; ⚠ the ℓ=0/1 stiffnesses `K0c`/`K_eta+2·T_Omega` + `Z0_ret`/`Z1_ret` are **023's**,
  NOT counted; `static=1`+`χ_Q=1` de-counted diagnostics (subsumed by `Λ_static` / by `v₅`); Part-II CALIB set UNCHANGED = 6; Codex `REGISTER_CLEAN`. **Consumption:**
  stage018's `{u₂,u₄,v₅}` the one CHECKABLE derive-vs-typed non-regression; 019/020/021 + the completed joint + 008's amplitudes + 009/010's bulk mode + `c_s` (R1) =
  PROVENANCE (no dual-site); z-space only (NO units-restored dim leg — 023's). ⭐⭐ **The `.wl` is RE-AUTHORED independent** via built-in `SphericalHankelH1`/`H2` +
  `SeriesCoefficient` (the source `.wl`'s transliterated `branchData` DISCARDED per the mirror-policy screen — UNLIKE 018/019/021's keep-native; 020-style
  genuine-authoring). ⭐⭐ **Directive review = Codex→Grok→Codex bookend:** Codex 5 BLOCKING folded (the `.wl` re-authoring; the de-rigged `Λ_static`+scanned-order;
  the ℓ=2 double-count + χ_Q subsumption + per-slot mutants; the Earned/Deferred SUBSET framing + explicit LOCAL/JOINT output; the checkable-consumption = stage018
  fingerprint) → **Codex confirm-pass caught a directive-level REASONING error** — the proposed `Λ_static` ablation "outgoing→incoming" is INERT (both branches give
  `−(ℓ+1)`), so the mutant must be a POLE-ORDER corruption `h_mut=z·h`, and `static=1` de-counted → **Grok-4.5 compute-verify `DIRECTIVE_CLEAN`** (SymPy-confirmed the
  cross-ℓ series, the pole-order-vs-inert mutant, the scan-vs-preselect counterexample, the χ_Q subsumption identity, the u₂/u₄ isolation) → final Codex confirm
  `DIRECTIVE_CLEAN`. Dual-engine SymPy 56 / Mathematica 65, CWD-independent. Tri-review: `FIDELITY_CLEAN` (independent re-derivation of the fingerprints, the
  pole-order-vs-inert `Λ_static` mutant, the scan-vs-preselect order, the χ_Q subsumption, the per-slot isolation — all COMPUTED, no 019/023 leak, genuine built-in
  `.wl`) + `ADVERSARIAL_ISSUES` (both confirm-pass catches HELD live-proven; all core physics teeth genuine; 2 LOW-severity redundancy teeth — F1 the `3e`
  `rerun_gate_logic` constant-`len(trace)==2`, F2 the read-set-excludes tooth subsumed by the exact-equality tooth) → Codex remediated F1 make-genuine (compare the
  two traced verdicts) + F2 honest de-count (a diagnostic; the exact-equality tooth retains the guard) + 2 `.wl` nits (a symbol typo + `IncomingLowerRealUnchanged`
  mutation-aware parity) → fresh-agent `REVERIFY_CLEAN` (coupling meta-test: F1 fires when the ablation is neutered + goes vacuous when the fix is reverted; F2's
  retained exact-equality tooth fires on a wired forbidden read; no regression, both engines exit 0 at repo root AND `/tmp`). Tallies 56/65 → **55/64** (net −1 per
  engine from the honest F2 de-count). Registration at count 22, PDF rebuilt (61pp). ⚠ **Lesson banked:** the Codex confirm-pass caught a directive-level reasoning
  error a design-review missed — an able-to-fail mutant "outgoing→incoming" that is INERT for the leading pole (both Hankel branches share the `z^{−(ℓ+1)}`
  singular behavior; only the imaginary sign flips), so `Λ_static` needs a POLE-ORDER mutant `h_mut=z·h`; the proof a normalization tooth is able-to-fail is a
  mutant that shifts the pole ORDER, not the radiation branch.
- **✅ II-G5b `ledger_stage023` DONE (2026-07-10) — pathA_34 `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` 2/2 (the native nullspace underdetermination departure — the
  COMPLETING, FAIL-DELIVERING leg); ⭐ pathA_34 fold COMPLETE (022∧023); ⭐⭐ CLOSES the Gate-1–5 gravity ladder.** SYMBOLIC/exact/float-free. EARNED: over **11
  genuine generator dofs** `[OmegaU,OmegaW,Rmix,gU,gW,D0,K0c,K_eta,T_Omega,Z0_ret,Z1_ret]` the collected Gate-5 named constraints `{P0_raw, K0c, K_eta+2·T_Omega}`
  have constraint-Jacobian rank **3** (a GENUINE `sp.Matrix(rows).rank()` on symbolic `diff` rows / constructive `NullSpace` — **NO zero-padding, NO hardcoded 8/2**;
  the v1 REJECTION locus repaired) → native nullspace dim **8** → return-augmented rank **5** → **return-moving nullity 2**: the return admittances `{Z0_ret, Z1_ret}`
  survive every constraint (with explicit unit-vector witnesses that preserve every constraint yet move `T0/T1`; `Z_is_premise=True`, pathA_29) →
  `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`. The **counterfactual** selector `{Z0_ret=K0c, Z1_ret=K_eta+2·T_Omega}` collapses the return-moving nullity 2→0 (native nullity
  8→6, rank 3→5) → `CROSS_L_RESIDUAL_PREDICTION` — the able-to-fail witness (the DEFAULT verdict is the predictive token, so the FAIL is EARNED from the computed
  return-moving nullity, NOT baked). The `A0/A1` scalar/dipole residuals FORWARD-consume 022's `{1,1/2}` (`A_ℓ=i·v_ℓ·(aω/c_s)^{2ℓ+1}·{M0|D1}·(1−T_ℓ)`, `ε_ℓ=Z_ℓ/K_ℓ`
  FORWARD not back-solved), checked vs the INDEPENDENT pathA_29 form (`expected_A1`'s `2c_s³` encodes the consumed `1/2` → corrupt `v1`→`A1_form` fires); the
  stage021-machinery dim gate (`[A0]=(0,1,−1)`, `[A1]=(1,1,−1)`, `[P0_physical]=(0,0,0)`; sourced-`[M0]`→`FAIL_DIMENSIONAL`, free-carrier `q_free`→`NO_FAIL`); the
  strengthened `FAIL_TAUTOLOGICAL` firewall (`class_matches_computed`; `ε_eff` magnitude classed `deferred_branch_data`; the NEW `emit_epsilon_magnitude_as_derived`
  mutation → `FAIL_TAUTOLOGICAL`) forbids the ε_eff/Z back-solve. The joint ladder DROPS the `quad_regression` rung (022 owns Gate-4, consumed `quad_regression=False`
  provenance — 023 does NOT rebuild the fingerprint core / probe 3e) + removes the inert `able_to_fail_bad`. Prints `AUDIT_STATUS=PASS` (script/teeth, exit 0) distinct
  from `PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE (2/2, COMPLETING)` (the earned characterized-FAIL — the stage003/pathA_36 pattern). ⭐ **Register: ZERO new
  counted CALIB knobs (set stays 6); `Z0_ret/Z1_ret` add zero new free dofs (aliases); `K0c/K1` add COUNTED `FREE-UNREDUCED` PENDING reduction-debt** — `Z0_ret/Z1_ret` = COORDINATE ALIASES of the existing `ε0/ε1` FREE-UNREDUCED debt (register L165; no double-count, no new dof);
  ⚠ `K0c` + the ℓ=1 sector `{K_eta,T_Omega}` (via `K1=K_eta+2·T_Omega`) = pathA_34-convention effective ℓ=0/1 stiffnesses, **`FREE-UNREDUCED` `PENDING` scalar-reduction — COUNTED as debt** (per the register rule pending debt stays counted until DERIVED; NOT `DERIVED`, NOT `CALIB`; dims `M T⁻²`≠013 `M L⁻¹T⁻²`/017
  `M L⁻³T⁻²` — the stage016 volume-vs-line convention trap); `q_free/eta_null/gain0/gain1` control-construction tracked-not-counted; new obligation edge **R42** (the
  cross-ℓ nullspace departure + the sharpened Gate-6 obligation: Gate 6 must supply 2 independent return equations — SHARPENS the `ε0/ε1` R24-family debt, adds no
  dofs); Part-II CALIB set UNCHANGED = 6; Codex `REGISTER_CLEAN` (pending Codex-verify at commit). ⚠ Consumption: stage022's `{1,1/2}` (CHECKABLE derive-vs-typed →
  `A0/A1`); the pathA_29 residual form (009/010); `Z_is_premise`; 008's `M0/D1/R0/R1`; 017's `P0_raw`; 013/017 (context for `K0c/K1`); `c_s`(R1)/`a`(CONV) = PROVENANCE.
  ⭐⭐ **The `.wl` is RE-AUTHORED independent** via a constructive `NullSpace` route (`Length[NullSpace[Jbase]]=8`, `MatrixRank[basis·Gᵀ]=2`,
  `Greturn·NullSpace[Jselector]=0`) — materially different from the `.py`'s `augRank−rank0` (mirror-policy screen; like 020/022, UNLIKE 018/019/021's keep-native).
  ⭐⭐ **Directive review = Codex→Grok→Codex bookend (⭐ the v1 REJECTION locus, the SHARPEST in Part II):** Codex design-review **7 BLOCKING** folded (the `.wl` decisive
  RE-AUTHOR; ⚠ the **"selector collapses the RETURN-MOVING nullity, not the nullity" math fix** [native nullity 8→6 not 0] + corrected isolated ablation teeth; the
  counterfactual-witness-not-proven-Gate-6-selector relabel; the provenance-cut fix [022's `{1,1/2}`=`cited_earned_input` + `assert_not_derive` rewired to the
  023-derived forward T0/T1 map + `gate4_prefactor` tag dropped]; the strengthened firewall [`emit_epsilon` tooth + de-counted `rerun_gate_logic` + fixed
  `able_to_fail_bad`]; the **K0c/K1 PENDING-not-DERIVED** register correction; the **Z0_ret/Z1_ret aliases-not-new-dofs** register correction) → Codex confirm-pass 5
  BLOCKING + 2 nits → final-confirm 2 BLOCKING + 1 nit (all consistency-sweep gaps) → **Grok-4.5 compute-verify `DIRECTIVE_CLEAN`** (independent SymPy confirmed rank
  3/nullity 8/return-moving 2, the selector flip with native nullity 8→6, `A1−expected_A1=0` iff `v1=1/2`, the dims, the `K0c/K1` dim-conflict + `Z_ret` alias
  conventions; validated the `.wl` constructive route incl. `Greturn·Nsel=0` a genuine identity; 1 honest-scope note folded — raw nullity 8 includes `K0c/K1`
  self-constraint bookkeeping, verdict rides return-moving 2) → closing Codex confirm. Dual-engine SymPy 116/Mathematica 123 → **111/117** (net −5/−6 from honest
  de-counts). Tri-review: `FIDELITY_CLEAN` (independent SymPy re-derivation of the rank audit, selector flip, `A0/A1`+`v1=1/2` consumption, dims — faithful, genuine
  rank not zero-padded, no ε back-solve, materially-different constructive `.wl`) + `ADVERSARIAL_CLEAN` (per-tooth ablation matrix, 15 mutations both engines —
  hardcoding-the-rank/zero-padding/faking-the-`.wl`-basis all FAIL, the 4 KEY anti-rig properties CONFIRMED; **4 non-blocking de-count nits** → Codex remediated 2
  make-genuine [witness-preservation recomputes each Jacobian-row dot product from the stored witness; neutralized-mutation uses a cache-distinct inert context +
  independence check] + 2 de-count [provenance-documentation → labeled PROVENANCE prints; the T/ε identity → a SELF-CONSISTENCY check] → fresh-agent `REVERIFY_CLEAN`
  coupling meta-test, no regression). Exit 0, CWD-independent (repo root AND `/tmp`), zero file I/O. Registration at count 23, PDF rebuilt. ⚠ **Lesson banked:** the
  bookend caught a directive-level MATH error (the selector collapses the *return-moving* nullity to 0, NOT the full native nullity — which goes 8→6) + a register
  over-claim (K0c/K1 "DERIVED" was unsupported — the pathA_34-convention dims `(0,1,−2)` don't match 013/017, the stage016 volume-vs-line trap recurring) + a
  double-count (Z0_ret/Z1_ret are the *same* ε0/ε1 freedoms, aliased) — three directive-level errors the Codex→Grok→Codex bookend caught before the build.
- **▶ NEXT = Cluster C `ledger_stage024`** (II-P1, pathA_43 `DENSITY_PORT_HOSTED` 1/4 — the ℓ=2 quadrupole radiative-port numerator `N0_den` as a DENSITY-NATIVE
  two-port over `(q2` wall, `Φ2` bulk-density`)` via Schur-elim / Green-DtN; consumes 029/010's bulk mode + 032/017's wall mode; the OLD EM `A_w`/`U,W` vector
  scaffold RETIRES). Source: `software/stage1_solver/tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}` (already contract-clean — no argparse/JSON/reads, the
  LIGHTEST reshape; the work is decomposition + preserving the anti-rig mechanisms = the COMPUTED taint/host-set guard + the lineage token-check, 024/025/026's — see
  the per-gate trip-ups L93). ⭐ **Cluster C = 024–029** (density port 024–027 + 2.5PN match-back 028 + PN DOI-cite 029); after Cluster C the Part-II gravity sector
  CLOSES → the scheduled **MIDWAY KNOB AUDIT** (parameter_register §"MIDWAY KNOB AUDIT" — the pathA_40 `Δr=2` codimension dry-run over Parts I–II + the held-out
  vs irreducible-route-less tally). ⭐ Author a running-start source map FIRST (per the calibrated pipeline). This stage is DENSITY-PORT (pathA_43), a DIFFERENT
  source family from pathA_30–34 (contract-clean, no scratch-YAML bridge) — see the reshape-cost map L78–79.

## Per-stage process (unchanged, calibrated)

Same as Part I: author reshape directive → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN` →
**⭐ FINAL Grok-4.5 headless design-review pass on the same prompt** (added 2026-07-08; no GLM on Parts I–VI, but Grok
DOES run now — it compute-verifies the math via SymPy and caught a kernel-preserving residual-tooth defect + 5 nits on
stage013 that Codex's xhigh pass missed) → assess/verify each Grok catch independently → fold → **Codex confirm-pass on
the Grok-folded directive** (final double-check; Codex bookends) → pre-exec USER GATE per
stage → Codex builds only the two scripts → dual-engine exit 0 → arbiter re-run → tri-review (arbiter + fidelity +
adversarial-scoped-to-reshape-integrity with ablation) → remediate → counts bump → parameter-register update +
Codex-verify → note/card/`\input`/registration → PDF → commit. Orchestrator authors notes/cards/LaTeX/registration;
Codex codes. Grok invocation + rationale: blueprint §6 + [[reference-grok-cli-review]].
