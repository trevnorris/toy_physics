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
- **▶ NEXT = Cluster A `ledger_stage011`** (II-G1a, pathA_30 frozen-reduction certificate: frozen wall → Helmholtz
  `ψ''+(ω/c_S)²ψ=0`, `c_S²=5Kρ*⁴/m`, projection measure + certificate).

## Per-stage process (unchanged, calibrated)

Same as Part I: author reshape directive → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN` (no GLM on
Parts I–VI) → pre-exec USER GATE per stage → Codex builds only the two scripts → dual-engine exit 0 → arbiter
re-run → tri-review (arbiter + fidelity + adversarial-scoped-to-reshape-integrity with ablation) → remediate →
counts bump → parameter-register update + Codex-verify → note/card/`\input`/registration → PDF → commit.
Orchestrator authors notes/cards/LaTeX/registration; Codex codes.
