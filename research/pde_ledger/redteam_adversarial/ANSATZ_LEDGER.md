# ANSATZ_LEDGER — the non-derived-value catalog & layer-3 branch-realization to-do list

**What this is.** Every value in the moving-throat superfluid PDE ledger that is *not* derived from first
principles, sorted by **why** it isn't derived. It is the layer-3 work-list: the values the toy model
currently *posits*, separated into "the real PDE would fix this" (the falsifiable to-do list) versus
"this is a literature fit / a self-calibration / a pure convention." The headline section (§1) is the
list the user asked for: **values we can't analytically derive but need, and that the right model running
would determine.**

**Generated from.** The Phase B genealogy corpus
`redteam_adversarial/provenance/_ansatz_corpus.yaml` (133 conceptual families = 103 free_choice + 30
published_target, consolidated from the Phase B genealogy, commit **2495589**), cross-checked against the
4C findings ledger `redteam_adversarial/provenance/_4c_findings_ledger.md` and the sourced external
benchmarks `redteam_adversarial/benchmarks.yaml`. Every (b) and (c) entry below carries a notes
citation (`path:line`) that was opened and verified, not taken from the excerpt cache.

## Four-category key

- **(b) BRANCH-DETERMINABLE** — the priority list. The real PDE *would* select a specific value; the
  ledger posits it pending the full solve. **Falsifiable**: run the right model and it either reproduces
  the value or it doesn't. Signals in the notes: "reference-branch numerical freeze", "carried/preferred",
  "selected/concrete branch", a geometry/branch pick, "not a theorem of the full PDE."
- **(c-lit) EXTERNAL-LITERATURE FIT** — back-solved to match a *published* constant. Exactly two anchors
  in this ledger (GR quadrupole 2G/5c⁵ chain; CODATA m_p/m_e). Cross-referenced to `benchmarks.yaml`.
- **(c-rec) INTERNAL RECORDED-BENCHMARK MATCH** — back-solved/pinned to the model's *own* recorded
  session output (Session-I…V readbacks), not to literature. Calibration-to-self, not a prediction.
- **(a) CONVENTION / GAUGE** — genuinely arbitrary; there is no "true" value to derive (normalizations,
  unblocked/limit illustrations, sign/branch labels, gauge picks, probe sample inputs). The model does
  **not** owe us these.

> **Tag-trust warning.** The corpus `published_target` tag is *unreliable* for the Family-1 canonical
> block. `r_F1` / `S_can` / `Pi_star` are tagged `published_target` by batches 16/17 but they rest on the
> **posited** L/a = 37/20 — no external number backs them. They are reclassified here as **derived-from-(b)**,
> NOT external fits. (Cross-batch split confirmed in `_4c_findings_ledger.md:30-31`.)

## Count summary

| category | families | note |
|---|---|---|
| **(b) branch-determinable** | **13** | the priority list; root posits + the canonical block that inherits them |
| **(c-lit) external-literature fit** | **9** | all trace to 2 anchors: GR 2G/5c⁵ (7) + CODATA 1836… (2) |
| **(c-rec) internal recorded-benchmark match** | **~24** | Session-I…V readbacks (λ_L, V_eff, Ξ_turn, helicity exports, γ_crit, …) |
| **(a) convention / gauge** | **~83** | wall-action coeffs, limit illustrations, free model coefficients, gauge picks |
| **§5 borderline / needs-judgment** | **4** | gamma_0, chi_Q, n=5 EOS, c0=3/4-vs-c0_min duplication |

(Family counts sum to the 133 corpus families; a handful of multi-parameter families straddle categories
and are placed by their load-bearing parameter, noted inline.)

---

## §1 (b) BRANCH-DETERMINABLE — branch-realization targets  ⭐ THE PRIORITY LIST

These are the values the user wants: *can't be analytically derived right now, are needed, and the right
PDE solve would determine them.* Lead row is the headline posit; the indented rows derive **from** it.

| # | concept | current posited value | stages | clean-rational? | what running the real model must reproduce | notes evidence |
|---|---|---|---|---|---|---|
| **1** | **L/a (aspect ratio Λ\*)** | **37/20 = 1.85** | 073, 121, 139, 146, 232 | **YES — and explicitly a placeholder**: notes call it "L/a ≈ 1.85" then freeze the clean 37/20; the real branch would select whatever (likely irrational) aspect ratio the throat geometry actually picks | the moving-throat geometry must select the throat axial-span / mouth-radius ratio; everything in the Family-1 block hangs off this one number | `notes/stages/moving_throat_pde_stage073_family1_geometry_map.md:56` ("**reference-branch numerical freeze**, not a new theorem of the unsolved moving-throat PDE"); `…stage121_geometric_r_selection.md:60` ("L/a = 37/20") |
| 1a | r_F1 (normalized hybridization) | √((12/π²)(37/20)²−1) = √(4107−100π²)/(10π) ≈ 1.77799353547498 | 122,123,124,161,162,165,168 | derived-irrational *from* 37/20 | NOT independent — fixed the instant L_W=L is imposed at the true L/a; reproduced iff 37/20 is | `…stage121_geometric_r_selection.md:60`; `…stage122_mouth_source_compensation_test.md:49` (r_F1 = √(4107−100π²)/(10π)) |
| 1b | r_c^F1 (= r_F1²) | ≈ 3.16126101219081 | 121, 162 | derived *from* 1a | inherits 37/20 via r_F1² | `…stage121_geometric_r_selection.md:~78` (boxed r_c^F1 = r_F1² ≈ 3.16126…) |
| 1c | S_can / Σ0_can / T_can / Π_can / Π_star / g_star (canonical Family-1 block) | S_can ≈ 4.651033550168…; the …867/…876 sub-1e-13 drift is a separate transposition note | 156, 158, 163 | derived *from* 1a | the whole canonical-point block is `internal_consistency` at the SOLVE level but bottoms out on 37/20 → reproduced iff the real branch reproduces L/a | corpus `fam_0228_s_can` & `fam_0209_pi_star` both list **origin_stage 121** → `…stage121…:60`; precision-drift note `_4c_findings_ledger.md:33-35` |
| **2** | **epsilon_r (Family-1 wall width)** | **1/20 = 0.05** | 073 | **YES — clean rational, sibling reference-branch freeze** (same geometry-freeze family as 37/20); a placeholder for the true wall/healing width | the real wall profile must set the thin-layer width ε_r (and hence ℓ = ε_r·a) | `…stage073_family1_geometry_map.md:18` ("balanced thin-layer-consistent Family-1 radial wall branch `epsilon_r = 0.05`"); `:30` ("epsilon_r = 0.05 = 1/20") |
| 2a | ell (support / healing scale ℓ) | ℓ = ε_r·a (controlled local closure) | 074 | derived *from* 2 | the parent GNLS medium must fix the support scale instead of leaving it free | `…stage074_family1_healing_lock.md:55` ("a **controlled local closure**, not yet a theorem of the full moving-throat PDE … first honest way to tie the support scale to the parent GNLS medium") |
| **3** | **c0 (minimal isotropic conservative module)** | **3/4** (c1 = 1/4) | 088, 090, 240 (as c0_min) | YES — "smallest viable" posit; a candidate, not forced | the real conservative quadrupole precursor must select the contact/pole split; ρ_α=4/3, ζ_req=1/3, Π_tr chain all derive from 3/4 | `…stage088_loading_ratio_from_minimal_module.md:81` ("the smallest viable isotropic conservative precursor to `c0 = 3/4, c1 = 1/4`"); `…stage090_updated_reduced_status.md:15` ("natural contact-plus-pole realization … selects ρ_α=4/3, ζ_req=1/3"); `…stage240…:249` ("The minimal isotropic conservative quadrupole module is c_0=3/4") |
| 3a | c_contact / c_0 / c0_min (= same 3/4 posit under different stage labels) | 3/4 | 088, 090, 240 | duplicate of 3 | same as 3 (see §5 dedup note) | `…stage088…:83` ("`c0 = 3/4,`"); `…stage240…:249` |
| **4** | **chi_Q (outgoing-normalization multiplier)** | **1** (alt branch value 9) | 100, 105 | YES (1) — a branch selection, not arbitrary | the actual retarded/outgoing l=2 branch must fix χ_Q (=1 was matched to the stage-104 outgoing DtN h₂⁽¹⁾ fingerprint at 105); "all genuinely new retarded uncertainty sits in one number only: χ_Q" | `…stage100_outgoing_normalization_factorization.md:79` ("all deviations from the exact compact outgoing `l=2` fingerprint are captured by the single multiplier `chi_Q`"); origin/fix at 105 per `_4c_findings_ledger.md:38` |
| **5** | **g_nat (normalized mouth-coupling ratio)** | branch ansatz (O(1), set by the geometric r_F1) | 122, 124 | branch-selected | the actual core must pick the normalized mouth-coupling 𝔤 at the geometrically fixed r_F1 ("what normalized mouth-coupling ratio does the actual core pick?") | `…stage122_mouth_source_compensation_test.md:30` ("a concrete branch ansatz, not as a theorem of the full PDE: the cleanest 'same mouth source, same normalized loading' core closure") |
| **6** | **K_m (local Robin wall closure)** | first natural Robin closure | 071, 073 | branch-selected | the parent PDE must fix the wall-shell Robin constant; "not a universal theorem of the parent PDE; the first natural local Robin closure" | `…stage071_tanh_wall_branch.md:95` ("This is not a universal theorem of the parent PDE; it is the first natural local Robin closure consistent with the same wall-shell support scale") |
| **7** | **L_W = L (mixed-tube ↔ throat identification)** | L_W := L | 121 | identification choice the solve would confirm | the real geometry must confirm the mixed D/N side-channel is the actual throat corridor (the identification that *forces* r_F1) | `…stage121_geometric_r_selection.md:26` ("\boxed{L_W=L.}") |
| **8** | **eps_2 / eps_4 (geometry-contamination numbers)** | posited 0 on the natural branch | 092, 093 | YES (0) — "or prove that they vanish" | the real branch must either vanish these or fix them; explicitly an open question, not a theorem | `…stage093_grouped_p2_status_update.md:55` ("determine the two geometry-contamination numbers `eps_2` and `eps_4`, or prove that they vanish on the natural branch") |

**§1 read-out (concept · value · what the model must give):**
1. **L/a = 37/20** — the throat aspect ratio the moving-throat geometry actually selects (clean-rational placeholder for a likely irrational). HEADLINE.
   - r_F1 ≈ 1.778, r_c^F1 ≈ 3.161, and the canonical block (S_can/Σ0_can/T_can/Π_can/Π_star/g_star) all inherit it.
2. **epsilon_r = 1/20** — the true Family-1 wall/healing width (and ℓ = ε_r·a).
3. **c0 = 3/4** — the conservative quadrupole contact/pole split the real precursor selects (ρ_α=4/3, ζ_req=1/3 follow).
4. **chi_Q = 1** — the outgoing l=2 normalization multiplier the actual retarded branch fixes.
5. **g_nat** — the normalized mouth-coupling ratio the actual core picks at the fixed r_F1.
6. **K_m** — the wall-shell Robin closure constant the parent PDE fixes.
7. **L_W = L** — the geometric identification (mixed tube = actual throat span) the solve must confirm.
8. **eps_2, eps_4 = 0** — the geometry-contamination numbers the natural branch must vanish (or fix).

---

## §2 (c-lit) external-literature fits

Two published anchors only. All families trace to one of them; all are **honestly disclosed** in the
notes/cards (StatusOpen / "target values" / "proton-proxy") — the audit confirmed disclosure adequacy,
no `paper_card_overclaim` on these.

| concept (family) | posited value | published constant matched | origin stage | benchmarks.yaml id | notes evidence |
|---|---|---|---|---|---|
| P0_target | 54·G·c_s⁵/(5·a⁵·c⁵) (= 54/5 ≡ 27/20 coeff) | GR mass-quadrupole 2G/5c⁵ (Peters 1964 / Maggiore Ch.3) | 022 | `bench_2e5805a359c4` | `notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md:256` ("The universal GR quadrupole target is") |
| N_Q | 54/5 | same GR 2G/5c⁵ | 022 | `bench_62250fff83d4` | `…stage030_selected_mode_normalization.md:191` ("N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)") |
| Gamma_5 | 2G/5c⁵ | same GR 2G/5c⁵ | 106 (target fit lives here) | `bench_a187d7865e33` | `…stage106_canonical_outgoing_reduced_closure.md:40` ("the canonical … quadrupole coefficients are fixed to their target values") |
| Gamma_5_normalized | GR/Burke–Thorne target | same GR 2G/5c⁵ | 106 | (chain of `bench_a187d7865e33`) | `…stage106…:53` ("the normalized odd coefficient is exactly the GR/Burke–Thorne target") |
| mhat_rad | back-solved to 54G/5c⁵ | same GR 2G/5c⁵ | 025 | (chain) | `…stage025_minimal_isotropic_normalization.md:101` ("… = 54 G c_s^5 / (5 a^5 c^5)") |
| K_req | 8 | back-solved so GR quad target holds | 026 | (chain) | `…stage026_concrete_axial_overlaps.md:218` ("on this branch the GR quadrupole target is equivalent to one concrete statement") |
| lambda_req / R_target / Lambda / alpha_star / Kbar_0 / T_quad / RQ_req / mathrm | back-solved to GR quad | same GR 2G/5c⁵ | 030/035/038/031/097/224/228/223 | (chain) | `…stage224…:79` ("T_{quad}:=54Gc_s^5/5a^5c^5"); `…stage030…:185` ("the target becomes a direct conservative spectral condition") |
| **m_s = μ_η** | 1836.15267343 | **CODATA 2018 m_p/m_e = 1836.152 673 43(11)** (exact to printed digits) | 250 | `bench_051666ee239b` | `…stage250…sympy_audit.md` §7.1 "Proton-proxy classical threshold speed" (corpus `fam_0072_e_max` excerpt, stage250 line 423) |
| **m_s_ratio_1836** | 1836.15267343 | same CODATA 2018 m_p/m_e | 250 | `bench_5d8d048b10f9` | same as m_s (corpus `fam_0347_chi_peak` carries m_s_ratio_1836) |

> Note: K_req/RQ_req/T_quad etc. each surface as their own corpus family but are all *carriers* of the
> single stage-022 GR anchor (`_4c_findings_ledger.md:16,59,97`). RQ_req additionally inherits the
> **illustrative** V_known barrier (see §3) — its "target" is the illustrative number, not literature.

---

## §3 (c-rec) internal recorded-benchmark matches

Back-solved or pinned to the model's **own** recorded Session-I…V outputs — calibration-to-self, not a
prediction and not literature. Disclosed in the notes (sometimes the *script* wording overclaims; flagged).

| concept (family) | posited value | matched to (internal session output) | stage | notes evidence |
|---|---|---|---|---|
| **lambda_L** | 0.26971918 | inverts the compiler so V_eff(r_soft) reproduces the **recorded Session-I softening point V_eff^sess = 1.74701126** | 247 | `…stage247…sympy_audit.md:497` ("which is fixed exactly by the recorded softening point") — script wording "independently derived" overclaims (tautological forward-check), `_4c_findings_ledger.md:81` |
| V_eff | 1.74701126 | recorded strongest stationary softening point | 247 | `…stage247…:367` ("Appendix B.2 … records the strongest stationary softening point at") |
| Vrebuild_soft | reproduces the Session-I softening point | recorded Session-I | 247 | `…stage247…:620` ("the Session-I strongest-softening point is reproduced by the concrete benchmark decomposition") |
| Xi_turn | 0.34437471 | carried Session-II hardcode (lowered-branch dynamic observable) | 248 | `…stage248…:419` ("The reported lowered-branch dynamic observables were") |
| v0_cross | above-threshold demo speed | Session-II report | 248 | `…stage248…:513` ("The Session-II report also gave an explicit above-threshold demonstration speed") |
| E_max / chi_peak / chi_raw / t_transit_min/max / g_UV | session readbacks | Session-III scan outputs | 248/250 | `…stage250…:483` ("The Session-III scan did not find an upper collapse edge…"); `:423` |
| G_W / Rmix / UVdrop_obs / Veff_obs / Wsess_obs / β_Q/β_U/β_W / η_leak / μ_w / ξ_R | Session-I one-port table | recorded Session-I table | 247 | `…stage247…:427` ("The Session-I table also records") |
| U_obs / V_obs / a_U / a_V / eps_eta_ref (=0.14313458) | representative relaxed-rigid-mouth point | Session write-up | 245 | `…stage245…:497` ("The session write-up records the representative relaxed-rigid-mouth point") |
| a_0 / b_0 / r_eval / r_sigma / a0 | Session-I source-side parameters | recorded Session-I | 246/247 | `…stage246…:472` ("Use the recorded Session-I source-side parameters") |
| gamma_crit / t_collapse_0 / t_cross | 6.94311167 (γ_crit) | **Session-IV envelope closure** (notes: NOT a microscopic theorem) | 251 | `…stage251…:562` ("is not a theorem of the minimal microscopic kernel. It is a property of the Session-IV envelope closure"); `:449` |
| K_turn / gamma_lattice_legacy / Upsilon_lat / Upsilon_lat_sess / gamma_lattice_red | 2.73855812 / 4.79562976 | Session-V reduced lattice-rate inputs (benchmark-only) | 252/253 | `notes/CHECKPOINT_CONSTANT_PROVENANCE.md:952` ("`4.79562976` is a carried forward legacy Session-V reduced lattice-rate input"); `:956` ("`…2.73855812` are benchmark-only") |
| H_int_aligned / H_sub_minus_final / hint_aligned / R_int / C_sigma | reported helicity-export outputs | Session helicity exports | 249 | `…stage249…:416` ("The reported helicity-export outputs were") |
| E_sub | 2.5 (fixed subbarrier energy) | session launch-energy convention | 248 | `…stage248…:470` ("At the fixed subbarrier energy E_sub=2.5, the launch speed is") |

> **Boundary call:** V_known ≈ 1.18190922 (222/223) is **not** here and **not** in §2 — the notes call it
> **illustrative**, not a recorded benchmark and not a published constant. Placed in §4(a) as an
> illustrative placeholder per `_4c_findings_ledger.md:55-56,102`.

---

## §4 (a) conventions / gauge — the model does NOT owe these

Genuinely arbitrary normalizations, limit illustrations, free model coefficients, gauge/coordinate
labels, and probe sample inputs. Grouped; brief. (~83 families.)

- **Normalizations / limits (no "true" value):** λ_μ = 1 (shell-weighted normalization) `…stage076…:66`;
  mhat_0 → 1 (point-particle source-map limit) `…stage195…:210` / `…stage022…:279`; ε_blk = 0
  (unblocked limit) `…stage085…:132`; λ_mu_1 (benchmark λ_μ=1) `…stage232…:149`; f'(0)=1 wall
  normalization `…stage065…:108`; mhat_0 the point-particle normalization (`fam_0469`).
- **Illustrative placeholders (explicitly "illustrative"):** V_known ≈ 1.18190922 `…stage222…:350`,
  `…stage223…:329`; ΔV_req(1) `…stage222…:436`; the stage-222 sample couplings K/M/Ω_U/Ω_W/a/c_s/lam_*
  /λ_B/λ_U/λ_W/λ_R/varpi (the `fam_0115_k` 20-param illustrative block, "The sample couplings are
  illustrative") `…stage222…:434`; λ_B alone `…stage222…:434`; chi_rm probe `…stage250…:388`.
- **Free model coefficients (a stage introduces a degree of freedom, no posited number):** K_s/K_q
  (shell/mixed stiffnesses) `…stage114…:97,98`; D_W_bare/D_0/D_A/D_01/D_a `…stage114/022/173`; B_0/B_star
  `…step_09…:50`, `…stage191…:253`; κ_0/κ_1/κ_q/M_q/M_minus/M `…stage133/134/178/228`; the wall-action
  constitutive coeffs K_η/T_w/T_Ω/μ_η ("fixed effective constitutive functions of the chosen [wall]",
  `…stage001_geometry_lift.md:214`); α/α_2/α_r/α_star-probe; Λ_mix/Λ_def/Λ_1/Lambda_mix; ε_2/ε_4/ε_β/
  ε_η/ε_blk; a_0/a0/A_0/A0_prime/A_q/A_w; V_0/V_0f; Ω_Q/Ω_A_l/OmU; Σ_source_density; T_U/T_quad-probe;
  Ξ_1/Ξ_prefactor_100; C_a/C_σ/C_tr_target; K_corr/K_geom/K_eta/K0_sym; N_n/N_Q_universal-label;
  S_EM/E_w/E_sub-energy; R_i/R_int; d_r/D_01; hat/mathrm/mathfrak_r (LaTeX-token / symbol families);
  Delta_norm/Delta_rm/Delta_def; β_path=2/beta_path; M_q_over_Sigma_m; Rratio_base; cref; KW_t; lambda_W;
  lambda_mu; Gamma_AB (open-system completion, "would enter after coupling", `…stage002…:112`).
- **Frozen-by-parent-action conventions (modeling choices upstream, not derived in-ledger):** n = 5 EOS
  exponent / n_eos / n_wall_eos / K_eos (`…phase0_spec.md`, `…stage062…:88`) — *see §5, borderline*;
  α_r = 10 wall depth; V_0 wall amplitude.
- **Sign/branch labels & gauge picks:** the gauge-invariant combos at `…stage021…:74` ("are exactly
  gauge invariant"); branch sign labels throughout.

---

## §5 borderline / needs-judgment — for human adjudication

1. **gamma_0 = (1+r_c)/9 — (b)/derived-from-(b) vs (a).** At stage 114 it is a *free* bare mixed
   low-frequency coefficient pair (κ_0, γ_0) `…stage114_concrete_core_schur.md:101`; at stage 115 the
   odd-preservation / coupling-balance condition *constrains* it (γ_c = 1/9 appears as a derived balance,
   `…stage115_core_balance_compensation.md:25`). **Competing reading:** (i) it is *derived* (an
   odd-preservation theorem of the posited core model) → not an ansatz at all; (ii) it inherits the
   posited core model (K_s,K_q,λ,r_c), so the *number* is branch-determinable → (b). The ansatz-value
   catalog lists it as the first postulated value, but the stage-115 notes present it as a derived
   condition. Recommend: **(b) derived-from-(b)** — the value is fixed once the posited core branch
   (which is itself §1-class) is fixed. Flagged in `_4c_findings_ledger.md:42`.

2. **chi_Q — placed in §1 (b) but tag is `free_choice` with values {1, 9}.** The corpus tags it
   free_choice (`fam_0345_chi_q`), and "9" appears as an alternative value, which *looks* like a free
   pick. But the 4C findings settled χ_Q = 1 as `internal_consistency`, **origin stage 105**, matched to
   the stage-104 outgoing DtN h₂⁽¹⁾ fingerprint (`_4c_findings_ledger.md:38`). I placed it in §1 (b)
   because the actual outgoing branch *determines* it. **Competing reading:** if "9" is a genuinely free
   alternate branch the user may want it noted as a branch *label* (a) rather than a value the solve
   fixes. Recommend keeping in (b); flag the {1,9} ambiguity.

3. **n = 5 frozen EOS exponent — (a) convention vs (b) parent-fixed.** Tagged free_choice; the notes say
   the parent 4D action *already* fixes the matter sector to the stiff EOS P(ρ)=K_EOS·ρ⁵
   (`…stage062_parent_action_gain.md:88`, `…phase0_spec.md`). **Competing reading:** (a) it is an upstream
   modeling convention the toy adopts (placed in §4); or (b) the "right model" (the parent GNLS) genuinely
   fixes n, making it branch/parent-determinable. Placed in §4(a) as a frozen-by-parent convention;
   adjudicate whether n is owed by the parent theory.

4. **c0 = 3/4 appears as three corpus families (c0, c_0/c_contact, c0_min) — dedup.** `fam_0333_c0`,
   `fam_0335_c_0`, `fam_0336_c_contact` (stages 088/090) and `fam_0334_c0_min` (stage 240) are the **same
   posit** (3/4) at different stages/labels. Counted once in §1 (#3/3a). Not a classification dispute —
   just confirm the user wants it as one entry, not four.

---

*End of ANSATZ_LEDGER. §1 is the falsifiable to-do list; §2/§3 are honest disclosures the audit
confirmed; §4 is what the model does not owe; §5 awaits human gate.*
