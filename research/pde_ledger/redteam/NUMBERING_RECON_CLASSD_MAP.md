# Numbering Reconciliation — Phase 2 (Class D) content-confirmed map

**Date:** 2026-06-03
**Status:** ✅ MAPPED + APPLIED + VERIFIED + COMMITTED 2026-06-03 (`fac17c7`, the 899-edit body).
The mapping (this doc) drove the application; see **§APPLICATION** and **§DEFERRED TAIL** at the bottom.
**Deferred tail:** ✅ VERIFIED + APPLIED 2026-06-03 (9 leads / 7 files, label-only — see **§DEFERRED TAIL**;
A1→118–119 & A2→183–185 corrected the original guesses by content). Only the **§OUT-OF-SCOPE broader
non-citing-file pass** remains, held for a separate decision.
**Scope:** the 121 flagged Class-D body-prose `Stage NNN` cross-refs from
`redteam/NUMBERING_RECON_SCAN.md` (§CLASS D) + the 3 special `.wl` stage227 cross-refs.

## Method (how each target was confirmed)
- Ground truth = the **filename / `\label{stage:NNN}`** number of the stage that actually
  provides the cited deliverable. Never an offset.
- For each row: read the FULL citing sentence in the citing file, extract the described
  deliverable, then OPEN the candidate target note and confirm its H1/boxed-result matches.
  Evidence = the matching line in the target note.
- The 121 rows were mapped by 7 parallel reading agents (one per citing-note cluster); the
  orchestrator independently spot-verified every medium-confidence, self-reference,
  unusual-offset, and same-stale-number-splits-by-content case against the actual files.

## Why offsets are useless here (and the proof)
Confirmed stale→canonical offsets in this set: **0, +17, +34, +51, +68, +85, +102, +103** —
mostly multiples of 17 (multiple realignment epochs), but **not reliably so**, and the *same*
stale number maps to *different* canonicals by content:
- "Stage 252" → **218** (closed the local mixed-ray search), → **201** (realization compiler),
  → **184** (three branch composites), → **235** (rigid-mouth packet split) — all in different notes.
- **stage134:70 "Stage 233 already fixed Π_*≈1.50882951349316" → 130** at offset **+103**
  (a NON-17-multiple). A sweep toward the dominant pattern would map +102 → **131**, which is
  WRONG: stage131 only "*Translate*[s] the explicit Family-1 compensation point" (L6) while
  **stage130** "gives the **unique explicit compensation point** Π_*" (L105). Content decides; the
  offset would have lied. This row is the canonical proof that offset-sweeping is unsafe.

## Confidence
All 124 rows are **HIGH** confidence (content-confirmed against the target file). Three rows
are **SELF-references** (the note cites its own deliverable with a stale number → maps to the
note's own canonical). No genuine forward references found among the flagged rows.

---

## The map (grouped by citing note; all under `notes/stages/` unless noted)
`offset = cited − canonical` shown only to expose non-uniformity; it was never used to decide.

| citing note | line | cited | **canonical** | offset | evidence (matching deliverable in target) |
| --- | --- | --- | --- | --- | --- |
| stage122_mouth_source_compensation_test | 34 | 221 | **119** | +102 | 119 boxes `1+r²=4(g−r)²` + closed form `g_comp^±(r)`; cited "compensated canonical branch requires 1+r²=4(g−r)²" |
| stage134_family1_mouth_fixedpoint | 70 | 233 | **130** | +103 | 130 §3 "gives the **unique explicit compensation point** Π_*≈1.50882951349316" (131 only *translates* it) |
| stage164_microscopic_log_channels | 5 | 248 | **163** | +85 | 163 reduces first off-family defect to single scalar `δ_⊥` |
| stage165_exact_branch_drifts | 5 | 249 | **164** | +85 | 164 = explicit microscopic log-imbalance channels |
| stage166_bundle_inversion_four_drifts | 5 | 250 | **165** | +85 | 165 = lower-branch drift laws `L_W, v_w0, T_m` |
| stage167_bundle_transport_tangent_compensation | 5 | 251 | **166** | +85 | 166 = bundle inversion of the four irreducible branch drifts |
| stage167_bundle_transport_tangent_compensation | 114 | 250 | **165** | +85 | 165 = lower-branch transport law `δln v_w0` (L116/L129) |
| stage168_off_bundle_slippage | 216 | 248 | **163** | +85 | 163 splits mouth-bias variation `δΠ=δΠ_tan+δΠ_nor` |
| stage168_off_bundle_slippage | 286 | 248 | **163** | +85 | 163 gives `δC=16σ_*/√(1+r_*²)·δ_⊥` |
| stage169_no_linear_p2_scalar_slippage | 5 | 253 | **168** | +85 | 168 = off-bundle slippage decomposition `(ε_L,ε_v,ε_T)` |
| stage171_microscopic_grouped_obstructions | 100 | 238 | **170** | +68 | 170 gives obstruction pair `(K_A,G_A)` |
| stage171_microscopic_grouped_obstructions | 476 | 238 | **170** | +68 | 170 collapses linear grouped-P2 outlet to `(K_A,G_A)` |
| stage172_physical_slope_collapse | 7 | 239 | **171** | +68 | 171 collapses the two outlet obstructions to `(𝔎_1,𝔊_1)` |
| stage172_physical_slope_collapse | 120 | 239 | **171** | +68 | 171 = "Linear Grouped Outlet Obstructions" |
| stage172_physical_slope_collapse | 265 | 238 | **170** | +68 | 170 = direct outlet map `δκ_W^(A)=…K_A` |
| stage173_axisymmetric_loading_mismatch | 5 | 240 | **172** | +68 | 172 collapse to physical slopes `u_2^(1), P_1` |
| stage173_axisymmetric_loading_mismatch | 191 | 240 | **172** | +68 | 172 hidden-even relation `u_4^(1)=8/9·u_2^(1)` |
| stage174_static_self_similarity | 5 | 241 | **173** | +68 | 173 one scalar `Ξ_load=δ_N−δ_D` |
| stage174_static_self_similarity | 367 | 241 | **173** | +68 | 173 even-preserving branch `D_21=−D_01/9` etc |
| stage175_wall_normalized_load_shape | 5 | 242 | **174** | +68 | 174 wall-referenced self-similarity fields `Σ_α^(B), Σ_r^(Z)` |
| stage175_wall_normalized_load_shape | 382 | 241 | **173** | +68 | 173 box "on the **even-preserving branch** `Δ_Q^(2m)=…Ξ_load`" (verbatim) |
| stage176_outgoing_load_factorization | 5 | 243 | **175** | +68 | 175 outgoing-load theorem `Λ_r=P_r/Δ_r` |
| stage176_outgoing_load_factorization | 178 | 243 | **175** | +68 | 175 conservative-shape outgoing-load theorem |
| stage176_outgoing_load_factorization | 274 | 243 | **175** | +68 | 175 outgoing-load theorem `2Σ_r ρ_r^(N)δlnΛ_r=δ_K` |
| stage177_weak_axisymmetric_outgoing_slippage | 149 | 244 | **176** | +68 | 176 boxed first-order defect field `Σ_r^(N)=…` |
| stage177_weak_axisymmetric_outgoing_slippage | 252 | 241 | **173** | +68 | 173 "**weak-axisymmetric** `P_1/P_0=N_01/N_0−D_01/D_0=Ξ_load`" |
| stage177_weak_axisymmetric_outgoing_slippage | 364 | 244 | **176** | +68 | 176 defines slippages `M_r, I_r, H_r` |
| stage178_outgoing_port_coloading_theorem | 5 | 245 | **177** | +68 | 177 collapse to one scalar `Ξ_1=P_1/P_0` |
| stage178_outgoing_port_coloading_theorem | 103 | 243 | **175** | +68 | 175 outgoing-load factor on conservative-shape branches |
| stage178_outgoing_port_coloading_theorem | 410 | 245 | **177** | +68 | 177 collapse to one scalar `Ξ_1` |
| stage179_transfer_shape_theorem | 5 | 246 | **178** | +68 | 178 `Ξ_1=ν̄_N−κ_1` portwise formula |
| stage179_transfer_shape_theorem | 238 | 246 | **178** | +68 | 178 `Ξ_1=Σ_r ρ_r^(N)(ν_r−κ_1)` |
| stage179_transfer_shape_theorem | 425 | 246 | **178** | +68 | 178 reduces defect to `Ξ_1=ν̄_N−κ_1` |
| stage180_effective_transfer_shape_collapse | 5 | 247 | **179** | +68 | 179 `Ξ_1=2Σ_r ρ_r^(N)τ_r` transfer-shape theorem |
| stage180_effective_transfer_shape_collapse | 61 | 247 | **179** | +68 | 179 boxed portwise factorization `N_A0^(r)=K_A·T_{A,r}²` |
| stage180_effective_transfer_shape_collapse | 362 | 241 | **172** | +69 | 172 "on the **weak axisymmetric branch** `Δ_Q^(2m)=…P_1/P_0`" lane pattern (≠ 173's even-preserving box) |
| stage180_effective_transfer_shape_collapse | 388 | 247 | **179** | +68 | 179 outgoing-weighted transfer-shape slope |
| stage181_coherent_tracking_defect | 6 | 248 | **180** | +68 | 180 effective transfer-shape collapse `T_eff,A²` |
| stage182_microscopic_coherent_slippage | 5 | 249 | **181** | +68 | 181 boxed `Ξ_1` on coherent tracking branch |
| stage182_microscopic_coherent_slippage | 102 | 249 | **181** | +68 | 181 support-blindness, `ζ` drops out (`∂_ζ Ξ_1=0`) |
| stage182_microscopic_coherent_slippage | 216 | 249 | **181** | +68 | 181 Support-Blindness Theorem |
| stage182_microscopic_coherent_slippage | 347 | 249 | **181** | +68 | 181 tracking factor `R_tr` |
| stage183_triangular_normal_form | 6 | 250 | **182** | +68 | 182 four-slippage + one dressing-slippage reduction |
| stage183_triangular_normal_form | 85 | 250 | **182** | +68 | 182 tracking/nontracking split |
| stage183_triangular_normal_form | 142 | 250 | **182** | +68 | 182 `Θ_1` formula (verbatim) |
| stage184_branch_invariant_coordinates | 188 | 251 | **183** | +68 | 183 inverse reconstruction `Σ_tr ∝ Θ_1` |
| stage184_branch_invariant_coordinates | 334 | 251 | **183** | +68 | 183 rigidity equivalence `Σ_tr=Σ_nt=Σ_η=0` |
| stage185_microscopic_monomials | 6 | 252 | **184** | +68 | 184 three branch-composites reduction |
| stage185_microscopic_monomials | 92 | 250 | **182** | +68 | 182 microscopic drifts `Σ_χ, Σ_δ` |
| stage185_microscopic_monomials | 322 | 251 | **183** | +68 | 183 triangular observable normal form |
| stage185_microscopic_monomials | 462 | 252 | **184** | +68 | 184 three branch composites |
| stage186_similarity_orbit_closure | 6 | 253 | **185** | +68 | 185 three-monomials reduction |
| stage186_similarity_orbit_closure | 60 | 253 | **185** | +68 | 185 `δln 𝔠_tr,*=Σ_tr` (verbatim) |
| stage186_similarity_orbit_closure | 351 | 251 | **183** | +68 | 183 rigidity equivalence |
| stage187_orbit_quotient_closure | 5 | 237 | **186** | +51 | 186 five-parameter similarity-orbit tangent problem |
| stage187_orbit_quotient_closure | 108 | 237 | **186** | +51 | 186 five-parameter family `𝒢_*` |
| stage188_branch_observables_completion | 216 | 238 | **187** | +51 | 187 finite tangent-quotient packet |
| stage188_branch_observables_completion | 457 | 238 | **187** | +51 | 187 exact finite quotient problem |
| stage189_transfer_shape_prefactor_compiler | 104 | 239 | **188** | +51 | 188 branch-observable completion `𝔑_*`; locked by 189's self-label (L96 "Stage 240"=189 ⇒ "239"=188) |
| stage192_orbit_quotient_projectors | 108 | 242 | **191** | +51 | 191 finite Packet B orbit-lock |
| stage195_source_map_reduction_of_canonical_outgoing_branch | 248 | 245 | **194** | +51 | 194 isotropic DtN deformation algebra |
| stage196_higher_odd_irrelevance_theorem | 284 | 246 | **195** | +51 | 195 factorization `m²χ_Q N_Q=1` |
| stage197_conditional_packetA_closure_theorem | 30 | 246 | **195** | +51 | 195 factorization `m²χ_Q N_Q=1` |
| stage197_conditional_packetA_closure_theorem | 227 | 245 | **194** | +51 | 194 isotropic DtN deformation algebra |
| stage197_conditional_packetA_closure_theorem | 327 | 247 | **196** | +51 | 196 higher-odd irrelevance theorem |
| stage200_reference_free_home_stretch_theorem | 411 | 242 | **191** | +51 | 191 reduced endgame to two packets |
| stage202_free_quintuple_target_graph | 14 | 252 | **201** | +51 | 201 realization compiler `x∈Z_* ⟺ …` |
| stage202_free_quintuple_target_graph | 497 | 252 | **201** | +51 | 201 canonical orbit projection |
| stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem | 533 | 253 | **202** | +51 | 202 free-quintuple target graph over 5 free coords |
| stage205_directional_hessian_and_quadratic_root_refinement | 14 | 238 | **204** | +34 | 204 scalar root predictor `Φ_s(τ)=1` |
| stage205_directional_hessian_and_quadratic_root_refinement | 70 | 238 | **204** | +34 | 204 finite graph invariance (on orbit `O_*` ∀τ) |
| stage205_directional_hessian_and_quadratic_root_refinement | 395 | 238 | **204** | +34 | 204 log-ray scalar-root problem |
| stage206_certified_ray_ranking_and_local_bracketing_theorem | 56 | 239 | **205** | +34 | 205 `h_s(τ):=ln Φ_s(τ)` |
| stage206_certified_ray_ranking_and_local_bracketing_theorem | 301 | 239 | **205** | +34 | 205 turning-point closure test |
| stage207_primitive_ray_hessian_envelopes_and_certified_ray_table | 14 | 238 | **204** | +34 | 204 primitive free-direction table `e_λ,e_c,…` |
| stage207_primitive_ray_hessian_envelopes_and_certified_ray_table | 206 | 238 | **204** | +34 | 204 dependent graph exponents |
| stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy | 55 | 241 | **207** | +34 | 207 certified ray-table bracket data `H_0, Γ_i` |
| stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy | 469 | 241 | **207** | +34 | 207 certified ray table + off-diagonal elimination |
| stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem | 453 | 242 | **208** | +34 | 208 pairwise mixed rays, two-ray audit |
| stage210_three_coordinate_mixed_simplex_audit | 14 | 243 | **209** | +34 | 209 pairwise ratio optimizer (pairwise problem closed) |
| stage210_three_coordinate_mixed_simplex_audit | 550 | 243 | **209** | +34 | 209 pairwise boundary problem finite/exact |
| stage213_..._support_cardinality_4_gate | 14 | 246 | **212** | +34 | 212 full primitive-triple closed simplex (support-≤3) |
| stage213_..._support_cardinality_4_gate | 136 | 246 | **212** | +34 | 212 per-triple closed simplex |
| stage213_..._support_cardinality_4_gate | 155 | 246 | **212** | +34 | 212 per-triple closed simplex |
| stage213_..._support_cardinality_4_gate | 174 | 246 | **212** | +34 | 212 per-triple closed simplex |
| stage213_..._support_cardinality_4_gate | 193 | 246 | **212** | +34 | 212 per-triple closed simplex |
| stage213_..._support_cardinality_4_gate | 688 | 246 | **212** | +34 | 212 support-≤3 finite/exact (budget 600) |
| stage214_..._finite_candidate_reduction | 561 | 247 | **213** | +34 | 213 quadruple-screen audit (4 faces + 2 interior) |
| stage214_..._finite_candidate_reduction | 584 | 247 | **213** | +34 | 213 support-cardinality-4 gate |
| stage215_full_primitive_quadruple_ranking_theorem | 416 | 248 | **214** | +34 | 214 quadruple interior candidate set (54) |
| stage215_full_primitive_quadruple_ranking_theorem | 433 | 246 | **212** | +34 | 212 support-≤3 budget 600 |
| stage216_..._support_cardinality_5_gate | 14 | 249 | **215** | +34 | 215 support-≤4 closed-simplex ledger |
| stage216_..._support_cardinality_5_gate | 91 | 249 | **215** | +34 | 215 imported closed-simplex intervals |
| stage216_..._support_cardinality_5_gate | 327 | 250 | **216 (SELF)** | 0 | own deferral gate ("…deferred to Stage 251"=217 forward; "Stage 250 keeps…"=this note) |
| stage216_..._support_cardinality_5_gate | 494 | 249 | **215** | +34 | 215 support-≤4 ledger |
| stage217_..._finite_candidate_reduction | 95 | 249 | **215** | +34 | 215 imported closed-simplex intervals |
| stage218_..._local_mixed_ray_search_closure | 14 | 249 | **215** | +34 | 215 global support-≤4 search |
| stage218_..._local_mixed_ray_search_closure | 287 | 251 | **217** | +34 | 217 support-5 interior finite candidate set |
| stage218_..._local_mixed_ray_search_closure | 367 | 251 | **217** | +34 | 217 unique 5-coord interior optimizer candidate set |
| stage219_..._square_law_suppression_test_sympy_audit | 14 | 252 | **218** | +34 | 218 final local mixed-ray closure |
| stage219_..._square_law_suppression_test_sympy_audit | 22 | 252 | **218** | +34 | 218 local mixed-ray search sieve closed |
| stage219_..._square_law_suppression_test_sympy_audit | 41 | 253 | **219 (SELF)** | 0 | own square-law suppression verdict (this note's output #6) |
| stage221_..._linear_survival_window_sympy_audit | 518 | 252 | **218** | +34 | 218 final local mixed-ray closure |
| stage222_..._residue_linewidth_survival_test_sympy_audit | 15 | 238 | **221** | +17 | 221 reduces linear dynamic survival to a residue/linewidth inequality |
| stage222_..._residue_linewidth_survival_test_sympy_audit | 36 | 239 | **222 (SELF)** | 0 | own outputs-list conclusion (static/dynamic tension) |
| stage222_..._residue_linewidth_survival_test_sympy_audit | 48 | 238 | **221** | +17 | 221 only linear corridor left = resonant dispersive |
| stage227_..._first_co_loading_no_go_sympy_audit | 22 | 243 | **226** | +17 | 226 pure-transfer subcorridor `D01=D21=D41=0` |
| stage228_..._first_actual_dynamic_window_test_sympy_audit | 24 | 244 | **227** | +17 | 227 load-factor corridor `Ξ_1=2δlnΛ, Λ=P/Δ` |
| stage228_..._first_actual_dynamic_window_test_sympy_audit | 100 | 244 | **227** | +17 | 227 load-factor pure-transfer corridor |
| stage229_..._softening_depth_crossover_theorem_sympy_audit | 21 | 245 | **228** | +17 | 228 `Ξ_1=2(π_1−δ_1)` + numerator/denominator-rigid subcorridors |
| stage231_..._dynamic_class_map_sympy_audit | 165 | 247 | **230** | +17 | 230 classifier-to-dynamic-window compiler, two global statements |
| stage231_..._dynamic_class_map_sympy_audit | 529 | 247 | **230** | +17 | 230 boxed `inf B_dyn>B_stat` (verbatim) |
| stage233_..._static_turbulence_gate_sympy_audit | 19 | 249 | **232** | +17 | 232 support/source numerically safe; bottleneck = orbit-lock placement |
| stage233_..._static_turbulence_gate_sympy_audit | 82 | 239 | **188** | +51 | 188 generic first-order observable compiler `(Θ_1,Ξ_1,ℛ_1)`, `𝔑_*` (≠ 234's specialized one) |
| stage233_..._static_turbulence_gate_sympy_audit | 183 | 241 | **224** | +17 | 224 weak-axisymmetric ceiling transport `P̄_0(1+|εΞ_1|)≤P_crit` |
| stage233_..._static_turbulence_gate_sympy_audit | 303 | 249 | **232** | +17 | 232 support/source overshoot factor ~7.4, not the bottleneck |
| stage235_..._codimension_two_orbit_lock_point_sympy_audit | 25 | 251 | **234** | +17 | 234 direct branch-observable compiler `(Θ_1,Ξ_1,ℛ_1,c_η)` |
| stage235_..._codimension_two_orbit_lock_point_sympy_audit | 356 | 251 | **234** | +17 | 234 strip `|εΞ_1|≤B` |
| stage236_..._static_only_restoration_gap_sympy_audit | 442 | 252 | **235** | +17 | 235 static strip codim-1 inside codim-2 orbit-lock |
| stage239_..._cartesian_orbit_lock_theorem_sympy_audit | 40 | 253 | **236** | +17 | 236 microscopic dependent-plane carriers of `q_η` |
| stage239_..._cartesian_orbit_lock_theorem_sympy_audit | 61 | 253 | **236** | +17 | 236 quotient-failure image `y_rm` on `Δ_T=0` plane |

### Special — `.wl` cross-refs (Class-D in `mathematica/`; flagged by P4-53)
| file | line | cited | **canonical** | offset | evidence |
| --- | --- | --- | --- | --- | --- |
| mathematica/moving_throat_pde_stage227_..._mathematica_audit.wl | 154 | Stage-242 | **225** | +17 | 225 = "First-Order Conservative Compensation Surface"; computes the compatibility wall stiffness (kCompat). Consistent w/ stage227's note-frame `243→226` |
| mathematica/moving_throat_pde_stage227_..._mathematica_audit.wl | 243 | Stage-243 | **226** | +17 | 226 boxes raw basis `t₁≈(−4.359…)`, `t₂≈(1.909…)` (the "transfer basis" vectors) |
| mathematica/moving_throat_pde_stage227_..._mathematica_audit.wl | 244 | Stage-243 | **226** | +17 | same as :243 |

---

## Self-references (3) — map to the note's OWN canonical, not upstream
These are stale numbers the note uses for ITS OWN deliverable. On application, rewrite to the
note's own canonical number (it's referring to itself):
- `stage216 :327` "Stage 250 keeps only the first exact gate…" → **216** (adjacent "deferred to
  Stage 251" is a forward ref to **217**, the interior optimizer — out of the 121, noted).
- `stage219 :41` "So Stage 253 keeps the same-charge mixed corridor alive…" → **219** (restates
  this note's own square-law suppression verdict).
- `stage222 :36` "So Stage 239 keeps the dynamic same-charge idea alive…" → **222** (closing line
  of this note's own outputs list).

## Adjacent residuals discovered while mapping (NOT in the 121; for the eventual cleanup)
The reading agents surfaced additional stale numbers in the SAME files that the scan's
high-signal filter did not flag (mostly stale **self-titles in `## N. …` section headers** and
"What Stage N changes" transition lines that name the note's OWN or the NEXT stage). Examples:
stage183 §5 "after Stage 251"(=183); stage185 §7 "What Stage 253 changes"(=185); stage188 §7
"after Stage 239"(=188); stage189 L96 "So Stage 240 is the exact compiler"(=189) + "Stage 239
direct branch observables"(=188); stage191 L311 "Stage 242 home-stretch theorem"(=191); stage192
L57 "So Stage 243 turns…"(=192); stage205 §9 / stage206 §5 / stage208 §8 / stage209 §11 /
stage210 §9 / stage216 §9 self-titles; stage203 L~534 "Stage 237 goes one step further"(=203);
stage236 L443 "Stage 253 now says…"(=236); stage239 L31 already-correct "Stage 238"→238 (offset 0).
These are real but out of scope for the 121-row Class-D pass — recommend folding into the same
content-keyed application sweep (each maps to its own/next canonical, content-confirmed).

## §APPLICATION (2026-06-03 — DONE)
Applied by 7 fresh cluster-agents (one per citing-note cluster), orchestrator-applied per WHO-APPLIES
(doc-only label edits, NOT Codex). Each agent applied its confirmed map rows + re-scanned its files for
**adjacent residuals** (repeated upstream cites, stale `## k.` section self-titles, "what remains for
Stage N" next-stage transitions) and content-mapped+fixed those too.

- **899 label-token line-edits across 58 files** (57 citing notes + the stage227 `.wl`); git diff = 899/899,
  perfectly balanced (no lines added/removed).
- Composition: the **124 content-confirmed rows** above (121 Class-D body + 3 `.wl`, incl. the 3 self-refs)
  + **~775 adjacent residuals** (same-target repeats, self-titles, next-stage transitions), each content-confirmed.
- **VERIFICATION (orchestrator, mechanical + spot):**
  - **0 label-only violations across all 899 line-pairs.** A normalizer collapsed every
    `Stage`/`STAGE`/`Stages N–M`/`post-Stage`/compound `A/B`/`X and Y` ref to a placeholder and confirmed
    each changed line is otherwise **byte-identical** old vs new → no math/value/prose/section-number byte
    changed anywhere. (`redteam/tmp_prompts/verify_labelonly.py`, gitignored.)
  - New targets all in 1–253. Offset histogram exposed the spread (−17,+17,+34,+51,+68,+69,+85,+102,+103);
    **every non-17-multiple / unusual class was content-checked** and confirmed: the lone `−17`
    (stage166 "60"→**77**, stage077 = "first concrete branch value of Θ_w ≈4.06863…"; stage060 has none),
    `+69` (stage180→172 lane pattern), `+102` (parent-core 118/119/121), `+103` (stage134 Π_*→130).
  - `.wl` confirmed landed (Stage-225 / Stage-226); sampled map rows confirmed (219:22→218, 134:70→130, 180:362→172);
    **no `Stage 24x/25x` stale tokens remain** in any citing file.

## §DEFERRED TAIL — ✅ VERIFIED + APPLIED 2026-06-03 (next session, this commit)
All 9 ambiguous leads were content-verified against the actual target files (read-only verification agent +
orchestrator review), then applied label-only (number-token swaps only; 9 lines / 7 files; git diff byte-identical
otherwise; en-dashes/hyphen preserved). **Two leads were CORRECTED by content** (the guesses below were wrong):
A1 resolved to **118–119** (not the guessed ~114–118), A2 to **183–185** (not the guessed 184–186). Every map
relied on the cited deliverable matching the target note's boxed result — never an offset sweep.

**A. Range refs (2) — ✅ APPLIED:**
- `stage164:28` — "…throat-core closure already fixed at **Stages 169–170**" → **118–119** (+51). CONFIRMED: stage164's
  OWN body cites this same closure canonically — `stage118:256` "the first explicit throat-core closure" (boxes K_s,K_q,λ,g_s,g_q),
  `stage164:147` "the same explicit overlap closure already used at Stage 118", `stage164:43` "The Stage 119 parent ratios are…".
  So "169–170" was an old-epoch duplicate of the 118/119 the note already names. (Lead's "~114–118" was wrong; offset+content agree at +51.)
- `stage187:364` — "**Stages 234–236** already identified the first-order observable bridge from the quotient coordinates" →
  **183–185** (+51, NOT the guessed 184–186). CONFIRMED equation-by-equation: `stage183:133–135` defines (Σ_tr,Σ_nt,Σ_η)
  (triangular normal form); `stage184:48–52` δln𝔗_*=Σ_tr / δln𝔑_*=Σ_nt / δlnε_η=Σ_η; `stage185:24–30` the **verbatim 𝔠-notation**
  δln𝔠_{tr,*}=Σ_tr,… (stage185's MAIN RESULT). stage186 only *carries* it ("Direct microscopic monomials from Stage 185") → bridge ends at 185.

**B. Compound-secondary numbers (5 lines; primary kept, trailing stale at the +51 epoch) — ✅ APPLIED:**
- `stage178:155` "Stage 173/**228**" → **177** (+51; stage177 collapses to one scalar Ξ_1=P_1/P_0, identified w/ the physical outgoing-prefactor slope).
- `stage179:266` "Stage 176/**228**/**229**" → **177**/**178** (+51; the three slippage languages 176/177/178).
- `stage179:286` "Stage 177/**229**" → **178** (+51; stage178 microscopic slopes 𝔪_r, Ξ_1=ν̄_N−κ_1).
- `stage179:319` "Stage 176/**228**/**229**" → **177**/**178** (+51).
- `stage182:313` "Stage 180/**232**" → **181** (+51; stage181 owns the selected-branch identity Ξ_1=−η_1/(1−ε_η)−ℛ_1).

**C. Low-band stale (2; invariant-stem old-number evidence) — ✅ APPLIED:**
- `stage171:69` "**Stage 6** … full grouped real P2 bundle" → **023** (`full_grouped_bundle`; stage023 H1 = "Full Grouped Real `P2` Bundle"
  verbatim; corroborated — stage024's own purpose uses the same stale "Stage 6" for canonical 023).
- `stage173:21` "**Stage-7** axisymmetric grouped signature (λ_2m=1,½,−1)" → **024** (`overlap_isotropy`; stage024 boxes
  `λ_(20)=1, λ_(21)=1/2, λ_(22)=-1` verbatim as its "First Axisymmetric Splitting Law"). Hyphen preserved → `Stage-024`.

**D. Reviewed, already-correct, left untouched (NO action):** `stage180:141/252` "Stage 21"(=021),
`stage182:83/138` "Stage 30"(=030), `stage166:145` "Stage-022/6"(=022/026), `stage229:5` "Stage 143/093",
`stage221` "Stage 220"(=220), `stage239` "Stage 238/237"(=238/237), `stage236:478` "Stage 237"(forward, =237).

## §OUT-OF-SCOPE discovery (broader deferred pool — separate future pass)
Stale `Stage NNN` refs also live in files **never in the Class-D citing set** (so untouched by this pass):
e.g. `stage121` body "before Stage 223" (→121), `stage193` `\mathcal S_{244}` subscript. These belong to the
scan's deferred ~1,200-mention project-wide pool. Recommend a separate pass over **all ~243 notes**
(content/stem-keyed, never-sweep), not just the 57 citing files.
