# ledger_stage017 — grouped-P2 lane isotropy (Check II-G3b)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 sector). The COMPLETING leg of the 2-way split of `pathA_32`: this stage
carries the **grouped-P2 lane-isotropy component (2/2) of the joint `ISOTROPY_CALIBRATED`** and **CONSUMES** stage 016's SO(3)
covariance theorem to **land the JOINT** and **export the ℓ=2 port kernel**. (Stage 016 = the EARNED-FIRST angular covariance theorem.)

**Verdict.** `ISOTROPY_CALIBRATED` (JOINT, 2-stage) — landed here as **COMPLETE** (016 EARNED-cited ∧ 017 EARNED-here). Ledger
earned-label: `GROUPED_P2_LANE_ISOTROPY_EARNED`.

**Status.** Mixed exact-symbolic + numeric-sample: the lane algebra, D-lanes, raw-D/normalized-u defects, and `cs_equal` degeneracy are
exact symbolic (`compact(expr)==0` residuals, `expect_zero`/`expect_bool`); the response-probe "moves under ε" tests + the
stability/denominator windows are numeric over a calibration sample/window (`f_eval` floats, `PROBE_TOL`/`SYMBOLIC_TOL`/`NUMERIC_TOL`
guards — as the source). Dual-engine SymPy **118 PASS** / Mathematica **127 PASS**, exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_32_grouped_p2_isotropy_{sympy.py,.wl}` +
> `reports/pathA_32_grouped_p2_isotropy.md` (the 017 slice = report :3, :5, :12–14, :18–25, :28–40, :58–63). The report is cited for
> provenance only; the derivation below is self-contained.

---

## 1. What this stage earns

The ℓ=2 wall response, assembled into the three grouped `{20,21,22}` lanes from stage 016's consumed covariance, responds
isotropically — and the response is genuinely able-to-fail. This is what makes the joint verdict `ISOTROPY_CALIBRATED` (the earned
structural content) rather than a bare assertion.

### 1.1 The grouped-lane assembly (rides 016's consumed covariance)
Each channel is assembled from the consumed eigenvalue `λ_m = 6` and the K₂ angular-stiffness form:
```
M₂ = c_self·M̃,   K₂ = c_self·(K̃ + λ_m·T̃_Ω),   with c_self = ∫Y_m² dΩ = 1 (the Gram diagonal, cited from 016);
D-lanes:  D0 = K₂ − (B̃0+Z̃0),   D2 = −(M̃ + B̃2+Z̃2),   D4 = −(B̃4+Z̃4).
```
The three grouped lanes are `20`, `21 = avg(21c,21s)`, `22 = avg(22c,22s)`. This grouped `{M₂, K₂, D0/D2/D4, B̃/Z̃}` **is the exported
ℓ=2 port kernel**. `λ_m = 6` and the K₂ form are CITED from stage 016 (not re-derived) via a genuine cross-stage dual-site (§4).

### 1.2 Raw-D defects = 0 — the PRIMARY isotropy test
With the P₂-projection `defect_pair(triple) = ( (2·D20 − D21 − D22)/10, (D21 − D22)/2 )`, the grouped D-lane defects vanish at every
order:
```
{a_D0, b_D0, a_D2, b_D2, a_D4, b_D4} = 0   (exact),   and the within-group c/s degeneracy D21c=D21s, D22c=D22s holds.
```
This is the isotropy statement: the grouped ℓ=2 response is angularly degenerate because the lanes ride 016's single eigenvalue
`λ_m = 6` (the covariance). `raw_defects_zero ∧ all(cs_equal)` is the `lane_collapse_ok` gate.

### 1.3 Normalized-u defects = 0 — the CROSS-CHECK
The DtN-observable normalized responses `u₂ = −D2/D0`, `u₄ = (D2² − D0·D4)/D0²` per lane also have vanishing P₂-projection defects:
`{a2, b2, a4, b4} = 0` (exact). ⚠ This is a CROSS-CHECK, **not** the primary test: a pure-prefactor anisotropy MOVES the raw-D
defects but leaves the normalized-u defects zero (the `pure_prefactor` probe, §2) — so normalized-u alone would MISS a real
anisotropy. Both are kept as distinct asserts; **raw-D is decisive**.

### 1.4 The joint landing (COMPLETES `ISOTROPY_CALIBRATED`)
The full `verdict_from_gates` runs all eight gates: `dimensional_ok`, `covariant`, `tautology_clear` (CITED from 016);
`dynamic_retained` (`D2` carries `M̃`, `D0` does not), `stability_ok` (`M̃>0 ∧ K₂>0` over the frozen window), `denominator_guard_ok`
(`D_A0>0`), `lane_collapse_ok` (raw-D=0 ∧ cs_equal), `able_to_fail_ok` (017's 6-probe aggregate ∧ 016's cited 3-probe aggregate — the
joint 9-probe battery reconstituted across the two stages). With `calibration_inputs` non-empty, the verdict lands
`ISOTROPY_CALIBRATED` (COMPLETE = 016 ∧ 017).

---

## 2. The able-to-fail battery (017-owned)

017's aggregate battery is its six response probes; each `computed_probe_gate_flag` is a COMPUTED conjunction, and **neutering any one
flips `able_to_fail_ok` false**.
| probe | mutation → verdict | integrity (load-bearing) |
|---|---|---|
| `pure_prefactor_anisotropy` | `w=1+ε·P₂` uniform → `FAIL_ANISOTROPIC_BRANCH` | `raw_D_moves` (numeric) ∧ `normalized_u_defects_stay_zero` (exact) — the raw-vs-normalized discriminator |
| `sector_selective_anisotropy` | `w=1+ε·P₂` on B̃/Z̃ only → `FAIL_ANISOTROPIC_BRANCH` | `raw_D_moves ∧ u_defects_move` |
| `m_dependent_profile` | `β_22=(1+δ)β_20` → `FAIL_ANISOTROPIC_BRANCH` | `raw_D_moves` |
| `degenerate_beta_zero` | `β₂=0 → M₂=K₂=0` → `FAIL_STABILITY` | `computed_fail_gate ∧ self_ablation.fail_suppressed` |
| `singular_denominator` | `B0+Z0=K₂ → D_A0=0` → `FAIL_SINGULAR_RESPONSE` | `computed_fail_gate ∧ self_ablation.fail_suppressed` |
| `static_drop_inertia` | drop `M̃` from `D_A2` → `FAIL_STATIC_RESPONSE` | `computed_fail_gate ∧ self_ablation.fail_suppressed` |

⚠ **The three anisotropy probes' `verdict` is FORCED via `case_verdict(lane_collapse_ok=False)`**, so their `expected_probe_verdicts_match`
is vacuous — their load-bearing integrity rides ONLY the computed move-flags (`raw_D_moves`/`normalized_u_defects_stay_zero`/
`u_defects_move`), not the forced verdict, and no narrative field (e.g. `normalized_u_defects_may_cancel_for_pure_profile_scale`) is a
gate flag. The three 016-owned probes (`wrong_eigenvalue`/`tautology_hash_collision`/`dimensional_corrupt_T_Omega`) are cited via 016's
3-probe aggregate, not re-run here.

---

## 3. Honest scope

- **CALIBRATED not PASS (the calibration partition).** 017 EARNS the grouped-lane isotropy (raw-D/normalized-u genuinely vanish, riding
  016's covariance), but the wall profile + the radial/support scalars — `β₂(w)`, `R0(w)`, `M̃/K̃/T̃_Ω`, `B̃/Z̃`, the calibration window —
  are FROZEN calibration inputs, **not derived from the Gate-1 `R0` support equation** (report :5). That frozen-ness is WHY the joint is
  `ISOTROPY_CALIBRATED` rather than `ISOTROPY_PASS` (the source verdict ladder makes `PASS` require `β₂` derived from `R0` AND `T_Ω`
  derived). Edge **R36** records that debt (the `ISOTROPY_PASS` target).
- **Deferred (Gate 4/5/6, sim-deferred).** The 54/5 quadrupole normalization, the outgoing odd-N coefficients, and the solved nonlinear
  branch data (`G = GENUINE_BLOCKED`) are downstream work, not 017's.

---

## 4. Consumed / exported / register

- **Consumed — 016's covariance theorem via a GENUINE cross-stage dual-site.** ⭐ Because 016↔017 share the SAME pathA_32 convention
  (VOLUME densities / dimensionless `β₂`) — unlike 016's provenance-only cite of 013 (line-vs-volume) — the cited `λ_m=6` + K₂ form are a
  CHECKABLE relation. The dual-site: **Site A** = the `λ_m` fed into the lane assembly equals 016's exported eigenvalue; **Site B** = the
  ASSEMBLED lane K₂'s `T̃_Ω`-coefficient equals the cited `λ_m` (guards the assembly formula — reads the assembled `ungrouped[name]["K2"]`,
  the fix from the adversarial finding). ⚠ **The λ-dimensionless trap:** a corruption `λ:6→4` is dimensionally silent, leaves the
  stability/denominator windows positive, AND leaves raw-D isotropy incidentally unbroken (all-lanes λ=4 keeps raw-D=0) — so an incidental
  gate cannot catch it. A one-site λ/K₂-form corruption fires an EXPLICIT integrity guard; the coordinated-both-sites corruption is caught
  by a single-harmonic `(−Δ_S²)Y₂₀ = λ·Y₂₀` echo (a genuine harmonic-shaped residual). 017 does NOT re-derive the covariance theorem (a
  single-`Y20` echo only). Also consumed as PROVENANCE (not checkable relations): `μ_η`, `T_w` (counted at 013; `K_η=T_w β²` R29
  non-transferable across the convention gap), `R0(w)`, the Gate-1 D/N boundary provenance (012 R28, 011 R26). `c_S` is NOT consumed.
- **Register.** ⭐ **Two new counted CALIB knobs `{T_Ω, β₂}`** — the THIRD calibration-adding Part-II stage. `T_Ω` (the ℓ=2 angular/
  rotational stiffness density, a genuinely new physical DOF absent from the ℓ=0 breathing packet) and `β₂(w)` (the frozen ℓ=2 radial
  profile) are **calibration inputs not derived** from the Gate-1 `R0` support equation (edge R36). The radial scalars `M̃/K̃/T̃_Ω` are
  DERIVED manifestations `= ∫density·β₂²` (edge **R35**, the ℓ=2 analogue of R29) — not counted afresh. The ℓ=2 port-kernel support/Maxwell
  scalars `{B̃0/2/4, Z̃0/2/4}` are frozen CALIB inputs first load-bearing here but **tracked/downstream-pinned** (the isotropy verdict is
  value-independent of them; the `Z̃` Maxwell couplings pin at 018–024 + the EM knit/Gate-4 — irreducibility a Part-VII adjudication). New
  edges **R35** (radial reduction, `DERIVED`) + **R36** (the `R0`-support debt, `PENDING`); **R34** backfilled into the edges table (016's
  covariance provenance, was rollup-prose only). **Part-II counted CALIB set → `{μ_η,T_w,β}`(013) + `{Vp0/ℓ_c}`(015) + `{T_Ω,β₂}`(017) = 6.**
- **Exported — the ℓ=2 PORT KERNEL (017's, not 016's).** The grouped `M₂`, the angular `K₂ = K̃ + 6·T̃_Ω`, the support scalars `B̃/Z̃`, and
  the D-lanes (`D0/D2/D4`) → stages **018–021** (pathA_33 quadrupole wall mode) + **022/023** (pathA_34 cross-ℓ ℓ=2 map) + **024** (pathA_43
  density-port wall mode).

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a genuinely independent NATIVE route (native
`Integrate`/`FullSimplify`/`D`/`N`, its own `verdictFromGates`/`caseVerdict`/`dimOf`) — no `Get`/`Import`/`Export`/YAML, no mirroring; the
source pair's scratch-YAML engine-agreement handoff is severed. Agreement is transcript-level: both engines emit the grouped lanes,
raw-D=0, normalized-u=0, the six probe verdicts + self-ablations, the six-probe aggregate, and the joint `ISOTROPY_CALIBRATED`. Arity
self-check + unevaluated-leakage scan carried.

**Directive review** used the Codex→Grok→Codex bookend: Codex `DIRECTIVE_CLEAN` first-pass; the Grok compute-verify pass returned 1
BLOCKING (the §5 EXCLUDE list wrongly stripped the harmonics/`intS2` infrastructure 017 needs for its anisotropy coefficients + the echo)
+ 2 nits (the export coefficient pinned to `K̃+6·T̃_Ω`; the three anisotropy probes' forced-verdict flagged vacuous so the load rides the
move-flags), all folded; a final Codex confirm → `DIRECTIVE_CLEAN`. Grok independently compute-confirmed the core math (raw-D=0, the
`defect_pair` projection, the pure_prefactor discriminator, and — decisively — that all-lanes λ=4 keeps raw-D=0, so the explicit
dual-site is provably necessary).

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (independent re-derivation of the anisotropy coefficients `{2/7,1/7,−2/7}`, the Y20 echo
residual, the pure_prefactor discriminator, the six probe verdicts, and the joint landing — no dropped check, no transliteration error,
no hardcoded-should-be-computed value) + `ADVERSARIAL_ISSUES` (42 per-tooth ablations, 41 firing at their own assert) with ONE genuine
finding: the `.wl` Site-B (K₂-form) residual reconstructed K₂ fresh from `lambdaByChannel` (algebraically identical to Site A) instead of
reading the assembled lane, so an assembly-formula K₂ corruption fired only downstream, not at the named `dual_site.K2_form` guard.
**Remediated** (Codex): the `.wl` Site-B now reads the assembled lane `D[laneAssoc[#]["K2"], TomegaTilde] − cSelf·lambdaRef` (after
`ungrouped` is assembled), mirroring the SymPy sibling, plus matching assembly-formula ablation teeth in both engines. A fresh-agent
`REVERIFY_CLEAN` confirmed the coupling meta-test: the fixed tooth fires at its own named assert on an assembly-formula corruption in
both engines (residual 1), neutering the fix makes it go vacuous (proving it load-bearing), and no regression (SymPy 118/0, WL 127/0,
the one-site λ / coordinated-echo / six-probe-neuter behaviors all preserved).
