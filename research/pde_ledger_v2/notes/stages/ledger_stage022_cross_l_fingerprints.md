# ledger_stage022 — the cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the Gate-4 non-regression (Check II-G5a)

**Part / anchor.** Part II — Gravity (the frozen-throat cross-ℓ radiative sector). The EARNED-FIRST leg of a 2-way split of
`pathA_34`: this stage carries the **cross-ℓ outgoing DtN `−(ℓ+1)/Λ_ℓ` fingerprints (ℓ=0,1,2) + the Gate-4 quadrupole
non-regression component (1/2) of the joint `FAIL_UNDERDETERMINED_NOT_PREDICTIVE`** — the EARNED exterior-wave signature that
generalizes stage 018's ℓ=2-only fingerprint across ℓ and checks it does not regress the completed quadrupole. The native
nullspace underdetermination departure — dim-8 / return-nullity-2, the residuals vs pathA_29, the `Z0_ret/Z1_ret` Gate-6
selector need — that DELIVERS the FAIL is **stage 023** (II-G5b).

**Verdict.** `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (JOINT, 2-stage) — landed here as an **EARNED-FIRST PARTIAL** component (022
EARNED; 023 PENDING, and it is 023 that delivers the FAIL). ⭐ **Two distinct verdicts, both printed:**
```
LOCAL_AUDIT_VERDICT: CROSS_L_FINGERPRINT_OK        (022's own gate — fingerprints + non-regression; the exit-0 gate)
JOINT_LANDING_LABEL (PARTIAL; FAIL not evaluated by 022): FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)
```
022's SCRIPT passes (exit 0) — the earned content is correct and every tooth fires; the FAIL token is the joint gate's physics
label, carried as a 1/2 partial. 022 is the earned half of a gate that ultimately fails (cf. stage 003 = earned transverse
photons + the characterized `FAIL_CAUCHY` departure).

**Status.** Symbolic / exact / float-free / **z-space** (dimensionless — no units-restoration leg, unlike 018): the outgoing
cross-ℓ spherical-Hankel logarithmic derivatives and their small-frequency series are exact `sympy.series` / native `Series`
expansions; `bool_zero`/`expect_bool`/`expect_fail` residual asserts, no `scipy`/`numpy`/floats/tolerances. Dual-engine SymPy
**55 PASS** / Mathematica **64 PASS**, exit 0, CWD-independent (post-remediation; the build's 56/65 became 55/64 after an
honest de-count of one subsumed tooth per engine — §5).

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_34_cross_l_unification_{sympy.py,.wl}` +
> `reports/pathA_34_cross_l_unification.md` (the 022 slice = report :3, :5, :40, :47) + the original directive
> `directives/pathA_34_cross_l_unification.md`. The report/directive are cited for provenance only; the derivation below is
> self-contained.

---

## 1. What this stage earns

The outgoing radiative response of EACH multipole channel ℓ=0,1,2, matched to a small-frequency expansion, has a rational
**fingerprint** — a DERIVED signature that generalizes stage 018's ℓ=2-only result across ℓ. The ℓ=2 leg reproduces the
completed pathA_33 quadrupole (the Gate-4 non-regression), so the cross-ℓ machinery does not regress the earned result.

### 1.1 The cross-ℓ outgoing DtN log-derivative (self-contained exterior spherical-Hankel algebra)
For each ℓ the exterior outgoing solution is the spherical Hankel function of the first kind `h_ℓ⁽¹⁾ = j_ℓ + i·y_ℓ`
(`z = aω/c_s`). The Dirichlet-to-Neumann (DtN) eigenvalue for channel ℓ is the dimensionless logarithmic derivative, and the
normalized radiative response is
```
Λ_ℓᵒᵘᵗ(z) = z · h_ℓ⁽¹⁾′(z) / h_ℓ⁽¹⁾(z),      Ŷ_ℓᵒᵘᵗ(z) = −(ℓ+1) / Λ_ℓᵒᵘᵗ(z).
```
The static-DtN slot `Λ_ℓ(z→0) = −(ℓ+1)` is **DERIVED** from the Hankel log-derivative (`h_ℓ ~ z^{−(ℓ+1)}` as z→0 ⟹
`z·h'/h → −(ℓ+1)`) — read off `lam_series.coeff(z,0)`, NOT the hand-set numerator (the source's `−(ℓ+1)` set-then-compare is
de-rigged). `h_ℓ` is built from explicit rational·sin/cos (`.py`) / the built-in `SphericalHankelH1` (`.wl`); no stored
artifact is read (022 is self-contained exterior-wave algebra; the interior port kernel is cited as provenance — §4).

### 1.2 The cross-ℓ fingerprint (the EARNED headline)
Series-expanding each `Ŷ_ℓᵒᵘᵗ` about `z=0`:
```
ℓ=0:  Ŷ₀ᵒᵘᵗ(z) = 1 + i·z + O(z²)                                  → radiative_coeff = 1     at order ω¹ (z¹)
ℓ=1:  Ŷ₁ᵒᵘᵗ(z) = 1 + ½z² + i·½z³ + O(z⁴)                          → radiative_coeff = 1/2   at order ω³ (z³)
ℓ=2:  Ŷ₂ᵒᵘᵗ(z) = 1 + (1/9)z² + (4/81)z⁴ + i·(1/27)z⁵ + O(z⁶)      → radiative_coeff = 1/27  at order ω⁵ (z⁵)
```
so the leading radiative coefficient is `{ℓ=0: 1, ℓ=1: 1/2, ℓ=2: 1/27}` at leading order `2ℓ+1 = {ω¹, ω³, ω⁵}` (via
`z = aω/c_s`). The `radiative_coeff` is the imaginary series coefficient at `z^(2ℓ+1)`, extracted as `coeff(z, 2ℓ+1)/i`, and
checked **derive-vs-typed** (`bool_zero(derived − expected)`) — the derived side is the ACTUAL series of `−(ℓ+1)/(z·h'/h)`,
EMITTED, not hardcoded. The **radiative order** is verified by SCANNING the imaginary series for its first nonzero power and
asserting all lower imaginary coefficients vanish (NOT the source's dodgeable "check the preselected `2ℓ+1` coeff is nonzero"
— a `+i·z` corruption of the ℓ=2 series leaves `z^5` nonzero yet the true leading imaginary order becomes 1; the scan catches
it). The static slot `static = Ŷ_ℓ(0) = 1` is a de-counted printed DIAGNOSTIC (algebraically `−(ℓ+1)/Λ_static = 1`, subsumed
by the derived `Λ_static`).

### 1.3 The incoming-sign-flip (the genuine sign tooth)
The incoming branch `h_ℓ⁽²⁾ = j_ℓ − i·y_ℓ` flips ONLY the radiative (imaginary) sign — `radiative_coeff → {−1, −1/2, −1/27}`,
with every lower real coefficient unchanged (`bool_zero(incoming + outgoing) = 0`). This sign-flip is the genuine, computed
sign content. ⚠ The `χ_Q = radiative_coeff / (canonical slot) = +1` equality is a de-counted printed DIAGNOSTIC — it is
algebraically subsumed by the `v₅ = 1/27` non-regression (`27·v₅ − 1 ≡ 27·(v₅ − 1/27)`), so once the ℓ=2 slot passes,
`χ_Q = 1` is automatic; it is not an independent tooth.

### 1.4 The Gate-4 quadrupole non-regression
The ℓ=2 leg (`−(ℓ+1)/Λ_ℓ|_{ℓ=2} = −3/Λ₂`) reproduces stage 018's earned quadrupole fingerprint:
```
u₂ = 1/9,   u₄ = 4/81,   v₅ = 1/27   (the z², z⁴ reactive coefficients + the z⁵ radiative coefficient).
```
This is a **derive-vs-typed** non-regression: 022 re-derives the ℓ=2 tuple from its own Hankel and checks it against stage
018's independently-earned TYPED literals `{1/9, 4/81, 1/27}` (NOT a same-machinery reconstruction of the typed side — no
subsumed X≡X). The cross-ℓ `−(ℓ+1)` normalization gives the RIGHT `−3` at ℓ=2, so the general machinery does not regress the
quadrupole. Each of `u₂, u₄, v₅` has its OWN coefficient-isolating mutant (the `3e` incoming mutation flips ONLY `v₅`, so
`u₂`/`u₄` need their own derivation-copy mutants — the stage 018 pattern); breaking the ℓ=2 branch (→ incoming) fires
`FAIL_QUAD_REGRESSION`.

---

## 2. The able-to-fail battery (022-owned)

The verdict runs a SCOPED gate chain — an **022-LOCAL** verdict over the fingerprints + the ℓ=2 non-regression ONLY (023's
nullspace/return/`base_verdict` is NOT built here; the local read-set is a computed `{cross_l_fingerprints,
ell2_non_regression}` that provably excludes them). The 022 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| cross-ℓ radiative coefficient (per ℓ) | corrupt the Hankel derivation → `radiative_coeff ≠ {1, 1/2, 1/27}` → `FAIL_FINGERPRINT` | derive-vs-typed; the derived side is independent of the typed target |
| derived `Λ_static = −(ℓ+1)` | **pole-order** corruption `h_mut = z·h` → `Λ_static: −(ℓ+1)→−ℓ` → fires | ⚠ NOT "outgoing→incoming" (inert — H1↔H2 both give `−(ℓ+1)`); reads the computed `lam_series.coeff(z,0)` |
| scanned radiative order | add `+i·z` to ℓ=2 → first nonzero imaginary power = 1 ≠ 5 → fires | ⚠ NOT the preselected-`2ℓ+1`-nonzero check (which the `+i·z` corruption dodges) |
| incoming-sign-flip | branch outgoing→incoming → radiative sign flips → fires | the genuine computed sign content |
| ℓ=2 non-regression per-slot (u₂, u₄, v₅) | an isolating derivation-copy mutant changes ONLY its slot → its own assert fires | vs stage 018's typed `{1/9, 4/81, 1/27}`; u₂/u₄ have their own mutants (3e flips only v₅) |
| `3e_break_gate4` + dynamic self-ablation | incoming ℓ=2 → `v₅` flips → `FAIL_QUAD_REGRESSION`; the self-ablation is a DYNAMIC two-verdict re-run (`with ≠ without`) re-scoped to the 022-LOCAL verdict | neutering the mutation is detected |
| 022-LOCAL read-set | wire a nullspace/`base_verdict` object into the verdict → the read-set assert fires | proves 022 builds no nullspace (the DEFAULT-verdict-is-PASS trip-up avoided) |

De-counted printed DIAGNOSTICS (NOT verdict teeth): `static = 1` (subsumed by `Λ_static`) and `χ_Q = 1` (subsumed by `v₅`).

<!-- §5 tri-review fills the per-tooth-ablation tally -->

---

## 3. Honest scope

- **EARNED cross-ℓ signature + non-regression / the FAIL is 023's.** 022 DERIVES the outgoing cross-ℓ fingerprints
  (`{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}`, with `Λ_static = −(ℓ+1)` derived + the scanned order) and the Gate-4 quadrupole
  non-regression — the report's `Earned:` FINGERPRINT / raw-order / Gate-4 subset (`:47`). The rest of the report's earned
  content — the residual form/sign/order (conditional on a positive bounded branch) — PLUS the FAIL-delivering native
  nullspace underdetermination (dim-8 / return-nullity-2, `ε_eff` left free at the linear Gate-5 level) and the deferred
  scalar/dipole magnitude/nonzero prediction belong to **stage 023**.
- **The joint FAIL is DEFERRED to 023, not evaluated here.** The `JOINT_LANDING_LABEL` is a printed string; 022 builds no
  nullspace and does NOT back-solve `ε_eff`/`Z` (the `FAIL_TAUTOLOGICAL` firewall, 023's, is not reintroduced).
- **z-space only.** The fingerprints are dimensionless z-space rationals (`_z` names); there is NO units-restored dim leg (the
  `z = aω/c_s` units realization + the residual dimensional gate are 023's, on the ℓ=0/1 return amplitudes).
- **Gate 6 sim-deferred.** The `Z0_ret/Z1_ret` selector that would fix the ℓ=0/1 return (the first concrete Gate-6 input) is
  the deferred nonlinear solve — 023's characterized departure, not 022's.

---

## 4. Consumed / exported

- **Consumed — the CHECKABLE non-regression (one dual-site) + PROVENANCE cites.**
  - **stage 018's fingerprint `{u₂=1/9, u₄=4/81, v₅=1/27}` — the CHECKABLE derive-vs-typed non-regression** (the one genuinely
    checked consumption): 022's re-derived ℓ=2 fingerprint reproduces stage 018's earned TYPED values, a `3e`/per-slot-tooth-
    backed check (mutate the ℓ=2 branch → the derived coeff changes → the match fails; the emitted-vs-checked test passes). The
    "second site" is stage 018's independently-earned LITERAL, NOT a same-machinery reconstruction.
  - **stage 019/020/021 + the completed pathA_33 `QUAD_CALIBRATED` joint — PROVENANCE / context** (the prefactor
    `P(ω)=D₀N/D^cons²`, the `54/5` magnitude, the μ̂₀-free dimensional closure are cited DONE, NOT re-derived; 022's
    non-regression is FINGERPRINT level only).
  - **008's raw ℓ=0/1/2 outgoing amplitudes** (the channel structure; 022 self-derives the DtN responses) + **009/010's bulk
    Helmholtz mode** — cited PROVENANCE (no dual-site).
  - **`c_s`** (the density sound speed, stage 005 R1 `c_s²=5Kρ⁴/m`) — the PROVENANCE of the units symbol in `z = aω/c_s`, NOT a
    consumed value (the earned rationals are `c_s`-free). ⚠ Distinct from the frozen-wall Helmholtz speed `c_S` (011–017).
- **Register.** ZERO new counted knobs (an EARNED/structural fingerprint + non-regression slice, like 016/018 / 011/012/014).
  The cross-ℓ DtN fingerprints introduce no calibration (pure exterior-wave rationals); `c_s` is R1-DERIVED (cited PROVENANCE,
  a units carrier); `a` is the `CONV` pin (R2-family). ⚠ The ℓ=0/1 stiffnesses `K0c`/`K_eta+2·T_Omega` and the return
  admittances `Z0_ret`/`Z1_ret` are **023's** (the transfers/nullspace), NOT counted at 022. New structural edge **R41** (the
  cross-ℓ outgoing-DtN `−(ℓ+1)/Λ_ℓ` fingerprint provenance — the ℓ=0/1/2 radiative signature `{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}`
  with `Λ_static = −(ℓ+1)` derived + the scanned order, plus the Gate-4 quadrupole non-regression reproducing the stage 018
  fingerprint component of the now-completed pathA_33 fold; discharges nothing — earned exterior-wave structure + a consistency
  non-regression, not a reduction). Part-II counted CALIB set unchanged at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) +
  `{T_Ω, β₂}` (017) = 6.
- **Exported.** The cross-ℓ `−(ℓ+1)/Λ_ℓ` fingerprints + the ℓ=0/1 radiative coefficients `{1, 1/2}` at orders `{ω¹, ω³}` (→
  **023**'s return residual amplitudes `A0/A1`, which realize `z^(2ℓ+1) → (aω/c_s)^(2ℓ+1)` × the return transmission `(1−T_ℓ)`)
  + the ℓ=2 non-regression result (the earned quadrupole survives the cross-ℓ generalization) → **023** (the nullspace
  departure builds on the earned fingerprints) + the Part-VII cross-ℓ consistency record.

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a **genuinely independent route,
RE-AUTHORED** (not a kept transliteration): it uses Wolfram's BUILT-IN `SphericalHankelH1`/`SphericalHankelH2` +
`SeriesCoefficient` — a materially different construction from the `.py`'s hand-built `sin`/`cos` `j_ℓ`/`y_ℓ` +
`Series`/`Coefficient` (the source `.wl`'s `branchData` was a near-transliteration of the `.py` `dtn_branch`, so it was
discarded, not kept — the mirror-policy transliteration screen). Agreement is transcript-level: both engines emit the derived
`radiative_coeff = {1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}`, `Λ_static = {−1, −2, −3}` (with `h_mut = z·h` firing each and H1↔H2 shown
inert), `static = 1` + `χ_Q = 1` as de-counted diagnostics, the scanned order, the ℓ=2 non-regression `{1/9, 4/81, 1/27}`, the
per-slot u₂/u₄/v₅ mutants, the incoming-sign-flip, and `3e` → `FAIL_QUAD_REGRESSION`. The stage-007 unevaluated-leakage failure
mode is guarded (arity self-check + transcript scan).

**Directive review** used the Codex→Grok→Codex bookend. Codex's design-review returned 5 BLOCKING + 1 nit — all genuine and
folded: the `.wl` must be re-authored independent (not the transliterated `branchData`); the normalization (`Λ_static`) and
raw-order teeth must be genuinely derived (a `Λ_static` from the log-derivative + a coefficient SCAN, replacing the source's
X≡X set-then-compare and preselected-order check); the ℓ=2 tuple was double-counted and χ_Q was algebraically subsumed (remove
the duplicate, de-count χ_Q, add per-slot u₂/u₄ mutants); the report's `Earned:` partition is a subset relation (022 = the
fingerprint/order/Gate-4 subset) with an explicit `LOCAL_AUDIT_VERDICT`/`JOINT_LANDING_LABEL` output; the checkable consumption
is stage 018's fingerprint (not the "complete fold"). Codex's confirm-pass then caught a further directive-level reasoning error
(the proposed `Λ_static` ablation "outgoing→incoming" is INERT — both branches give `−(ℓ+1)` — so the mutant must be a
pole-order corruption `h_mut = z·h`; and `static = 1` is subsumed → de-counted). A Grok-4.5 compute-verify of the folded
directive returned `DIRECTIVE_CLEAN`, independently confirming all five checks (the cross-ℓ series, the pole-order-vs-inert
mutant, the scan-vs-preselect counterexample, the χ_Q subsumption identity, the u₂/u₄ isolation). A final Codex confirm closed
the bookend `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents (arbiter re-run reproduced the build SymPy 56 / Mathematica 65, exit 0, CWD-independent):
`FIDELITY_CLEAN` (an independent read hand-re-derived the cross-ℓ fingerprints `{1, 1/2, 1/27}` at `{ω¹, ω³, ω⁵}`, the derived
`Λ_static = −(ℓ+1)` with the pole-order `h_mut = z·h` mutant firing while `H1↔H2` is inert, the scan-vs-preselect order
counterexample, the `χ_Q` subsumption identity, and the per-slot u₂/u₄/v₅ isolation — everything COMPUTED not stamped, no 019
prefactor / 023 nullspace leakage, the `.wl` a genuine built-in-`SphericalHankelH1` route) + `ADVERSARIAL_ISSUES` (per-tooth
ablation: both confirm-pass reasoning catches HELD live-proven — the pole-order mutant is load-bearing [removing the `z` pole
makes it not fire], `static=1`/`χ_Q=1` are print-only in both engines — and all core physics teeth are genuine and
live-verified; **2 LOW-severity redundancy teeth flagged**: F1 the `3e` `rerun_gate_logic` checked a constant `len(trace)==2`
[vacuous], F2 the "read-set excludes forbidden" tooth was subsumed by the adjacent exact-equality "reads exactly {the two}"
tooth). Codex remediated: F1 → make-genuine (compare the two traced verdicts, so a neutered ablation fires it), F2 → honest
de-count (a printed diagnostic; the exact-equality tooth retains the exclusion guard), plus two `.wl` nits (a symbol typo +
`IncomingLowerRealUnchanged` mutation-aware parity with the `.py`). Fresh-agent `REVERIFY_CLEAN` (coupling meta-test: F1 fires
when the ablation is neutered and goes vacuous when the fix is reverted; F2's retained exact-equality tooth fires on a wired
forbidden read; the parity fix fires on a real mutation; no regression, both engines exit 0 at repo root AND `/tmp`). Tallies
56/65 → **55/64** (net −1 per engine from the honest F2 de-count). Symbolic per-tooth ablation, mutations on copies.
