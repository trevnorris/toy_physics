# ledger_stage019 — the squared-denominator prefactor algebra `P(ω)=D₀N/D^cons²` (Check II-G4b)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative sector). The SECOND leg of a 4-way split of
`pathA_33`: this stage carries the **squared-denominator prefactor-algebra component (2/4) of the joint `QUAD_CALIBRATED`**.
The outgoing ℓ=2 DtN Hankel fingerprint + χ_Q sign was **stage 018** (II-G4a, done); the `54/5=2·27/5` provenance partition
+ the CALIBRATED verdict label is **stage 020** (II-G4c); the μ̂₀-free `[P₀^phys]=1` dimensional closure is **stage 021**
(II-G4d).

**Verdict.** `QUAD_CALIBRATED` (JOINT, 4-stage) — landed here as a **PARTIAL** component (018 done, 019 EARNED; 020/021
PENDING). Ledger earned-label: `PREFACTOR_ALGEBRA_EARNED`.

**Status.** Exact closed-form / symbolic (float-free): the prefactor coefficients are exact rational functions of the
abstract port scalars, read off a symbolic `sympy.series` / native `Series` of `D₀N/D^cons²`; `expect_zero`/`expect_bool`/
`expect_fail` residual asserts, no `scipy`/`numpy`/floats/tolerances (the light integer sample substitution is a
transcript-level cross-check only). Dual-engine SymPy **18 PASS** / Mathematica **24 PASS** (the +6 = the `.wl` arity +
unevaluated-leakage self-check block), exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}` +
> `reports/pathA_33_quadrupole_normalization.md` (the 019 slice = report :15–18, :41) + the original directive
> `directives/pathA_33_quadrupole_normalization.md` (§2.2 the prefactor algebra, §3g the wrong-prefactor-object probe, §4
> the anti-hardcode firewall). The report/directive are cited for provenance only; the derivation below is self-contained.

---

## 1. What this stage earns

The ℓ=2 radiative port's normalization prefactor — the frequency-dependent factor multiplying the fingerprint — has a
DERIVED squared-denominator algebra: the coefficients `P₀/P₂/P₄` of `P(ω)=D₀N/D^cons²` in terms of the abstract port
scalars, with the `−2D₂N₀` factor-of-2 the provable signature of the squared (not plain) denominator. This is where stage
018's deferred literal `N_n/D_n` port-kernel consumption lands — but at the **provenance** level (§4).

### 1.1 The squared-denominator prefactor object
The port numerator and the conservative ("consumption") denominator are
```
N(ω)      = N₀ + N₂ ω² + N₄ ω⁴,          (the port numerator, N-moments)
D^cons(ω) = D₀ + D₂ ω² + D₄ ω⁴,          (the D-lanes / conservative scalars)
P(ω)      = D₀ · N(ω) / D^cons(ω)²        ⭐ the SQUARED denominator (NOT plain N/D).
```
The `D_n` are stage 017's exported ℓ=2 port-kernel D-lanes and the `N_n` are the concrete port N-moments of
`build_port_moments`; here they are carried as **abstract, port-agnostic symbols** — the algebra below holds for ANY nonzero
`D₀..N₄`, so no port value is consumed (the numerical `(D_n,N_n)` are deferred Gate-6 branch data).

### 1.2 The prefactor algebra (the EARNED headline)
Expanding `P(ω)` about `ω=0` (via `(1+x)⁻²=1−2x+3x²−…` with `x=(D₂/D₀)ω²+(D₄/D₀)ω⁴`) and reading the coefficients off the
ACTUAL series:
```
P₀ = N₀/D₀
P₂ = (D₀N₂ − 2·D₂N₀) / D₀²
P₄ = (D₀²N₄ − 2·D₀(D₂N₂ + D₄N₀) + 3·D₂²N₀) / D₀³
```
⭐ **The `−2D₂N₀` term is the signature of the squared denominator** — it arises only from the `−2x` term of `(1+x)⁻²`. The
coefficients are `coeff(ω,n)` off the ACTUAL symbolic series of `D₀N/D^cons²` (`.py` `series_no_o` / `.wl` `serW`), with the
object + series EMITTED — not hardcoded and compared. The derive-then-check `coeffs[n] − expected[n] == 0` has an
INDEPENDENT typed reference side: mutating the object (to plain `N/D`, or perturbing a `D_n`/`N_n` power) changes the
series-extracted coefficient and fires exactly that coefficient's assert; a swapped-in **correct** object flips the
discriminator to True — a genuine positive control (the made-genuine tooth 16, §2).

### 1.3 The N/D self-check (the sharp discriminator)
The plain (non-squared) object exposes the factor of two directly. `N/D^cons` series-expands (via `(1+x)⁻¹=1−x+…`) to
```
P₂^plain = (D₀N₂ − D₂N₀) / D₀²      →      P₂^plain − P₂ = D₂N₀/D₀² ≠ 0,
```
i.e. `−D₂N₀` versus the correct `−2D₂N₀`. The gap `D₂N₀/D₀²` is **computed** (not asserted). The counterfactual probe
`3g_wrong_prefactor_object` replaces `D₀N/D^cons²` with plain `N/D`, so `plain_equals_correct_P2 = False` → the probe fires
`FAIL_PREFACTOR_ALGEBRA`, while the correct object does NOT fire. Its self-ablation is a genuine **DYNAMIC** re-run of an
**019-local** verdict gate (the prefactor match ∧ the N/D self-check; `with_mutation ≠ without_mutation`, never a constant,
and never the joint `base_verdict` — the pathA_33 v1 constant-`self_ablation` trip-up avoided).

### 1.4 Units-free scope (no dimensional leg)
019 introduces NO `c_s`/`a`/`G`/`μ̂₀` and has NO dimensional leg — the algebra is over the abstract `{ω, D₀..N₄}` only (the
μ̂₀-free `[P₀^phys]=1` closure that assigns `[N₀]/[D₀]` is stage 021's; the `54/5`/`G` magnitude is stage 020's). A runtime
free-symbol guard asserts the earned expressions' symbols are exactly `{ω, D₀,D₂,D₄,N₀,N₂,N₄}`.

---

## 2. The able-to-fail battery (019-owned)

The verdict runs a SCOPED gate chain (019's gates only — the prefactor match ∧ the N/D self-check; 018/020/021's
fingerprint/partition/dimensional-closure gates are NOT computed here). The 019 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| series-extracted P₀/P₂/P₄ (per-coefficient) | corrupt the object (plain `N/D`, drop the squaring, perturb a `D_n`/`N_n` power) → a coefficient ≠ its expected form → assert fires | derived side is `coeff(ω,n)` off the actual series, independent of the typed reference |
| N/D self-check | `plain_equals_correct_P2` = computed `bool(D₂N₀/D₀² == 0)` = False; force True → `3g` stops firing | the gap `D₂N₀/D₀²` is computed, not stamped |
| `3g_wrong_prefactor_object` | plain `N/D` → missing factor of 2 → `FAIL_PREFACTOR_ALGEBRA`; correct object → `NO_FAIL` | DYNAMIC 019-local self-ablation (`with ≠ without`), not a constant, not the joint `base_verdict` |
| swapped-in correct-object positive control (tooth 16) | put `correct_obj` in the discriminator's plain slot, series-extract → discriminator flips to True | a genuine positive control (made-genuine in remediation; not a mirror of the P₂ tooth) |
| units-free live-symbol guard | inject a stray symbol beyond the series order, or shrink the allowed set → guard fires | asserts free symbols == `{ω, D₀..N₄}` (the units-free enforcement — not a source grep) |

Adversarial per-tooth ablation (`ADVERSARIAL_ISSUES` → remediated → `REVERIFY_CLEAN`): 21 live `.py` mutations + static
`.wl` analysis. The central firewall (series-extraction), the N/D factor-of-2 self-check, and the DYNAMIC 019-local
ablation are all genuine. Two subsumed/mirror teeth were flagged and fixed: **tooth 16** (was a mirror of the P₂ tooth) →
**made-genuine** as the swapped-in correct-object positive control (coupling meta-test decisive: mutate its object → fires;
neuter the fix → vacuous); **tooth 11** (an assert on the constant `rerun_gate_logic:True`) → honestly **de-counted** (the
dynamic property is guaranteed by the surviving `with==FAIL ∧ without==EARNED ∧ fail_suppressed` teeth). No vacuous tooth
remains; the port-moments block is a labeled non-check (asserts nothing).

---

## 3. Honest scope

- **EARNED prefactor structure / CALIBRATED magnitude.** 019 DERIVES the port-agnostic algebra (`P₀/P₂/P₄` + the N/D
  factor-of-two self-check). The `54/5=2·27/5` normalization magnitude and `G` are CALIBRATED (`G=GENUINE_BLOCKED`), landed
  at stage 020; the DtN fingerprint is stage 018; the μ̂₀-free `[P₀^phys]=1` dimensional closure is stage 021.
- **The numerical port scalars are deferred.** The abstract `D_n`/`N_n` are 017's D-lanes + `build_port_moments`' concrete
  N-moments carried symbolically; the numerical `(D_n,N_n)` from the solved nonlinear branch remain Gate-6 sim-deferred
  work (report :49).
- **Units-free by construction.** No `c_s`/`a`/`G`/`μ̂₀` appears as a live symbol; the units restoration + dim closure are
  021's (this is the key difference from 018, which carried `c_s` as a units carrier).

---

## 4. Consumed / exported

- **Consumed — PROVENANCE ONLY (NOT a checkable cross-stage relation; NO guard/dual-site).** ⚠ This is where 018's deferred
  literal `N_n/D_n` port-kernel consumption lands — at the **provenance** level, not a value-consumption dual-site. Because
  the algebra is PORT-AGNOSTIC (holds for any nonzero `D₀..N₄`) and `build_port_moments` is **emitted-but-never-checked**
  (its output flows only to the payload/serialization, read by no residual/match/verdict — grep-confirmed by both the
  directive Codex→Grok→Codex pass and the tri-review), there is NO checkable relation to guard (the same provenance-only
  landing as stage 018; a guard on the unused/deferred moments would be a vacuous tooth). Cited as PROVENANCE (narrative):
  - **017's ℓ=2 port kernel** — the abstract `D₀,D₂,D₄` are 017's exported D-lanes (the conservative scalars, + support/
    Maxwell scalars `B̃/Z̃`), carried as generic symbols.
  - **`build_port_moments`' concrete N-moments** — the `N₀,N₂,N₄` are the port N-moments (`N_A0_r=P_port²/Δ²`, etc.),
    DEFERRED Gate-6 branch data, asserted against nothing.
  - **018's fingerprint context** — 019's `P(ω)` is the prefactor the 018 fingerprint's `v₅` slot feeds; narrative cite only
    (019 does NOT re-derive 018's `1/9,4/81,1/27`/`χ_Q`).
- **Register.** ZERO new counted knobs (an EARNED/structural prefactor-algebra slice, like 018 / 016 / 011/012/014). New
  structural edge **R38** (the squared-denominator prefactor-algebra provenance: `P₀=N₀/D₀`, `P₂=(D₀N₂−2D₂N₀)/D₀²`,
  `P₄=(D₀²N₄−2D₀(D₂N₂+D₄N₀)+3D₂²N₀)/D₀³` SERIES-EXTRACTED, the `−2D₂N₀` factor-of-2 the squared-denominator signature;
  discharges nothing — earned port-agnostic abstract algebra, not a reduction). Part-II counted CALIB set unchanged at
  `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6.
- **Exported.** The prefactor algebra — `P₀=N₀/D₀` (→ stage 021's μ̂₀-free `[P₀^phys]=1` dim closure builds on it),
  `P₂/P₄` + the N/D self-check (→ stage 020's `54/5` partition context) → 020/021 + 027 (pathA_43 closure).

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a **genuinely independent native route**
(native `serW`/`Coefficient`/`FullSimplify` on its own `prefObj`/`plainObj`, typing its own expected `P₀/P₂/P₄`) — no
`Get`/`Import`/`Export`/YAML, no mirroring of the `.py`; the source pair's scratch-YAML engine-agreement handoff is severed.
Agreement is transcript-level: both engines emit the series-extracted `P₀=N₀/D₀`/`P₂=(D₀N₂−2D₂N₀)/D₀²`/`P₄=…`, the plain
`P₂=(D₀N₂−D₂N₀)/D₀²`, the gap `D₂N₀/D₀²`, and `plainEqualsCorrectP2=False`. The stage-007 unevaluated-leakage failure mode
is actively guarded (arity + `SeriesData`/`Coefficient`/`Series`-leakage self-check).

**Directive review** used the Codex→Grok→Codex bookend: Codex `DIRECTIVE_ISSUES` (3 nits, no BLOCKING — a line-ref fix
L512→L401, the D-lane/N-moment provenance split, the sibling-symbol wording; all folded); Grok-4.5 compute-verify
`DIRECTIVE_CLEAN` (independently re-derived `P₀/P₂/P₄` + the plain-`N/D` gap `D₂N₀/D₀²`, grep-confirmed `build_port_moments`
emitted-but-never-checked) + 1 genuine nit (D0..N4-only sample subs, dropping the source `.wl`'s `a/cs/c/G` carriers) folded;
a final Codex confirm-pass (folds 1/3/4 clean + 2 straggler wording spots fixed) → clean.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read hand-re-derived `P₀/P₂/P₄` via the `(1+x)⁻²`
expansion + the `N/D` gap `D₂N₀/D₀²` — no dropped check, no transliteration error, the units-freeness enforced by the
runtime symbol guard) + `ADVERSARIAL_ISSUES` (per-tooth ablation: the central firewall / N/D self-check / dynamic 019-local
ablation all genuine; two subsumed/mirror teeth — tooth 16 made-genuine as a swapped-in correct-object positive control,
tooth 11 de-counted — plus a de-obfuscation of narrative tokens that had been string-concatenated to dodge a source grep
(the pathA_41 anti-pattern) and an `assert_no_float` legibility fix). Fresh-agent re-verify `REVERIFY_CLEAN` (the coupling
meta-test on tooth 16 decisive: mutate its object → fires; neuter the fix → vacuous; no live-symbol leak from the
de-obfuscation; no regression). Tallies 19/25 → **18/24** (net −1 per engine from the honest tooth-11 de-count). Symbolic
per-tooth ablation, mutations on copies.
