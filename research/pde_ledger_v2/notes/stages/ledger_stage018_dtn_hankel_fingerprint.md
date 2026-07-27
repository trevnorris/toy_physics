# ledger_stage018 — the outgoing ℓ=2 DtN Hankel fingerprint + χ_Q sign (Check II-G4a)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative sector). The EARNED-FIRST leg of a 4-way split of
`pathA_33`: this stage carries the **outgoing ℓ=2 DtN spherical-Hankel fingerprint + the χ_Q sign component (1/4) of the
joint `QUAD_CALIBRATED`** — the EARNED exterior-wave signature. The prefactor algebra is **stage 019** (II-G4b); the
`54/5=2·27/5` provenance partition + the CALIBRATED verdict label is **stage 020** (II-G4c); the μ̂₀-free `[P₀^phys]=1`
dimensional closure is **stage 021** (II-G4d).

**Verdict.** `QUAD_CALIBRATED` (JOINT, 4-stage) — landed here as a **PARTIAL** component (018 EARNED; 019/020/021 PENDING).
Ledger earned-label: `DTN_HANKEL_FINGERPRINT_EARNED`.

**Status.** Exact closed-form / symbolic (float-free): the outgoing ℓ=2 spherical-Hankel logarithmic derivative and its
small-frequency series are exact `sympy.series` / native `Series` expansions; `expect_zero`/`expect_bool`/`expect_fail`
residual asserts, no `scipy`/`numpy`/floats/tolerances. Dual-engine SymPy **59 PASS** / Mathematica **65 PASS** (the +6 =
the `.wl` arity + unevaluated-leakage self-check block), exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}` +
> `reports/pathA_33_quadrupole_normalization.md` (the 018 slice = report :3, :9, :11) + the original directive
> `directives/pathA_33_quadrupole_normalization.md` (§2.1 the fingerprint derivation, §2.4 the χ_Q sign, §4 the
> anti-hardcode firewall). The report/directive are cited for provenance only; the derivation below is self-contained.

---

## 1. What this stage earns

The ℓ=2 sector's outgoing radiative response, matched to a small-frequency expansion, has a rational **fingerprint** — a
DERIVED signature the quadrupole normalization (019–021) then rides. Together with the outgoing-vs-incoming **sign** it is
the earned exterior-wave content of Gate 4.

### 1.1 The outgoing DtN log-derivative (self-contained exterior spherical-Hankel algebra)
The exterior ℓ=2 outgoing solution is the spherical Hankel function of the first kind,
```
h₂⁽¹⁾(z) = j₂(z) + i·y₂(z),   z = a ω / c_s,
j₂(z) = (3/z³ − 1/z) sin z − (3/z²) cos z,
y₂(z) = (1/z − 3/z³) cos z − (3/z²) sin z.
```
The Dirichlet-to-Neumann (DtN) eigenvalue for the ℓ=2 channel is the dimensionless logarithmic derivative, and the
normalized radiative response is
```
Λ₂ᵒᵘᵗ(z) = z · h₂⁽¹⁾′(z) / h₂⁽¹⁾(z),      Ŷ₂ᵒᵘᵗ(z) = −3 / Λ₂ᵒᵘᵗ(z),
```
where the `−3 = −(ℓ+1)` for ℓ=2 is the outgoing static DtN slot (self-derived, not consumed). `j₂`/`y₂` are built from
explicit rational·sin/cos — no stored artifact is read (018 is self-contained exterior-wave algebra; it does NOT consume
the interior port kernel — see §4).

### 1.2 The fingerprint (the EARNED headline)
Series-expanding `Ŷ₂ᵒᵘᵗ` about `z=0`:
```
Ŷ₂ᵒᵘᵗ(z) = 1 + (1/9) z² + (4/81) z⁴ + i·(1/27) z⁵ + O(z⁶),
```
so the dimensionless coefficients are `u₂ᶻ = 1/9`, `u₄ᶻ = 4/81`, `v₅ᶻ = 1/27` (the odd z¹, z³ terms vanish; the z⁵ term is
purely imaginary — the RADIATING part, extracted as `coeff(z,5)/i`). Restoring units through `z = aω/c_s`:
```
u₂ = a²/(9 c_s²),   u₄ = 4a⁴/(81 c_s⁴),   v₅ = a⁵/(27 c_s⁵).
```
⭐ **The `27` is EARNED here** as `1/v₅ᶻ` — the radiative-slot denominator, series-derived, NOT a typed literal. This is the
`derived_in_gate` "27" that stage 020's `54/5 = 2·27/5` provenance partition rides. The even `u₂, u₄` are the reactive
(non-radiating) static-response coefficients; the imaginary `v₅ω⁵` is the leading radiation.

The coefficients are read off the ACTUAL symbolic series of `−3/(z·h₂⁽¹⁾′/h₂⁽¹⁾)` (`coeff(z,n)` / `Coefficient[…,z,n]`),
with `h₂⁽¹⁾`/`Λ₂ᵒᵘᵗ`/`Ŷ₂ᵒᵘᵗ` EMITTED — not hardcoded and compared. The derive-then-check `out[k] − expected[k] == 0` has
an INDEPENDENT derived side: perturbing the typed target (e.g. `u₂ᶻ` expected `1/9 → 1/8`) fires the assert (residual
`−1/72`), and mutating the DERIVATION itself (multiply `h` by `(1+z⁴)` → `u₄ᶻ` becomes `112/81`, leaving `u₂ᶻ`, `v₅ᶻ`
untouched) fires exactly that coefficient's `expect_fail` — a genuine per-coefficient firewall (directive §4).

### 1.3 The χ_Q sign (EARNED — outgoing vs incoming)
Against the canonical radiative slot `a⁵/(27 c_s⁵)`,
```
χ_Q = v₅ / (a⁵/(27 c_s⁵)) = +1   (outgoing h₂⁽¹⁾),      χ_Q = −1   (incoming h₂⁽²⁾ = j₂ − i·y₂).
```
The sign is COMPUTED: it EMERGES from `j₂ ± i·y₂` (the z⁵ coefficient flips sign between the two branches), never a typed
`χ_Q = +1` (a typed value would be a tautology, directive §4). The standing branch `j₂` alone gives `Λ_stand(0) = +2` (=+ℓ,
NOT the outgoing `−(ℓ+1) = −3`), `Ŷ_stand(0) = −3/2 ≠ 1`, and NO imaginary radiating term — proving the `+1/27` radiative
rational is outgoing-BC-selected, not universal.

### 1.4 Passivity (not-static)
The radiating imaginary `v₅ω⁵` arises from the genuine outgoing-wave BC, not an inserted sink: the `3b_imposed_dissipation`
probe replaces the outgoing BC with a phenomenological damping term and is REJECTED (`FAIL_NOT_OUTGOING`) by source-label
provenance, and its self-ablation is a genuine DYNAMIC re-run of the scoped verdict (`with_mutation ≠ without_mutation`, not
a constant — the pathA_33 v1 trip-up avoided). The static↔dynamic consistency is the fingerprint's own static slot
`Ŷ₂ᵒᵘᵗ(ω→0) = 1` (the z⁰ coefficient — the DC limit of the SAME outgoing response; ⚠ NOT the prefactor `P₀ = N₀/D₀`, which
is 019's).

### 1.5 The units-restored dim leg
With `[a] = L`, `[c_s] = L T⁻¹`, the physical coefficients carry `[u₂] = [a²/c_s²] = T²`, `[u₄] = T⁴`, `[v₅] = T⁵`. The dim
function runs on the REAL physical coefficients (not a typed tuple); corrupting a coefficient's `a`/`c_s` power makes it
inhomogeneous → `FAIL_DIMENSIONAL`. ⚠ This is NOT the μ̂₀-free `[P₀^phys] = 1` gate (that is 021's — 018 does not touch
`P₀`, `μ̂₀`, `N₀`, `D₀`; see §3).

The Mathematica source has exactly ten dimension-valued objects. Its artifact coverage is explicit here because the
cross-engine comparator can enforce artifact name-set symmetry but cannot discover source objects that neither artifact
emits:

| `.wl` object and definition locus | artifact status | coverage reason / read locus |
|---|---|---|
| `zeroDim` (`.wl:84`) | not emitted | Private neutral element returned by `dimOf` for constants and empty symbol sets (`.wl:104,116`), not a named physical quantity or an expression-walker result. |
| `dimL` (`.wl:85`) | emitted as `a` | Read through the live `dimRules[a]` binding at the `DIM` print site (`.wl:388`). |
| `dimSpeed` (`.wl:86`) | emitted as `c_s0_dim` | Read through live `dimRules[cs]` (`.wl:389`); this is the register's bulk-density `c_s0`, not stage012's frozen-wall `c_S`. |
| `dimT2` (`.wl:87`) | not emitted | Declared expected-target literal, not a walker result; a divergent target fires `expectZero` against live `u2Dim` (`.wl:395`). |
| `dimT4` (`.wl:88`) | not emitted | Declared expected-target literal, not a walker result; a divergent target fires `expectZero` against live `u4Dim` (`.wl:396`). |
| `dimT5` (`.wl:89`) | not emitted | Declared expected-target literal, not a walker result; a divergent target fires `expectZero` against live `v5Dim` (`.wl:397`). |
| `u2Dim` (`.wl:229`) | emitted as `u2` | Computed by `dimOf[out["u2"], dimRules]` and printed directly (`.wl:390`). |
| `u4Dim` (`.wl:230`) | emitted as `u4` | Computed by `dimOf[out["u4"], dimRules]` and printed directly (`.wl:391`). |
| `v5Dim` (`.wl:231`) | emitted as `v5` | Computed by `dimOf[out["v5"], dimRules]` and printed directly (`.wl:392`). |
| `corruptedU2Dim` (`.wl:233`) | emitted as `corrupted_u2_dim` | Computed by the same walker from the live corruption expression and printed directly (`.wl:393`), so cross-engine comparison covers the ablation tooth's nonzero `L` channel. |

---

## 2. The able-to-fail battery (018-owned)

The verdict runs a SCOPED gate chain (018's gates only — the fingerprint match, the χ_Q/passivity gate, and the units-dim
leg; 019/020/021's prefactor/partition/dimensional-closure gates are NOT computed here). The 018 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| fingerprint derivation (per-coefficient) | mutate the derivation (`h·(1+z⁴)`, `−3→−2`) → a coefficient ≠ `1/9,4/81,1/27` → `FAIL_FINGERPRINT` | isolates per-coefficient; the derived side is independent of the typed target |
| χ_Q sign | branch outgoing→incoming → `v₅` sign flips → `χ_Q = −1 ≠ +1` → `FAIL_FINGERPRINT` | the sign is the computed ratio, never typed |
| standing contrast | `j₂` → `Λ_stand(0) = +2`, no radiating term | proves `+1/27` is outgoing-BC-selected |
| passivity (`3b_imposed_dissipation`) | inserted phenomenological sink → `FAIL_NOT_OUTGOING` | DYNAMIC self-ablation (`with ≠ without`), not a constant |
| units-restored dim leg | corrupt an `a`/`c_s` power → `[u₂] ≠ T²` → `FAIL_DIMENSIONAL` | reads the real physical coefficient (not a dim-tautology) |
| μ̂₀-free forbidden-token guard | inject a `mtilde0`/`N0`/prefactor token → the string-scan guard fires | keeps the v1 rig out of 018 |

Adversarial per-tooth ablation: 30/32 targeted mutations fired at their own assert (the 2 non-fires were harness artifacts
— a SymPy `trigsimp` timeout on non-Hankel input and a Mathematica auto-eval — re-covered by corrected runs + the
twin-engine cross-check). No vacuous tooth (`X≡X` / stamped literal / dim-tautology / subsumed guard) was found.

---

## 3. Honest scope

- **EARNED exterior-wave signature / CALIBRATED magnitude.** 018 DERIVES the outgoing fingerprint shape (`1/9, 4/81,
  1/27`) + the `χ_Q = +1` sign — the exterior radiative structure. The `54/5` normalization magnitude and `G` are
  CALIBRATED (`G = GENUINE_BLOCKED`), landed at stage 020; the prefactor algebra is stage 019; the μ̂₀-free dimensional
  closure is stage 021.
- **`c_s` is a units carrier, not a live consumed value.** The earned rationals `1/9, 4/81, 1/27` and `χ_Q` are `c_s`-FREE
  (dimensionless z-space; `χ_Q`'s `c_s⁵` cancels). `c_s` (the density sound speed, R1) appears only as the units-restoring
  symbol in the physical coefficients — the fingerprint is value-independent of it.
- **Deferred (Gate 4/5/6, sim-deferred).** The actual-branch `a`-scaling of `P₀ = N₀/D₀`, the numerical port scalars
  `N_n/D_n`, and the solved nonlinear branch data remain downstream work, not 018's (report :49).
- **⚠ The χ_Q reconciliation (a tracked Part-VII item — NOT merged here).** pathA_33 gives `χ_Q = +1` (the outgoing-DtN
  Hankel context); pathA_22b gave `χ_Q ≈ 0.712` (an older minimal-combination context) — same name, DIFFERENT
  computations (blueprint §8). 018 lands `χ_Q = +1` in ITS context and records the reconciliation as a Part-VII open item.

---

## 4. Consumed / exported

- **Consumed — PROVENANCE ONLY (NOT a checkable cross-stage relation; NO guard/dual-site).** ⚠ 018 is SELF-CONTAINED
  exterior spherical-Hankel algebra: the fingerprint is built from explicit `j₂`/`y₂` and never references the interior
  port-kernel D-lanes or a stored bulk-mode object. So — UNLIKE stage 017's genuine cross-stage dual-site on 016's
  λ/K₂-form — there is NO checkable relation to guard here (a guard on an unused object would be a vacuous tooth). Cited as
  PROVENANCE (narrative):
  - **017's ℓ=2 port kernel** (the grouped `M₂`, angular `K₂ = K̃ + 6·T̃_Ω`, support scalars `B̃/Z̃`, D-lanes) — the
    INTERIOR wall mode whose EXTERIOR radiative response 018 computes; the literal `N_n/D_n` consumption is **stage 019's**
    (the prefactor `P = D₀N/D^cons²` matches the exterior/interior at the port).
  - **009/010's bulk Helmholtz mode** — the exterior outgoing solution's bulk companion (018 reconstructs `h₂⁽¹⁾`
    self-contained).
  - **`c_s`** (the density sound speed, stage 005 R1 `c_s² = 5Kρ⁴/m`) — the PROVENANCE of the units symbol, NOT a consumed
    value. ⚠ Distinct from the frozen-wall Helmholtz speed `c_S` (011–017; same R1 form, evaluated at the wall density ρ*
    vs the bulk ρ0). `c_S` is not a live symbol here.
- **Register.** ZERO new counted knobs (an EARNED/structural fingerprint slice, like 016 / 011/012/014). `c_s` is
  R1-DERIVED (cited PROVENANCE, a units carrier — not a new knob); `a` is the `CONV` pin (R2-family); the port scalars
  `N_n/D_n` are 019's deferred Gate-6 branch data. New structural edge **R37** (the outgoing-DtN ℓ=2 Hankel-fingerprint
  provenance: the exterior radiative signature `u₂ = a²/9c_s², u₄ = 4a⁴/81c_s⁴, v₅ = a⁵/27c_s⁵` with the `27` computed from
  `v₅ᶻ`, and `χ_Q = +1` outgoing / `−1` incoming; discharges nothing — earned exterior-wave structure, not a reduction).
  Part-II counted CALIB set unchanged at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6. The χ_Q
  reconciliation (`+1` vs pathA_22b `≈0.712`) is flagged as a Part-VII item.
- **Exported.** The outgoing ℓ=2 fingerprint — `v₅ = a⁵/27c_s⁵` (→ the `27` that stage 020's `54/5 = 2·27/5` partition
  rides), `u₂ = a²/9c_s²`, `u₄ = 4a⁴/81c_s⁴` (→ stage 019's prefactor context), `χ_Q = +1` (→ stage 020's `Γ₅ = 2χ_Q·G/5c⁵`
  equivalence), and the outgoing/incoming/standing branch classification → 019/020 + 022 (pathA_34 non-regression) + 027
  (pathA_43 closure).

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a **genuinely independent native route**
(native `j2`/`y2` from explicit rational·sin/cos, native `serZ`/`Series`/`Coefficient`/`FullSimplify`, its own `dimOf`) —
no `Get`/`Import`/`Export`/YAML, no mirroring of the `.py`; the source pair's scratch-YAML engine-agreement handoff is
severed. Agreement is transcript-level: both engines emit `Λ₂ᵒᵘᵗ`/`Ŷ₂ᵒᵘᵗ`, the derived `u₂ᶻ=1/9`/`u₄ᶻ=4/81`/`v₅ᶻ=1/27`,
`χ_Q=+1`/incoming `−1`, the standing static slot, and the 018 PARTIAL landing. The stage-007 unevaluated-leakage failure
mode is actively guarded AND able-to-fail (the `SeriesData`/`Derivative`-leakage guards fire under injection); arity
self-check carried.

**Directive review** used the Codex→Grok→Codex bookend: Codex `DIRECTIVE_CLEAN` first pass (compute-verifying
`−3/(z·h₂⁽¹⁾′/h₂⁽¹⁾) = 1 + z²/9 + 4z⁴/81 + i·z⁵/27`, incoming `χ_Q=−1`, standing `Λ_stand(0)=+2`, the c_s-free-ness);
Grok-4.5 compute-verify `DIRECTIVE_CLEAN` (independently re-confirming the same series + branch signs) + three non-blocking
clarity nits folded (the `P₀` name collision between the prefactor `N₀/D₀` and the fingerprint static slot; the
`passivity_from_source` scope; the 3a/3b self-ablation scoped to 018's local gates); a final Codex confirm →
`DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read hand-re-derived `j₂`/`y₂`, `Λ₂ᵒᵘᵗ`, `Ŷ₂ᵒᵘᵗ`, the full
series `1 + z²/9 + 4z⁴/81 + i·z⁵/27`, `χ_Q=±1`, the standing `Λ_stand(0)=+2` — no dropped check, no transliteration error,
the μ̂₀-free cut grep-confirmed) + `ADVERSARIAL_CLEAN` (per-tooth ablation: 30/32 targeted mutations fired at their own
assert, the 2 non-fires harness artifacts re-covered by corrected runs; the fingerprint firewall, the χ_Q-sign branch
flip, the DYNAMIC passivity self-ablation, the μ̂₀-free forbidden-token guard, the no-theatrical-dual-site, and the
independent native `.wl` all confirmed; no vacuous tooth). No remediation required (the two documented caveats — a SymPy
`trigsimp` timeout on non-Hankel mutations re-covered by the twin engine, and the verdict-router-scoped self-ablation with
the physics separately ablated by the derivation-copy teeth — are non-defects). Symbolic per-tooth ablation, mutations on
copies.
