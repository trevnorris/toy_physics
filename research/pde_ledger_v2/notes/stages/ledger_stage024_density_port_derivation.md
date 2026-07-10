# ledger_stage024 — the density-native ℓ=2 quadrupole two-port derivation (Check II-P1)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative-port sector, Cluster C). The EARNED-FIRST derivation
leg of a 4-way split of `pathA_43`: this stage carries the **density-native ℓ=2 quadrupole two-port derivation (1/4) of the
joint `DENSITY_PORT_HOSTED`** — the DERIVED coupling STRUCTURE `N0_den` built as a two-port over `(q₂` wall, `Φ₂`
bulk-density`)`, retiring the old EM `A_w`/`U,W` vector scaffold. The COMPUTED vector-freedom taint is **stage 025** (II-P2);
the continuity-moment lineage token-check is **stage 026** (II-P3); the able-to-fail port checks + closure overlay (the
COMPLETING leg that LANDS the joint) is **stage 027** (II-P4).

**Verdict.** LOCAL `DENSITY_TWO_PORT_DERIVED` (the derivation-integrity token, exit 0) + a printed JOINT PARTIAL
`DENSITY_PORT_HOSTED (1/4)`. ⚠ 024 does NOT own the joint verdict — that LANDS at 027 (the COMPLETING leg). 024 is the
EARNED-FIRST leg (the 018/016/022 pattern), NOT the FAIL-delivering completing pattern of 023.

**Status.** Exact closed-form / symbolic (float-free): the 2×2 static-operator inverse and the resulting rational port
numerator are exact `sympy.Matrix.inv()` / native `Solve` algebra; `expect_zero`/`expect_bool` residual asserts, no
`scipy`/`numpy`/floats/tolerances. Dual-engine SymPy **7 PASS** / Mathematica **10 PASS** (the +3 = the `.wl` arity + the
Module-result-shape + the unevaluated-leakage self-checks), both exit 0, CWD-independent (repo root + `/tmp`), zero file I/O.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}`
> (the 024 slice = `schur_density_expression` / `densityExpression`) + `reports/pathA_43_density_quadrupole_port.md` (the
> JOINT `DENSITY_PORT_HOSTED`, report `:5–14`, `:57`) + the original directive
> `directives/pathA_43_density_quadrupole_port.md` (§1/§2/§3A the physics). ⚠ pathA_43 is CONTRACT-CLEAN — both engines were
> already standalone print-only with zero file I/O, so there is NO bridge to sever; the reshape was DECOMPOSITION + a
> derivation-genuineness UPGRADE (§5). The report/directive are cited for provenance only; the derivation below is
> self-contained.

---

## 1. What this stage earns

pathA_43 answered the Phase-A A1 question — is the ℓ=2 quadrupole radiative-port numerator `N₀` DENSITY-native, or was the
old EM vector (`A_w`/`U,W`) channel load-bearing? The answer is a genuine density two-port on `(q₂` wall, `Φ₂` bulk-density`)`,
and the EM vector scaffold RETIRES. **024 builds the DERIVATION** — the 2×2 static coupling operator over `(q₂,Φ₂)`, its
inverse, and the resulting two-port numerator `N0_den` — and EXPORTS `N0_den` for the sibling legs (025 vector-freedom taint,
026 continuity lineage, 027 port checks + closure) to guard.

**Physical picture** (`docs/conceptual_foundation.md` §3/§5): gravitational quadrupole radiation rides the medium's OWN
density/sound (`c_s`) ripple — `c_s` is the analog "speed of light" of the GRAVITATIONAL sector — NOT an EM vector field. The
moving throat's ℓ=2 coupling to the outgoing quadrupole `c_s` wave is a density two-port: the wall mode `q₂` coupled to the
bulk-density Helmholtz mode `Φ₂` through a projected-continuity interface coupling. (This is NOT "route (c)" — that term is
the unrelated χ_B wall-provenance fork.)

### 1.1 The static two-port operator (self-contained reduced-Lagrangian algebra)
The reduced, mass-normalized static Lagrangian of the ℓ=2 sector is a 2×2 quadratic form over the coordinate pair
`(q₂, Φ₂)` (the wall angular amplitude and the bulk-density amplitude), driven by a continuity/interface source:
```
M = [[ ϖ_q2 , −λ_c   ],          (g_q, g_φ) = g_base·(η_q, η_φ),
     [ −λ_c , ϖ_Φ2  ]]           g_base = √rho_eff · c_s² · I25 · Ξ_Q / a^(7/2).
```
`M` is symmetric with a NEGATIVE off-diagonal `−λ_c` (the wall↔bulk mixing); `g_base` is the continuity/interface source
amplitude carrying the `a^(−7/2)` of the structured ℓ=2 projection; `(η_q, η_φ)` split it into the two coordinate channels.
The operator entries are ABSTRACT symbols whose literal throat VALUES are SIM_DEFERRED (§3) — the two-port STRUCTURE holds for
any admissible `ϖ_q2, ϖ_Φ2, λ_c`.

### 1.2 The Φ₂ response and the port numerator (the EARNED headline)
The observed radiative amplitude is the Φ₂ component of the inverted response `M⁻¹·(g_q, g_φ)`. Writing the determinant
`Δ` and the adjugate·source numerator `P_den`,
```
Δ = ϖ_q2·ϖ_Φ2 − λ_c²,        P_den = ϖ_q2·g_φ + λ_c·g_q,
Φ₂-response = (M⁻¹·source)[Φ₂] = P_den / Δ,
```
so the ℓ=2 quadrupole radiative-port numerator is the SQUARED response:
```
N0_den = (Φ₂-response)² = P_den² / Δ²
       = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)² / ( a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)² ).
```
⭐ This is the EARNED content: a genuine `(q₂, Φ₂)` two-port with a continuity-DERIVED coupling, its host-set
`⊂ {q₂, Φ₂, c_s, a, rho_eff, I25, Ξ_Q, η_q, η_φ, ϖ_q2, ϖ_Φ2, λ_c}` and NONE of the vector symbols
`{A_w, F_μw, J^w, U, W, Ω_U, Ω_W, R_mix, g_U, g_W}`. `N0_den` is COMPUTED via `Matrix.inv()` (SymPy) / `Solve`+eliminate (WL),
NOT typed. The `(λ_c² − ϖ_Φ2·ϖ_q2)` in the denominator is the sign-equal print of `−Δ`; the canonical factored `N0_den` is
identical between the two engines.

### 1.3 The two independent routes (the dual-derivation)
The two engines derive the SAME `N0_den` by materially different algorithms — this is what makes the agreement a REAL
dual-derivation, not two typings of one closed form:
- **SymPy — the full inverse.** Form the 2×2 `static_operator`, take `static_operator.inv()`, read the Φ₂ component
  `response = (M⁻¹·source)[Φ₂]`, then `N0_den = response²`.
- **Mathematica — the Green-DtN elimination (KEEP-NATIVE).** A genuinely different elimination order: solve the WALL row for
  `q₂` first (`qRule = q₂ → (λ_c·Φ₂ + g_q)/ϖ_q2`), substitute into the BULK row
  (`phiEq = ϖ_Φ2·Φ₂ − λ_c·q₂ == g_φ /. qRule`), `Solve` the reduced scalar interface equation for Φ₂
  (`phiResponse = Φ₂ /. First[Solve[phiEq, Φ₂]]`), then `n0 = phiResponse²`.

`Solve`+substitution (`qRule`/`phiEq`) is a DtN/Green-function interface matching that SymPy's symmetric `Matrix.inv()` never
invokes → **keep-native is DEFENSIBLE** (like 018/019/021, unlike 022/023's re-author). Agreement is transcript-level and
verified by the ARBITER: `simplify(N0_den_py − N0_den_wl) == 0` (a canonical/algebraic comparison, not raw-string — the
leading-minus over `(λ_c² − ϖ_q2·ϖ_Φ2)` would false-flag a string diff). Neither engine reads the other.

### 1.4 The genuineness upgrade (the load-bearing inverse)
In the CURRENT pathA_43 SymPy source the full-inverse `response` was DECORATIVE — `N0_den` was assembled directly from the
hand-written `P_den²/Δ²`, so a corrupted inverse would NOT move `N0_den`. The reshaped 024 makes the inverse LOAD-BEARING:
`N0_den = make_N0(response)` with `make_N0(r) := r²`, so the DATAFLOW is `response → N0_den`. `P_den`/`Δ` are retained ONLY as
the independent adjugate/determinant trace form, guarded by the factorization cross-check `compact(response − P_den/Δ) == 0`.
The WL side (`n0 = phiResponse²`) was already load-bearing. ⭐ Because a typed `P_den²/Δ²` passes the cross-check AND prints an
identical transcript, a code-read alone cannot catch a decorative inverse; the enforcement is the RUNTIME dataflow probe
`make_N0(response + delta_probe) ≠ make_N0(response)` — a decorative `make_N0(r) := P_den²/Δ²` (ignoring `r`) makes the two
sides equal and FIRES the probe at its own named assert. (This is the analog of stage023's "genuine `Matrix.rank()`, not
zero-padded / hardcoded".)

### 1.5 The nonsingular domain and the coupling boundary
The inverse/elimination is defined ONLY where the determinant `Δ = ϖ_q2·ϖ_Φ2 − λ_c² ≠ 0`; positivity of the three operator
entries does NOT imply it (they can satisfy `λ_c² = ϖ_q2·ϖ_Φ2`). 024 makes `Δ≠0` an EXPLICIT domain condition + a NAMED
runtime guard BEFORE the inversion/`Solve` in both engines. The stronger `Δ>0` is the physical-stability reading — NOT gated
here, noted as the SIM_DEFERRED stability condition. The derivation boundary is the coupling: `g_base→0 ⟹ N0_den=0`, COMPUTED
(the zero-control recomputes the whole derivation with `g_base=0` and yields exactly 0), not asserted — while baseline
`N0_den ≠ 0` for nonzero symbolic coupling. (The standalone `FAIL_PORT_VANISHES` VERDICT is 027's.)

### 1.6 The provenance (EXHIBITED, not yet PROVEN)
The three operator entries carry their physical origins as a printed `physical_relations` trace:
```
ϖ_q2  ← pathA_32 wall :  K₂/M₂ = (c_s/a)²·κ_q   (the wall angular ℓ=2 operator; K₂ = ∫[T_w β₂'² + (K_η+6·T_Ω)β₂²])
ϖ_Φ2  ← pathA_29 bulk :  (c_s/a)²·(6 + (m·a)²)   (the bulk Helmholtz ℓ=2 mode; 6 = ℓ(ℓ+1) covariance)
λ_c   ← projected continuity/interface :  (c_s/a)²·λ̂_Q
```
024 EXHIBITS this provenance (the `physical_relations` print); the COMPUTED taint-set PROOF that the tags are genuine (not
self-asserted) is 025's. `I25` is a typed ℓ=2 continuity-moment input (its `∫Y₂*S_leak` lineage is validated in stage026, a
forward reference — like 018 used the deferred port kernel `N_n/D_n`).

---

## 2. The able-to-fail battery (024-owned)

The verdict runs a SCOPED gate chain (024's derivation-integrity teeth only — 025's taint graph, 026's lineage token-check,
and 027's dim/scaling/sign/closure gates are NOT computed here). Per-tooth ablation: every tooth fired at its own named
assert. The 024 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| A inverse/factorization + dataflow | corrupt `P_den` (drop `λ_c·g_q`) / `Δ` (wrong `λ_c²` sign) / the `M⁻¹` → `response ≠ P_den/Δ`; a decorative `make_N0(r):=P_den²/Δ²` → probe `make_N0(r+δ)=make_N0(r)` fires | 024's CORE; the RUNTIME `delta_probe` catches the decorative inverse a code-read cannot |
| B nonsingular `Δ≠0` guard | singular control `λ_c² → ϖ_q2·ϖ_Φ2` (`Δ→0`) → the named guard FAILS before the inversion, both engines | positivity of entries does NOT imply `Δ≠0` |
| C coupling-vanishes | baseline `N0_den ≠ 0`; zero-control bypasses `g_base=0` while retaining the expectation → recomputed numerator nonzero → `N0_den\|g_base=0 = 0` assert fires | TWO named asserts (baseline + zero-control) |
| D density-only host-set | inject a vector symbol (`·Ω_U/Ω_W`) → `∩ {A_w,…,g_W} ≠ ∅` → the membership assert fires | MANIFEST-construction check, NOT the 025 taint graph |
| E `.wl` arity / leakage | plant a call-site arity mismatch or leak an unevaluated `Solve`/helper into a transcript object → the named arity/leakage assert fires | the stage007 lesson (a mismatched call silently skips at exit 0) |

⚠ Two checks are DELIBERATELY NOT in-script teeth: the **two-route agreement** is the ARBITER's oracle (ACCEPTANCE-level —
neither engine reads the other, so it cannot fire at its own in-script assert; a one-engine mutation makes the two transcripts
DISAGREE → the comparator fails); the **`.wl` algorithm independence** (`Solve`+eliminate ≠ `Matrix.inv()`) is a
TRANSLITERATION-REVIEW property (replacing `Solve` with `Inverse` can preserve every value, so no algebraic assert can reject
it) — kept in the Codex→Grok→Codex screen, not the tooth list.

---

## 3. Honest scope

- **EARNED structure / CALIBRATED–SIM_DEFERRED magnitude.** 024 DERIVES the reduced coupling STRUCTURE (the two-port form,
  its host-variable set, the inverse/elimination). The literal magnitudes — `I25` value, `Ξ_Q`, `η_q`, `η_φ`, the `λ_c` throat
  value — are SIM_DEFERRED; `G`, `2/5`, `54/5` are CALIBRATED (`G = GENUINE_BLOCKED`). The 54/5 partition + closure overlay +
  the PN match-back are the SIBLING legs (027/028), NOT 024. 024 is the FORM derivation, not a magnitude derivation.
- **`c_s`/`a` are units carriers.** `c_s` is the density/sound (ripple) speed — the analog "speed of light" of the
  GRAVITATIONAL sector (stage005 R1 `c_s² = 5Kρ⁴/m`); `a` is the `CONV` pin. The structure holds symbolically in them.
- **⚠ `rho_eff` is NOT stage005's `ρ0`.** `rho_eff` is the effective reduced-3D MASS density `[M L⁻³]` that mass-normalizes
  the `(q₂, Φ₂)` coordinates (the `√rho_eff` in `g_base`). Stage005's registered `ρ0` is a 4D NUMBER density `[L⁻⁴]` — a
  DIFFERENT quantity; pathA_29 supplies NO reduction from `ρ0` to `rho_eff`. Provenance = the pathA_29 bulk-mode
  mass-normalization (STRUCTURAL); the literal reduction `rho_eff ← {ρ0, m, geometry}` is SIM_DEFERRED/GAP.
- **Deferred (sim-deferred / Gate-6).** The `ϖ_q2`/`ϖ_Φ2`/`λ_c` throat values, the `I25`/`Ξ_Q`/`η` magnitudes, and the `Δ>0`
  physical-stability reading remain downstream work, not 024's.

---

## 4. Consumed / exported

- **Consumed — PROVENANCE ONLY (NO cross-stage dual-site fires).** 024's consumptions are the operator entries' physical
  origins; the load-bearing genuineness check is INTRA-stage (the two routes + the inverse-is-load-bearing factorization
  cross-check). Cited as PROVENANCE (narrative):
  - **009/010's bulk ℓ=2 Helmholtz mode + 016's `ℓ(ℓ+1)=6` covariance** → `ϖ_Φ2 = (c_s/a)²·(6 + (ma)²)` (abstract entry,
    tag `pathA_29_bulk`).
  - **016/017's wall ℓ=2 angular operator** → `ϖ_q2 = K₂/M₂ = (c_s/a)²·κ_q` (abstract entry, tag `pathA_32_wall`;
    ⭐ a downstream pin of 017's tracked port scalars).
  - **The projected-continuity/interface operator + the `I25` moment** → `λ_c = (c_s/a)²·λ̂_Q` (lineage VALIDATED in stage026,
    forward ref; `λ_c` carries all three tags). PROVENANCE.
  - **stage005 (`c_s² = 5Kρ⁴/m`) + `a` CONV** — units carriers. ⚠ `rho_eff` = a NEW tracked SIM_DEFERRED/GAP mass density,
    NOT stage005's `ρ0` (§3). ⚠ NOT 024: 018's `dtn_sign`/`χ_Q=1`, 021's dim gate, the `54/5` closure — those are 027's.
- **Register.** ZERO new counted CALIB knobs (a DERIVATION/structure slice, like the EARNED-first legs 018/016/022). The
  operator entries `ϖ_q2`/`ϖ_Φ2`/`λ_c` are ABSTRACT symbols with SIM_DEFERRED values → tracked as `GAP`/deferred (Gate-6
  throat data), NOT counted CALIB; `I25`/`Ξ_Q`/`η_q`/`η_φ` = structural symbols with SIM_DEFERRED magnitudes →
  tracked-not-counted; `rho_eff` (`[M L⁻³]`) = a NEW tracked `GAP`/SIM_DEFERRED quantity (NOT `ρ0`, NOT a counted CALIB).
  Part-II counted CALIB set unchanged at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = **6**. New
  structural edge **R43** — the density-native ℓ=2 two-port ALGEBRAIC derivation provenance (the `N0_den` structure earned via
  the inverse/elimination; the density construction EXHIBITED; SCOPED to 024's algebra only), discharging nothing (a
  structure/provenance edge, like R37/R39/R41). ⚠ R43 does NOT record "the vector scaffold retired" — that is the JOINT result
  PROVEN by 025's vector-taint proof + the completing 027. Symbol dims recorded as PROVENANCE (`[I25]=L^(5/2)`, `[c_s]=LT⁻¹`,
  `[rho_eff]=ML⁻³`, `[a]=L`, `[ϖ_q2]=[ϖ_Φ2]=[λ_c]=T⁻²`, `[η_q]=[η_φ]=[Ξ_Q]=1` ⟹ `[N0_den]=L⁻¹M`); the units-restored dim
  CHECK itself is 027's, not run here.
- **Exported.** `N0_den` (the CENTRAL export, factored) → **025** (its `free_symbols` feed the taint graph + the ablation),
  **026** (the `I25` moment must appear in `N0_den.free_symbols`), **027** (the dimension `[N0_den]=L⁻¹M`, the a-scaling
  `P0_physical = (c_s/a)²·N0_den/D0` a-power −5, the radiative sign, and the K̄ closure) + Part VII. Also the static operator +
  the `physical_relations` provenance + the density-only host-set. ⭐ 024 EXHIBITS the retirement of the EM `A_w`/`U,W` vector
  scaffold (blueprint §3 Part II A1); 025 PROVES it vector-free.

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O — and pathA_43 was ALREADY contract-clean (no
`Export`/`Put`/`Write`/`>>`/`Import`; the `results.yaml` was written by an external orchestrator, not the engines), so the
reshape "sever the bridge" step is a NO-OP. This is the LIGHTEST reshape in Part II; the work was DECOMPOSITION + the
genuineness UPGRADE (§1.4) + the LOCAL/PARTIAL verdict framing. The `.wl` is a **genuinely independent native Green-DtN route**
(native `qRule`/`phiEq`/`Solve` eliminate-q2, its own `clean`/`fmt`/`$Assumptions`, no `Get`/`Import`/`Export`, no mirror of
the `.py`). Agreement is transcript-level: both emit `N0_den` (canonical factored), `port_picture: ii two-port(q2,Phi2)`, the
Φ₂-response `= P_den/Δ`, the `Δ≠0` guard status, the `physical_relations` provenance, the density-only host-set, the LOCAL
`DENSITY_TWO_PORT_DERIVED`, and the JOINT PARTIAL `DENSITY_PORT_HOSTED (1/4)`. SymPy 7 PASS / Mathematica 10 PASS (the +3 = the
`.wl` arity scan + the Module-result-shape check + the unevaluated-leakage transcript scan), both exit 0, CWD-independent.

**Directive review** used the Codex→Grok→Codex bookend: a Codex xhigh design-review folded 4 BLOCKING (the `Δ≠0`
nonsingular-guard + singular control; the `ρ`→`rho_eff` mis-provenance correction; the per-tooth independent-ORACLE
restructure; R43 scoped to 024's algebra only) + a confirm-pass folded 2 further BLOCKING (tooth A's RUNTIME dataflow probe
that FIRES on a decorative typed `P_den²/Δ²`; tooth C's two named asserts) — Codex compute-verified the `N0_den` formula,
`response ≡ P_den/Δ`, `[N0_den]=L⁻¹M`, and the a-powers −3→−5; a Grok-4.5 compute-verify pass re-confirmed them.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read re-derived `N0_den` by a THIRD `linsolve` route,
diff = 0 — no dropped check, no transliteration error) + `ADVERSARIAL_CLEAN` (per-tooth ablation: every tooth A–E fired at its
own named assert; the decorative-inverse rig was caught at RUNTIME by the `delta_probe` in BOTH engines; the arbiter two-route
agreement confirmed by `simplify(N0_den_py − N0_den_wl) == 0`; no vacuous tooth). ZERO remediation required — the
contract-clean lightest-reshape stage, like stage018/012.
