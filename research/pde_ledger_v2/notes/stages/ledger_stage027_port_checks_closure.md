# ledger_stage027 — the density-port checks + closure (Check II-P4)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative-port sector, Cluster C). The COMPLETING leg of a 4-way
split of `pathA_43`: this stage carries the **port checks + closure (4/4) and OWNS the joint `DENSITY_PORT_HOSTED`** — it runs
the 6 able-to-fail port checks + the K̄ closure over the density-port numerator `N0_den` (derived at stage 024, proven
vector-free at stage 025, continuity-lineaged at stage 026), ASSEMBLES the joint `origin_ok`/`DENSITY_PORT_HOSTED` from the
025 (vector-freedom) + 026 (continuity-lineage) certificates + its own checks, and RECORDS the EM `A_w`/`U,W` vector-scaffold
RETIREMENT with the completed joint. The DERIVATION is **stage 024** (II-P1); the vector-freedom taint proof is **stage 025**
(II-P2); the continuity-lineage token-check is **stage 026** (II-P3).

**Verdict.** The JOINT `DENSITY_PORT_HOSTED` (the exit-0 gate; 027 OWNS it — the stage023 pattern, UNLIKE the PARTIAL legs
024/025/026 which printed a 1/4·2/4·3/4 partial). ⚠⚠ **COMPLETE ≠ PASS: the joint stays CALIBRATED not PASS** —
`G=GENUINE_BLOCKED` (020's provenance; `54/5`/`Γ̄₅`=`external_bridge_input`, the `27`=018's `derived_in_gate`); the MAGNITUDE
is SIM_DEFERRED, only the STRUCTURE is hosted.

**Status.** Exact symbolic (dimension-vector `(L,M,T)` algebra, spherical-Hankel series, exact rationals; `expect_bool`/
`expect_zero` residuals, no `scipy`/`numpy`/floats/tolerances). Dual-engine SymPy **55 PASS** / Mathematica **55 PASS**, both
exit 0, CWD-independent (repo root + `/tmp`), zero file I/O.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}` (the 027 slice
> = the dim engine `dim_of`/`scale_power`/`BASE_DIMS`/`BASE_A_POWERS` L82–268, `dtn_sign` L521–534, `closure_overlay`
> L701–722, `evaluate`/the verdict chain L537–698, `assert_gate` L773–803, `controls` L725–770). ⚠ pathA_43 is
> CONTRACT-CLEAN — both engines were already standalone print-only with zero file I/O, so there is NO bridge to sever; the
> reshape was DECOMPOSITION + the joint-owning verdict framing + the `.wl` RE-AUTHOR. The report/directive are cited for
> provenance only; the proof below is self-contained.

---

## 1. What this stage earns

Stage 024 DERIVED the ℓ=2 density-port numerator `N0_den`; stage 025 PROVED it vector-free; stage 026 EARNED the continuity
lineage of its ℓ=2 moment `I25`. **027 verifies the port is well-formed as a radiative port** — it has the right dimension,
a-scaling, radiative sign, is nonzero, and its magnitude closes the 2.5PN Burke–Thorne moment consistency — then ASSEMBLES
the joint verdict and records the scaffold retirement.

**Physical picture** (`docs/conceptual_foundation.md` §3/§5): gravity's quadrupole radiation rides the medium's OWN density
ripple at `c_s` (the "speed gravitational changes propagate"), NOT an EM vector field. The moving throat's ℓ=2 coupling to
the outgoing quadrupole `c_s` wave is the density two-port `N0_den` (024). 027 confirms it is a genuine radiative port on the
density mode — which is what RETIRES the old EM `A_w`/`U,W` vector scaffold (the pathA_43 diagnostic sliver).

### 1.1 The subject: 024's `N0_den`, cited not re-derived
027 CITES stage 024's factored export (it does not re-run the 2×2 inverse):
```
N0_den = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)² / ( a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)² ).
```
The **consumption-integrity oracle** computes `N0_den.free_symbols` off the actual cited expression and asserts it equals the
exported **10-symbol host contract** `{a, c_s, rho_eff, I25, Ξ_Q, η_q, η_φ, ϖ_q2, ϖ_Φ2, λ_c}` (`{q2, Phi2}` are coordinate
metadata, NOT in the numerator). ⚠ `rho = rho_eff` (reduced-3D mass density, not stage005's `ρ0`).

### 1.2 The 6 able-to-fail port checks
1. **DIM** — `[N0_den] == (−1,1,0) = L⁻¹M` in the `(L,M,T)` tuple convention (via the `dim_of`/`BASE_DIMS` dim engine over the
   sourced port dims `[I25]=(5/2,0,0)`, `[c_s]=(1,0,−1)`, `[rho_eff]=(−3,1,0)`, `[a]=(1,0,0)`, `[ϖ_q2]=[ϖ_Φ2]=[λ_c]=(0,0,−2)`).
   ⚠ The `(L,M,T)` order differs from the register header's `[L,T,M]` — a cross-convention transfer trap (the stage016
   lesson); the check is self-contained in its convention. The μ̂₀-free discipline (021): the gate reads the SOURCED port
   dims, not a back-solved carrier — corrupt-`[N₀]` fires `FAIL_PORT_MALFORMED(dimensional)`; corrupt-`[G]` is NO_FAIL
   (`G ∉ N0_den.free_symbols`, a genuine scope diagnostic wired as a real mutant re-run).
2. **SCALE** — `P0_physical=(c_s/a)²·N0_den/D0` has a-power `−5` (021's μ̂₀-free `[P0_phys]=1` machinery). ⭐ The `scaling`
   control (a-power −3) puts `I_wrong2` (dim `(2,0,0)`) + `g_base ∝ a⁻³` into the port, which DELIBERATELY compensates in the
   L-slot (`2−3 = 5/2−7/2 = −1`), so `dim_ok` STAYS True and ONLY `scaling_wrong` fires → **`FAIL_PORT_MALFORMED(scaling)`
   single-tag** (not `(dimensional,scaling)`).
3. **SIGN** — the outgoing DtN `+i z⁵/27`: `coeff(z⁵)/i == 1/27` (χ_Q=+1), COMPUTED from the spherical-Hankel `ŷ₂ = −3/(z·h₂'/h₂)`
   with `h₂=j₂+i·y₂`; the `incoming_sign` control (`h₂=j₂−i·y₂`) → `−1/27` (χ_Q=−1) → `FAIL_PORT_MALFORMED(sign)`. Not a typed χ_Q.
4. **NONZERO** — `N0_den ≠ 0`; the `zero_coupling` control (`g_base=0 → N0_den=0`) → `FAIL_PORT_VANISHES` (via the 027-LOCAL
   subject-binding OR-arm keeping `origin_ok` True — §1.4).
5. **DEFERRED** — the honest SIM-inconclusive branch: `deferred_uncertified` (`Xi_deferred` a-power undecidable) →
   `PORT_INCONCLUSIVE_SIM_DEFERRED` (neither FAIL nor spurious PASS); the converse `proven_deferred` → `DENSITY_PORT_HOSTED`.
6. **CLOSURE** — the K̄ overlay: `K̄₄ − 4K̄₂²/K̄₀ == 0` ∧ `Γ̄₅ − 2G/(5c⁵) == 0` (both STANDALONE asserts), with
   `K̄₀=54Gc_s⁵/(5a⁵c⁵)`, `K̄₂=6Gc_s³/(5a³c⁵)`, `K̄₄=8Gc_s/(15ac⁵)`, `Γ̄₅=K̄₀·a⁵/(27c_s⁵)`. Both residuals are identically 0
   (`4K̄₂²/K̄₀ = (8/15)Gc_s/(ac⁵) = K̄₄`; `Γ̄₅ = 54G/(5·27·c⁵) = 2G/(5c⁵)` = Burke–Thorne). ⚠ This is the A3 2.5PN CONSISTENCY,
   SHARED with stage 028's full match-back (marked shared, not double-counted); it is a consistency over the CALIBRATED
   moments (`G=GENUINE_BLOCKED`), NOT a first-principles `Γ̄₅`/`G` derivation.

### 1.3 The a-power `−7/2` seam with 026 (027 owns the SCALING, 026 owns the EARNING)
026 owns the `I25`-vs-`I_wrong2` earning gate (which SYMBOL a lineage earns); 027 owns the SCALING VERDICT (`p0_power == −5`).
The `scaling` control's port (a-power −3) is built by SPECIALIZING the cited factored `N0_den` export (substituting
`g_base ∝ a⁻³`, moment `I_wrong2`) — NOT by re-running 024's inverse. 027's subject-binding membership uses the CURRENT
fixture moment (`port_moment = I25` at the `−7/2` baseline, `I_wrong2` under the scaling fixture) — it does NOT re-run 026's
earning gate.

### 1.4 The joint assembly (027 OWNS `DENSITY_PORT_HOSTED`) — from atomic sibling certificates + a 027-local subject binding
027 CONSUMES the ATOMIC sibling facts (NOT a sibling's already-de-counted composite predicate) and assembles the joint
LOCALLY:
- **025's atomic vector-freedom facts:** the taint set (→ `continuity_interface ∈ tags` ∧ `vector_port ∉ tags`),
  `¬vector_host_symbols`, `source_map_complete`, the `vector_free` (P2) certificate. ⚠ The tags are **025's** (026 does not
  export tags); consumed over the BASELINE `N0_den`, not live-retainted (under `zero_coupling` the live taint of `N0_den=0`
  is empty).
- **026's atomic continuity-lineage facts:** `moment_valid=True`, the validated `I25` symbol, `lineage_valid`
  (↔ 026's `lineage_certificate=PASS`). ⚠ 027 does NOT consume 026's flat `continuity_dependency_ok`.
- **027's LOCAL subject binding:** `subject_binding = (port_moment ∈ N0_den.free_symbols) OR (compact(N0_den)==0 ∧
  coupling_zero)`. ⭐ The vanished-coupling OR-arm is 027-OWNED (it SUPERSEDES the source's outer
  `vanished_continuity_coupling`) — this is what makes `zero_coupling` fail at `nonzero_ok` (`FAIL_PORT_VANISHES`) NOT at
  `origin_ok`; the meta-test (disable the OR-arm → `zero_coupling` flips to `FAIL_NOT_DENSITY_DERIVED`) proves it load-bearing.
- `origin_ok = (continuity_interface ∈ tags) ∧ (vector_port ∉ tags) ∧ (¬vector_host_symbols) ∧ source_map_complete ∧
  vector_free ∧ (lineage_valid ∧ moment_valid) ∧ subject_binding`; the joint verdict via the chain: `¬origin_ok →
  FAIL_NOT_DENSITY_DERIVED`; `¬nonzero_ok → FAIL_PORT_VANISHES`; malformed → `FAIL_PORT_MALFORMED(<bad list>)`; all-true incl
  `closure_ok → DENSITY_PORT_HOSTED`; else → `PORT_INCONCLUSIVE_SIM_DEFERRED`.

### 1.5 The EM-scaffold retirement (CONDITIONAL on the joint landing)
With the joint `DENSITY_PORT_HOSTED`, 027 RECORDS the EM `A_w`/`U,W` vector-scaffold RETIREMENT + the pathA_43 diagnostic
sliver CLOSED. ⚠ CONDITIONAL: R43/R44 explicitly DEFERRED the retirement to this joint — a `FAIL_*` /
`PORT_INCONCLUSIVE_SIM_DEFERRED` does NOT record it (the port has not displaced the scaffold; the sliver reopens). The
retirement tooth proves this (force a FAIL → not recorded).

---

## 2. The able-to-fail battery (027-owned; per-tooth ablation — each fired at its OWN named assert, then went vacuous when neutered)

| tooth | rig / mutation → routed assert | notes |
|---|---|---|
| DIM | corrupt `[N0_den]` → the DIM `[N0_den]=L⁻¹M` assert → `FAIL_PORT_MALFORMED(dimensional)` | μ̂₀-free; corrupt-`[G]` NO_FAIL is a genuine dim re-run (G∉free_symbols) |
| SCALE | a-power −3 (moment `I_wrong2`, `g_base∝a⁻³`) → the SCALE `p0_power=−5` assert → **single-tag** `FAIL_PORT_MALFORMED(scaling)` | the L-slot compensation keeps `dim_ok` True |
| SIGN | incoming `j₂−i·y₂` → the SIGN `coeff/i=1/27` assert → `FAIL_PORT_MALFORMED(sign)` | χ_Q DERIVED, no typed-χ_Q tautology |
| NONZERO | `coupling_zero` → the NONZERO assert → `FAIL_PORT_VANISHES` | via the 027-local OR-arm keeping `origin_ok` True |
| OR-arm meta | disable the subject-binding 2nd arm → `zero_coupling` flips `FAIL_PORT_VANISHES → FAIL_NOT_DENSITY_DERIVED` | proves the OR-arm 027-local + load-bearing |
| CLOSURE K̄₄ | `K̄₄+δ₄` (post-construction, NOT `K̄₀`) → the STANDALONE K̄₄ residual assert | isolated (Γ̄₅ residual stays 0) |
| CLOSURE Γ̄₅ | `Γ̄₅+δΓ` → the STANDALONE Γ̄₅ residual assert | the added standalone assert (source gated Γ̄₅ only via `closure_ok`) |
| DEFERRED | `deferred_uncertified` → INCONCLUSIVE; `proven_deferred` → HOSTED | both fire — honest, not baked toward hosted |
| ASSEMBLY (×7 atomic facts) | corrupt each 025/026 atomic fact one-at-a-time → its own `origin_ok` assert → `FAIL_NOT_DENSITY_DERIVED` | every conjunct individually load-bearing; `∧ vector_free` genuine |
| SUBJECT-INTEGRITY | replace a live `N0_den` symbol with a same-dim/a-power foreign symbol → ONLY the host-contract assert fires | dim/scaling stay valid (preservation) |
| POSITIVE | baseline → `DENSITY_PORT_HOSTED`; corrupt ANY of the 6 checks → FLIPS | able-to-PASS + able-to-flip |
| RETIREMENT-CONDITIONAL | force a FAIL → the retirement is NOT recorded | the retirement rides the landed joint |
| WL arity / leakage | a wrong `.wl` `Module` call arity / a leaked unevaluated head → the runtime scanners fire | the stage007 lesson |

Each rig is exercised by a coupling meta-test (fires at its OWN named assert AND, when neutered, stops firing). The
adversarial re-review ran the full per-tooth ablation matrix on BOTH engines and found no vacuous / subsumed / tautological
construct.

---

## 3. Honest scope

- **027 verifies the STRUCTURE, not the magnitude.** The joint `DENSITY_PORT_HOSTED` lands the density-hosted, vector-free,
  continuity-lineaged, dimensionally/radiatively/scaling well-formed, closure-consistent STRUCTURE; the MAGNITUDE stays
  CALIBRATED (`G=GENUINE_BLOCKED`; `Γ̄₅`/`54/5`=`external_bridge_input` = the GR Burke–Thorne bridge; the `27`=018's earned
  fingerprint) / SIM_DEFERRED (`Ξ_Q`/`λ̂_Q`/`rho_eff` throat values, Gate-6). **COMPLETE ≠ PASS.**
- **The K̄ closure is the A3 2.5PN boundary, SHARED with 028.** 027's closure = the two-residual CONSISTENCY
  (`Γ̄₅=2G/(5c⁵)`, `K̄₄=4K̄₂²/K̄₀`); stage 028 does the full INV1–INV5 match-back + the 11-mutation coherent-rescale matrix.
  Marked shared, not double-counted — 027 does NOT rebuild 028's INV5 anchors.
- **The EM-scaffold retirement is CONDITIONAL on the joint landing** (R43/R44 deferred it here). A failed port does NOT
  displace the scaffold.

---

## 4. Consumed / exported

- **Consumed — PROVENANCE + the checkable subject contract.**
  - **stage 024's `N0_den`** — the SUBJECT of the dim/scaling/nonzero checks (cited factored; the oracle asserts
    `free_symbols == HOST_CONTRACT`; `rho=rho_eff`).
  - **stage 025's atomic vector-freedom facts** (the taint set / `vector_port∉tags` / `vector_host_symbols=∅` /
    `source_map_complete` / `vector_free`).
  - **stage 026's atomic continuity-lineage facts** (`moment_valid=True` / validated `I25` / `lineage_valid`).
  - **stage 018's DtN fingerprint** (`+i z⁵/27`, `χ_Q=+1`/`−1`) — the sign gate. **stage 021's μ̂₀-free dim machinery**
    (`[P0_phys]=1`, a-power `−5`) — cited. **stage 020's `54/5`/`G=GENUINE_BLOCKED`** — the closure magnitude, cited.
  - **stage 005 (`c_s`) + `a` (CONV) + `c` (GR bridge) + `G` (`GENUINE_BLOCKED`) + `D0`** — units/closure carriers.
- **Register.** ZERO new counted CALIB knobs (a checks+closure/proof slice — the K̄ moments are calibrated functions of
  `{G,c_s,a,c}`; `c`/`D0`/`Xi_deferred`/`I_wrong2` not new knobs). Part-II counted CALIB set unchanged at `{μ_η, T_w, β}`
  (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = **6**. New structural edge **R46** — the density-port checks + closure +
  the JOINT landing + the scaffold-retirement record; discharges nothing at the knob level (a proof/provenance edge, like
  R43/R44/R45). ⭐ **Upgrades the `N0_den` row — the `[N0_den]=L⁻¹M` dim CHECK now LANDS here** (was PROVENANCE-only at 024).
- **Exported.** ⭐ the completed joint `DENSITY_PORT_HOSTED` (CALIBRATED) + the EM-scaffold RETIREMENT record (the diagnostic
  sliver CLOSED) + the K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` → **stage 028** (the 2.5PN match-back, the A3 boundary). ⭐ **027 CLOSES the
  pathA_43 density-port fold (024∧025∧026∧027).**

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O — and pathA_43 was ALREADY contract-clean, so the reshape
"sever the bridge" step is a NO-OP; the work was DECOMPOSITION + the joint-owning verdict framing + the `.wl` RE-AUTHOR. ⚠
**Like stages 025/026 (and unlike stage 024's keep-native derivation lane), stage 027 RE-AUTHORS its `.wl`** — the source
`.wl` 027 slice (`dimOf`/`scalePower`/`dtnSign`/`closureOverlay`/`evaluate`) is a one-for-one transliteration of the `.py`, so
the `MATHEMATICA_MIRROR_POLICY` required genuinely independent routes: the DtN via native `SphericalHankelH1`/`H2` +
`SeriesCoefficient` (vs the hand-built `j₂±i·y₂` + `Series`); the dimension/scaling via unit-monomial / weight-rescaling +
`Exponent` (vs the recursive `dim_of`/`scale_power` walker); the closure via dimensionless ratio invariants
(`K̄₄K̄₀/(4K̄₂²)−1`, `Γ̄₅/(2G/5c⁵)−1`) with `Together`/`Cancel` (vs the additive residual `ok`); the verdict via an ordered
`SelectFirst`/`Association` dispatch (vs the `if/elif` chain). Agreement is transcript-level (both emit the same 6 checks +
closure residuals + joint verdict + per-control verdicts); neither engine reads the other. SymPy 55 PASS / Mathematica 55
PASS, both exit 0, CWD-independent.

**Directive review** used the Codex→Grok→Codex bookend, which caught issues at EVERY leg. The Codex xhigh design-review folded
**4 BLOCKING** (the certificate mis-scoping + the masked vanished-coupling OR-arm → 027-local subject binding; a dedicated
SUBJECT-INTEGRITY tooth for the host-contract oracle; the `.wl` re-author of ALL verdict-bearing computations — the source
`.wl` is a transliteration, the 018/021 keep-native precedent does not transfer; the isolated post-construction K̄₄/Γ̄₅
closure mutations). A **Grok-4.5 compute-verify** pass then caught **1 BLOCKING** (the hardcoded `validated_I25` subject
binding broke the `scaling` control — under a-power −3 the port's moment is `I_wrong2`; the fix binds the current-fixture
`port_moment`) + 2 nits (the certificate tags are over the baseline, not live-retainted; `vector_independence_ok` is not a
025 export). A Codex confirm-pass caught **3 BLOCKING** consistency-sweep gaps (a stale `validated_I25` in §3; §2(d)
contradicting the OR-arm fold; the explicit `origin_ok` formulas omitting `∧ vector_free`) + 2 nits; a Codex final-confirm
caught 1 prose nit → `DIRECTIVE_CLEAN`. Grok/Codex independently compute-verified the a-power L-slot compensation
(single-tag scaling), `[N0_den]=(−1,1,0)`, both K̄ residuals `==0` (+ the isolation), and the DtN `±1/27`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read re-derived `N0_den` from stage024's 2×2 inverse [diff=0],
the DtN `±1/27` via SymPy's own `jn`/`yn`, both K̄ residuals `==0`, and the scaling single-tag; confirmed no dropped check
[coverage diff], the `.wl` genuinely independent function-by-function, and the consumption matching 024/025/026's actual
exports) + `ADVERSARIAL_CLEAN` (the full per-tooth ablation matrix on BOTH engines — every tooth fired at its OWN named assert
and went vacuous when neutered; the OR-arm meta-test, the single-tag scaling, both isolated closure residuals, each atomic
certificate corruption, the SUBJECT-INTEGRITY tooth, and the conditional retirement all confirmed load-bearing; no
vacuous / subsumed / tautological / stamped construct; the `.wl` genuinely independent, not a mirror). **ZERO remediation** —
the clean stage came through both legs clean (like stage 024). Symbolic per-tooth ablation, mutations on copies.
