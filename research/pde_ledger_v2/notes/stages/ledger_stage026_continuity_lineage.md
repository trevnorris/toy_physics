# ledger_stage026 — the continuity-lineage token-check (Check II-P3)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative-port sector, Cluster C). The CONTINUITY-LINEAGE leg of
a 4-way split of `pathA_43`: this stage carries the **continuity-lineage validation (3/4) of the joint `DENSITY_PORT_HOSTED`**
— it proves the ℓ=2 moment `I25` that sources stage 024's density-port numerator `N0_den` is a genuine `∫Y₂*·S_leak`
continuity moment DESCENDED from pathA_29's projected-continuity operator (the SAME operator that produced pathA_28's ℓ0
`M0` / ℓ1 `D1`), not a self-asserted / relabeled tag. This is the leg that **DISCHARGES two forward references**: it EARNS
`moment_valid=True` (that stage 025 consumed as a typed input) and validates `I25`'s lineage (that stage 024 cited as a typed
input). The DERIVATION of `N0_den` is **stage 024** (II-P1); the vector-freedom taint proof is **stage 025** (II-P2); the
able-to-fail port checks + closure overlay (the COMPLETING leg that LANDS the joint) is **stage 027** (II-P4).

**Verdict.** LOCAL `CONTINUITY_LINEAGE_EARNED` (the lineage-earning token = G1 ∧ G2, exit 0) + a printed JOINT PARTIAL
`DENSITY_PORT_HOSTED (3/4)`. ⚠ 026 does NOT own the joint verdict — that LANDS at 027 (the COMPLETING leg). 026 is a PARTIAL
leg (the 018/016/022/024/025 pattern).

**Status.** Exact symbolic / lexical-token set-algebra (float-free): the lineage gate is an EXACT-token subset over the
ℓ0→ℓ1→ℓ2 moment strings; the earning gate is a symbol-identity distinction; `expect_bool` residual asserts, no
`scipy`/`numpy`/floats/tolerances. Dual-engine SymPy **27 PASS** / Mathematica **27 PASS**, both exit 0, CWD-independent (repo
root + `/tmp`), zero file I/O.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}` (the 026 slice
> = the `CONTINUITY_*` tokens L192–197, `contains_all` L305–307, `continuity_lineage_valid` L310–322,
> `continuity_moment_symbol` [the `I25`-vs-`I_wrong2` gate] L325–332, `lineage_for` L392–417, the `continuity_dependency_ok`
> conjunct L585–592) + `reports/pathA_43_density_quadrupole_port.md` (the JOINT `DENSITY_PORT_HOSTED`; the two-rig history) +
> the original directive (§10 the two-rig changelog, NIT-1 = the self-asserted-tag Rig 2). ⚠ pathA_43 is CONTRACT-CLEAN —
> both engines were already standalone print-only with zero file I/O, so there is NO bridge to sever; the reshape was
> DECOMPOSITION + the EXACT-TOKEN genuineness UPGRADE + the `.wl` RE-AUTHOR (§5). The report/directive are cited for
> provenance only; the proof below is self-contained.

---

## 1. What this stage earns

stage 024 DERIVED the ℓ=2 density-port numerator `N0_den` and CITED `I25` (the ℓ=2 continuity moment symbol) as a typed
input. stage 025 PROVED `N0_den` vector-free but CONDITIONAL on a typed `moment_valid=True`. **026 EARNS both** — it runs the
COMPUTED lineage token-check that makes `moment_valid=True` genuine (025's forward ref) and validates `I25`'s `∫Y₂*·S_leak`
ℓ=2 lineage (024's forward ref). 024 is the EARNED-FIRST derivation, 025 is the anti-rig proof over `free_symbols`, **026 is
the ANCESTRY validation that ties the port's continuity moment to pathA_29's projected-continuity operator.**

**Physical picture** (`docs/conceptual_foundation.md` §3/§5): the density port's ℓ=2 source coupling `g_base =
√rho_eff·c_s²·I25·Ξ_Q/a^(7/2)` carries the ℓ=2 continuity moment `I25`. pathA_29's projected-continuity operator produces the
mass-rate ℓ0 moment `M0 = ∫S_leak d³x` and the dipole/momentum-rate ℓ1 moment `D1_i = ∫x_i S_leak d³x + ∫j_i d³x` (pathA_28's
exports; the return-cancellation targets `R0=−M0`, `R1=−D1`). The SAME operator, projected onto the ℓ=2 spherical harmonic
`Y₂`, gives the ℓ=2 moment `Q2_m = ∫Y₂*·S_leak d³x` — whose reduced value is `I25`. **026 proves the port's ℓ=2 moment is the
structured continuation of that ℓ0→ℓ1→ℓ2 ancestry.**

### 1.1 The subject: 024's `N0_den`, cited not re-derived
026 CITES stage 024's factored export (it does NOT re-derive the 2×2 inverse — that would blur the 024/026 cut):
```
N0_den = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)² / ( a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)² ).
```
The **consumption-integrity oracle** (tooth I) computes `N0_den.free_symbols` off the actual cited expression and asserts it
equals the exported **10-symbol host contract** `{a, c_s, rho_eff, I25, Ξ_Q, η_q, η_φ, ϖ_q2, ϖ_Φ2, λ_c}`. ⭐ `I25` enters as
the EXTERNAL `I25²` factor — dropping/renaming it removes `I25` (the drop-`I25` mutation routes here); dropping `η_q·λ_c`
removes only `η_q` (`λ_c` survives in the denominator). The earned moment `I25` must appear in `N0_den.free_symbols` — the
tie between "the validated lineage" and "the actual port's moment."

### 1.2 The ℓ0→ℓ1→ℓ2 continuity lineage
The lineage is the structured record that `I25` descends from pathA_29's projected-continuity operator:
```
operator_id : pathA_29_projected_continuity
  ℓ0  M0    = Integral( S_leak d3x )                            (pathA_28 mass rate;   the return target R0=−M0)
  ℓ1  D1_i  = Integral( x_i*S_leak d3x ) + Integral( j_i d3x )  (pathA_28 dipole rate; the return target R1=−D1)
  ℓ2  Q2_m  = Integral( Y2_m_star*S_leak d3x )                  (the ℓ=2 continuation;  reduced value = I25)
  ℓ2 kernel : Y2_m_star*S_leak
```
The ℓ0/ℓ1 tokens ARE pathA_28's exported moments (stage 008, PROVENANCE); the operator is pathA_29's (stage 009/010,
PROVENANCE). 026 does not reconstruct 008/009/010's exports (that would over-engineer the consumption) — it verifies the
lineage EXHIBITS the ancestor tokens.

### 1.3 The DECISIVE gate G1 — the EXACT lexical-token check (the genuineness UPGRADE)
`continuity_lineage_valid` COMPUTES `operator_id == "pathA_29_projected_continuity"` AND (each ancestry level's required token
SET ⊆ the string's actual token set): l0 ⊇ {`M0`,`Integral`,`S_leak`,`d3x`} ∧ l1 ⊇ {`D1_i`,`Integral`,`x_i`,`S_leak`,`j_i`,
`d3x`} ∧ l2 ⊇ {`Q2_m`,`Integral`,`Y2_m_star`,`S_leak`,`d3x`} ∧ l2_kernel ⊇ {`Y2_m_star`,`S_leak`}. ⭐⭐ **The reshape UPGRADES
the source's raw SUBSTRING `contains_all` to an EXACT LEXICAL-TOKEN check** — the source's `token in text` is defeated by a
token-STUFFED forgery (`NOT_M0`⊃`M0`, `FakeIntegral`⊃`Integral`, `S_leakage`⊃`S_leak`, `d3xyz`⊃`d3x`). ⚠ **The token
alphabet is LOAD-BEARING: identifier tokens `[A-Za-z0-9_]+`** (Python `re.findall(r'[A-Za-z0-9_]+', s)`; Mathematica
`StringCases[s, (WordCharacter | "_")..]`, NOT bare `WordCharacter..` — which EXCLUDES `_` and would FRAGMENT the genuine
tokens `S_leak`/`Y2_m_star`/`D1_i`/`Q2_m`/`x_i`/`j_i` and false-fail the genuine lineage). The genuine `CONTINUITY_L*` strings
yield the required tokens INTACT (baseline PASSES); the stuffed forgery does NOT (the token-stuffing control FAILS). ⭐ The
check reads tokens OUT of the actual moment strings — **NEVER a self-asserted `valid` flag** (§1.6).

### 1.4 The DECISIVE gate G2 — the `I25`-vs-`I_wrong2` earning gate (makes the lineage LOAD-BEARING)
Source-faithful: `moment_valid ≡ lineage_valid`; `earned_moment = I25 iff (lineage_valid ∧ a_power == −7/2) else I_wrong2`.
On the consumed `−7/2` baseline the genuine lineage earns `(I25, moment_valid=True)`; an invalid/stuffed lineage earns
`(I_wrong2, moment_valid=False)`. ⭐ **An INVALID lineage substitutes the DIFFERENT symbol `I_wrong2` (dim `(2,0,0)` ≠ `I25`'s
`L^(5/2)`) into the port**, so a mis-lineaged port is a different object — this is what makes the token-check MATTER (not a
string decoration). G1 ∧ G2 EARN the two forward-referenced outputs: `moment_valid` (025's) and the validated `I25` (024's).
⚠ The a-power `−7/2` is a seam with 027: a valid lineage with a-power ≠ −7/2 yields `(I_wrong2, moment_valid=True)` and lands
`FAIL_PORT_MALFORMED(scaling)` at **027** — 026 uses the `−7/2` baseline only and does NOT run the wrong-power case against
the fixed baseline `N0_den`.

### 1.5 The dependency binding — first arm de-counted, OR-arm the mandatory tooth G
`continuity_dependency_ok = lineage_valid ∧ moment_valid ∧ (earned_moment ∈ N0_den.free_symbols  OR  (compact(expr)==0 ∧
coupling_zero))`. ⚠ The FIRST-ARM membership witness `earned_moment ∈ N0_den.free_symbols` is a **DE-COUNTED redundant
witness** (printed, not tallied): because the port is consumed from 024 with a FIXED contract that includes `I25` and the
baseline lineage earns `I25`, any mutation making `I25 ∉ N0_den.free_symbols` also breaks the consumption-integrity contract
oracle (tooth I) and any lineage corruption is caught by G1 first — so it cannot independently fail (subsumed by G1 + tooth
I; the drop-`I25` mutation routes to tooth I). The SECOND OR-arm `(compact(expr)==0 ∧ coupling_zero)` is a genuine mandatory
tooth **G**, exercised by 4 ISOLATED probes on synthetic `(expr, coupling_zero)` inputs (so tooth I does not subsume it).

### 1.6 Why the check must be COMPUTED over tokens, never a self-asserted flag (Rig 2)
The pathA_43 gate history holds **Rig 2 — the self-asserted `continuity_interface` tag** (the "material" rig, caught by the
GLM tertiary at the source): a devious executor writes a NEW equation ≡ the retired vector port `P_A = Ω_U²·g_W + R_mix·g_U`,
tags it `continuity_interface`, and wires it as `N0`'s ancestor — the origin check passed by fiat. The earned mechanism 026
preserves: the validator checks structural TOKENS (the ℓ0→ℓ1→ℓ2 ancestry), NOT a shared `valid` flag. ⭐ **The rig battery
SELF-ENFORCES this**: since the source lineage dicts have no `valid` field, every NEGATIVE lineage here carries a DECOY
`"valid": True` + a self-asserted `continuity_interface` claim, and is passed through the SAME validator (no config-name /
expected-verdict lookup is verdict-bearing). The coupling meta-test — substitute a flag-only validator (reads
`lineage["valid"]`) — makes the decoy-flagged negatives wrongly PASS, so the audit FAILS: the exact-token check is proven
load-bearing. ⚠ Do NOT enforce the lineage with a source grep for tokens ([[feedback-grep-acceptance-dodgeable]]) — the real
property is the COMPUTED exact-token check over the actual moment strings.

---

## 2. The able-to-fail battery (026-owned; per-tooth ablation — each fired at its OWN named assert, then went vacuous when neutered)

| tooth | rig / mutation → routed assert | notes |
|---|---|---|
| A fake_continuity | operator_id `mis_tagged_vector_formula` + decoy `valid=True` + vector-port expr → the **G1 lineage-gate** assert (operator_id + tokens absent) fires → `(I_wrong2, False)` | ⭐ Rig 2; a flag-read PASSES it |
| B attack2 | correct operator_id + decoy `valid=True` + BAD ℓ2 (`GARBAGE…`) → the **G1 ℓ2-token** assert fires | the harder self-asserted case |
| C operator_id / l0 / l1 / l2_kernel token ablations | neuter each ancestry conjunct (each with a decoy `valid=True`) → the corresponding **G1 per-level token** assert fires | ⭐ every conjunct of the ℓ0→ℓ1→ℓ2 check independently load-bearing |
| C token-STUFFING forgery | `NOT_M0`/`FakeIntegral`/`S_leakage`/`d3xyz` (PASSES a substring check) → the **exact-token** assert fires; reverting to substring → the tooth stops firing | ⭐ the genuineness-upgrade meta-test |
| D earning gate (G2) | neuter the substitution (always `I25`) → the **G2 distinguishes** assert fires (a corrupted lineage wrongly earns `I25`) | the load-bearing "makes the lineage matter" tooth |
| F baseline-valid POSITIVE (able-to-PASS) | the genuine lineage PASSES; corrupt ANY ancestry token → FLIPS to FAIL | the anti-over-rejection / reversibility control |
| G(i–iv) vanished-coupling OR-arm (mandatory) | 4 ISOLATED synthetic probes: `expr==0∧cz=True`→PASS via arm 2; `expr==0∧cz=False`→FAIL; nonzero missing `earned_moment`+`cz=True`→FAIL; neuter arm 2 → probe (i) FAILs | isolated from tooth I's fixed `N0_den` |
| I subject_integrity | corrupt the cited `N0_den` — drop the external `I25²` factor (removes `I25`, routes HERE) / drop `η_q·λ_c` (removes only `η_q`) / replace `rho_eff` → `free_symbols ≠ HOST_CONTRACT` → the **exact host-contract** assert fires | the consumption oracle |
| H′ arity | a wrong call arity → the def/call **arity scanner** fires | the stage007 lesson |
| H′ leakage | a leaked unevaluated authored-helper head → the **unevaluated-leakage scanner** fires | — |

Each rig is exercised by a coupling meta-test: it must fire at its OWN named assert AND, when neutered, stop firing. ⚠ Two
checks are DELIBERATELY NOT counted teeth: the **first-arm dependency witness** `earned_moment ∈ N0_den.free_symbols` (§1.5, a
de-counted subsumed witness) and the **`.wl` algorithm independence** (a TRANSLITERATION-REVIEW property — a copied
`StringContainsQ`+`If` `.wl` can preserve every value, so no runtime assert rejects it; kept in the Codex→Grok→Codex screen).
The tri-review de-counted one further subsumed tally (`baseline continuity_dependency_ok is True`, §5).

---

## 3. Honest scope

- **026 validates the LINEAGE, not the magnitude.** It proves `I25` is a genuine `∫Y₂*·S_leak` ℓ=2 continuity moment descended
  from pathA_29's operator; the literal moment VALUE (`I25`, `Ξ_Q`, `λ̂_Q`, `rho_eff` throat values) remains SIM_DEFERRED
  (024/027, Gate-6), and `G`/`2/5`/`54/5` are CALIBRATED (027).
- **026 discharges two forward references.** It EARNS `moment_valid=True` (stage 025's LOCAL `DENSITY_PORT_VECTOR_FREE` was
  CONDITIONAL on this) and validates `I25`'s lineage (stage 024 CITED `I25` as a typed input). The register `I25` row's
  "lineage → 026" forward ref is now DISCHARGED (the LINEAGE is earned; the magnitude stays `GAP`/SIM_DEFERRED).
- **The joint + the retirement are 027's, not 026's.** 026 proves the continuity-lineage CONJUNCT; the joint
  `DENSITY_PORT_HOSTED` (024 derivation ∧ 025 vector-freedom ∧ 026 lineage ∧ 027 dim/scaling/sign/closure) and the EM
  `A_w`/`U,W` scaffold retirement land at 027 (which also ASSEMBLES the full `origin_ok` from the 025 + 026 certificates).
  026 asserts ONLY the lineage predicates + the baseline dependency witness — it does NOT compute `tags` /
  `source_map_complete` (025's) or a full `origin_ok` (027's).

---

## 4. Consumed / exported

- **Consumed — PROVENANCE + the one checkable subject contract.**
  - **stage 024's `N0_den`** — the SUBJECT (cited factored export; the consumption-integrity oracle asserts its computed
    `free_symbols == HOST_CONTRACT`, the one checkable cross-stage constraint; `rho`=`rho_eff`).
  - **stage 008's `M0`/`D1`** — the ℓ0/ℓ1 ancestor moments (the `CONTINUITY_L0`/`L1` token strings). PROVENANCE.
  - **stage 009/010's projected-continuity OPERATOR** (`CONTINUITY_OPERATOR_ID`). PROVENANCE.
  - **stage 005 (`c_s`) + `a` CONV** — units carriers (in the consumed `N0_den`; 026 runs no dim gate).
- **Register.** ZERO new counted CALIB knobs (a LINEAGE/proof slice — it consumes 024's `N0_den` + 008/009/010 provenance,
  introduces NO new physical symbols; `I25` is already tracked [SIM_DEFERRED magnitude], `I_wrong2` is a control-scar,
  tracked-not-counted like stage023's `q_free`). Part-II counted CALIB set unchanged at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}`
  (015) + `{T_Ω, β₂}` (017) = **6**. New structural edge **R45** — the density-port continuity-lineage validation (`I25` is a
  genuine `∫Y₂*·S_leak` ℓ=2 continuity moment descended from pathA_29's operator, COMPUTED over the ℓ0→ℓ1→ℓ2 exact lexical
  tokens + the `I25`-vs-`I_wrong2` earning gate), scoped as VALIDATION of the cited lineage certificate (NOT a re-derivation
  of the upstream 008/009/010 operator); discharges nothing at the KNOB level (a proof/provenance edge, like R43/R44). No new
  dims (026 runs no dim gate — `[N0_den]=L⁻¹M` / `[I25]=L^(5/2)` are 027's / provenance-recorded at 024).
- **Exported.** ⭐ `moment_valid=True` (DISCHARGES the 025 forward ref) + the validated `I25` (DISCHARGES the 024 forward ref)
  + the lineage certificate → **027** (as one conjunct of the ASSEMBLED joint `DENSITY_PORT_HOSTED`). 026 does NOT export
  `tags` / `source_map_complete` / a full `origin_ok`.

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O — and pathA_43 was ALREADY contract-clean, so the reshape
"sever the bridge" step is a NO-OP; the work was DECOMPOSITION + the EXACT-TOKEN genuineness UPGRADE + the `.wl` RE-AUTHOR +
the LOCAL/PARTIAL verdict framing. ⚠ **Like stage 025 (and unlike stage 024's keep-native derivation lane), stage 026
RE-AUTHORS its `.wl`** — the source `.wl` lineage machinery (`containsAll`/`lineageValidQ`/`continuityMomentSymbol`/
`lineageFor`) was a near line-for-line transliteration of the `.py`, so the `MATHEMATICA_MIRROR_POLICY` required a genuinely
independent computation: the decisive gates use a native `Graph` + `FindPath` + `FoldList` ℓ0→ℓ1→ℓ2 ancestry walk whose
terminal ℓ2 node PRODUCES the earned moment `I25` only after the full walk validates (with `StringCases[s,(WordCharacter|"_")
..]` + `SubsetQ` exact-token subsets and a `Cases`-over-`HoldComplete` free-symbol extractor) — materially different from the
`.py`'s flat token-set conjunction. Agreement is transcript-level (both emit the cited `N0_den` + `free_symbols == contract`,
the lineage dict, `continuity_lineage_valid=True`, the earned `(I25, moment_valid=True)`, the per-rig outcomes, the LOCAL
verdict + JOINT PARTIAL); neither engine reads the other. SymPy 27 PASS / Mathematica 27 PASS, both exit 0, CWD-independent.

**Directive review** used the Codex→Grok→Codex bookend, which caught a directive-level error at EVERY leg. The Codex xhigh
design-review folded **5 BLOCKING** — most importantly (1) the source's `contains_all` is raw SUBSTRING, defeated by a
token-stuffed forgery → the reshape UPGRADES to EXACT LEXICAL-TOKEN semantics + a required stuffing control; (2) the negatives
must carry a DECOY `valid=True` + pass through the SAME validator (the source dicts have no `valid` field) + a coupling
meta-test; (3) the dependency binding is LOGICALLY SUBSUMED (no non-subsumed mutation exists) → de-counted witness, drop-`I25`
→ tooth I; (4) source-faithful moment semantics `moment_valid ≡ lineage_valid` (a valid lineage + wrong a-power gives
`(I_wrong2, moment_valid=True)`, 027's scaling); (5) the `continuity_interface ∈ tags` conjunct is 025's tag machinery (026 is
tag-free; the tag is also from `Ξ_Q`/`η_q`/`η_φ`/`λ_c`) → removed from 026, `origin_ok` assembly → 027. A **Grok-4.5
compute-verify** pass then caught **1 BLOCKING** — the exact-token ALPHABET is load-bearing: a naïve Mathematica
`WordCharacter` split EXCLUDES `_` and would fragment the genuine tokens and false-fail the lineage → pinned
`[A-Za-z0-9_]+` (Python `\w`; Mathematica `(WordCharacter|"_")..`). A Codex confirm-pass caught **1 BLOCKING** (tooth G left
optional → made mandatory with 4 isolated probes) + 2 nits; a Codex final-confirm caught **1 BLOCKING** (the whole-predicate
de-count contradicted the now-mandatory tooth G → scoped the de-count to the first-arm witness only); a Codex final-confirm
v5 returned `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read re-derived `N0_den` from stage024's 2×2 inverse, DIFF=0;
verified the exact-token tokenizer keeps the genuine tokens intact while rejecting the stuffed forgery, the source-faithful
G2 semantics, no 024/025/027 leakage, the de-counted first-arm witness, and the genuine independent `Graph`/`FindPath`/
`FoldList` `.wl` route) + `ADVERSARIAL_CLEAN` (per-tooth ablation on BOTH engines — 22 mutations, every one fired at its OWN
named assert and went vacuous when neutered; the three decisive proofs held: reverting to substring makes the stuffing
forgery wrongly pass, a flag-only validator makes the decoy-flagged negatives wrongly pass, and bare-`WordCharacter`
false-fails the genuine lineage; the H′ arity/leakage scanners fire natively; the `.wl` is genuinely independent). The
adversarial leg flagged 1 [minor] non-blocking subsumed tally (`baseline continuity_dependency_ok is True` — for the baseline
its content reduces to conjuncts already tallied at G1/G2 plus the de-counted membership witness; no mutation makes only it
fail) → DE-COUNTED to a labeled diagnostic print (tally 28→27 per engine, no other check touched, Codex's mutation proof:
forcing it False keeps exit 0 + 27/27), then fresh-agent `REVERIFY_CLEAN` (the diff is exactly the de-count; the earned teeth
unchanged and still able-to-fail; both engines 27/27).
