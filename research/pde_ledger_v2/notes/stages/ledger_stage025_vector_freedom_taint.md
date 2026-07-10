# ledger_stage025 — the density-port vector-freedom taint proof (Check II-P2)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative-port sector, Cluster C). The ANTI-RIG PROOF leg of a
4-way split of `pathA_43`: this stage carries the **vector-freedom taint proof (2/4) of the joint `DENSITY_PORT_HOSTED`** —
it consumes stage 024's exported density-port numerator `N0_den` and PROVES it is computationally **vector-free** (its
COMPUTED provenance taint is the density ancestry only; the retired EM `A_w`/`U,W` vector scaffold cannot be hiding in it).
The DERIVATION of `N0_den` is **stage 024** (II-P1); the continuity-moment lineage token-check is **stage 026** (II-P3); the
able-to-fail port checks + closure overlay (the COMPLETING leg that LANDS the joint) is **stage 027** (II-P4).

**Verdict.** LOCAL `DENSITY_PORT_VECTOR_FREE` (the vector-freedom-proof token, exit 0; ⚠ CONDITIONAL on the cited
`moment_valid=True`, a typed forward reference stage 026 earns) + a printed JOINT PARTIAL `DENSITY_PORT_HOSTED (2/4)`. ⚠ 025
does NOT own the joint verdict — that LANDS at 027 (the COMPLETING leg). 025 is a PARTIAL leg (the 018/016/022/024 pattern).

**Status.** Exact symbolic / set-algebra (float-free): the taint is a UNION of provenance tags over the actual
`N0_den.free_symbols`; `expect_bool` residual asserts, no `scipy`/`numpy`/floats/tolerances. Dual-engine SymPy **18 PASS** /
Mathematica **18 PASS**, both exit 0, CWD-independent (repo root + `/tmp`), zero file I/O.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_43_density_quadrupole_port_{sympy.py,.wl}` (the 025 slice
> = `VECTOR_SYMBOLS` L179–190, `BASE_SOURCE_TAGS` L271–302, `taint_for_expr`/`source_graph_for_expr`/`vector_ablated_expr`
> L335–389, the verdict booleans L580–614) + `reports/pathA_43_density_quadrupole_port.md` (the JOINT `DENSITY_PORT_HOSTED`;
> the two-rig history `:16–19`, `:45–54`, `:99–106`) + the original directive (§0/§3A/§3B the vector-independence machinery,
> §10 the two-rig changelog). ⚠ pathA_43 is CONTRACT-CLEAN — both engines were already standalone print-only with zero file
> I/O, so there is NO bridge to sever; the reshape was DECOMPOSITION + the `.wl` RE-AUTHOR (§5). The report/directive are
> cited for provenance only; the proof below is self-contained.

---

## 1. What this stage earns

stage 024 DERIVED the ℓ=2 density-port numerator `N0_den` and EXHIBITED a density-only construction. **025 PROVES that
construction is genuinely vector-free** — that `N0_den` cannot be a disguised copy of the retired EM vector port
`P_A = Ω_U²·g_W + R_mix·g_U`. This is the computational half of the EM-vector-scaffold retirement (the retirement itself is
the JOINT 025∧027 result, recorded at 027).

**Physical picture** (`docs/conceptual_foundation.md` §3/§5): gravitational quadrupole radiation rides the medium's OWN
density/sound (`c_s`) ripple — NOT an EM vector field. The old `A_w`/`U,W` vector scaffold was a heavier bookkeeping device the
ledger does not carry; 024 re-hosted the ℓ=2 port on the density mode, and 025 proves nothing vector survives in `N0_den`.

### 1.1 The subject: 024's `N0_den`, cited not re-derived
025 CITES stage 024's factored export (it does NOT re-derive the 2×2 inverse — that would blur the 024/025 cut):
```
N0_den = I25²·Ξ_Q²·c_s⁴·rho_eff·(η_φ·ϖ_q2 + η_q·λ_c)² / ( a⁷·(λ_c² − ϖ_Φ2·ϖ_q2)² ).
```
The **consumption-integrity oracle** computes `N0_den.free_symbols` off the actual cited expression and asserts it equals the
exported **10-symbol host contract** `{a, c_s, rho_eff, I25, Ξ_Q, η_q, η_φ, ϖ_q2, ϖ_Φ2, λ_c}` (the coordinate hosts `q₂, Φ₂`
are projected out of the numerator — they extend the *allowable* density-host universe to 12 but are NOT in
`N0_den.free_symbols`). A dropped/renamed symbol in the cited subject breaks the contract (tooth I).

### 1.2 The provenance ledger (the tag map) and the COMPUTED taint
Each symbol carries a fixed provenance tag set (`BASE_SOURCE_TAGS`, faithful to the pathA_43 source with the 024 renames):
```
c_s, rho_eff, ϖ_Φ2   → {pathA_29_bulk}          ϖ_q2                → {pathA_32_wall}
Ξ_Q, η_q, η_φ        → {continuity_interface}    a                   → {pathA_29_bulk, pathA_32_wall}
λ_c                  → {continuity_interface, pathA_29_bulk, pathA_32_wall}   (the interface object joining the modes)
I25                  → {continuity_interface, pathA_32_wall}   (installed from the cited moment_valid input)
[vector + relabel + σ_hidden] → {vector_port}     free_carrier → ∅ (empty)
```
The **taint** is COMPUTED — the UNION of `tag_map[s]` over `s ∈ N0_den.free_symbols` — not typed. ⭐ The whole point: the
check is over PROVENANCE TAGS, never a `set(names) & FORBIDDEN` name-check.

### 1.3 The two separated predicates (the DECISIVE gate)
- **P1 `baseline_ancestry_ok`** — the genuine `N0_den`'s computed taint equals EXACTLY `{continuity_interface, pathA_29_bulk,
  pathA_32_wall}` (the baseline density ancestry). Asserted for the genuine subject ONLY.
- **P2 `vector_free`** — `source_map_complete` (no symbol missing from the tag map, no symbol with an empty tag set) AND
  `vector_host_symbols = N0_den.free_symbols ∩ VECTOR_SYMBOLS = ∅` AND `vector_port ∉ taint`. The GENERAL vector-freedom
  predicate every passing case satisfies.

`N0_den` satisfies BOTH → LOCAL `DENSITY_PORT_VECTOR_FREE = P1 ∧ P2 ∧ moment_valid`.

### 1.4 The expression ablation is a RETAINED WITNESS (de-counted, NOT the decisive gate)
The expression-level vector-ablation (substitute every vector-tainted symbol → 0, compare) is retained as a printed
consistency witness but is **NOT a counted gate**: once P2 holds (`vector_port ∉ taint` ∧ `∩ VECTOR_SYMBOLS = ∅`), the
ablation set is empty so `ablated ≡ N0_den` trivially — it CANNOT independently fail (a subsumed guard). It is also `nan` on
the singular `P²/Δ²` vector-port fixtures. So the DECISIVE catch is the COMPUTED taint-set (P1/P2) + `source_map_complete`;
the ablation is honestly de-counted. (The SymPy witness uses substitution; the `.wl` witness uses `D[expr, v]` partial
derivatives — a materially different algorithm proving the same fact.)

### 1.5 Why the taint must be computed over TAGS (the two caught rigs)
The pathA_43 gate history holds two rigs a NAME-check would pass "by fiat":
- **Rig 1 — the relabel rig.** A density-LOOKING two-port built from `{omega_wall, omega_rho, r_mix, g_rho, g_qold}` — symbols
  NOT in `VECTOR_SYMBOLS` but tagged `vector_port`. A `free_symbols ∩ VECTOR_SYMBOLS` name-check PASSES it; the COMPUTED taint
  catches it (`vector_port ∈ taint`). (Caught by Codex xhigh at the source's v0.)
- **The hidden-vector rig.** `N0_den · σ_hidden`, `σ_hidden` tagged `vector_port` but NOT in `VECTOR_SYMBOLS` — a vector
  intermediate surviving into the numerator even though the rest of the host set "looks" density-only; caught by the taint.
The adversarial re-review PROVED this is load-bearing: downgrading the machinery to a name-check (stripping `vector_port` from
`σ_hidden` + the relabel symbols) makes the relabel/hidden META asserts FIRE in BOTH engines; hardcoding a self-asserted
taint or neutering the assert also fires. (The self-asserted-`continuity_interface`-tag rig, "Rig 2", is 026's lineage
token-check — 025 CITES `moment_valid`, it does not re-earn the lineage.)

---

## 2. The able-to-fail battery (025-owned; per-tooth ablation — each fired at its OWN named assert, then went vacuous when neutered)

| tooth | rig / mutation → routed assert | notes |
|---|---|---|
| A relabel_rig | a density-looking port over `vector_port`-tagged symbols → the **computed-taint** assert (`vector_port ∈ taint`) fires | ⭐ the anti-fiat core; a name-check MISSES it |
| B hidden_vector | `N0_den · σ_hidden` (σ_hidden tagged `vector_port`, not in VECTOR_SYMBOLS) → the **computed-taint** assert fires | a 2nd fixture shape of the same tag detection |
| C vector_injection | `N0_den · Ω_U/Ω_W` (genuine `VECTOR_SYMBOLS`) → the **vector-host** assert (`∩ VECTOR_SYMBOLS ≠ ∅`) fires | deliberately caught by the name set |
| D provenance_less_rider | `N0_den · free_carrier` (empty tags) → the **source_map_complete** assert (empty-tag half) fires | a symbol with no provenance cannot ride along |
| E missing_symbol | `N0_den · missing_rider` (absent from the tag map) → the **source_map_complete** assert (missing half) fires | the complementary half of the guard |
| F tagged_carrier (able-to-PASS) | `N0_den · free_carrier` tagged `pathA_34_dimensionless_free_carrier` → PASSES P2 with a 4-tag taint, NOT P1; neuter the tag → `source_map_complete` FLIPS to false | the anti-over-rejection / reversibility control (stage020 MIXED-control discipline) |
| G raw_vector_port | the bare old vector port `Ω_U²·g_W + R_mix·g_U` → the **vector-host** assert fires | the deliberately-obvious negative |
| I subject_integrity | corrupt the cited `N0_den` (`rho_eff → foreign_subject`) → `free_symbols ≠ HOST_CONTRACT` → the **exact host-contract** assert fires | the consumption oracle |
| H′ arity | a wrong call arity → the def/call **arity scanner** fires (compares the real held-call arity to the real `DownValues` def arity) | the stage007 lesson |
| H′ leakage | a leaked unevaluated authored-helper head in a transcript object → the **unevaluated-leakage scanner** fires | — |

Each rig is exercised by a coupling meta-test (`exercise_rig`): it must fire at its OWN named assert AND, when the rig is
neutered (the tag repaired / the mutation reverted), stop firing — proving the tooth is coupled to the thing it claims to
check. ⚠ Two checks are DELIBERATELY NOT counted teeth: the **expression ablation** (§1.4, a de-counted witness) and the
**`.wl` algorithm independence** (a TRANSLITERATION-REVIEW property — replacing the native `Graph` traversal with a copied
`Variables`+`Lookup`+`Union` can preserve every value, so no runtime assert rejects it; kept in the Codex→Grok→Codex screen).

---

## 3. Honest scope

- **025 proves the vector-freedom STRUCTURE, CONDITIONAL on 026.** The LOCAL verdict is conditional on the cited
  `moment_valid=True` — a typed forward reference to stage 026 (which earns the `∫Y₂*S_leak` continuity-moment lineage). 025
  owns the tag-map + the computed taint; 026 owns the lineage VALIDATION. If 026 fails to earn `moment_valid`, 025's LOCAL
  result is void — hence "conditional."
- **The retirement is 027's, not 025's.** 025 proves the vector-freedom CONJUNCT; the EM `A_w`/`U,W` scaffold RETIREMENT is
  the JOINT `DENSITY_PORT_HOSTED` (025 vector-freedom ∧ 026 lineage ∧ 027 dim/scaling/sign/closure). If 027 later fails, the
  density port has NOT displaced the scaffold. 027 records the retirement with the completed joint.
- **Magnitudes SIM_DEFERRED.** 025 is a proof over the STRUCTURE of `N0_den`; the literal `Ξ_Q`/`λ_c`/`rho_eff` throat values,
  the `G`/`2/5`/`54/5` calibration, and the port closure are downstream (024 SIM_DEFERRED, 027 CALIBRATED).

---

## 4. Consumed / exported

- **Consumed — PROVENANCE + the one checkable subject contract.**
  - **stage 024's `N0_den`** — the SUBJECT (cited factored export; the consumption-integrity oracle asserts its computed
    `free_symbols == HOST_CONTRACT`, the one checkable cross-stage constraint).
  - **stage 026's `moment_valid`** — the `I25` continuity-moment validity boolean, cited as a typed input (forward ref; the
    LOCAL verdict is conditional on it). PROVENANCE.
  - the density symbols' PROVENANCE tags (`BASE_SOURCE_TAGS`) — sourced from 024's `physical_relations` (ϖ_q2←pathA_32 wall,
    ϖ_Φ2←pathA_29 bulk, λ_c←continuity) + stage005 (`c_s`) + `a` CONV.
- **Register.** ZERO new counted CALIB knobs (a PROOF/taint slice — it consumes 024's `N0_den`, introduces NO new physical
  symbols; the vector/relabel/free-carrier symbols are CONTROL fixtures, tracked-not-counted, like stage023's
  `q_free`/`eta_null`). Part-II counted CALIB set unchanged at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017)
  = **6**. New structural edge **R44** — the density-port vector-freedom PROOF (the computed taint = the density ancestry,
  `source_map_complete`, `vector_host_symbols=∅`), records ONLY the vector-freedom CONJUNCT (the retirement is the JOINT
  025∧027 result, 027 records it), discharging nothing (a proof/provenance edge, like R43/R41/R39). No new dims (025 runs no
  dim gate — `[N0_den]=L⁻¹M` is 027's).
- **Exported.** The COMPUTED taint set `{continuity_interface, pathA_29_bulk, pathA_32_wall}` + the LOCAL verdict token
  `DENSITY_PORT_VECTOR_FREE` + the `source_map_complete` certificate → **026** (which earns `moment_valid`) and **027** (which
  consumes the vector-freedom verdict as one conjunct of the joint `DENSITY_PORT_HOSTED`).

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O — and pathA_43 was ALREADY contract-clean, so the reshape
"sever the bridge" step is a NO-OP; the work was DECOMPOSITION + the `.wl` RE-AUTHOR + the LOCAL/PARTIAL verdict framing. ⚠
**Unlike stage 024 (which KEPT its `.wl` native for the derivation lane), stage 025 RE-AUTHORS its `.wl`** — the source `.wl`
taint machinery was a near line-for-line transliteration of the `.py`, so the `MATHEMATICA_MIRROR_POLICY` required a
genuinely independent computation: the decisive taint + `source_map_complete` use a native directed `Graph` of
`DirectedEdge["sym:"→"tag:"]` + `VertexOutComponent` reachability (NOT `Variables`+`Lookup`+`Union`), and the ablation
witness uses `D[expr, v]` partial-derivative independence. Agreement is transcript-level (both emit the cited `N0_den`, the
computed taint set, `source_map_complete`, the per-rig outcomes, the LOCAL verdict + JOINT PARTIAL); neither engine reads the
other. SymPy 18 PASS / Mathematica 18 PASS, both exit 0, CWD-independent.

**Directive review** used the Codex→Grok→Codex bookend. The Codex xhigh design-review folded 5 BLOCKING — most importantly
that the expression-level ablation is LOGICALLY SUBSUMED by the taint-set gate (reframed to a retained de-counted witness;
the COMPUTED taint-set identity + `source_map_complete` made the decisive gate) + the two-predicate split
(`baseline_ancestry_ok` P1 vs `vector_free` P2, so the properly-tagged 4-tag carrier passes the general predicate not the
exact identity) + R44 records only the vector-freedom conjunct + the simplified cite-the-factored-`N0_den` consumption + tooth
H split into a review-acceptance and two runtime scanners. A Grok-4.5 compute-verify pass returned `DIRECTIVE_CLEAN` (SymPy
re-confirmed: the subsumption, the relabel/hidden rigs' `vector_port` taint vs the name-check ∅, the 4-tag reversibility
control, the `source_map_complete` halves) + 3 nits folded; a Codex confirm-pass folded 4 consistency-sweep gaps; a Codex
final re-confirm returned `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read re-derived `N0_den` from stage024's 2×2 inverse, DIFF=0;
confirmed the tag-map faithfulness with the 024 renames, the P1/P2 separation, the conditional verdict, the de-counted
ablation witness, and the genuine independent `.wl` route) + `ADVERSARIAL_CLEAN` (per-tooth ablation on BOTH engines: every
rig fired at its OWN named META assert and went vacuous when neutered; the decisive name-check-dodge test proved the
tag-taint is load-bearing — downgrading to a name-check fires the relabel/hidden asserts; the ablation is non-verdict-bearing;
the `.wl` arity scanner compares real held-call vs real `DownValues` arity). The adversarial leg flagged 2 [minor]
non-blocking checks that could not independently fail — a subsumed `<= DENSITY_HOST_UNIVERSE` assert (subsumed by the exact
`== HOST_CONTRACT`) and a masked `moment_valid=True` assert (masked by `P2 source_map_complete` firing first) — both
DE-COUNTED to labeled diagnostic prints (tally 20→18 per engine, no other check touched), then fresh-agent `REVERIFY_CLEAN`
(the diff is exactly the 2 de-counts; the anti-rig core unchanged and still able-to-fail; both engines 18/18).
