# ledger_stage028 — the 2.5PN Burke–Thorne match-back (Check II-P5)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative-port sector, Cluster C). The CONSISTENCY cap of the
pathA_43 density-port fold (024 derivation ∧ 025 vector-freedom ∧ 026 continuity-lineage ∧ 027 port-checks+closure): the
2.5PN Burke–Thorne match-back over the port's CALIBRATED closure moments. Source gate: `pathA_2_5pn_matchback` (the Phase-A
A3 consolidation artifact). The last COMPUTED stage of the Part-II gravity sector (029 = a thin PN DOI-cite; then the
scheduled MIDWAY KNOB AUDIT).

**Verdict.** `MATCHBACK_CONSISTENT` (the LOCAL exit-0 consistency gate = baseline all-zero ∧ every one of the 11 mutation
caught-by rows == expected ∧ the exact-rational no-float discipline). ⚠ **There is NO `FAIL_*`/`PASS` PHYSICS verdict token
— 028 is a CONSISTENCY artifact.** ⚠⚠ **CALIBRATED consistency, NOT a first-principles `Γ̄₅`/`G` derivation:** the K̄ moments
are hardcoded CALIBRATED closure inputs (`G=GENUINE_BLOCKED`, 020's provenance; `54/5`/`Γ̄₅`=`external_bridge_input`, the
`27`=018's `derived_in_gate`); the full 1PN→4PN from-throat re-derivation stays SIM-DEFERRED (Gate 6).

**Status.** Exact symbolic (exact-rational residuals in moment/z-space; `expect_zero`/`expect_bool`; no `scipy`/`numpy`/
floats/tolerances; a no-float guard over raw + reduced residuals). Adds **zero** new counted knobs (a consistency slice).
New structural edge **R47** (the full match-back consistency, SHARED with 027's R46 port-closure slot — not double-counted).

---

## 1. What this stage earns

The separate audited PN corpus `research/4d_2_5pn` proved a **conditional** 2.5PN theorem whose SINGLE open item is the
scalar normalization `m̂₀²Γ₅ = 2G/(5c⁵)` on the moving-throat branch (`4d_2_5pn.tex:57–60`; boxed target `:819–824`;
equivalently the moment-invariant pair `K̄₄K̄₀ = 4K̄₂²` ∧ `Γ̄₅ = 9K̄₂^{5/2}/K̄₀^{3/2}`, `:469–473`). pathA_43 built the ℓ=2
quadrupole radiative port natively on the density/`c_s` mode (024–027); its CALIBRATED closure moments

```
K̄₀ = 54 G c_s⁵ / (5 a⁵ c⁵)      K̄₂ = 6 G c_s³ / (5 a³ c⁵)      K̄₄ = 8 G c_s / (15 a c⁵)
```

**MEET that open item at reduced-closure.** Stage 028 is the self-contained, dual-engine, able-to-fail artifact that
verifies it, via five invariants — each an exact-rational residual reduced to `0` in BOTH engines:

- **INV1** `K̄₄·K̄₀ − 4·K̄₂² = 0` — the corpus moment invariant. (VERIFIED: `4K̄₂²/K̄₀ = (8/15)Gc_s/(ac⁵) = K̄₄`.)
- **INV2** `K̄₀·a⁵/(27·c_s⁵) − 2G/(5c⁵) = 0` — the **pathA_43 form → Burke–Thorne** (`Γ̄₅ = K̄₀·a⁵/(27c_s⁵) = 54G/(5·27·c⁵)
  = 2G/(5c⁵)`). The `/27` is the outgoing ℓ=2 density-Hankel fingerprint EARNED at 018 (the `+i z⁵/27`), NOT calibrated.
- **INV3** `9·K̄₂^{5/2}/K̄₀^{3/2} − 2G/(5c⁵) = 0` — the corpus's **OWN** form → Burke–Thorne. (VERIFIED: `9·(6/5)^{5/2}/
  (54/5)^{3/2} = 2/5`, using `54 = 9·6` so the `√6` cancels: `9·6^{5/2}/(5·27·6^{3/2}) = 9·6/135 = 2/5`.) ⭐ This tie-out is
  028's genuinely-new content — 027's closure only checks the pathA_43 form.
- **INV4** `path_form − corpus_form = 0` — the cross-form agreement. ⚠ **INV4 ≡ INV2 − INV3 IDENTICALLY** (`(path−bt) −
  (corpus−bt) = path − corpus`), so it is a REDUNDANT cross-form DIAGNOSTIC (retained for transcript readability), NOT an
  independent tooth.
- **INV5** the 8 **INDEPENDENT literal anchors** pinning `{K̄₀,K̄₂,K̄₄,BT}` to `{54/5, 6/5, 8/15, 2/5}` + the structural
  literals `{27, 9, 5/2, 3/2}`. ⭐⭐ **The anti-rig teeth** (see §2).

**The `54/5 = 2·27/5` split.** The magnitude entering the normalization is `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵)` — the `27` is EARNED
(the density-Hankel tail, 018/pathA_43), the `2/5` Burke–Thorne factor and `G` are CALIBRATED/external. So `c⁵` (light) is a
GR-matching / `λγ` units bridge (`P₀ ∝ c_s⁵/c⁵`), NOT EM propagation.

**Physical picture** (`docs/conceptual_foundation.md` §3/§5): gravity's quadrupole radiation rides the medium's OWN density
ripple at `c_s`. This stage is the *cheapest decisive falsifier* for the density-mode port: had the calibrated moments been
mutually inconsistent, or inconsistent with Burke–Thorne, the residuals would be nonzero and the gate would FAIL.

### 1.1 The A3 boundary (028 ↔ 027) — SHARED, not double-counted
027's `closure_overlay` checked the **two** port-closure residuals `K̄₄−4K̄₂²/K̄₀=0` ∧ `Γ̄₅−2G/(5c⁵)=0` (the A3 SLOT it
OWNS, R46). **028's INV1 and INV2 ARE those two residuals** re-expressed: INV1 = (027-residual-1) × K̄₀ (same content,
`K̄₀≠0`); INV2 = 027-residual-2 exactly (`Γ̄₅ = K̄₀·a⁵/(27c_s⁵)`). **028 ADDS** INV3 (the corpus's own form), INV4 (cross-form),
INV5 (the anchors), and the 11-mutation coherent-rescale matrix — its genuinely-new content. 028 **CONSUMES** 027's exported
K̄ moments as PROVENANCE (restated locally, self-contained; the port-closure consistency is counted ONCE at R46), and does
NOT rebuild the port derivation (024) or the port checks (027). ⭐ An authoring-time **A3 fidelity comparison** confirms
028's restated K̄ moments equal 027's actual `closure_overlay` exports (a stale-citation guard; NOT a runtime cross-file
tie).

---

## 2. The able-to-fail battery (per-tooth ablation; each fired at its OWN named assert)

⭐⭐ **The anti-rig invariant (per-gate trip-up L95): the INV5 literal anchors are 028's OWN independently-restated
constants, computed NOWHERE from the moment/config coefficients — never deduped against the closure moments.** They are the
ONLY residuals that catch the **coherent-rescale** mutation (below); dedupe them and the rig PASSES silently.

**The GENUINE independent EXIT-1 teeth (the 11 caught-by rows + the no-float guard — 12 counted):**
- **⭐⭐ Row 1 — the coherent-rescale anti-rig.** `{K̄₀,K̄₂,K̄₄,BT}×2` passes INV1–INV4 (scale-covariant: INV1 ~ λ²;
  path/corpus/bt ~ λ with corpus degree `5/2−3/2 = 1`), caught **ONLY** by the 4 INV5 anchors `{INV5_Kbar0, INV5_Kbar2,
  INV5_Kbar4, INV5_BT}` (`K̄₀·a⁵c⁵/(Gc_s⁵) = 2·54/5 ≠ 54/5`). The decisive ablation: DEDUPE the anchor RHS against the
  rescaled moments → the row goes EMPTY (the rig PASSES) → the row's `actual == expected` assert EXITS 1.
- **Row 2 — the coupled mutation.** `{K̄₀,K̄₂,K̄₄}×2`, BT fixed → INV4 stays 0 (cross-form still agrees) while
  `{INV2,INV3,INV5_Kbar0,INV5_Kbar2,INV5_Kbar4}` fire — shows INV4 alone does not close the rescale gap.
- **Rows 3–11 — the 9 single-parameter mutations** (K̄₄→8/14, K̄₄ sign-flip, K̄₂→7/5, K̄₀→55/5, denom 27→26, corpus 9→8,
  exp 5/2→3/2, exp 3/2→1, BT 2/5→3/5), each firing its own INV subset (its own caught-by row assert).
- **The no-float guard** — a `sp.Float`/`_Real` atom in any residual (raw or reduced) fires.

Each row's assertion is `actual == EXPECTED_CAUGHT_BY[name]` with the oracle table held IMMUTABLE; the in-script 11-mutation
matrix runs at EXIT 0 (the diagnostic). The adversarial tri-review leg ablates each tooth to EXIT 1 (actual-side mutation,
oracle immutable) + a coupling meta-test (neuter the planted drift → the expected failure disappears → the harness fails).

**Review / acceptance GATES (verified by the tri-review, NOT runtime EXIT-1 residual teeth):** the `.wl` per-function
independence + arity/leakage scan (G1), and the honest-scope print (G2).

**DE-COUNTED (retained for fidelity/transcript, NOT independent teeth):** baseline-all-zero (`all(rᵢ==0)`, the aggregate
positive landing — subsumed by the per-residual checks); the `actual`-non-empty check (entailed by `actual == expected`);
and INV4 (≡ INV2−INV3).

---

## 3. Honest scope

**IS:** a dual-engine consistency check that the density-mode ℓ=2 port's CALIBRATED closure reproduces the 2.5PN
Burke–Thorne normalization `research/4d_2_5pn` left open (`Γ̄₅=2G/(5c⁵)`) + the moment invariant `K̄₄=4K̄₂²/K̄₀`. The cheapest
decisive falsifier for the density-mode port.

**IS NOT** a first-principles derivation of `Γ̄₅`/`G`. The moment coefficients `K̄₀,K̄₂,K̄₄` are hardcoded CALIBRATED
closure inputs, not solved from the port numerator `N0_den`; `G = GENUINE_BLOCKED`; `Γ̄₅`/`54/5` = `external_bridge_input`
(the GR Burke–Thorne bridge); only the `27` is `derived_in_gate` (018). A full first-principles 1PN→4PN re-derivation from
the throat interior stays SIM-DEFERRED (Gate 6) — out of reach (solver tractability, not hardware). 028 is the reachable
consistency check, consistent with the calibrate-predict discipline and the sim-deferred guardrail (completing the SPEC, not
proving the theory).

---

## 4. Consumed / exported

**Consumed (PROVENANCE):** 027's exported K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` (the A3 SLOT; INV1/INV2 SHARED with R46, checked by
the A3 fidelity comparison) + 018's `27` density-Hankel fingerprint (INV2's `/27`) + 020's `54/5=2·27/5` /
`G=GENUINE_BLOCKED` (the calibrated magnitude) + the corpus's OWN form `9K̄₂^{5/2}/K̄₀^{3/2}` (`4d_2_5pn.tex:469`, INV3,
imported) + `c` (GR-units bridge, benchmark) + `a` (CONV). No runtime cross-consumption of any peer/report/source/note file
(zero file I/O).

**Exported (→ 029 + Part VII):** the 2.5PN consistency landing — the `research/4d_2_5pn` single open item MET at
reduced-closure. ⚠ Cluster C continues **029** (PN corpus DOI-cite — a THIN cite-only stage note, NO scripts) → then the
scheduled **MIDWAY KNOB AUDIT** (Part-II gravity sector CLOSES — the pathA_40 `Δr=2` codimension dry-run over Parts I–II +
the held-out vs irreducible-route-less tally). 028 is the CONSISTENCY cap of the pathA_43 density-port fold — the last
COMPUTED stage of the Part-II gravity sector.

---

## 5. Dual-engine and verification

**Contract-clean, no bridge to sever** (like 024): both source engines were already standalone print-only, zero file I/O.
The reshape was rename + LOCAL-ledger idioms + the A3 PROVENANCE citation + the honest-scope framing + the `.wl` RE-AUTHOR.

**⭐ The `.wl` is RE-AUTHORED independent** (like 020/022/023/025/026/027, and unlike a keep-native lane): the source `.wl`
028 artifact was a one-for-one transliteration of the `.py` (identical config Association, `mutations`, `expectedCaughtBy`,
residual forms, and `FullSimplify[Cancel[Together]]` ↔ `factor∘cancel∘simplify` reduction), so the mirror policy required a
genuinely independent route — INV1 decided via a `GroebnerBasis` polynomial-identity test; INV2–INV5 via `PossibleZeroQ`
over positive-domain `Refine`/`PowerExpand` expressions; the mutations replace complete moment EXPRESSIONS directly (no
Python-style coefficient/config pipeline); the caught-by rows via `Pick` from the independent zero flags; machine-real
detection via `FreeQ[…, _Real]`. Reduced expressions are used only for transcript display, not verdicts. Transcript-level
agreement, neither engine reading the other.

**Dual-engine: SymPy 12 / Mathematica 12** counted EXIT-1 teeth (the 11 caught-by rows + the no-float guard; the baseline
zero-tests + arity + leakage checks reported separately as 14 diagnostics per the de-count rules), both exit 0,
CWD-independent (repo root + `/tmp`, byte-identical transcripts). The transcript prints the A3 fidelity line
(`A3_FIDELITY: restated literals match stage027 closure_overlay`, an authoring-time comparison — the fidelity leg
independently confirmed 028's restated K̄ moments equal 027's actual `closure_overlay` exports).

**Directive review — the Codex→Grok→Codex bookend caught issues at EVERY leg:** the Codex xhigh design-review folded 2
BLOCKING (subsumed teeth — baseline-all-zero/`actual`-non-empty/INV4≡INV2−INV3 cannot all be independent own-assert teeth;
and §4 needed the adversarial per-row EXIT-1 ablation protocol + coupling meta-test, since the in-script matrix runs at exit
0) + 2 nits (the consumption rationale reword + the post-build A3 fidelity comparison; the consolidation-note path + file
lengths); the Grok-4.5 compute-verify pass returned `DIRECTIVE_CLEAN` (it independently SymPy-verified the coherent-rescale
caught only by the 4 INV5 anchors, the forbidden-dedupe rig passing, INV1=027-res-1×K̄₀ + INV2=027-res-2, all 12 BASE
residuals = 0 incl. the corpus fractional-power reduction, INV4≡INV2−INV3); the Codex confirm-pass folded 2 BLOCKING +
1 nit (teeth-list consistency — the `.wl`-independence + honest-scope are review GATES not EXIT-1 teeth, coherent-rescale is
row 1 not a separate tooth; the source map had to be kept in sync); a Codex final-confirm folded 1 BLOCKING (2 stale
source-map headline spots) → `DIRECTIVE_CLEAN`.

**Tri-review** (fresh agents, both CLEAN → ZERO remediation, like the contract-clean 024/027): arbiter re-run reproduced
both engines (exit 0, 12/12, `MATCHBACK_CONSISTENT`, A3 fidelity match); `FIDELITY_CLEAN` (INV1–INV5 faithful, all 12
baseline residuals → 0 incl. the INV3 fractional-power tie-out; the caught-by matrix faithful with no mutation dropped; the
A3 fidelity comparison genuine not a self-compare; the `.wl` genuinely independent; honest scope + honest de-count; no
math-flags); `ADVERSARIAL_CLEAN` (⭐ the decisive coherent-rescale anti-rig LOAD-BEARING — deduping the `ANCHOR_*` literals
[independent module-level, `.py` L163–170 / `.wl` L105–112, never read from the moment config] makes the coherent-rescale
row compute empty → EXIT 1 in both engines; all 11 caught-by rows fire at their own assert + the coupling meta-test holds
universally; the caught-by matrix genuinely runtime-computed; the `.wl` not a mirror; the de-counted set genuinely
non-load-bearing — INV4 has no isolated knob; all 12 counted teeth load-bearing; no reshape rig).
