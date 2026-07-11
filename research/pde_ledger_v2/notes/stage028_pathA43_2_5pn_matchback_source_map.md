# II-P5 (ledger_stage028) source map — the 2.5PN Burke–Thorne match-back (`pathA_2_5pn_matchback`; consistency over CALIBRATED moments)

> Running-start prep captured 2026-07-10 (post stage027 = pathA_43 II-P4 port CHECKS + CLOSURE DONE, commit
> `8f131184`; the pathA_43 density-port fold COMPLETE 024∧025∧026∧027), before authoring the II-P5 reshape directive, so
> the directive can be written without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-10**
> by a full read of BOTH engines (`pathA_2_5pn_matchback_sympy.py` = **249 lines**, `pathA_2_5pn_matchback.wl` = **215
> lines**) PLUS the consolidation note `notes/pathA_43_2_5pn_matchback.md` (177 lines) that these scripts back. No line
> drift (the sources untouched since the 2026-07-06 A3 commit).
> Companion: `part2_gravity_atomic_split.md` (Cluster C row **028** L53 = "INV1–INV5 (both Γ̄₅ forms → Burke–Thorne,
> cross-form, K̄₄=4K̄₂²/K̄₀, the INV5 literal anchors) + the 11-mutation matrix; `G=GENUINE_BLOCKED` carried"; the ▶ NEXT
> stage028 entry **L680–690**; the **per-gate trip-up L95** [⚠ "pathA_2_5pn_matchback (028): the INV5 literal anchors are
> the anti-rig teeth vs the coherent-rescale mutation — never dedupe them against the closure values"]; the reshape-cost
> map **L78–79** [pathA_43 + matchback = the lightest reshape class, contract-clean; the work is decomposition + preserving
> the anti-rig mechanisms]; the cross-stage flows **L109** ["024–027 export N0_den + K̄ moments → 028"]) + the built
> **stage027 pair + source map + directive** (the A3 boundary 028 shares — 027's K̄ CLOSURE overlay = the two-residual
> port-closure CONSISTENCY, 028 = the FULL INV1–INV5 match-back + the coherent-rescale matrix; SHARED, not double-counted)
> + the **stage024 pair** (the CONTRACT-CLEAN lightest-reshape exemplar — decomposition, zero bridge, the anti-rig teeth
> preserved). Build-order id **028**, Part II, Cluster C (the 2.5PN consistency cap).
> Source top-line: there is **NO `FAIL_*`/`PASS` physics verdict token** — the artifact emits `SYMPY_MATCHBACK: PASS` /
> `WOLFRAM_MATCHBACK: PASS` (the exit-0 consistency gate). The ledger headline (split table L53) is **"(consistency over
> CALIBRATED moments)"**. Proposed target stem: `ledger_stage028_2_5pn_matchback` (confirm slug at directive authoring).
> Proposed LOCAL verdict token: `MATCHBACK_CONSISTENT` (the exit-0 gate; the honest scope printed alongside — CALIBRATED
> consistency, `G=GENUINE_BLOCKED`, NOT a `Γ̄₅`/`G` derivation).

## ⭐ The FIVE headline points (READ FIRST)

1. **⚠ 028 is a CONSISTENCY artifact, NOT a derivation gate — the honest scope is load-bearing.** It checks that the
   density-mode ℓ=2 port's CALIBRATED closure moments `{K̄₀,K̄₂,K̄₄}` (derived-structure/calibrated-magnitude at 024/027,
   the `54/5` calibrated + the `27` earned per 020/018) reproduce the 2.5PN Burke–Thorne normalization
   `Γ̄₅ = 2G/(5c⁵)` + the moment invariant `K̄₄ = 4K̄₂²/K̄₀` that the separate audited corpus `research/4d_2_5pn` left as
   its **single open item** (`4d_2_5pn.tex:57–60`, boxed target `:819–824`; equivalently `Γ̄₅ = 9K̄₂^{5/2}/K̄₀^{3/2}` +
   `K̄₄K̄₀=4K̄₂²`, `:469–473`). ⚠⚠ **It is NOT a first-principles derivation of `Γ̄₅`/`G`** — the moments are hardcoded
   CALIBRATED closure inputs, `G = GENUINE_BLOCKED`; the full 1PN→4PN from-throat re-derivation stays **SIM-DEFERRED
   (Gate 6)**. The note must carry this scope verbatim-class (consolidation note §4: "IS a dual-engine consistency check …
   IS NOT a first-principles derivation of `Γ̄₅`"). The able-to-fail content is: *had the calibrated moments been mutually
   inconsistent, or inconsistent with Burke–Thorne, the residuals would be nonzero and the gate would FAIL* — proven by the
   11-mutation matrix.

2. **⭐⭐ THE EARNED CONTENT = INV1–INV5 (each an exact-rational residual reduced to 0, dual-engine) + the 11-mutation
   able-to-fail matrix with an asserted caught-by matrix.** The five invariants (`.py` `build_residuals` L157–188, raw
   forms **L166–179**; `.wl` `buildResiduals` L136–164, raw **L145–158**):
   - **INV1** `K̄₄·K̄₀ − 4·K̄₂² = 0` — the corpus moment invariant (`.py` L167).
   - **INV2** `K̄₀·a⁵/(27·c_s⁵) − 2G/(5c⁵) = 0` — the **pathA_43 form → Burke–Thorne** (`.py` L168; `path_form` L163).
     The `/27` is the outgoing ℓ=2 density-Hankel fingerprint EARNED at 018 (`+i z⁵/27`), NOT calibrated.
   - **INV3** `9·K̄₂^{5/2}/K̄₀^{3/2} − 2G/(5c⁵) = 0` — the corpus's **OWN** form → Burke–Thorne (`.py` L169;
     `corpus_form` L164). ⭐ **This is the tie-out pathA_43/027 does NOT carry** — 028's genuinely-new content (027's
     closure only checks the pathA_43 form).
   - **INV4** `path_form − corpus_form = 0` — the redundant cross-form agreement, built DIRECTLY (`.py` L170). Genuinely-new.
   - **INV5** the **literal anchors** pinning `{K̄₀,K̄₂,K̄₄,BT}` to `{54/5, 6/5, 8/15, 2/5}` and the structural literals
     `{27, 9, 5/2, 3/2}` (`.py` L171–178; `RESIDUAL_ORDER` slots 5–12, L26–33). ⭐⭐ **THE ANTI-RIG TEETH — INDEPENDENTLY
     RESTATED literals, NEVER deduped against the closure moments** (per-gate trip-up L95): they are the ONLY residuals
     that catch the coherent-rescale mutation (below).

3. **⭐⭐ THE COHERENT-RESCALE ANTI-RIG (the #1 design invariant to preserve).** The 11-mutation probe (`.py` MUTATIONS
   **L64–82**, `EXPECTED_CAUGHT_BY` **L84–144**; `.wl` `mutations` **L46–64**, `expectedCaughtBy` **L66–126**) feeds
   data-only corruptions through the SAME residual path and asserts an expected caught-by matrix. The nine single-parameter
   mutations (K̄₄→8/14, K̄₄ sign-flip, K̄₂→7/5, K̄₀→55/5, denom 27→26, corpus 9→8, exp 5/2→3/2, exp 3/2→1, BT 2/5→3/5)
   each fire an obvious INV subset. **The two DECISIVE mutations:**
   - **`coherent_scale_all_moments_and_BT_x2`** (`k0,k2,k4,BT × 2`, L74–77): **passes INV1–INV4** (the moment invariant +
     both BT forms + the cross-form are all scale-covariant), caught **ONLY by the INV5 anchors** (L131–136). ⭐ **This is
     WHY the INV5 anchors exist and must stay independent** — dedupe them against the (rescaled) moments and this rig
     PASSES silently.
   - **`coupled_moments_x2_BT_fixed`** (`k0,k2,k4 × 2`, BT fixed, L78–81): leaves INV4 zero (cross-form still agrees) while
     INV2/INV3/anchors fire (L137–143). Shows INV4 alone is not sufficient.
   ⚠ **The directive MUST pin: the INV5 anchor RHS literals are 028's OWN independently-restated constants, computed
   NOWHERE from the `k*_coeff`/`k*_scale` config; the coherent-rescale caught-by row is `{INV5_Kbar0, INV5_Kbar2,
   INV5_Kbar4, INV5_BT}` and NOTHING else.** The EXECUTABLE tooth is ROW 1 of the 11 caught-by rows — the coherent-rescale
   row's immutable-oracle `actual == expected` assertion (the adversarial leg ablates the ACTUAL side: dedupe the anchors →
   the row goes EMPTY → EXIT 1; §4-1), NOT separately-tallied anchor asserts.

4. **⭐ THE A3 BOUNDARY (028↔027) — SHARED, not double-counted.** 027's `closure_overlay` checked the **two** port-closure
   residuals `K̄₄−4K̄₂²/K̄₀=0` ∧ `Γ̄₅−2G/(5c⁵)=0` (the A3 SLOT it owns; R46). **028's INV1 and INV2 ARE those two
   residuals** re-expressed: INV1 = `K̄₄·K̄₀ − 4K̄₂²` = (027-residual-1) × K̄₀ (same content, `K̄₀≠0`); INV2 =
   `K̄₀·a⁵/(27c_s⁵) − 2G/(5c⁵)` = 027-residual-2 exactly (`Γ̄₅ = K̄₀·a⁵/(27c_s⁵)`). **028 ADDS** INV3 (the corpus's own
   form), INV4 (cross-form agreement), INV5 (the anchors), and the 11-mutation coherent-rescale matrix — its genuinely-new
   content. ⚠ **The directive states the cut: 028 CONSUMES 027's exported K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` as PROVENANCE (the
   baseline of INV1–INV4), does NOT rebuild the port derivation (024) nor the port checks (027); the port-closure
   consistency is SHARED with R46, counted ONCE.** At the register: ZERO new knobs, ZERO new dims; edge R47 = the full
   match-back consistency, SHARED with R46's port-closure slot (not double-counted).

5. **⭐ CONTRACT-CLEAN — there is NO bridge to sever (the lightest reshape class, shared with 024).** BOTH engines are
   ALREADY standalone print-only, `raise SystemExit(1)` / `Exit[1]` on failure, **zero file I/O** (`.py` imports only
   `sympy`/`dataclasses`; `.wl` grep for `Export|Put|Write|OpenWrite|>>|Save|Import|Get` = ZERO hits; both explicitly
   "self-contained by design: restates the calibrated literals locally, does not read or write repo artifacts, and does
   not consume the peer engine", `.py:2–7`/`.wl:1–6`). So the reshape "sever the bridge" step is a **no-op**. The 028 work
   is: **(a)** rename to the ledger stem + LOCAL-ledger idioms (`banner`/`subbanner`/`expect_zero`/tally/`OVERALL PASS`/
   nonzero exit for `.py`; `heading`/`fail[]`/`Exit[1]` for `.wl`) replacing the bespoke `main()`/`Print` loop; **(b)** the
   A3 PROVENANCE citation of 027's K̄ moments (headline 4) + the honest CALIBRATED-consistency scope (headline 1); **(c)**
   ⚠⚠ **RE-AUTHOR the `.wl`** — the source `.wl` is a **one-for-one TRANSLITERATION** of the `.py` (identical config
   Association, identical `mutations` list, identical `expectedCaughtBy` matrix, identical residual forms; the same
   `factor∘cancel∘simplify` ↔ `FullSimplify[Cancel[Together]]` reduction) → the `MATHEMATICA_MIRROR_POLICY` screen
   REQUIRES a genuinely independent route (the stage027 lesson: "the source `.wl` being one-for-one ports means the
   018/021 keep-native precedent does NOT transfer"); **(d)** the per-tooth-ablation battery over the GENUINE EXIT-1 teeth =
   the 11 caught-by row assertions (coherent-rescale = ROW 1) + the no-float guard (the `.wl` independence + honest-scope are
   review/acceptance GATES, not teeth; baseline-all-zero + `actual`-non-empty + INV4≡INV2−INV3 DE-COUNTED — §4); **(e)** the
   clean LOCAL verdict `MATCHBACK_CONSISTENT` + the honest
   scope. ⚠ **NOTHING to DROP** (no `.py`-only vestigial constructs; the scripts are lean).

## §1 The 028 slice (the WHOLE artifact — both engines are single-purpose)

Unlike 024–027 (slices of the ONE pathA_43 monolith), the matchback scripts are DEDICATED single-purpose artifacts —
028 owns the WHOLE of each. No cut, no siblings to carve away.

### §1a The SymPy artifact (`pathA_2_5pn_matchback_sympy.py`, 249 ln)
- **Symbols** `G, c_s, a, c` (positive) **L19**; `R=sp.Rational`, `Z=sp.Integer` **L16–17**.
- **`RESIDUAL_ORDER`** **L21–34** — the 12 residual labels (INV1–4 + the 8 INV5 anchor slots). The stable print/assert order.
- **`MatchConfig`** dataclass **L37–51** — the calibrated coefficients (`k0/k2/k4/bt_coeff`, `path_denominator=27`,
  `corpus_factor=9`, `exp_k2=5/2`, `exp_k0=3/2`) + the four `*_scale=1` mutation hooks.
- **`BASE`** **L53–62** — the calibrated closure values `{54/5, 6/5, 8/15, 2/5, 27, 9, 5/2, 3/2}` (⭐ the A3-consumed
  moments; INV5 anchors reference the SAME literals on the RHS but as 028's OWN independent constants, headline 3).
- **`MUTATIONS`** **L64–82** + **`EXPECTED_CAUGHT_BY`** **L84–144** — the 11-mutation probe + the asserted caught-by matrix
  (⭐ the coherent-rescale L74–77 → INV5-only L131–136; the coupled L78–81 → INV2/INV3/INV5 L137–143).
- **`compact`** **L147–148** = `factor(cancel(simplify(...)))`; **`assert_no_float`** **L151–154** (the no-float guard,
  `atoms(sp.Float)`).
- **`build_residuals(cfg)`** **L157–188** — constructs `k0/k2/k4/bt` (with scales) L158–161, `path_form` L163,
  `corpus_form` L164, the 12 raw residuals **L166–179**, asserts no-float on raw + reduced, returns the reduced dict.
- **`fired_labels`** **L191–192**, **`vector_text`** **L195–196**, **`print_provenance`** **L199–205** (the CALIBRATED /
  `G=GENUINE_BLOCKED` / upstream-27 / imported-corpus-form / runtime-isolation provenance block).
- **`main`** **L208–241** — provenance, baseline residual vector (assert all-zero L214–215), the mutated vectors +
  caught_by L217–225, the EXPECTED caught-by assertion L227–235, the PASS prints L237–240; `SystemExit(main())` +
  `AssertionError → SystemExit(1)` **L244–249**.

### §1b The Mathematica artifact (`pathA_2_5pn_matchback.wl`, 215 ln) — ⚠ A TRANSLITERATION → RE-AUTHOR
- `fail[]`→`Exit[1]` **L10**; `$Assumptions` **L12**; `residualOrder` **L14–27**; `baseConfig` Association **L29–42**;
  `mergeConfig` **L44**; `mutations` **L46–64**; `expectedCaughtBy` **L66–126**; `stripConditionalZero`/`compact`
  **L128–129** (`FullSimplify[Cancel[Together[...]]]`); `assertNoFloat` **L131–134**; `buildResiduals` **L136–164** (raw
  **L145–158**); `firedLabels` **L168**; `printProvenance` **L170–177**; the main body **L179–214**; `Exit[0]` **L215**.
- ⚠⚠ **This `.wl` is a near line-for-line MIRROR of the `.py`:** same `baseConfig` keys ↔ `MatchConfig` fields, same
  `mutations`/`expectedCaughtBy` structure, same raw residual forms, same `FullSimplify[Cancel[Together]]` ↔
  `factor(cancel(simplify))` reduction, same print sequence. **Native SYNTAX ≠ an independent route.** The mirror-policy
  screen (blueprint §5.3) REQUIRES 028's `.wl` to be RE-AUTHORED — like stage020/022/023/025/026/027's re-authored `.wl`.
  ⭐ Candidate independent routes (a SUGGESTION menu for the directive, NOT a pre-design — Codex owns the code): compute
  each invariant's zero-test by a materially-different mechanism — e.g. `PossibleZeroQ`/`Simplify[resid]===0`/
  `Reduce[resid==0, {G,c_s,a,c}]` truth-values rather than `!=0` on a `FullSimplify∘Cancel∘Together` residual; build the
  moment invariant via `Resultant`/`GroebnerBasis`/`Eliminate` (INV1 as a polynomial identity); build the caught-by matrix
  by directly perturbing the moment EXPRESSIONS and testing `PossibleZeroQ` rather than mutating a coefficient config. ⚠
  **The DECISIVE computations (the 12 residual zero-tests + the caught-by matrix, esp. the coherent-rescale INV5-only row)
  are what must be verdict-bearing-AND-independent** (the 025/026/027 lesson: a copied computation retained for the real
  decision remains a mirror). Arity discipline (def/call scan + unevaluated-leakage transcript scan — the stage007 lesson).

## §2 The A3 consumption resolution (028 CONSUMES 027's K̄ moments — PROVENANCE, not a live dual-site tie)

⭐ **The consumption is a PROVENANCE citation, NOT a live dual-site integrity tie (the load-bearing design decision).**
Rationale: **(1)** zero file I/O — 028 cannot read 027's output; **(2)** ⚠⚠ **a live dual-site tie of the moments would
work AGAINST the anti-rig** — the whole point of INV5 (headline 3) is that the anchors are INDEPENDENT of the moments; if
028 "consumed" 027's moments AND re-derived the INV5 anchors from the SAME consumed values, the coherent-rescale rig would
scale both and PASS (the exact L95 trip-up). So 028 RESTATES 027's K̄ moments locally (self-contained, as the source
already is), CITES 027 as the PROVENANCE of those values (the density-port closure), and the INV5 anchors stay 028's OWN
independent literals. **(3)** A heavy set-then-compare "the baseline equals a cited 027 literal" would be exactly the
stage009 set-then-compare-to-self rig (a single-source corruption exits 0) — AVOID it. The genuine anti-rig is the
internal 11-mutation matrix + the independent INV5 anchors, already present.

- **Consumes (PROVENANCE):** 027's exported K̄ moments `{K̄₀=54Gc_s⁵/(5a⁵c⁵), K̄₂=6Gc_s³/(5a³c⁵), K̄₄=8Gc_s/(15ac⁵),
  Γ̄₅=2G/(5c⁵)}` (the A3 SLOT; INV1/INV2 = 027's two closure residuals re-expressed, SHARED) + 018's `27`
  density-Hankel fingerprint (in INV2's `/27`) + 020's `54/5=2·27/5` / `G=GENUINE_BLOCKED` (the calibrated magnitude) +
  the corpus's OWN form `9K̄₂^{5/2}/K̄₀^{3/2}` imported from `4d_2_5pn.tex:469` (INV3, NOT re-derived) + `c` (GR-units
  bridge) + `a` (CONV). ⚠ **NO runtime cross-consumption** of any peer/report/source/note/`_scratch` file (`.py` L205 /
  `.wl` L176 print the runtime-isolation guarantee — preserve it).
- **Exports (→ Part VII):** the completed 2.5PN consistency landing (the `research/4d_2_5pn` single open item MET at
  reduced-closure) — the last consistency cap before the PN DOI-cite (029) closes Cluster C.

## §3 Reshape cost (NO bridge to sever) + the 028 verdict (a CONSISTENCY landing, no FAIL/PASS physics token)

**⭐ CONTRACT-CLEAN (shared with 024; part2 L78–79).** NO argparse/JSON/scratch-YAML/`Export`/reads to strip. The reshape =
**rename + LOCAL-ledger idioms + the A3 PROVENANCE citation + the `.wl` RE-AUTHOR + the honest-scope framing**, preserving
the INV1–INV5 residual battery + the 11-mutation matrix + the caught-by assertion verbatim-content (never dedupe the INV5
anchors).

**⭐ THE 028 VERDICT (a CONSISTENCY landing).** 028 has **no `FAIL_*`/`PASS` physics verdict** (unlike 022/023's
characterized-FAIL or 024–027's `DENSITY_PORT_HOSTED`). It emits a LOCAL consistency verdict `MATCHBACK_CONSISTENT` (the
exit-0 gate = baseline all-zero ∧ every mutation caught-by matches expected ∧ no-float), with the honest scope printed:
**CALIBRATED consistency, `G=GENUINE_BLOCKED`, NOT a first-principles `Γ̄₅`/`G` derivation** (the full 1PN→4PN from-throat
= SIM-DEFERRED Gate 6). The note frames it as the *cheapest decisive falsifier* for the density-mode port (consolidation
note §4): had the calibrated moments been mutually inconsistent, or inconsistent with Burke–Thorne, the residuals would be
nonzero and the gate would FAIL — the 11-mutation matrix proves each residual able-to-fire.

**Acceptance (dual-engine, both exit 0, CWD-independent):**
- Each engine runs from repo root AND a foreign CWD (`/tmp`), print-only, exit 0, no files written.
- Both transcripts print: the provenance block (CALIBRATED / `G=GENUINE_BLOCKED` / upstream-`27` / imported corpus form /
  runtime-isolation); the baseline residual vector (all 12 = 0, the aggregate positive); the 11 mutated residual vectors +
  their caught_by; the EXPECTED caught-by matrix (each mutation's `actual == expected`, the load-bearing per-row assertion —
  the in-script matrix runs at **EXIT 0**); the no-float / mutation-probe PASS lines; the LOCAL `MATCHBACK_CONSISTENT`
  verdict + the honest CALIBRATED-consistency scope.
- Dual-engine agreement is transcript-level (both print the same 12 residual zero-tests + the same caught-by matrix);
  neither reads the other; the `.wl` passes the per-function transliteration screen (§1b — genuinely re-authored).
- ⭐ **The GENUINE independent EXIT-1 teeth (what the adversarial leg ablates, actual-side + coupling meta-test):** the **11
  per-residual `actual == expected` caught-by rows** (the coherent-rescale INV5-only anti-rig is ROW 1, not a separate tooth)
  + the **no-float guard**. ⚠ **The in-script 11-mutation matrix runs at exit 0 — NOT a substitute for the adversarial
  tri-review leg's EXIT-1 per-tooth ablation** (blueprint §6; §4 protocol). The `.wl` independence + honest-scope are
  REVIEW/ACCEPTANCE gates (not EXIT-1 residual teeth). **DE-COUNTED (retained for fidelity/transcript, NOT teeth):**
  baseline-all-zero, the `actual`-non-empty check, and INV4 (≡ INV2−INV3).
- **⭐ The A3 FIDELITY comparison (authoring / tri-review, NOT runtime; §2):** 028's restated baseline K̄ moments match 027's
  ACTUAL exported `closure_overlay` K̄ moments (a stale-citation guard).

## §4 Teeth candidates (028-specific, per-tooth ablation MANDATORY — mutate the named object, confirm EXIT-1 AT its own assert)

⭐⭐ **The per-gate trip-up L95 + [[feedback-per-tooth-ablation]]: the INV5 anchors are independently-restated (never
deduped); each GENUINE independent tooth ablatable to EXIT 1 at its own named assert.**

⚠⚠ **THE PROTOCOL (aligned with the directive v2/v3 — read before the tooth list).** The in-script 11-mutation matrix runs at
**EXIT 0**: `main` computes each mutation's caught-by row and asserts `actual == EXPECTED_CAUGHT_BY[name]`. **That in-script
matrix is a diagnostic that PASSES on the normal run — it is NOT a substitute for the adversarial tri-review leg's per-tooth
EXIT-1 ablation** (blueprint §6). The adversarial leg (a fresh agent, on the built scripts) must, for EVERY genuine tooth:
**(i)** perturb the ACTUAL side (the row's mutation FIXTURE, or a residual's construction — with the `EXPECTED_CAUGHT_BY`
oracle table held IMMUTABLE) so the computed caught-by row drifts from expected → the row-specific `actual == expected`
assert EXITS 1; then **(ii)** the coupling meta-test — restore/neuter the planted drift → the expected failure DISAPPEARS →
the adversarial harness itself must FAIL. ⚠ Do NOT ablate by mutating the ORACLE (`EXPECTED_CAUGHT_BY`) — ablate the
ACTUAL-side computation. For residual-binding checks, use ISOLATED POST-CONSTRUCTION sentinel perturbations.

**The GENUINE independent EXIT-1 teeth (the 11 caught-by rows + no-float; each ablated actual-side + coupling meta-test):**

1. **⭐⭐ The coherent-rescale anti-rig row (ROW 1 of the 11 — the load-bearing test).** The `coherent_scale_..._x2` mutation
   → passes INV1–INV4, caught ONLY by `{INV5_Kbar0, INV5_Kbar2, INV5_Kbar4, INV5_BT}`. ⭐ **The decisive ablation
   (actual-side, oracle immutable): DEDUPE the INV5 anchor RHS against the (rescaled) moment coefficients → the computed
   coherent-rescale row becomes EMPTY (the mutation now PASSES) → the row's `actual == expected` assert EXITS 1** (proves the
   anchor independence is load-bearing — the exact L95 trip-up); coupling meta-test restores it. The row asserts the
   coherent-rescale caught-by is EXACTLY the 4 INV5 anchor labels.
2. **The coupled-moments row (`coupled_..._BT_fixed` → `{INV2,INV3,INV5_Kbar0,INV5_Kbar2,INV5_Kbar4}`).** ⚠ **Ablate the
   ACTUAL side, NOT the oracle:** perturb the coupled-mutation FIXTURE (e.g. also scale BT, or drop one moment scale) so its
   COMPUTED caught-by row drifts from the immutable expected → the row's `actual == expected` assert EXITS 1; coupling
   meta-test restores it. (This row also documents that INV4 stays 0 under coupled scaling — a redundant-diagnostic
   observation, not an independent tooth.)
3. **Each single-parameter mutation (9 rows) — its own EXIT-1 ablation.** K̄₄→8/14 (INV1+INV5_Kbar4); K̄₄ sign-flip
   (INV1+INV5_Kbar4); K̄₂→7/5 (INV1+INV3+INV4+INV5_Kbar2); K̄₀→55/5 (INV1+INV2+INV3+INV4+INV5_Kbar0); denom 27→26
   (INV2+INV4+INV5_denom); corpus 9→8 (INV3+INV4+INV5_corpus); exp 5/2→3/2 (INV3+INV4+INV5_exp_k2); exp 3/2→1
   (INV3+INV4+INV5_exp_k0); BT 2/5→3/5 (INV2+INV3+INV5_BT). For each: perturb the ACTUAL side (the mutation fixture or the
   residual construction it fires) with the expected oracle immutable → the row's `actual == expected` assert EXITS 1;
   coupling meta-test restores it. ⚠ Note the caught-by rows list INV4 alongside INV2/INV3 (the redundant diagnostic ≡
   INV2−INV3) — this is source-faithful transcript content, NOT INV4 as an independent tooth.
4. **The no-float guard tooth (independent).** Inject a `sp.Float`/machine-real into a residual → `assert_no_float` EXITS 1
   (both engines, raw + reduced); coupling meta-test: remove the float → it stops firing. The exact-rational discipline.

**Review / acceptance GATES (verified by the tri-review agent — NOT runtime EXIT-1 residual teeth; do NOT tally as teeth):**

- **G1. The `.wl` per-function INDEPENDENCE + arity integrity (§1b).** Confirm the re-authored `.wl` computes the 12
  zero-tests + the caught-by matrix by a materially-different route (not the `.py` config-dict/`FullSimplify` mirror);
  def/call arity scan + unevaluated-leakage transcript scan. (A review/authoring gate, not an ablatable residual assert.)
- **G2. ⭐ The honest-scope gate (no over-claim).** Confirm the transcript prints the CALIBRATED-consistency scope
  (`G=GENUINE_BLOCKED`, NOT a `Γ̄₅`/`G` derivation) — the 020 CALIBRATED-not-PASS discipline; do NOT dress the consistency as
  an earned magnitude. (A review print-content gate, not an EXIT-1 residual tooth.)

**⚠ DE-COUNTED (retained for fidelity/transcript, NOT independent teeth):** baseline-all-zero (`all(rᵢ==0)`, the aggregate
positive — subsumed by the per-residual checks, a corrupted BASE coefficient trips an earlier residual assert first); the
`actual`-non-empty check (entailed by `actual == expected` since every expected row is non-empty); and INV4 (≡ INV2−INV3
identically — a redundant cross-form diagnostic, headline 3). These are NOT tallied as independent able-to-fail teeth.

⚠ **NOT 028 (do not build):** any re-derivation of the port numerator `N0_den` (024's), the port checks (027's), the
DtN fingerprint series (018's — 028 CITES the `27`), or the moments themselves (calibrated closure inputs; `G` blocked).
028 is the consistency artifact ONLY.

## §5 Consumed / exported

- **Consumes (PROVENANCE):** 027's exported K̄ moments `{K̄₀,K̄₂,K̄₄,Γ̄₅}` (the A3 SLOT; INV1/INV2 SHARED with R46, not
  double-counted) + 018's `27` fingerprint (INV2's `/27`) + 020's `54/5` / `G=GENUINE_BLOCKED` (the calibrated magnitude) +
  the corpus's OWN form `9K̄₂^{5/2}/K̄₀^{3/2}` (`4d_2_5pn.tex:469`, INV3, imported) + `c` (GR bridge) + `a` (CONV).
- **Exports (→ 029 + Part VII):** the 2.5PN consistency landing (the `research/4d_2_5pn` single open item MET at
  reduced-closure). ⚠ **Cluster C continues 029** (PN corpus DOI-cite — a THIN cite-only stage note, NO scripts) → then
  the scheduled **MIDWAY KNOB AUDIT** (Part-II gravity sector CLOSES — the pathA_40 `Δr=2` codimension dry-run over Parts
  I–II + the held-out vs irreducible-route-less tally).

## §6 Register expectation (orchestrator authors; CONFIRM at register + Codex-verify)

- **⭐ ZERO new counted CALIB knobs** (028 is a consistency slice — it CONSUMES 027's K̄ moments + 018/020's provenance,
  introduces NO new physical symbols). The K̄ moments are calibrated functions of `{G, c_s, a, c}` — `G=GENUINE_BLOCKED`
  (registered), `c_s`=R1, `a`=`CONV`, `c`=the GR-units bridge (benchmark, registered at 020). Part-II CALIB set stays
  **6** (`{μ_η,T_w,β}`(013) + `{Vp0/ℓ_c}`(015) + `{T_Ω,β₂}`(017)). Like the proof/consistency legs 020/021/024/025/026/027
  — ZERO new knobs.
- **⭐ ZERO new dims** (the K̄ moment dims are 020's/027's; 028 checks exact-rational residuals in z/moment space, no new
  units-restored quantity).
- **New structural edge (propose R47):** the full 2.5PN Burke–Thorne match-back consistency over the consumed K̄ moments —
  INV1 (`K̄₄K̄₀=4K̄₂²`) ∧ INV2 (pathA_43 form → `2G/(5c⁵)`) ∧ INV3 (corpus form `9K̄₂^{5/2}/K̄₀^{3/2}` → BT) ∧ INV4
  (cross-form agreement) ∧ INV5 (independent literal anchors defeating the coherent-rescale rig), dual-engine + able-to-fail
  (the 11-mutation matrix). ⚠ **SHARED with R46** (027's port-closure slot = INV1+INV2 re-expressed) — 028 adds INV3/INV4/
  INV5 + the coherent-rescale matrix; the port-closure consistency is counted ONCE (R46). Structural (a consistency/
  provenance edge, like R43/R44/R45/R46); discharges NOTHING at the knob level. ⚠ **CALIBRATED-not-derivation:** the match
  is a CONSISTENCY over the calibrated moments (`G=GENUINE_BLOCKED`), NOT a first-principles `Γ̄₅`/`G` derivation (Gate-6/
  sim-deferred).
- **⚠ CONFIRM at register:** the corpus factor `9`/exponents `5/2,3/2` are the corpus's OWN imported form literals (INV3),
  NOT new knobs; the INV5 anchors are 028's independently-restated constants (the anti-rig), NOT new knobs. No double-count
  of 018/020/027's registered symbols.

## §7 Deliverables (per the calibrated pipeline)

`ledger_stage028_2_5pn_matchback_{sympy_audit.py,mathematica_audit.wl}` (Codex builds) + self-contained note
`notes/stages/ledger_stage028_2_5pn_matchback.md` (orchestrator; full inline — INV1–INV5, the corpus open item, the A3
boundary with 027, the coherent-rescale anti-rig + INV5 independence, the CALIBRATED-consistency scope `G=GENUINE_BLOCKED`)
+ paper card `paper/stages/stage_028.tex` + `\input` (in `paper/appendices/stage_appendix_part02.tex`) + registration
(provenance index / coverage / manifest, count → 28) + parameter register update (R47, ZERO new knobs, SHARED with R46) +
Codex-verify + PDF + commit + docs/memory sync. ⚠ **Cluster C continues 029 (PN DOI-cite) → then the scheduled MIDWAY
KNOB AUDIT (Part-II gravity sector CLOSES).** ⭐ **028 is the CONSISTENCY cap of the pathA_43 density-port fold — the last
computed stage of the Part-II gravity sector (029 is a thin cite-only note; the MIDWAY KNOB AUDIT is the sector-close
checkpoint).**
