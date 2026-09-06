# S11c-c2 self-energy fold — physics review adjudication (gate record)

**Artifact reviewed:** `research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (1053 lines;
astra/`gpt-6-astra` build). Emitted output `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out` (499 MB), navigated by
`_measurements/S11c_c2_sympy_stdout_index.json` (154-tag byte-offset index).
**Role:** script / physics-bearing build → **review-until-clear** (G2/G4), not one-pass.
**Authorship → legs (G1):** script is Codex(astra)-written ⇒ two legs = **fresh general-purpose Claude agent +
Grok**, orchestrator-launched, identical prompt. My own verification (below) is required (G4/rule 13) but is
**not** a leg.

**Identical prompt:** `directives/_legs/S11c_c2_physics_review_prompt.md` (13,386 bytes).
**Blindness by absence (S1):** for the review's duration, astra's builder report + provenance narrative were
relocated to `ext_logs/c2_review_hold/`; the legs saw only the mechanical byte-offset index. astra's ~36
quarantined self-review files stayed out of tree (`ext_logs/c2_builder_overstep/`) and were disregarded. Held
files restored before adjudication.

**Handed inputs:** the script; the physics authority `directives/S11c_c2_SHARED_PHYSICS.md`; the real imports
`scripts/{S11c_b_exports.py,S11c_c1_exports.py,ledger_fold.py}`; predecessor records
`steps/S11c_c1_curved_bulk_closure.md`, `steps/S11c_b_variable_coefficient_operator.md`; the emitted `.out` + its
index. **Withheld (by not handing):** the build directive (its frozen symbol map would make the review a fidelity
check), the builder report, the provenance narrative, astra's self-reviews.

**Leg reports (separately attributed):**
- fresh Claude: `_measurements/S11c_c2_physics_claude_leg.md` (verdict: 0 confirmed physics defects; E unsettled/
  deferred, F not-a-defect + caveat).
- Grok (`grok-4.6` high): `_measurements/S11c_c2_physics_grok_leg.txt` (verdict: fold *wiring* right, but B/E/F/G
  flagged CONFIRMED; "the increment as defined is not the S11b-decoupled self-energy the spec claims").

---

## The three disagreements — adjudicated by MY OWN computation (rule 13; legs wrong in both directions)

Both legs independently agree A, C, D1–D6 and the fold **wiring** (δp not J_s, operator inverse, kernel bridge,
V_s→face_velocity, computed jets, ε-strip) are SOUND, each with shown CAS evidence. The legs SPLIT on B/E/F/G. I
resolved each myself; scripts + literal stdout committed beside this record.

### F (uniform limit) — Grok CONFIRMED "genuinely nonzero coupling"; Claude "genuine coupling decouples"
**Verified: Claude is right; Grok's F is a FALSE POSITIVE.**
Command: `_measurements/S11c_c2_adjudication_verify_F.py` → `…_verify_F.out`. Method (independent of the leg's
`Trial`-label heuristic): the increment = extract(close) − extract(open); extract(close) substitutes ALL δp
slots ⇒ carries **no** bare `delta_p_±`; every bare-`delta_p` term therefore comes only from `−extract(open)`.
Zeroing the bare slots leaves `extract(close)|_uniform`; `.doit()` the residual Integrals; test == 0.
**Literal result (all 4 cases, both nonzero blocks THETA & E_W):** the surviving `closed_part` is
`coeff·Test·Integral(0, …)` — integrand **literally 0** — `closed_part_zero(doit)=True`. So the genuine
closure-induced coupling is **identically zero in the uniform limit**: the transverse↔thickness coupling
decouples exactly as required. Grok saw `Integral(...)` present and reported a surviving `Z₀·μ_θ` coupling
**without evaluating that its integrand is 0** (`expand` leaves `Integral(0,…)` unevaluated; only `.doit()`
collapses it). (Grok's c1 "tangential freeze" was likewise a false positive — a pattern.)
**Residue (real, not a defect):** the nonzero part of the emitted uniform-limit object is the
`coeff·(δp_-+δp_+)·Test` bare open-slot bookkeeping from `−extract(open)` (δp is a free symbol in the open
operator). **⇒ the increment does not literally vanish; the genuine coupling decouples.**

### B (ε-grade) — Grok CONFIRMED "increment not O(εη)-only"
**Same object as F; not a physics defect.** Both legs agree the strip-one is correct and there is no O(1) and no
O(ε²) term. The O(ε) (η-free) grade Grok flagged **is** the `−extract(open)` open-slot bookkeeping (the F
residue), i.e. a property of the §3c increment representation, not a fold error. The genuine induced self-energy
is O(εη) (it decouples at η→0, verified in F).

### E (N6 rep-invariance) — Grok CONFIRMED "does not vanish (σ channel)"; Claude "leading order holds, σ_W deferred"
**⚠ CORRECTED 2026-09-06 (Codex-sol compact-prep verify, rule 13): E is OPEN — my earlier "false positive / cleanly
deferred" was an OVER-REACH; Grok's E is NOT disproven.**
Command: `_measurements/S11c_c2_adjudication_verify_EG.py` → `…_verify_EG.out`. The test proved only:
`REP_INVARIANCE_RESIDUAL.subs(sigma_W→0) == 0` for every block, both densities (`full_zero=False`,
`sigmaW->0_zero=True`) — every non-invariance term carries a σ_W factor, so leading O(ε)/O(εη) rep-invariance holds
**in the σ_W→0 projection**. ⛔ **That is confinement, NOT full N6.** The retained remnant sits at grades
`(1,0,1)=O(εσ_W)` and `(1,1,1)=O(εησ_W)` — **LINEAR in σ_W**; and σ_W is an **independent first-shape-order
bookkeeper** alongside η (`S11c_c1_SHARED_PHYSICS.md:223`), so `O(εσ_W)` is **retained physics at the same order as
O(εη)**, NOT a deferrable higher term. My association with the **O(σ_W²) drain-projection** carry-forward was WRONG
(a linear-σ_W residual is not that second-order object), and the only explicit c2 ≥64 GB deferral is the full
cross-engine / four-giant-family work (`S11c_c2_SHARED_PHYSICS.md:39`), NOT this already-emitted single-engine N6
residual. **⇒ DISPOSITION: N6 rep-invariance is established only after projecting to σ_W=0; the retained
`O(εσ_W)`/`O(εησ_W)` residual is UNRESOLVED and may be a `representation_pullback` build defect.** It must be
resolved — as a pullback defect fix, or a separately-reviewed scope/spec decision — **BEFORE the WL build**, ⛔ not
filed as a cleared/deferred note. (The fresh-Claude leg's "UNSETTLED / possible pullback incompleteness" was closer
to correct than my "false positive.")

### G (adjointness) — Grok CONFIRMED "blocks not adjoint; suppressed check hides asymmetry"; Claude "honest omission"
**Verified: SOUND; Grok's G is a FALSE POSITIVE (refuted a claim the builder never made).**
Command: same `…_verify_EG.py`/`.out` (block directionality) + builder-report line 70.
**Literal result:** in `SELF_ENERGY_INCREMENT`, the THICKNESS→TRANSVERSE reverse blocks are **identically zero**;
only TRANSVERSE→{θ,e_W} is nonzero — the induced self-energy is **one-way**, physically expected since δp lives
only in the θ/mechanical rows and depends on the transverse velocity (closing induces coupling *into* those rows,
not back). `CLOSED_COUPLING_KERNEL` carries **both** directions nonzero (the pre-existing open coupling), and
**both blocks are emitted** ⇒ the asymmetry is visible, not hidden. Builder report line 70 states the actual
reason: *"There is no independent adjointness construction in this implementation, so no
SELF_ENERGY_ADJOINTNESS_RESIDUAL is emitted. Both off-diagonal blocks are emitted."* — exactly spec §3b ("emit
both blocks and say there is no independent second route rather than dressing a structural zero as a check"). The
builder never claimed "adjoint by construction"; Grok attributed that and called it dishonest. Directional ≠
defect; omission is correct.

---

## VERDICT (⚠ CORRECTED 2026-09-06) — the c2 fold WIRING + A/C/D + F are SOUND; E/N6 is OPEN

The fold **wiring** (δp not J_s, operator-inverse, kernel bridge, V_s→face_velocity, computed jets, ε-strip) and
the A/C/D1–D6 constructions are supported (both legs + my verification). **F** resolves narrowly (the genuine
closed coupling decouples; the literal §3c increment retains the open slot → a wording fix). **G** is a measured
**directional/one-way** SymPy increment (both blocks emitted; ⛔ NO dissipativity/passivity claim is established —
only directionality; a physics interpretation awaits the WL comparison). **B** is F's residue, not a defect.
⛔ **E/N6 is OPEN** (corrected E section above): rep-invariance holds only in the σ_W→0 projection; the retained
`O(εσ_W)`/`O(εησ_W)` residual is UNRESOLVED (possible `representation_pullback` defect) — Grok's E is NOT disproven.
**⇒ this is NOT "0 confirmed defects": the core fold is sound, but N6 has an open retained-order finding to resolve
before the WL build.** [The earlier verdict "physics SOUND / F/G/E all false positives" was an over-reach.]

> ⚠⚠ **SUPERSEDED 2026-09-06 (E) — see `_measurements/S11c_c2_N6_spec_adjudication_record.md`.** Re-grounding E the
> corrected way (Codex-vetted question → 2-leg spec adjudication) found the deeper truth: **c2 §5c MIS-SPECIFIED N6**
> — it compared the two *anchorings* (distinct physics), but N6 (parent S11c-a §5a / sibling c1) is Eulerian-vs-
> material-coordinate **within a fixed anchoring**. The current `REP_INVARIANCE_RESIDUAL` is the wrong object (a
> nonzero value is *expected*, not a defect); **c2's real N6 was never run — it is UNTESTED**. The "σ_W→0 residual /
> possible `representation_pullback` defect" framing above is retired. Spec §5c corrected (2-leg review-until-clear);
> a build correction (the material-coordinate route) + the F/G re-grounding remain owed. The biased `verify_EG` E
> conclusion is retired.

**Carry-forwards / caveats:**
1. **F** (wording, not a defect): the emitted uniform-limit object is non-vanishing purely from the §3c
   `−extract(open)` bare-slot representation. Interpret as **"the genuine closure-induced coupling decouples,"**
   ⛔ not "the increment vanishes." ⇒ **light spec clarification owed** to `S11c_c2_SHARED_PHYSICS.md` §5e (the
   "increment must vanish" wording) + §3c (the increment carries the open-slot O(ε) piece by construction) — a
   wording fix that does not change what is computed; folds under spec review-until-clear. ⚠ After the spec edit,
   the committed export's `BUILD_INPUT_DIGESTS` (which pins the spec) goes stale ⇒ regenerate or lawfully repin +
   reverify the export before advancing.
2. **E** (⛔ OPEN — NOT a cleared/deferred note): resolve the retained `O(εσ_W)`/`O(εησ_W)` N6 residual as a
   `representation_pullback` defect fix or a separately-reviewed scope/spec decision BEFORE the WL build. The record
   must NOT claim full rep-invariance (only σ_W→0-projected).
3. **G:** the induced self-energy is **directional** (transverse→{θ,e_W}; reverse block identically zero) —
   expected from δp residing only in the θ/mechanical rows; both blocks emitted; no adjointness residual (correct
   per §3b). Record it as a physics feature.
4. The two S11c-b sign conventions (face-generalized-force, #90 closure-fold) are cross-engine-UNVALIDATED and
   multiply the increment ⇒ the blind-Wolfram cross-engine comparator surfaces them (spec §3c/§7); the §3d.4
   mechanical-power pairing adjudicates the face-force sign. (Deferred to the Wolfram engine + comparator.)

**Next (unchanged by this adjudication):**
- Commit this reviewed baseline (the script + review artifacts) **before** any repair overwrites it (this commit).
- **STEP B — publication-only export repair** (still viable; physics sound): `publish`@807 / call@952 only —
  drop `SELF_ENERGY_INCREMENT` (+ operands, term_origins, parity, §3d, §5 controls) to EMIT-only; export **both**
  `S11CC2_CLOSED_SLAB_OPERATOR` + `S11CC2_CLOSED_COUPLING_KERNEL` (all 4 cases) transparent-factored + the casewise
  `canonicalize(expanded − decode(compact))==0` check. ⛔ do not touch build_case/extract/kernel/payloads.
- **STEP C** — re-review the export repair (publication-only) until clear → commit → then the light §5e/§3c spec
  clarification (review-until-clear) → the blind Wolfram engine → its 2 legs → c2 T7 comparator + reconcile → c2
  step record.
