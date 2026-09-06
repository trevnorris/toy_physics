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
**Verified: leading-order rep-invariance HOLDS; the remnant is σ_W-sector only.**
Command: `_measurements/S11c_c2_adjudication_verify_EG.py` → `…_verify_EG.out`. Independent test: does
`REP_INVARIANCE_RESIDUAL.subs(sigma_W→0) == 0`? **Literal result (both densities, all 6 blocks):**
`full_zero=False` but `sigmaW->0_zero=True` — setting σ_W→0 annihilates the entire residual ⇒ **every**
non-invariance term carries a σ_W factor. So the leading O(ε) and O(εη) representation-invariance holds exactly;
the non-invariance is confined to the **σ_W drain sector** (the known c1 drain-projection O(σ_W²) carry-forward;
full evaluation matches the ≥64 GB deferred rep-invariance family). **⇒ not a load-bearing c2 defect;** surfaced +
carried forward. The step record must claim rep-invariance only at leading order, with the σ_W remnant open.

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

## VERDICT — the c2 fold physics is SOUND (zero confirmed physics defects)

The fresh-Claude leg's verdict holds. Grok's four "CONFIRMED" flags are all false positives (F, G) or a label
dispute where leading order is clean (B is F's residue; E's leading order cancels). Both legs + my verification
agree the load-bearing fold construction is physically correct.

**Carry-forwards / caveats (step-record interpretation, ⛔ NOT build defects):**
1. **F:** the emitted uniform-limit object is non-vanishing purely from the §3c `−extract(open)` bare-slot
   representation. Interpret as **"the genuine closure-induced coupling decouples,"** ⛔ not "the increment
   vanishes." ⇒ **light spec clarification owed** to `S11c_c2_SHARED_PHYSICS.md` §5e (the "increment must vanish"
   wording) + §3c (the increment carries the open-slot O(ε) piece by construction) — a wording fix that does not
   change what is computed; folds under spec review-until-clear.
2. **E:** leading-order (O(ε), O(εη)) representation-invariance ESTABLISHED; the σ_W-sector remnant is SURFACED +
   deferred (drain-projection O(σ_W²); ≥64 GB rep-invariance family). Record must not claim full rep-invariance.
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
