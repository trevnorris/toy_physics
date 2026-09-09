I'll ground this in the reconcile question, comparator result, committed N6, and the collapse-directive saga before recommending a path.The comparator picture is clear: matched zeros are trivial, and map 4 is the recurring failure. Next I’ll check the engines, S11c-d, and the c1 precedent for a cleaner instrument.I have the comparator and map-4 history. Next I’ll inspect the thickness convention in both engines and what S11c-d actually consumes.The map-4 commutator is structural, not a documentation bug. I’m checking one last point: whether Φ even sees `R_W`, and what S11c-d actually binds.**Recommend B — stop the collapse instrument, surface the graded operands as UNDECIDED, write the step record.** Do not spend v4 on this dictionary. Path A has a worse correctness trap than B’s incompleteness. Path C spends the hours on the broken piece and drops the maps that actually commute.

This is not “we got tired.” The third map-4 finding is a structural commutator. A coefficientwise-per-grade production rewrite of `R_W` cannot exist under the retained-order contract you already froze.

---

## 1. Path: B. Why.

What N6 actually certified is already in hand, dual-engine, and is the object S11c-d cares about:

- Per-engine `R_cov = 0` (Reading B, committed, `S11c_c2_N6_RESOLVED.md`).
- Blind WL independently reproduced the same vanishing, with knives (`STATUS.md`, WL rebuild).
- Comparator: `N6COV_R_COV` 160 ZERO / 0 NONZERO on matched keys; carrier bridge 320/0; support 400/0 (`run_tally.txt`, `run_data.md`).

Those matched zeros are `(0)−(0)` and do **not** prove operand agreement — the question doc is right about that (`reconcile_question.md` §2). They **do** prove both engines found the same invariant statement of covariance. That is N6’s load-bearing claim (`S11c_c2_SHARED_PHYSICS.md` §5c; `S11c_decisions.md` N4/N6: agree after the field redefinition, not identical raw components).

S11c-d binds the **SymPy closed operators** (`CLOSED_SLAB_OPERATOR`, `CLOSED_COUPLING_KERNEL` — `S11c_c2_sympy_build_directive.md:199–205`). It does not bind N6 residuals, source leaves, or Φ. WL N6 is a control engine, not a second production export. Operand collapse is corroboration of channel (b), not a d input.

The collapse instrument was supposed to turn that corroboration from “weak (trivial zeros)” into “strong (operands match under a frozen bridge).” Four directive rounds have shown the proposed bridge is the wrong object for the emitted leaves. Continuing is not “one more fix”; it is insisting on a rewrite that does not live at the stage the comparator emitted.

c1 is the right precedent, and you already used it for the stronger object: DtN **kernel** earned AGREE; raw `DTN_OPERATOR`, ENERGY, `t_s`, seal-5 density stayed UNDECIDED (`S11c_c1_comparator_reconcile.md:17–22, 244–247`). Here the kernel-analogue (`R_cov`, carrier bridge) already vanished. The 40/76/18 are the whole-form analogue.

---

## 2. What was missed

### The graded coefficient table is the wrong instrument for channel (b)

Map 4 as an **ungraded constructor identity** is correct. Opus independently re-derived all 19 rows, including E14’s triangular mixing (`collapse_directive_gate.md` v3; Grok’s v3 review § Map 4). Grok then found, and you verified, that those `R_W` rows do not commute with grading.

Worked example, already in Grok’s report, grounded in both engines:

- WL density term: `cCoupling · WBg · θ · eW` (`.wl:217–222`).
- SymPy: `C · W_bg · θ · E` with `E = W_0 e_W / W_bg` (`brane_operator_sympy_audit.py:1589–1590, 1631`), which **cancels to** `C · W_0 · θ · eW` **before** grading.
- Ungraded map `cCoupling → (W_0/W_bg)·C` is right.
- Both engines then expand `WBg → W_0(1+η·profile)` **before** grade extraction (WL `profileRules` + `finish`, `.wl:127–141`; SymPy `grades()` `xreplace(inputs.profiles)`, `diagnostic_sympy.py:245–248`).
- Graded η⁰ leaves: WL `cCoupling · W_0 · θ · eW`, SymPy `C · W_0 · θ · eW`. Applying live `R_W` reinserts `W_bg` that neither leaf contains.

That is not a documentation omission. **`R_W · W_bg = W_0` is a product of two η-dependent quantities.** After Series, the identity mixes grades:

\[
(R_W)_0 (W_{\mathrm{bg}})_0 = W_0, \qquad (R_W)_0 (W_{\mathrm{bg}})_1 + (R_W)_1 (W_{\mathrm{bg}})_0 = 0.
\]

WL’s η¹ piece of E02 is the measure-expansion `cCoupling · W_0 · w1 · θ · eW`. SymPy never had that term: `W_bg` cancelled in the constructor. The ungraded identity removes it (`R_W · WBg = W_0`, constant). A coefficientwise-per-grade rewrite **cannot**, unless it mixes grades.

The retained-order contract forbids that mix (`reconcile_question.md` §5; directive grade-combining tripwire; the L-CAS `σ_W→0` cousin you already stripped from §4). So the v4 options on the table are:

| v4 idea | What actually happens |
|---|---|
| Per-grade images, no η-expansion of `R_W` | η⁰ can map `cCoupling → C`. η¹ E02 **does not collapse** (the leftover you already flagged). Same for every `R_W` row (E02, E03, E04, E05, E13, E14). |
| Expand `W_bg` on both operands, then compare | Zero-jet freeze `W_bg → W_0` is `R_W → 1`. Same as the row above. Full `W_bg → W_0(1+η w1)` **is** expanding `R_W` in η. Tripwire. |
| Apply map-4 then grade, on ungraded atoms | The comparator did not emit ungraded atoms. Reconstructing `∑ η^i σ^j coeff_{ij}` from **already-expanded** grades, then multiplying by `R_W`, is the same stage error. Doing it properly means a new engine emit path. |

There is no sound production rewrite of `R_W` that is simultaneously (i) coefficientwise on the four independent grades, (ii) free of live `W_bg` reinsertion, (iii) able to cancel the measure-expansion terms. That is why map 4 has now drawn three *distinct* findings (not frozen → needs scales → graded stage). Round 4 will be a 19×4 table that is a disguised grade-mix, or an honest UNDECIDED you could have declared today.

The recurring difficulty is a mis-decomposition signal. You are trying to compare **coordinate components** of μ after a field redefinition `E = (W_0/W_bg) e_W` has been pushed through EL, source-solve, *and* the grade projector. The invariant objects were already compared: `R_cov` and `C_E − C_M`. They vanished. The leftover is the E-vs-`e_W` frame, which is a second representation axis sitting on top of Eulerian-vs-material. Reading B already says components need not match across a field-variable frame (`RESOLVED.md:10–16`). Applying that to thickness as well is the same physics, not a new gap.

### Cleaner constructions (and why not to take most of them)

**Compare ungraded constructor coefficients.** That is what the 19-row table *is*. Opus already did it as an independent re-derivation. A script that reprints `cCoupling` vs `R_W·C` from those sites is a table-identity check, not a source-operand test. Useful as a record citation; not worth a G1 cycle.

**Re-emit pre-expansion operands, then apply ungraded map 4, then grade the residual.** This is the unique *sound* collapse of channel (b). It is an engine-scope change (new emit path in both constructors), not a v4 dictionary tweak. Cost ≥ Path A, plus two engine reviews. Do not do this unless you later decide channel-(b) operand AGREE is load-bearing for a specific d claim. It is not, today.

**Numeric cross-engine PIT on constructors** (same numeric jets, live `W_bg`, compare source numbers). Sidesteps the symbolic table, still needs a reviewed dual-engine probe that *calls constructors*, not graded leaves. Similar cost, same “nice-to-have” status.

**Absorb the thickness convention upstream.** Do not. That is folding engines toward agreement (M1). The disagreement is the measurement.

**Convention-invariant object.** Already have it: `R_cov` and the carrier bridge. A weak-form pairing `μ_e δe` vs `μ_E δE` reintroduces `δE = R_W δe + e δR_W` and the same commutator.

**Φ is a different, actually well-posed slice.** Map 4 is scoped *off* Φ (`directive` typed-scope table). Both engines map `e_W` / `eW`, not `E` (`covariance_sympy.py:75–76`; `.wl:236–239`). The 18 `FROZEN_PHI` residuals do not need `R_W`. A Φ-only spelling collapse (maps 1/3/5) could corroborate caveat 1’s “same declared map, independently constructed” half. It would **not** close “is Φ physically correct from the motion” (`RESOLVED.md:39–40`). I would not spend a G1 cycle on it unless you specifically want that corroboration; caveat 1 stays carried either way. If you do anything extra, this is the only cheap slice that is not map-4.

**Census (§3d) does not need the instrument.** Eight leaves, frozen 6-row crosswalk already in v3, WL-only `THICKNESS` / `MATERIAL_NORMAL` stay one-sided. That is a fact-lookup in the step record.

### Is per-engine covariance enough for c2 → d?

Yes, for what d binds. Carry the three caveats open, as `RESOLVED.md:47–55` already required regardless of disposition. Do not pretend matched `(0)−(0)` is operand AGREE — the step record should say the opposite, the way c1 said “kernel closed, whole-form UNDECIDED.”

N6’s job is: the increment is the same operator in two representations. Both engines, independently, with biting knives, found `R_cov = 0` at retained order (SymPy δ ≈ 2.6e-22; WL own primes/bound). Cross-engine, every *matched* covariance and carrier-bridge key vanished. That is sufficient to hand d the closed operators with N6 marked covariance-satisfied, operands representational-unadjudicated.

### Correctness traps

**The trap in A is false AGREE.** Applying `R_W` at the wrong stage, or mixing grades to kill the η¹ measure terms, manufactures collapse of a stage artifact. That is the same class as L-CAS (a scale that does not live at the retained-grade stage, used as if it did). Four rounds of map-4 findings is the architecture telling you the dictionary is becoming the answer. Opus OBS-2 already warned: do not weight rename-only collapses (`collapse_directive_gate.md`). A v4 that “lands mostly collapsed” by freezing `W_bg → W_0` on η⁰ and quietly dropping η¹ E02 is an over-clear.

**The trap in B is not “hiding a disagreement A would call DISAGREE.”** Your own question forbids that inference: non-collapse is UNDECIDED unless bridge completeness is independently proven (`reconcile_question.md` §3b, §6; question-vet finding 2). Even a perfect Path A run that left velocity / source-solve / live-density in the residual — all **excluded** from the dictionary (directive Excluded block; map 6 removed) — would still be UNDECIDED on those pieces. A does not close two-errors-cancel on V or source-solve. `RESOLVED.md` already keeps V open (`V_E ≡ V_M` is builder agreement, not a derived transform).

**What B actually forgoes:** leftover *shape*. A spelling-only pass (maps 1/2/3/5/7/8, no map 4) would print whether the 76 source residuals *are* the predicted `R_W`/measure-expansion pattern or something else. That is the only two-errors-cancel catch still available from comparator leaves, and Path C throws it away to chase `R_W`. I would not build that pass now; I would record that the leftover was not inspected for shape, so a later reader must not upgrade UNDECIDED to “we know it’s only thickness.”

**Two-errors-cancel on constitutive physics is already bounded.** If the energy bases differed by more than the verified 19-row ungraded table, Opus’s re-derivation would have failed. Shared-Φ error is the remaining caveat-1 risk; WL-blind vs SymPy-N4-prolonged is independent construction of Φ, but without collapsing Φ images you have not shown those maps agree as *expressions*. That is a real gap. It is the same gap `RESOLVED.md` already told you to carry, and it is not closed by map 4.

**Path C is the worst of the three.** It still needs the impossible stage fix, tests only the non-commuting rows, and drops the maps that *do* commute. If you were going to build anything, you would drop map 4 from production, not keep only map 4.

---

## 3. If you ignore this and do A or C anyway

The map-4 stage “fix” as specified (per-grade images *or* consistent expansion, never expanding `R_W` in η) is **internally consistent and will not collapse the load-bearing `R_W` rows at η¹.** That is not a review miss; it is the identity. E14 is worse (`R_W · (γ_{14} − κ_{θW}/W_bg)`: two `W_bg`s, chain rule, `σ_W` from `gW`).

The only deeper fix is ungraded-then-grade, which requires pre-expansion emits you do not have. Do not author a 19×4 graded image table by hand. That is the route-2 failure mode you already invoked rule 15 to escape.

If a v4 review comes back SOUND, treat that as a red flag: the legs probably checked the ungraded table again (Opus’s false-negative) or accepted a grade mix.

---

## How to write B so it is honest, not a cop-out

In the c2 step record, state all of the following; do not collapse them into “cross-engine UNDECIDED”:

1. **Earned:** per-engine operator covariance, both engines; matched-key `R_cov` and `C_E−C_M` vanish (trivial as residual-of-zeros, still dual-engine confirmation of the *vanishing*); support agrees where matched; Reading B stands (`d21c8ff5`).
2. **Not earned:** operand AGREE on `C_E`/`C_M` (40), `SOURCE_*` (76), Φ (18). Schema-unmatched `R_N6` / channels / RC sources stay UNDECIDED as schema, not heaviness (`reconcile_question.md` §1.iii).
3. **Why the collapse stopped:** the load-bearing representational identity (`R_W`, thickness `E` vs `e_W`) is an ungraded constructor fact that does not commute with both engines’ grade extraction. Continuing the graded dictionary would either reinsert `W_bg`, mix grades, or leave the same residuals UNDECIDED. Cite Grok v3 MUST + your G4 verification (`collapse_directive_gate.md`).
4. **What *was* checked about constitutive physics:** the 19-row ungraded energy-basis table is independently re-derived and correct as a constructor identity (Opus v3). That is not operand collapse. Do not upgrade it to “sources agree.”
5. **Carry open, as already required:** Φ physical correctness; V transform; extracted-block leakage; c1 ENERGY; 2 S11c-b signs; 6 §3d items; fix `I_{M→E}` terminology. No per-substep card.
6. **Census:** do the 8-leaf semantic walk in the record (junk spelling, WL-only `THICKNESS` / `MATERIAL_NORMAL`). Do not call it metadata.

Do not call this “weak N6.” Call it: **covariance dual-engine confirmed; constitutive *components* not put in a common thickness coordinate at graded order.** That is accurate. “Weak because we only have 0−0” undersells `R_cov` and oversells what a v4 collapse would have added.

---

**Bottom line:** Path B. Kill the collapse-instrument series at v3 NOT-SOUND. The hours buy either a stage artifact, a tripwire violation, or the same UNDECIDED you can write this afternoon. S11c-d does not need the operands collapsed. The two-errors-cancel story that makes A feel mandatory is mostly already answered by the ungraded table plus dual-engine `R_cov`; what remains is excluded from A’s dictionary too. The real A-risk is an over-clear, not a missed disagreement.
