I'll review the folded S11c-c2 N6 ablation-harness directive against the engines and the round-2 diff, starting with the review-leg method and the handed files.I'll verify the fold against the live engines and the failed-build digest, starting with the WL and SymPy sites named in the findings.Next I’ll check whether any SymPy knife is case-gated, and whether the digest/PIT contract is actually implementable on the engines.I need to confirm whether K_junk’s DEAD set still holds once the knife is live, and whether `n.sha` can actually hash the emitted PIT payloads.Checking one last liveness detail: whether K_rank’s rank≥2 branch can fire on the pinned SymPy case, and whether `n.sha` is defined on parsed JSON rather than on `Node`.# Independent review — S11c-c2 N6 ablation-harness directive, FOLD round 2

**Artifact:** `research/pde_ledger_v3/directives/S11c_c2_N6_ablation_harness_directive.md`  
**Method:** source-first against the four engines + failed-build WL digest; fold diff for frozen-knife / wording checks. No execution.

**Overall: SOUND.** The fold resolves F1/F2/F3, introduces no cert-wrong defect, and leaves the ten-knife design byte-identical.

---

## Q1 — F1 fully resolved, nothing else broken by the repin

**SOUND.**

**WL pin = `actualJunkCase`, and it is the only `case===` gate.**

- `audit.wl:21`: `actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"}` — identical to the new iterator at directive:122–124 and clause 5:53–54.
- `audit.wl:845`: `activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0]` — the only `case ===` in the file (also written through at :978, :985 as metadata, not a second gate).
- Under the new pin, `activeJunk = actualJunkCoefficient`. K_junk FORM `actualJunkCoefficient 0→1` therefore feeds `materialAmplitude[..., jk]` at :247 (`el[pulled] + jk junkMu eW`) and `muM` at :847.

**The other two WL knives are not case-gated and still bite.**

- **K_carrier:** `materialNormalKnife` is a top-level assignment (`audit.wl:17`). Used at `:860` `materialGeometry[anchor, s, materialNormalKnife]` with no `If[case===...]`. The FORM site `:289` `invT[[1, 4]] += knife jet["WBg", {2}]` is inside that builder. Source velocity is a separate call with `0` at `:861`, as the directive states.
- **K_split_route:** `sourceChannel` at `:892` is a straight assignment in `buildCase`; the affine flag is not case-gated.

So `{"MATERIAL_ADVECTED","RHO4_CONSTANT"}` is exactly the unique case in which all three WL knives are simultaneously live. `K_carrier` / `K_split_route` were already live on the old pin; the repin does not mute them.

**SymPy pins unchanged, and no SymPy knife is case-gated the way WL K_junk was.**

Directive:234, :331, :434 still pin `--anchoring LAB_HELD --density RHOBR_CONSTANT`.

| Knife | Site | Case gate? |
|---|---|---|
| cov K_junk | `ACTUAL_JUNK` :35 → `actual_amplitudes(..., ACTUAL_JUNK)` :187–188; `actual.append(kappa_j * junk * b.e_W)` :145 | None. Module constant, always appended. |
| cov K_circular | `run` amplitude block :186–188 | None. |
| cov K_rank | `prolonged_phi.image_of` :92–97 | None. `h` is LAB_HELD-only (`:74`), which is the pinned anchoring; rank collapse is on path length, not case. |
| diag K_EW_rowdrop | `face_factory` rows :311 | None. |
| diag K_slotdrop | `run` `slots` :788 | None. |
| rec K_normal | `build_material_carrier` :91 | None. |
| rec K_source_route | `run` `source` :260 | None. |
| rec K_operand_swap | `ds` :211–212 | None. |

Diagnostic `N4_ADVECTION` controls (`:860`) emit only for `RHOBR_CONSTANT`. The unchanged SymPy pin is the case that makes those certified control groups exist.

**No knife — WL or SymPy — is inert on its harness’s pinned case.**

(The covariance docstring at :14–15 *mentions* `ACTUAL_JUNK` for `MATERIAL_ADVECTED.RHO4_CONSTANT`; that is usage prose, not a `case===` gate.)

---

## Q2 — F2 fully resolved, no residual leak

**SOUND.**

**Carve-out is gone.** `grep` of the folded directive finds **zero** hits for `small scalar`, `direct payload`, or `already-compact`. Clause 1:31–35 now states `R_N6` / `R_cov`/`R_COV` / `SPLIT_CHECK` / slot-closure guard residuals are certified objects: compact PIT fingerprint (tally + bitmap) + harness digest, never the raw tag/table.

That matches the engines:

- WL: `addNumeric["RC","R_N6", rawResidual]` and `SPLIT_CHECK` at `:897–899`; `addNumeric["COV","R_COV", covDelta]` at `:883`; guards via `numericObjects`. Every numeric payload is wrapped at `:804–809` with `"ARITHMETIC" -> (circuitExpression /@ values)` and `PROBE_NUMERATORS` SparseArray — a heavy PIT object, not a scalar.
- SymPy: `add('R_N6', R)` / `add('SPLIT_CHECK', ...)` (reconcile :270, :273) and `('R_COV', rcov)` (covariance :201) are tables passed to `n.pit` (`:668`, emit at `:769–772` with `numerator_denominator` + `nonzero_modular_numerator`).

**Surviving “direct payload” allowance cannot be satisfied by those objects.** Clause 1:34–35: only a genuine scalar integer/rational with **no** `ARITHMETIC` or table payload. `R_N6`/`R_COV`/`SPLIT_CHECK`/guards all have both.

**Remaining paths checked:**

| Location | Leak? |
|---|---|
| Clause 1 | No. Fingerprint = column keys + bitmap/tally only; `numerator_denominator` and DAG forbidden. Full payload is hashed then discarded. |
| Clause 4 | No. Spool stays in `/tmp`; error path may not echo a full symbolic object, DAG, or raw engine transcript. |
| Clause 5 | Budget/pin only. |
| Clause 8 | Failure diagnostics bounded: exit, stderr tail/size, spool SHA-256, missing tags — “never the raw object-bearing stdout or a symbolic payload.” |
| Deliverables :98–101 | Explicitly forbids full symbolic object, symbolic difference, arithmetic DAG, PIT sample matrix, **raw residual table**, raw engine transcript. |
| Per-harness certified/DEAD lists | Name tags only; printing is governed by clause 1. `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` are the named non-PIT objects (coverage note :578–579): digest + genuine scalars, not the substitution map / coverage table. |

The failed-build WL driver at `ablation_harness_wl.py:180–181` still dumps full stdout on subprocess failure; the **new** contract (clause 4+8) forbids reproducing that. Same for the old Python worker `print(raw.read_text())` at covariance harness :156–159.

The `e.g. R_N6 -> 18/288` in clause 1:33 is a fingerprint-shape illustration, not an acceptance test and not a path for the raw object.

---

## Q3 — F3 faithful + implementable

**SOUND.**

**Well-defined on both engines.**

- **WL.** Production emit is `WriteString[$Output, name, " = ", ToString[..., InputForm], "\n"]` (`audit.wl:52–53`) or a single streamed `name = <|...|>` line (`:59–73`). The failed-build driver already does `name,sep,payload=line.partition(' = '); emitted_hashes[name]=sha(payload.strip().encode())` (`ablation_harness_wl.py:188–190`). Hashing that assignment text (directive:24–26 includes `ARITHMETIC` and `PROBE_NUMERATORS`) is implementable. `emitted_object_sha256` / `RUN_PROVENANCE` have **zero** hits in `audit.wl`; they are harness records, as the fold now states.
- **SymPy.** `n.sha` is `diagnostic_sympy.py:91–92`: `hashlib.sha256(sp.srepr(c.cas(e)).encode())`. `c.cas` (`selfenergy_fold_sympy_audit.py:63–72`) walks `str` / `bool` / `Mapping` / `list`/`tuple` and otherwise `sp.sympify`. A parsed JSON PIT payload (`columns`, `numerator_denominator`, `nonzero_modular_numerator` from `pit` `:769–772`) is JSON-native and cas-able. An in-memory `Node` (`:107–112`) is not — `sp.sympify(Node)` is exactly the failure the fold names. `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` are ordinary JSONL `data` dicts (`covariance_sympy.py:105–124`); hashing the parsed emit (not a `Node`) is the same path. Clause 1:26–28 (“full pre-projection emitted/parsed payload is hashed and then discarded”) covers the non-PIT pair; they are not left without a digest.

**Same-support witness.**

- WL full assignment text includes `ARITHMETIC` (`:804`) **and** `PROBE_NUMERATORS` (`:808`). A FORM/COEFFICIENT change that keeps the nonzero support but retargets the circuit or the modular values changes that text → digest moves.
- SymPy hashed PIT payload includes `numerator_denominator`, which clause 1 then forbids printing. Same-support with different sample values changes that matrix → `n.sha` moves. The compact fingerprint (bitmap/tally only) would not.

**True blind spot against the frozen knife list?** Not one that makes a live knife invisible.

A digest-only object is `FROZEN_PHI` / `PHI_DOMAIN_CENSUS`. The live knife that must move them is cov **K_rank**, which rewrites `images[atom]` (`:92–97`) and therefore `substitution_map` inside the emitted `FROZEN_PHI` payload (`:120`). That string changes; the digest must move. K_junk / K_circular list those objects as DEAD and do not touch `prolonged_phi`.

A PIT object whose only *additional* witness is the digest (same-support bite) still has `numerator_denominator` (SymPy) or `ARITHMETIC`+`PROBE_NUMERATORS` (WL) inside the hashed payload. Every frozen FORM/COEFFICIENT on a live certified object changes one of those (coefficient ×2 is the textbook same-support case). A DAG change that is also PIT-value-identical would be invisible to SymPy’s chosen hash, but that is the engine’s own PIT witness, not a harness hole, and it is not a specific frozen knife.

---

## Q4 — no new defect; knives frozen

**SOUND.**

**Knives frozen.** The round-2 diff touches only: clause 1 (output contract), clause 5 (WL pin + independence sentence), Deliverables sentence, Harness 1 iterator + explanation, CHANGE-LOG. All eleven `### K_*` blocks (directive:167–535) — scope/site, Old, FORM, NO-OP, COEFFICIENT ×2, DEAD — are absent from the diff.

**No new inconsistency.**

- Clause 1, coverage note :573–582, and CHANGE-LOG :593–602 all say: harness-computed full-payload digest; `R_N6`/`R_COV`/`SPLIT_CHECK`/guards through compact PIT; WL pin = `actualJunkCase`; SymPy pins unchanged.
- “Unique all-three-live case” is correct: K_carrier and K_split_route are live on every case; K_junk is live only on `actualJunkCase`.
- Budget is still one case, `drawCount = 4`, `timeout --kill-after=5 600`. Not a regression. `caseOrdinal` seeds will differ from the four-case production ordinal of this case; directive:149–151 already requires printing actual seeds and comparing only within the restricted schedule.

**Opus Q3 sampler caveat is correctly scoped.** Joint-rejection sampler `audit.wl:774–791` can change the accepted sample set when a knife changes singularity structure, so a DEAD object’s **digest** can move with an unchanged object (false positive on one-sidedness, never a false negative on a live bite; nonzero tally is the robust channel). CHANGE-LOG :603–607 labels this a **later executable BUILD review** note, not a directive defect. Agreed.

---

## Overall

**SOUND.**

F1: WL K_junk is live on the repinned case; K_carrier and K_split_route remain live; every SymPy knife is live on `LAB_HELD`/`RHOBR_CONSTANT`; no knife is inert on its pin.

F2: the heavy-residual direct-payload carve-out is gone; the leftover scalar allowance cannot admit `R_N6`/`R_COV`/`SPLIT_CHECK`/guard tables; clause 1/4/8 + Deliverables close success and error transcripts.

F3: digest is harness-computed, well-defined on both engines, and is a same-support witness; no frozen knife is left with both fingerprint and digest blind.

Frozen knife design intact; no new cert-wrong ambiguity.
