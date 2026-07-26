# Integration-test thread — HANDOFF

⛔⛔ **SUPERSEDED IN PART, 2026-07-26. READ `manifests/PIVOT_TO_REWRITE.md` FIRST.**
The survey-first approach this document describes was **abandoned** after two measurements showed
(a) `notes/parameter_register.md` covers only 46.5 % of the scripts' dimension-bearing quantities and
cannot seed a module, and (b) there is no single idiom to consolidate — 6 basis conventions over 4 axis
sets. **We now REWRITE each script's dimension handling onto a designed module** rather than surveying
43 scripts to discover idioms to accommodate.

**Still valid below:** the process lessons (§8), the checker asymmetry (§7), the paused threads (§9),
and the decisions in §2 that concern the schema/validator (which are committed and parked, not deleted).
**No longer the plan:** §1's "run the survey", §3's next steps, §5's pilot numbers.
Repair round 8 (`35dc6aa0`) IS committed and did fix real defects — `Dim()` was unrecordable in 13
scripts, `basis_locus` applicability contradicted itself, and two claims about Python scoping were false.

---

⭐ **READ ORDER:** this file → `manifests/LEDGER_WIDE_PLAN.md` (PLAN OF RECORD) →
`docs/development_pipeline.md` (process).
📋 **`manifests/RESTART_PROMPT.md`** — paste-ready session-restart prompt with the ordered next steps.
⛔ `manifests/FANOUT_PLAN.md` is SUPERSEDED — the parallel fanout is CANCELLED.
Do NOT read vision docs; they re-confuse.

**Operating model: Claude is a THIN CONDUCTOR.** Agents return ≤12-line checkable verdicts; full
reports go to disk and are read only on a negative verdict; the orchestrator independently re-runs
acceptance commands.
⚠ **Violated early this session** — the orchestrator read three full review files and hand-folded
them, and *every one of those folds introduced a defect*. The fix is
`_scratch/pass1_dim_survey/FOLD_PROTOCOL.md`: a fork agent reads the review, verifies each finding
against source, folds, self-checks, returns a receipt. **Use it.** It has since caught defects the
orchestrator shipped, defects in its own drafts, and a circular justification between two directive
sections.

Branch `ledger-v2-rebuild`. `schemas/` is UNTRACKED; only `HANDOFF_NEXT_SESSION.md` is modified.

---
## 1. WHERE WE ARE

**Pass 1 §2.1 (the dimension survey) has NOT RUN.** This session built the apparatus to run it
safely, and hardened it through **eight adversarial verification passes (~1,000 break attempts)**.

**Built and verified:**
- ✅ **SymPy refactor baseline** — 43 scripts, all exit 0, **reproducible** (two runs identical).
  `_scratch/pass1_dim_survey/baseline/` + `BASELINE_NOTE.md`.
- ✅ **Mathematica `.out` determinism** — all 44 `.wl` re-run: 44/44 exit 0, **43/44 byte-identical,
  44/44 after normalising Mathematica's `$NNNNN` symbol counter**, zero non-counter differences, 442
  script-seconds total. `out_probe/PROBE_RESULT.md`.
- ✅ **Survey directive at r17** — 5 Codex review rounds (9→6→5→4→1), a Grok gate, and 9 verification
  folds. `_scratch/pass1_dim_survey/directive/SURVEY_DIRECTIVE.md`.
- ✅ **Two schemas + validator + fixtures** (`schemas/`, untracked) — arbiter re-run after round 7:
  **exit 0, 122 PASS / 0 FAIL**. Independently verified: **43/43 rules load-bearing under per-rule
  ablation**, **105/105 rejects minimal violations**, 17/17 accepts, 82/91 call sites individually
  protected, 91/124 exercised, 116/116 real register rows recordable.
- ✅ **Pilot directive** — `_scratch/pass1_dim_survey/directive/PILOT_DIRECTIVE.md`.
- ✅ **Rounds 5 (deletion) → 5b → 6 → 7 executed and independently verified.**
  ⚠ **Round 5 printed its done-marker with its own suite RED** — the orchestrator's arbiter re-run
  caught `exit=1, 45 FAIL`. **This is why the arbiter re-run is not ceremonial.** Round 7's verifier
  returned `ACCEPT_WITH_FINDINGS` with 0 invented semantics and 0 undisclosed removals — the first
  round clean on both.

---
## 2. DECISIONS MADE (do not re-open)

1. **Schema-first split (USER).** Directive owns **semantics**; committed schemas + validator own
   **shape**. Prose cannot pin a structured record for 44 parallel agents.
2. **TWO schema artifacts** — a record schema and a verification-output schema.
3. **`ownership` and `uses` are orthogonal** (§3.11); `declaration_required` derives from `ownership`
   ALONE. An earlier single `role` enum made `ASSERTED_TARGET` imply a declaration, which would have
   minted 8 spurious canonical declarations from stage038's `EXPECTED_UNIT_STATE` alone.
4. **§3.11 narrowed to dimension-VALUED bindings**; producers (`__mul__`, `dim_add`, `dim_of`…) belong
   to §3.10 contract features. Makes the recovery↔bindings correspondence a *stated* 1:1.
5. **`orders[].axes` carries INTERPRETED names**, not source spellings.
6. **Anchor applicability is CONDITIONAL** — ~40 `PRESENT` stages have no §4.4 anchor and declare
   `NONE` positively. CORPUS anchors and OWNERSHIP anchors are disjoint sets.
7. **§4.8 case (ii) gets FULL ownership coverage, not a sample** — this also defeats *shrinking* the
   registry, not just zeroing it.
8. **Grok is the DEFAULT tertiary; GLM only for a further opinion** (USER).
9. ⭐ **The shared dimension module is SymPy-ONLY (USER).** Mathematica stays an **independent
   verification engine**. Consequence: `.py`-vs-`.wl` dimension comparison is a **permanent standing
   cross-check**, not a one-time pre-unification sweep, and a genuine difference is a dual-engine
   finding the unit audits missed.
10. ⭐ **The parameter register is an evidence source (USER caught the omission).**
    `notes/parameter_register.md` is the ledger's running cross-stage dimensional ledger, per-stage,
    keyed by quantity identity, dual-engine-verified. It is DERIVED — scripts/`.wl` are SOURCE; on
    disagreement the source wins and the disagreement is a FINDING.
11. ⭐ **Semantic dimension adjudication is DESCOPED from the survey** (orchestrator, after two
    REJECTs). Both channels now record verbatim evidence + a **string-decidable** status; §4.9 owns
    monomial normalisation and adjudication later. See §4 below.
12. **Pin `model: opus` on every agent dispatch** (USER). `fork` inherits by contract;
    `general-purpose` does not.

---
## 3. IMMEDIATE NEXT STEPS

0. ✅ **ROUND 7 IS DONE AND VERIFIED — `ACCEPT_WITH_FINDINGS`, pilot READY.**
   `_scratch/pass1_dim_survey/reviews/REPAIR7_VERIFY.md`. Suite green (122 PASS / 0 FAIL, exit 0).
   **Lexical coincidence is dead:** renaming every axis to alien tokens with evidence untouched left
   17/17 honest records accepting — under r16 that flipped to REJECT. 10/10 previously-false-rejected
   lines accept, plus 14 more real lines the verifier drew itself. 43/43 rules load-bearing,
   105/105 minimal violations, 17/17 accepts, **0 invented semantics, 0 undisclosed removals**.

   ⚠ **TWO ITEMS MUST BE SETTLED BEFORE THE 44-WAY DISPATCH** (neither blocks the 4-script pilot —
   021/038/042 use tuple literals, 024 is `REAL+ABSENT`). Both are PRE-EXISTING (r16 rejected them
   too), not round-7 regressions:
   - ⭐ **`DIMENSIONLESS = Dim()` has NO honest escape** — 13 lines including `stage004:143`, a line
     §3.5 itself names. Zero components vs a declared 3-axis basis fails the arity test with no
     honest way out. **This is a fabrication-forcing rule and must be fixed in the DIRECTIVE**, not
     patched in the schema.
   - **`LENGTH = Dim(1, 0, 0)` rejects across 72 lines in 13 of 43 scripts.** An honest narrow-span
     escape (`(1, 0, 0)`) exists but is UNDOCUMENTED — document it or widen the rule.
   - Minor: the schema requires `basis_locus` on `NAMED_AXIS` while r17 §3.0 says it arises only on
     `POSITIONAL`; §3.5 reads the other way. The builder resolved the conflict silently instead of
     reporting it. Non-fabricating, but settle the §3.0-vs-§3.5 wording.

   The verifier's recommendation, adopted: **run the pilot and fold these two questions into the same
   pass** — real records will show whether the escapes are adequate.

1. **(historical, now closed) round-5 verdict** — `_scratch/pass1_dim_survey/reviews/REPAIR5_VERIFY.md`
   (dispatched at end of session; may have completed after the handoff was written). The decisive
   question it answers is **whether the suite was greened HONESTLY** — no relaxed error-code sets, no
   deleted examples, no weakened over-firing rule, no less-specific assertions — and whether the
   deletion round actually deleted rather than renaming the verdict machinery.
   **If it returns REJECT, do NOT run the pilot; fix and re-verify first.**
   ⭐ Always re-run the arbiter command yourself; do not trust a build report:
   `python3 research/pde_ledger_v2/schemas/validate_dimension_survey.py --test-examples`
2. **Then the 4-script pilot** — `PILOT_DIRECTIVE.md`: stage021, stage038, stage042, stage024.
   ⭐ Includes the **planted-defect positive control**: mislabel one ownership value in a *copy* of an
   accepted record and confirm a fresh verifier flags it. **If it does not fire, the verification leg
   is decorative and the 44-way dispatch must NOT proceed.**
4. **Capture the three decision numbers** (§5) → bring to the user with GO/NO-GO.
5. Then the remaining 40, module design (extremes first: stage042 fractional non-LMT, stage038 four
   axes), refactor, digest cascade.

---
## 4. ⭐ WHY ADJUDICATION WAS DESCOPED (read before re-adding it)

Four schema rounds tried to make an asserted `AGREE`/`DISAGREE` dimension verdict trustworthy and
failed the same way each time: **the verdict was always computed from fields the artifact itself
supplies.** Round 4's `MATCH` compared two verifier-authored fields; `named_axis_derivation: "x"`
passed; false `AGREE` was wholly unguarded; a false `DISAGREE` landed between `L⁻¹M` and `M L⁻¹` —
the same dimension the register spells two ways.

**Root cause: a document-shape validator was being asked to do symbolic mathematics.** Monomial
normalisation (superscripts, axis order `[L,T,M]` register vs `(L,M,T)` scripts, fractional exponents,
29 multi-parameter/`[X]=`-prefixed cells) is real work and belongs in a computed pass.

**What the survey does now:** records both dimensions verbatim with loci + a status decidable by
string comparison alone, on BOTH channels (register and `.wl`), under ONE pinned whitespace
convention. **§4.9** must later normalise monomials, bind each pair to its quantity via
`quantity_ref`, and partition every text-difference into same-dimension-different-spelling /
documented-non-transferable-identity / **genuine drift** — only the last being a finding.
**Tracked as task #13. Deferred, not dropped.**

---
## 5. NUMBERS THE PILOT MUST PRODUCE (they decide whether Pass 1 stalls)
1. `UNDETERMINED` leaf count per record.
2. ⭐ **`CONSTRAINED_BUT_NOT_STATED` count among `REAL+ABSENT` quantities** — the **Pass-1 stall
   signal**. If high, the notes/register do not carry the dimensions, §4.7 `REFACTOR READINESS`
   blocks, and the module needs a derivation phase first. *(Early evidence is encouraging: stage024's
   `I25` IS registered at `L^(5/2)`, dual-engine-verified.)*
3. Wall-clock + token cost per record, survey and verification legs separately.

---
## 6. OPEN FINDINGS — catalogued, NOT fixed (standing user instruction)
- **`RANGE_ENDPOINT_DRIFT`** — census `(39,49,10)` vs stage043's `(40,49,9)`; stage032 declares `Q_E`
  `{low:0,high:1}` while stage043 counts it at both endpoints. stage043 probably right; **NOT
  adjudicated.** The schema has no `reclassified`/`reconciled` lifecycle action.
- **The refactor gate's MEASURE is undefined** — §2.3 says "identical PASS count" but **16 of 43
  scripts emit no `PASS tally:` line** (001–003, 016–028). Options in `baseline/BASELINE_NOTE.md`;
  strongest is the full ordered list of emitted `PASS` markers, since a bare count is satisfiable by a
  script that stopped checking. **Pin it in the refactor directive.**
- **Verification honesty is unprovable from artifacts** — nothing distinguishes a verifier that
  re-derived every binding from one that asserted it; ~40 stages have no known-answer anchor. The
  pilot's planted-defect control is the only mitigation.
- **The §3.3↔§3.11 hard 1:1 is verified only against a sampled corpus** (~6 of 43).
- **stage044** is the only script with a filesystem side effect (`_scratch/stage044/verdict_py.json`,
  `:1343-1346`).

---
## 7. ⚠ THE ASYMMETRY (raised by the user, unresolved)
**The feeder is now specified far more rigorously than the detector.** `composite_build.py` — which
actually finds cross-stage desync — still has **C7 edges 1/18**, ~42 unfixtured issue codes, and the
`infer_closed_slice` replacement unbuilt. Everything hardened this session feeds a checker that has
not had the same treatment. **Give Pass 1b (task #7) the same adversarial break-attempt treatment
BEFORE 44 manifests are written against it.** A "44 stages in sync" claim rests on the checker.

---
## 8. PROCESS LESSONS BANKED (memories written this session)
- `feedback-no-fabrication-forcing-rules` — ⭐ a rule whose only honest outcome is an INVENTED value
  is worse than no rule; ask "can this be satisfied honestly in every legitimate case?" not only "can
  it fail?". Watch **vacuously-true conditions on mandated-empty collections**.
- `feedback-parsers-need-real-input` — ⭐ fixtures synthetic ≠ INPUT synthetic. An evidence file a
  parser READS must be exercised against the real committed file, and against rows the builder did
  not pick. This produced the round-3 REJECT.
- `feedback-gate-before-launch-argv-snapshot` — Codex snapshots its prompt into argv; editing the file
  after launch does nothing. Gate must PRECEDE launch. `pgrep -f "codex exec"` self-matches your shell.
- **Grok's empirical gate earns its keep** — asked to *construct a passing-but-useless survey*, it
  succeeded; an all-`CHECK_LOCAL` record satisfied every criterion and would have handed the module
  designer an empty registry. That produced §4.8.
- **Disclosure of removals must be demanded explicitly** — an undisclosed `minItems: 1` removal cost a
  guard. When demanded, the next round disclosed correctly and the removal was judged sound.
- **Gate comparisons on a completion log, not on "the file exists"** — comparing mid-write produced a
  false provenance alarm.
- ⭐ **THE ARBITER RE-RUN IS NOT CEREMONIAL.** Round 5 printed `___CODEX_BUILD_DONE___` with 45 of 110
  fixtures failing. Reading its report would have moved the project to the pilot on a broken
  validator, where four records would have failed for reasons unrelated to the survey. **Re-run the
  published acceptance command yourself, every time, and read its exit code — not the build report.**
- **When ordering a red suite greened, name the forbidden fixes.** Relaxing an expected error set,
  deleting the failing example, or weakening the over-firing rule all produce `exit=0` and destroy the
  property that makes fixtures worth anything (every reject fails for its OWN named rule, and is a
  minimal violation). Forbid each by name.
- **A deletion round must be measured as one** — demand rule count and validator size before/after. A
  deletion round that grows has usually re-implemented the removed machinery under a new name.

---
## 9. ⏸ PAUSED (do not lose)
- **044-v2** — redo stage 044 with a dynamical-Σ sleeve (un-freeze `S_hold`);
  `notes/stage044_v2_unfreeze_prep.md`. ⚠ its v1 reference manifest was destroyed by an agent in an
  earlier session (gitignored, unrecoverable) — re-extract from sources.
- **045** — VII-2b non-variational drain/return block; `notes/stage045_nonvariational_block_prep.md`.
- **Long-term (user's goal):** generalize this apparatus into a portable pipeline and test it on
  published arXiv papers. Generalize AFTER the ledger is done. Validate first against papers with
  KNOWN errata as positive controls — a pipeline that cannot rediscover a known error cannot be
  trusted to report a novel one.
