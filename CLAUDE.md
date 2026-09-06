# How we work

Seventeen controls. Every one exists because ignoring it cost a session — none is removed here. This
**2026-09-05 restructure** groups them so the rule is scannable and the evidence that earned it is kept but
out of the way; it changes presentation, **never authority** (design + preservation map:
`CLAUDE_streamline_proposal_2026-09-05.md`; the M/E/G/S controls were preservation-checked against the
prior version by two independent legs before commit).
Full process: `docs/development_pipeline.md`. What we're building: `docs/development_plan.md`.
Where we are: `STATUS.md`.

⭐ **Correctness override (R11):** cost is **never** a reason to drop a control, narrow a check, or skip a
leg. Only the user scales work down. Nothing below — no grouping, no budget, no "administrative" label —
weakens a gate.

---

## At a glance

1. **Identify the object and the artifact's role.** Supply verified equations; withhold the target-answer
   criterion (the expected-value acceptance test, only); construct independently; preserve disagreement. (M1, M2)
2. **Measure every claim.** Compute and emit *before* guards; interpret in records; attach the command and
   its literal stdout in a `_measurements/` file — **this grounding duty binds the orchestrator too, CAS claims
   included**. But I personally author/run only **mechanical fact-lookups**; a CAS computation that adjudicates a
   finding is a **Codex-written script reviewed under G1**, ⛔ never one I author or privately clear. (E1)
3. **Pick two non-author legs.** Orchestrator-written → Codex + Grok; Codex-written → fresh Claude + Grok.
   Launch on sight, before inspecting the result. (G1, G3)
4. **Use the right exit.** *Decision list:* two legs, verify, fold **once**, go. *Physics spec / script /
   record:* **review-until-clear** — iterate leg→fold→leg until nothing outstanding changes what is
   computed or may be claimed. The filename confers no exception. (G2, G4)
5. **Verify and preserve.** Adjudicate findings by evidence (a finding is not a mandate); ablate physics
   checks by FORM; keep varying quantities live; change a defect-breeding author; no commit before both
   reports; commit the reviewed baseline before a repair overwrites it. (E2, G4, M3)
6. **Enforce, don't merely promise.** Bound handoffs and leave visible evidence; never revive quarantine or
   drop a control for cost; follow the annex/GIN policy for v3 `.out`. (S1, R11)

**Scope precedence:** classify an artifact by **function, then authorship** — the filename confers nothing.
A packet serving several roles inherits **all** applicable controls. The decision-list one-pass exception
(G2) **never** caps review of a physics spec, even one embedded in a decision list. "Review-until-clear" =
G4's substantive exit (nothing outstanding changes what is computed or may be claimed), never consensus or
green labels.

### Artifact → review discipline

Reviewer key: **O** = orchestrator-written → Codex + Grok; **Cx** = Codex-written → fresh Claude agent +
Grok. Mixed/changed authorship needs an explicit valid non-author pairing — never a silent reuse or a
reduced count.

| Artifact / role | Legs (by author) | Method | Pass / fold discipline | Advance / commit gate |
|---|---|---|---|---|
| **Pre-builder decision list** | Two; O (always orchestrator-written) | Check the requested decisions + their supporting evidence; ⛔ no fictional-script ablation | **One two-leg pass**, verify findings, fold **once** — ⛔ no iterate-to-green | No builder before both reports + fold. An unresolved issue that changes what is computed or may be claimed routes to the applicable spec/build gate (R10 scope); physics-bearing content also meets the spec row. No commit before both reports. |
| **Physics spec / shared spec / physics-bearing directive** | Two; O or Cx by actual author | Review the requested physics, complete premises, recipe/answer leakage; substantiate physics claims independently. Defer executable script-control tests to the build. | **Review-until-clear**: repair/re-review any change to computation or claims | Clear before using as governing physics. Directive review never pays the build tax. Both reports before any commit; reviewed baseline preserved before overwrite. |
| **Script / physics-bearing build** | Two; O or Cx; normally Cx for `/build` | Independent derivation scripts + literal stdout; every load-bearing check ablated; **mandatory FORM ablation**; one-sided corruption for independence; emit-before-guard + output checks | **Review-until-clear**; change author when repairs breed defects (⛔ never a fourth fold) | Launch legs before inspecting results; both usable reports; preserve reviewed baseline before repair; accept only on substantive clearance. Serialize dual Mathematica ablations. |
| **Step record / `.tex` card / physics-bearing prose** | Two; O or Cx; ⛔ never chosen by extension | Source-first fidelity; quote both sides; ⛔ no build directive in the packet; a measured physics claim still needs its command + stdout; for cards check suppressed macro fields | **Review-until-clear** about what may be claimed | Own artifact/version review required; script review alone is not record review. Both reports before commit. |
| **Other / claimed non-physics** | Two if physics-bearing (G1); if you claim it is not, record why | First record whether it changes computation, premises, checks, or claims — if so it is physics-bearing, route to a row above | ⛔ Do not infer a one-pass or zero-review exemption from the suffix or an "administrative" label | The two-report commit gate has no non-physics exception. |

**Observable gate record** *(new 2026-09-05, from the approved proposal; implements R12's observable-artifact
principle — illustrative, not a newly mandated gate)*: artifact path/version + role(s);
author(s); handed-input manifest; the identical rendered review prompt; two separately-attributed reports
with script/stdout paths where required; verified finding dispositions; iteration/author-change record; the
reviewed-baseline commit and the accepted version. It documents controls; it does not substitute for
carrying them out, and it is not an answer quarantine.

---

## M · Specify and construct

**M1 — Two engines exist so they can disagree, and the disagreement is the measurement.** *(was R1, R6)*
Independent construction, not hidden answers. A disagreement is a **finding** — ⛔ never try to make
divergence impossible with more careful prose; that defeats the reason there are two engines. ⛔ Never treat
a disagreement as a builder target to eliminate; it is a finding on the orchestrator's side.

**M2 — Name the object; supply verified equations; withhold only the target answer.** *(was R3, R5)*
Name the **object**, ⛔ not the recipe: if a review argues *how* to derive something — is this quotient
well-defined, is this weight unique — the question was manufactured by specifying a derivation path. Ask for
the object; let the engine hand over what it built. The spec says what to **compute** — ⛔ never what
anything equals, is expected, or was measured. Withhold exactly one thing: an acceptance criterion
referencing an expected value (a builder iterating to exit 0 converges on any target it can see). Supply
everything else, as **equations**, and label supplied objects as supplied and unfalsifiable in that build.
Naming an object does not waive first-principles derivation, form ablation, or evidence. Prior art's result
for our object is never a supplied premise (M3).

**M3 — Prior art is an oracle, not a premise; keep every varying quantity live.** *(was R16, R17)*
Prior art is an **oracle**, never a premise: check our computed result against it; ⛔ never assume its result
for our object — its conditions may not be ours. A freeze is a red flag; a **required freeze is the finding**
— treating a quantity that VARIES as if it were constant manufactures a wrong-but-consistent answer both
engines share, which the comparator misreads as agreement (the measured account: L-R17). Keep every varying
quantity **LIVE** and differentiate it; when a method seems to require holding one fixed to proceed, that
requirement is the measurement, not a step. Caught only by a ground-truth anchor (verify *its* conditions —
never replace the derivation with its answer) or a variable-coefficient/form ablation (E2). *(evidence: L-R17)*

---

## E · Produce evidence

**E1 — A script prints computed objects; a record interprets them; every claim carries its command.**
*(was R2, R4)*
A script **prints** computed objects — ⛔ it never states conclusions. Emit both operands and the residual,
**then** guard (a residual asserted zero always prints `0` and carries no information). Interpretation
belongs to the step record. ⛔⛔ **This binds the orchestrator too** — the half I keep exempting: every
review-leg prompt I write says *"a prose derivation is worth nothing; show the script and its literal
stdout, or the claim is discarded,"* and the same standard applies to anything I claim. **A claim about an
artifact carries the command that produced it**, and the commands + their literal output go in a
`_measurements/` file beside the document — ⛔ **this grounding duty is unconditional; a CAS-backed claim is NOT
exempt from it.** ⛔⛔ **But grounding is NOT a licence to AUTHOR the instrument.** The only thing I author and
run **myself** is a **mechanical fact-lookup**: existence, verbatim retrieval, literal-match counting, or the
cardinality/shape of an explicitly named stored object — ⛔ with **no** derived predicate and **no** algebraic,
numeric, symbolic, unit/dimensional, `subs`, simplify, solve, differentiate, limit, tolerance, or truncation
step. ⭐ **Anything else used to adjudicate a review question is an INSTRUMENT** — regardless of language, length,
filename, or library — because whoever writes it fixes the framing, method, and truncation and then judges its
output (the single-engine bias this architecture exists to remove). ⇒ ⛔⛔ **I NEVER author a CAS
analysis/diagnostic/verify script (⛔ not in `_measurements/`, ⛔ not under "rule 13"), and I NEVER adjudicate
from an instrument I ran privately.** I own the **question** — and run it by Codex first, to agree it is the
right one, at the **retained order**, ⛔ not a convenient proxy — and the **adjudication** (G4). The instrument is
**Codex-written and reviewed under G1** (Codex-written → fresh Claude + Grok; the review must confirm both that
the script answers the question AND that the question is asked at the retained order, not a convenient proxy)
before its output is trusted. If I am deciding in prose what the engines should compute, I have **inverted the
method** — fix the instrument by
**routing it to Codex**, ⛔ don't reason around it, and ⛔ don't **become** the instrument. *(evidence: L-R2,
L-R4; and L-CAS — the 2026-09-06 E over-clear, resolved with my own biased `verify_EG` on a σ_W→0 proxy, is the
measured relapse.)*

**E2 — Ablate to test; FORM tests physics, COEFFICIENT tests arithmetic; a check against a nonexistent
script goes to the build review.** *(was R14, R17)*
**Ablate** to test; ⛔ don't read. A **FORM** control tests physics (it leaves the family); a **COEFFICIENT**
control tests arithmetic (it does not). Demand a script and its literal stdout — a prose re-derivation is the
same defect relocated into the review. The variable-coefficient/form ablation is what catches the M3 freeze.
Controls written against a script that **does not exist** cannot be ablated, so they get reviewed by reading
(the weakest instrument): put them in the **build** review, ⛔ not in a document about a build. Computed means
*dependent on the derivation*, not merely CAS-shaped: a hand-typed algebraic payload, an answer-bearing tag,
or a suppressed-identical payload is not computed — the build-skill's action/ansatz-only construction,
object-only tag names, independent routes, and value-independent emission still bind. *(evidence: L-R14)*

---

## G · Review and advance

**G1 — Whatever writes does not review; pick the two legs by authorship.** *(was R7)*
**Whatever writes does not review.** Two independent legs on anything **physics-bearing** — and a spec both
engines read is physics-bearing, because an error there makes both engines agree on the same wrong thing.
**Orchestrator-written → Codex + Grok. Codex-written → fresh Claude agent + Grok.** Choose legs by
authorship, ⛔ never by file type; record mixed/changed authorship and establish a valid non-author pairing
before proceeding (⛔ never silently reuse a contributor or reduce the count). My own verification of findings
is required (G4) but is **never** a leg. *(evidence: L-R7)*

**G2 — The decision-list gate — and ONLY the decision list gets the one-pass exception.** *(was R7 trigger;
disambiguated 2026-09-05)*
**TRIGGER — no builder launches until its decision list has had two legs.** The decision list is
orchestrator-written and is the one artifact the *builder* trusts: everything downstream is checked twice,
the list is checked zero times. **One pass, then fold and go — never iterated to green.** ⛔⛔ **The one-pass
exception is the decision list's alone.** A physics **spec** (`SHARED_PHYSICS`, or anything both engines
read), a **script**, and a **physics-bearing record** are **reviewed until clear** under G4 — the filename
confers no exception, and a decision-list pass never clears a spec embedded in it. If a packet serves both
roles, it satisfies **both** gates. *(evidence: L-R7)*

**G3 — Launch legs on sight, before inspecting the result; no commit before both report.** *(was R8, R9)*
Launch legs **on sight**, before I look at the result — a self-check discharges the felt need for an
independent one, and it is most convincing when it finds things. Before launch, verify only that the
deliverable exists, is non-empty, and ran plausibly — ⛔ do not inspect its result content, and never feed
one leg the other's output. **No commit before both legs report.** The commit is the last step. Reviewing
the directive does not pay the tax for the build. *(evidence: L-R8)*

**G4 — Verify each finding, iterate to clearance, change the author if folds keep breeding defects.**
*(was R10, R13, R15)*
**Stop when nothing outstanding changes what is computed or what may be claimed** — ⛔ not when both legs are
green. "A leg that finds nothing is weak evidence" is my prior; ⛔ put it in a leg's prompt and it becomes a
quota (keep it in rationale, out of rendered instructions). A finding is not a mandate — **verify it
myself** (legs have been wrong in both directions); obtain **both** reports before adjudicating or editing.
⛔⛔ **"Verify it myself" = I OWN the judgment, ⛔ NOT that I author OR privately run the CAS instrument** — a
physics-CAS verification is Codex-written and reviewed under G1 (E1); I adjudicate from the reviewed script + its
filed stdout + the legs, and author/run only mechanical fact-lookups myself.
If successive revisions keep breeding defects in the material just changed, **change the author** — ⛔ don't
fold a fourth time (author change waives neither review nor clearance). Commit the exact **reviewed
baseline** before a repair overwrites it (recording unresolved findings); that preservation commit is not
acceptance — accept the repaired result only after its own review and G4 clearance. *(evidence: L-R10, L-R13)*

---

## S · Enforce through observable artifacts

**S1 — A prohibition is not a control; blindness is enforced by absence.** *(was R12)*
**A prohibition is not a control.** Blindness is enforced by **absence** — by bounding what the builder is
handed — ⛔ never by a sentence forbidding a read. A do-not-read list is a denylist, and a denylist means the
architecture is wrong. The measured failure is absence of computation, not anchoring; **quarantine is cut
and E1 replaced it** — ⛔ do not revive any cut quarantine mechanism under a new name (raw logs live outside
the tree for hygiene, not as a blindness claim). ⛔⛔ **This applies to these controls too:** every one is
prose I drift from under load, so the ones that hold are the ones that leave an artifact whose absence you
can **see** (the at-a-glance gate record is one illustration, not a new mandate). *(evidence: L-R12)*

Operational runbooks stay authoritative in the skills — `.claude/skills/build/SKILL.md` and
`.claude/skills/review-legs/SKILL.md`. A root summary cannot replace their identical-prompt,
first-principles-script/stdout, mandatory-form-ablation, source-staging + directive-exclusion,
authorship-pairing, fixed-baseline, kernel-serialization/budget/sandbox, leak-probe, orphan-memory
diagnosis, and real-deliverable-verification obligations; ⛔ do not infer a weaker reading from a summary.

---

## Repository infrastructure — the `.out` transcripts live in git-annex + GIN (set up 2026-09-01)

The v3 CAS audit transcripts (`research/pde_ledger_v3/scripts/out/*.out` and
`research/pde_ledger_v3/mathematica/out/*.out`, ~370 MB) are **git-annex content backed by GIN**
(`gin.g-node.org/trevnorris/toy_physics`, public), NOT plain git blobs — one exceeded GitHub's 100 MB/file cap.
The policy is the root `.gitattributes` (last-match-wins): everything is `annex.largefiles=nothing` (plain git)
**except** those two `out/*.out` globs, which are `anything` (annex). GitHub is `annex-ignore` by design (git +
tiny pointers only); GIN holds the bytes.

- **Generating/updating a v3 `.out`**: just `datalad save -m "…"` — the policy annexes `out/*.out` automatically
  and keeps everything else in git. Then publish with **both** `datalad push --to gin` (content → GIN) **and**
  `git push origin <branch>` (git + pointers → GitHub). ⛔ Never `git add -f` a big `.out`. ⛔ Never annex an
  `*_exports.py` — they are hash-chained plain-git inputs the next step imports.
- **After any fresh checkout/clone the `.out` are annex symlinks with NO content** until `datalad get <path>`.
  `grep`-by-line still resolves once content is present (the symlink is followed), so a script/directive/leg
  that reads a committed `.out` must `datalad get` it first — otherwise it reads a dangling link.
- The GIN token is stored in the datalad keyring under credential name **`gin`**; use `--credential gin` with any
  `create-sibling-gin`/publish. Full record, exact commands, and open follow-ups: the
  `project-datalad-gin-out-storage` auto-memory.

---

## Evidence ledger — why each control exists

The measured incident or core rationale behind each control, quoted from the original rule — or **excerpted**
where the rest of that rule's content now lives in its canonical M/E/G/S entry above (so a few entries are the
war-story extract, not the rule's full text). ⛔ Nothing invented. A control whose original text carried no
separate war-story has **no entry here** — its rule above is self-contained (R3, R5, R6, R9, R11, R15, R16
are such). Each control that has an entry references it below as `L-R#`.

- **L-R2 (E1) — the orchestrator's own typed conclusions.** Measured 2026-08-12: four export-chain designs
  and eight legs died on a question one `len(LEDGER)` answered; *"the imported action is in the LEDGER"*
  (there is no action row), *"no change to the publish guard"* (false under execution) and a class citation
  pointing at two sections that declare no classes each cost a full round. **All four were typed conclusions
  with no computation behind them — the exact defect this rebuild exists to remove, relocated into the
  orchestrator.**
- **L-R4 (E1) — the broken instrument.** The tell: many turns reasoning toward an answer a script would
  settle in one. The cause, measured: the instrument was broken, so no measurement was available.
- **L-CAS (E1) — the orchestrator AS the instrument.** Measured 2026-09-06: reading E1 ("this binds the
  orchestrator too"; commands go in `_measurements/`) + G4 ("verify it myself") as licence, I hand-authored the
  CAS scripts `verify_F.py`/`verify_EG.py`, ran them, and used them to overrule a leg's F/E/G findings — then
  **over-cleared E/N6** (rep-invariance holds only in the σ_W→0 projection, not the retained order I claimed).
  A Codex-sol pass caught the over-reach. The defect: **when I write the instrument I pick the framing and
  truncation and then judge my own output — one engine as builder and judge**, the bias the two-engine
  architecture removes. ⇒ a CAS computation that answers a review question is Codex-written and reviewed under
  G1 (fresh Claude + Grok) before its output is trusted; only mechanical fact-lookups are the orchestrator's to write.
- **L-R7 (G1, G2) — six spec defects.** Measured 2026-08-09: six spec defects, each costing a build round
  *plus* two legs, when two legs before the build would have caught them — three "level-above" misses, one
  exception named instead of the property (which bred a regression), measured counts stated four lines above
  "do not state the counts", and **an acceptance test that would have passed with the defect still in
  place.**
- **L-R8 (G3) — the self-check trap.** A self-check discharges the felt need for an independent one, and it
  is most convincing when it finds things.
- **L-R10 (G4) — the quota mechanism.** "A leg that finds nothing is weak evidence" is my prior; put it in a
  leg's prompt and it becomes a quota.
- **L-R12 (S1) — the failed quarantine.** The measured failure is absence of computation, not anchoring;
  quarantine is cut and rule 2 replaced it. This applies to these rules as well: every one is prose I drift
  from under load, so the ones that hold are the ones that leave an artifact whose absence you can see.
- **L-R13 (G4) — legs wrong in both directions.** Legs have been wrong in both directions.
- **L-R14 (E2) — the nonexistent-script problem.** Controls written against a script that does not exist
  cannot be ablated, so they get reviewed by reading, which is the weakest instrument we have.
- **L-R17 (M3) — the frozen mechanisms.** The one move that has cost this program a round, nearly every time,
  is treating a quantity that VARIES — a background field, a field's argument dependence, a jet/derivative
  order, a rate, a shape — as if it were constant. It manufactures a wrong-but-consistent answer both engines
  can share, so the comparator reads agreement and the defect hides (26=26 was a coincidence of two frozen
  mechanisms). Caught only by a ground-truth anchor or a variable-coefficient/form ablation.
