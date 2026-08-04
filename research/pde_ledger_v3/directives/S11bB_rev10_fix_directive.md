# S11b-B — WRITE REVISION 10 OF THE BUILD DIRECTIVE

⛔⛔ **DOCUMENT TASK. ⛔ Do NOT write, run, or modify any `.wl` or `.py` script.**

You authored revisions 7, 8 and 9. ⭐ **You are the author again for rev 10.**

## Artifacts

**Edit:** `directives/S11bB_SHARED_PHYSICS.md`, and `S11bB_py_header.md` (one asymmetry only — see finding
Q5-sym). Reassemble both directives, verify byte-identity of the shared block, report the new `sha256`.

## The review to fold

`/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11bB_dir_rev9_agent.txt`

⭐ Fix every finding that survives your judgement; say so if one does not. ⛔ Derive remedies — do not
transcribe them. A second reviewer found nothing; ⚠ that leg has now cleared three revisions that contained
real defects, so ⛔ **do not treat its silence as corroboration.**

## ⛔⛔⛔ CONTAMINATION — you have now leaked twice, and so have I

The review's **Q1** lists three residues (**L1, L2, L3**). ⚠⚠ **All three entered through well-intentioned
edits, and two are of one specific shape:**

⭐⭐⭐ **THE RULE THAT WOULD HAVE PREVENTED THEM: ⛔ NEVER JUSTIFY A PARAMETERIZATION, A SPECIAL CASE, OR A
CONVENTION BY REFERENCE TO A PREVIOUS REVISION, OR BY NAMING WHICH SPECIAL CASE IS THE "BASELINE".**
The engines cannot read prior revisions — §0b bars them — so such a statement is simultaneously
**uncheckable** and an **outcome hint**. ⚠ `"τ_A = τ_V = τ_X recovers revision 8 exactly"` is exactly this
and must go. ⛔ Do not replace it with another sentence that privileges the equal-time case; ⭐ if the
parameters are free, **say nothing further about which values are special.**

⛔ Also remove the answer-shape presupposition in **L3** — asking *which coefficient vanishes* presumes the
condition is of vanishing type. ⭐ Ask for the necessary-and-sufficient condition **in whatever form it
takes**. And fix **L2**: §1b asserts the roots are complex while §0 and B5 ask whether they are real,
propagate, decay, grow, or fail to exist. ⭐ §1b must be conditional throughout — ⛔ do not touch the branch
prescription itself, only the framing wrapped around it.

⚠ **This will be checked against your diff, and the last two declarations of compliance were both wrong.**

## The remaining findings

- **Q2 — the two-port check has three readings and only one discriminates.** At real `ω` the literal
  time-average of the slab's energy rate vanishes identically, so a literal reading **false-FAILs a correct
  derivation** and an on-shell reading **passes vacuously**. ⭐⭐ **State unambiguously which computation is
  required**, so that neither degenerate reading is available. ⚠ Also record its verified blind spot: it is
  insensitive to the face response, the closure and the affinity normalization, because those cancel
  between the two sides. ⭐ State that limit explicitly rather than leaving the check looking stronger than
  it is.
- **Q3 — the general unequal-relaxation-time face response has no external standard.** The one supplied
  acceptance value lives in the equal-time slice. ⭐ Either give the general response a standard that does
  not privilege equality, or **state plainly that it is unchecked** and why. ⛔ Do not paper over it.
- **Q4 — two cautions.** (i) The `+i/τ_I` marker is exact only for two of the three channels; the third's
  pole is displaced by the interface coupling, so a literal test can both miss a real defect and alarm on a
  correct derivation. ⭐ Restate the diagnostic in terms of the property that is actually invariant — the
  reviewer verified it is mechanizable by **pole location relative to the real axis**. (ii) An engine that
  rationalizes a complex denominator by its conjugate puts a conjugated kernel *literally* into an in-scope
  object. ⭐ Make the check about the object's **analytic structure**, ⛔ not about symbol appearance.
- **Q5 — the read-bar is still incomplete, and the largest hole is structural.**
  ⛔⛔ **The git object store is not barred**: every prior revision, every `_scratch/` file and every step
  record is recoverable through `git log -p` / `git show`, and commit subject lines alone are
  answer-adjacent. ⚠⚠ **It reaches both engines identically, so it produces CORRELATED contamination that
  dual-engine cannot detect by construction.** ⭐ Bar it explicitly and comprehensively — history, object
  store, and any command that reads them.
  ⛔ Also bar `research/pde_ledger_v3/paper/` (a registry relation is stated in its appendices, precisely
  what `reduction/` is barred to withhold, and the predecessor step's card is there), and the v3-root
  process documents and `_review_prompt*.md` files the review names.
  **Q5-sym:** the **py header alone** tells its engine that the registry contains a relation this step is
  asked to reinterpret. ⭐ Remove that hint — the other engine never sees it. ⛔ Change nothing else there.
- **Finding 1 — the `Z_perm` check cannot make one of its advertised catches**, because the comparison
  happens after the equal-time specialization. ⭐ Narrow the claim to what it actually catches. ⛔ A check
  credited with a catch it cannot make is worse than one that claims less.

## ⛔⛔ CONSTRAINTS

1. ⛔ **MINIMAL AND SURGICAL.**
2. ⛔ **NEVER SUPPLY AN EXPECTED RESULT, ITS REASON, OR ITS SHAPE.**
3. ⛔⛔ **KEEP THE FALSIFICATION CHANNEL OPEN.** A growing root remains first-class; any diagnostic that
   could reject one must be explicitly scope-bounded.
4. ⛔ **EVERY CHECK MUST BE ABLE TO FAIL**, with a one-line statement of what wrong derivation it catches —
   ⚠ **and that statement must be TRUE**, see Finding 1.
5. ⛔ **DO NOT RE-OPEN** what independent legs cleared: scope boundary · `B1`'s mass balance · the `A/B/C`
   split · §1b's **branch prescription** including its upper-half-plane extension · §3b's
   virtual-displacement rule · the **supplied affinity `𝒜`** · `B8` controls B/C/D · the closure count.
6. ⛔ **TWO ENGINES MUST NOT BE ABLE TO DIVERGE.**
7. ⛔ **No new scope.**

## Output

The edited files, plus a report **under 60 lines**: each finding, what you changed, the `sha256`, and for
every check added/kept/modified a **true** one-line statement of what wrong derivation it catches.
⭐ State explicitly where you were tempted to reference a previous revision or name a baseline case, and
confirm you did not.

⛔ Then stop and exit.
