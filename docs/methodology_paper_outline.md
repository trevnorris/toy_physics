# Methodology Paper — Outline and Theses

**Status:** Outline only. Deferred until at least one Stage-1 (branch realization) result is in hand — pass or trustworthy miss. The paper's claim sharpens qualitatively once there is a concrete falsification outcome to point at.

**Purpose of this doc:** Catalog the methodology elements that distinguish this project from "AI helped me build something cool" so they aren't forgotten. The user named the methodology piece as a separate planned artifact and wanted the structure written down now while it was fresh.

**Audience for the eventual paper:** Other solo researchers (programmers, independent scientists, hobbyists with formal background) who want to do credible speculative research without an academic appointment. Secondarily: methodologists studying AI-assisted research workflows.

**One-line thesis:** *AI-assisted research produces credible falsifiable results only when wrapped in an explicit discipline structure that pre-commits to falsification before exploration. The discipline is the thing; the AI is the multiplier.*

---

## 1. The problem this paper addresses

The current default for AI-assisted speculative research is: someone has an idea, asks an LLM to flesh it out, gets back impressive-looking derivations, publishes a preprint, claims a breakthrough. The output looks like physics, often *is* internally consistent, and has no scientific weight because none of the degrees of freedom were committed before seeing the output. The reader cannot tell whether the derivation is a discovery or a confabulation.

The trap is structural, not motivational. Honest researchers fall into it because:

- LLMs are happy to produce plausible derivations for almost any premise.
- Exploration without falsification structure produces "interesting numbers" that admit no clean pass/fail.
- Without a peer-review or PI-supervision check, the researcher's own motivated reasoning has no counter-force.
- The compute costs of speculative simulation are low enough that running until something looks right is tempting.

This paper documents one solo-programmer-plus-AI workflow that produces falsifiable, error-budgeted, pre-registered results despite operating outside the academic credentialing system. The workflow is described as a set of disciplines, not as a software tool. It is intended to be reproducible by another solo researcher with similar tools.

---

## 2. The core disciplines (one section per discipline in the paper)

Each of these is an independently defensible practice. Together they form the structure that lets a solo-programmer-plus-AI program produce results with scientific weight.

### 2.1 Pre-registration of test design

Before running any test, write down:

- which branch is being tested,
- which targets count as a pass and which as a fail,
- which observables will be extracted and how,
- which error budget the result must clear,
- what would constitute a trustworthy miss.

Write this as a brief that the researcher commits to before seeing data. After writing, do not edit it in response to results. The result has scientific weight in proportion to how strictly this rule is held.

Concrete artifact in this project: `docs/branch_realization_brief.md` and the work-package decomposition in `notes/5pn/5pn_stage354_355_computational_handoff.md`. Both written in the open, both unmodified after results begin to land.

### 2.2 Target-blind extraction

The numerical solver does not know the targets. The extraction pipeline does not know the targets. Targets are compared to extracted values only at a designated comparison step, *after* extraction is complete and frozen.

The solver's `freeze.target_blind = true` flag and `freeze.no_post_residual_refit = true` flag (see `research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md`) make this operational: the freeze hash covers the candidate parameters, solver tolerances, mesh, and source revision; any change after a residual is computed invalidates the packet.

Why this matters: without blindness, the unconscious tuning loop is undetectable from inside the researcher's own head. With it, "this matches the target" is a real signal.

### 2.3 Sequential audit gates with explicit human approval

Work proceeds in numbered batches. Each batch is reviewed (by a fresh agent, see §2.7) before the next batch can start. The human explicitly approves rolling forward; the system does not autonomously chain batches.

This is what catches errors before they propagate. In the present project, the red-team audit pipeline has run ~60 of 253 stages across multiple sequential batches, each gated by a user-approval point. Errors found at stage `N` can be fixed before they contaminate stage `N+1`'s derivation.

Memory anchors: [[feedback-sequential-audit-chunks]], [[feedback-review-agents]].

### 2.4 Validation suites before physics claims

No physics observable is reported until the solver has demonstrated trustworthiness on independent benchmarks:

- a known-analytic limit (free / linear regime),
- manufactured-solution tests per operator,
- a published benchmark with known answer (e.g. GPE vortex / soliton),
- mesh-refinement at 3+ levels showing convergence order,
- conservation diagnostics over the run,
- an explicit stated noise floor.

A result without this stack is not a result. The brief makes this non-negotiable; the methodology paper makes the case for *why* it is non-negotiable.

This is what distinguished the present plan from the prior `4d_1pn_sim` attempt: the prior project produced "weird numbers" that could not be certified as scientifically valid because the validation stack did not exist. The brief was structured exactly to prevent that recurrence.

### 2.5 Framing splits between exploration and claim

Private/exploratory conversation engages with speculative ideas directly, including "unification" language and conjectural ontology. Public-facing artifacts (papers, preprints, posted code) hold strict toy-analog framing only.

This is not duplicity; it is the standard separation between thinking aloud and asserting. The discipline is: the *artifact* is what's claimed, the *conversation* is what's explored, and the two have different epistemic statuses by design.

In this project: see [[feedback-framing-split]] and [[project-analog-framework-goal]] in auto-memory. The user's stated project goal is a working *analog* (mathematical bridge), not an ontological claim — the public framing reflects that even when private conversation goes further.

### 2.6 No mutation of tests after seeing results

A test is a question asked of nature (or, in a formal program, of the framework). Once results are in, the test is over. Mutating the test — adjusting the target, redefining the branch, adding a "we forgot this factor" correction, switching to a different observable — turns a falsification into an unbounded search that no negative result can end.

This is the rule most often broken in speculative work, and it is the single highest-value rule in the methodology. The brief states it explicitly (§3.4). The freeze-hash machinery in `NONLINEAR_PROTOCOL_V2.md` enforces it mechanically.

If the realized branch misses, the result is the miss. Find the next branch on its own pre-registered terms; do not let the current branch's miss send you searching for a pass.

### 2.7 Clean-agent reviews and context discipline

Each review or audit is performed by a freshly spawned agent with no prior context from the conversation that produced the work being reviewed. Reviewers get the artifacts plus the standards; they do not get the social pressure of having watched the work be built.

This is the AI-mediated version of the principle that you cannot review your own code. It works because LLM agents do not carry contamination across sessions; a fresh agent with the same instructions will catch errors a contaminated one would normalize.

Memory anchors: [[feedback-review-agents]], [[feedback-script-review-depth]], [[feedback-no-fake-scripts]].

### 2.8 Tool-of-record discipline

Different tools for different jobs, each authoritative for its own domain:

- **SymPy** for primary symbolic math.
- **Mathematica** for cross-check (single-seat — only one `math -script` invocation at a time; see [[feedback-mathematica-single-seat]]).
- **Markdown or YAML** for any file an LLM reads or writes (see [[feedback-no-json-for-llm-io]]).
- **JSON** only for pure machine-to-machine data.
- **Numerical solvers** (PETSc / Dedalus / etc.) for the final compute step, with explicit freeze hashes on configs.

The point is not these specific tools — it's that each artifact has a canonical tool, and answers across tools must cross-validate. This catches errors that any single tool would miss.

### 2.9 Human / AI division of labor

What the AI is for:

- drafting derivations, exploring framing, generating candidate approaches,
- mechanical verification (running scripts, checking algebra, cross-validating across tools),
- review and audit (as a fresh agent with no contamination),
- documentation and traceability.

What the human is for:

- scope decisions and pre-registration commitments,
- gate approvals between batches,
- framing judgment (private exploration vs. public claim),
- recognizing when a result is too good to be true,
- final scientific verdict.

Neither side alone produces credible results. The AI without the human produces plausible-but-uncommitted derivations. The human without the AI produces work too slowly to cover the surface area of a serious program. The combination, *with the discipline structure of §2.1–2.8 in place*, produces falsifiable artifacts.

---

## 3. Failure modes the discipline prevents

Each of these is a way speculative AI-assisted work commonly goes wrong. Each is addressed by specific elements of the discipline above. The paper can use these as motivating examples.

| Failure mode | Disciplines that prevent it |
|---|---|
| Confirmation-bias tuning ("just one more factor") | §2.1 pre-registration, §2.2 target-blind extraction, §2.6 no mutation |
| Drift from speculation to claim | §2.5 framing splits, §2.1 pre-registration |
| Noise-dominated "results" mistaken for signals | §2.4 validation suites, error budget statement |
| AI-generated plausible-but-wrong derivations | §2.3 audit gates, §2.7 clean-agent reviews, §2.8 tool cross-check |
| Context contamination across sessions | §2.7 clean-agent reviews |
| Unbounded "let's try another branch" searching | §2.6 no mutation, brief §3.4 |
| Loss of provenance / unreproducibility | §2.8 tool discipline, freeze hashes, written brief |

---

## 4. What this method costs

The paper should not undersell the cost. The discipline is expensive in three currencies:

1. **Time.** Months of patient gate-passing, where shortcut-tempting paths exist but are refused. The present project is at 1500 hours of work and Stage 1 hasn't run yet.
2. **Compute.** The §2.4 validation stack costs roughly 3–5× the headline production solve. A researcher used to "run until it looks right" will balk at spending the validation budget.
3. **Conceptual.** Refusing to claim ontology, refusing to mutate after seeing data, refusing to roll forward without audit — these are uncomfortable disciplines that feel like they slow the work down. They do; that is what they are for.

The methodology is not low-cost. It is low-cost *relative to the alternative of producing uncitable garbage*, which is what the unstructured version produces.

---

## 5. What this method enables

- **Pre-registered tests with publishable outcomes.** A pass is a real pass; a trustworthy miss is real science. Neither is "interesting numbers."
- **Reproducibility by another solo researcher.** The discipline structure is the artifact, not the specific physics. Another programmer with similar tools can apply the same workflow to a different speculative program.
- **AI tools used for credibility instead of productivity theater.** The current default for AI in research is "moved faster." The discipline structure shifts AI from accelerator to *witness* — it documents, audits, cross-checks, and reviews the work.
- **Falsifiable speculation outside the academic system.** The biggest claim of this paper. Speculative research has historically required either an institutional credential or a long apprenticeship. The discipline structure provides a substitute: institutional credentialing replaced by procedural rigor.

---

## 6. Concrete artifacts from this run (fill in after Stage 1)

Section to be written after Stage 1 lands. Will include:

- the PDE ledger (algebra-side framework),
- the red-team pipeline (252 stages of audits with verification gates),
- the Stage-1 validated solver and converged branch,
- the Stage-1 verdict (pass or trustworthy miss) with full error budget,
- the methodology disciplines as actually-executed (not as planned).

The paper's credibility depends on this section being concrete and citable. It is what distinguishes the methodology paper from a methodology-papers-everywhere preprint that nobody can ground.

---

## 7. Limitations and honest counterexamples

The paper should explicitly acknowledge:

- **The discipline does not solve the research-creativity problem.** It produces credible results from a candidate framework; it does not generate the framework. The framework here came from years of the user's prior conjectural work, not from the disciplines.
- **The discipline does not substitute for domain expertise.** It augments it. A researcher with no formal background in the physics being modeled would still need to learn enough to write a sensible brief and to recognize when a result is suspicious. The AI does not supply that judgment.
- **Not all speculative programs admit clean falsification structure.** Some questions are intrinsically open-ended or require empirical data the researcher cannot access. The methodology applies cleanly to formal programs (mathematics, theoretical physics, computational modeling); less cleanly to social-science or historical claims.
- **The time investment is substantial.** 1500 hours over multiple years to reach the gate of the *first* falsifiable test is a real number. The paper should report it honestly rather than implying the workflow is fast.
- **The methodology is conditional on tools that may not persist.** Specific LLM availability, specific compute access, specific scripting infrastructure — these change. The disciplines themselves are tool-agnostic; the specific instantiation is not.

---

## 8. Suggested structure for the eventual paper

Approximate section order (~30–40 pages):

1. **Introduction** — the trap of unstructured AI-assisted speculative research, and the case for procedural rigor as a substitute for institutional credentialing.
2. **Related work** — pre-registration in social sciences, replication crisis, AI-for-science literature, falsification methodology (Popper, Lakatos), the specific failure modes documented in retracted speculative-physics preprints.
3. **The discipline structure** — §2.1 through §2.9 above, one subsection each.
4. **Failure modes addressed** — §3 above, as motivating examples.
5. **Case study: the toy-superfluid analog program** — concrete artifacts §6 above, with the Stage-1 outcome reported either way.
6. **Costs and limitations** — §4 and §7 above, honestly.
7. **What this enables** — §5 above.
8. **Reproducibility appendix** — exact tool list, exact discipline checklist, exact filesystem / memory layout, exact failure-mode test cases. The thing another solo researcher would need to run the same workflow.

---

## 9. What this paper is not

- Not an AI-hype paper.
- Not a "we discovered new physics" paper.
- Not a "Claude is great" paper. The disciplines work with any sufficiently capable LLM; the paper should be careful to frame it that way.
- Not a methodology paper that exists in advance of a falsifiable result. Premature methodology papers in this space are abundant and underwhelming. The Stage-1 result is what makes this version different.

---

## 10. Sources to cite (preliminary)

- Brief and execution plan: `docs/branch_realization_brief.md`, `docs/branch_realization_execution_plan.md`
- Red-team discipline: `docs/redteam_thoroughness.md`, `docs/adversarial_audit_directive.md`
- Tool discipline: `feedback_*` memories under `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/`
- Failure-mode example: `/var/projects/4d_1pn_sim/docs/results_journal.md` (the unstructured-attempt case study; serves as the cautionary counter-example the methodology was developed to avoid)
- Pre-registration literature: existing methodology work in clinical trials and social-science replication initiatives
- AI-for-science context: relevant existing work on LLM-assisted theorem proving, automated mathematics, etc.
