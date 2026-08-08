# Independent review — prior-art findings document

## Artifact

`/var/projects/toy_physics/docs/method_prior_art_findings.md`

It is **orchestrator-written prose**: a research-findings document surveying prior art on how to verify
CAS-derived symbolic physics derivations. It was commissioned to answer whether the project's current
"dual-engine" architecture (two independently-built CAS scripts — one Wolfram, one SymPy — written from
one shared specification, run, and their emitted objects compared by a harness) is the right architecture.

It is **not** a physics derivation. Its claims are of three kinds, and each has a different failure mode:

1. **Repository measurements** — statements about what this codebase contains.
2. **Literature claims** — statements about what published work says, with numbers.
3. **Recommendations** — inferences drawn from 1 and 2.

## What to check

⭐ **The governing question: would acting on this document lead the project somewhere wrong?**

The document will be used to decide whether to change the verification method for ~20 remaining derivation
steps, and whether to redo two completed ones. A wrong number, a misquoted paper, or an inference that does
not follow is therefore expensive. Check, in this order of priority:

### ① The two repository measurements in §0 — ⭐ open the files, quote the lines

The document claims:
- `research/pde_ledger_v3/reduction/engine_output_checks.py:728` is `sp.simplify(lhs - rhs) == 0`, and
  that it is reached **only after** a naming-negotiation layer at `:799`, `:815`, `:867`, `:1157`.
- `research/pde_ledger_v3/reduction/` totals **8,081 lines**; `reduction/checks_S10.yaml` is **3,121**
  lines and its own header concedes it compares a subset.
- S11 engine 1 emits **3,750 tags** (sourced to `research/pde_ledger_v3/REBUILD_HANDOFF.md`).
- Of engine 1's nine defects, **none** was visible by reading and the ablating leg found **8 of 9**
  (same source).

⭐ **Verify each by reading the file and quoting the literal line.** ⛔ Do not accept a number because it is
plausible. ⚠ The stronger check is the **inference** drawn from them, not the digits: §0 concludes that the
comparator "compares NAMES, then SIMPLIFIES." Read enough of `engine_output_checks.py` to say whether that
characterisation is fair, unfair, or incomplete — in particular whether the naming layer is genuinely a
**precondition** of the equality test or merely an adjacent feature.

### ② The load-bearing literature numbers

These carry the document's central argument. For each: **can you verify it, and does the source say what
the document says it says?**

- DLMF/LaCASt (arXiv 2109.08899): **660 of 2,405** formulae verified symbolically = **27.4%**; **418**
  more verified numerically; the paper's own phrase *"symbolic simplification verification is not extremely
  effective"*; and ~**892 cases (51.1%)** above threshold, attributed largely to missing constraint
  semantics and branch cuts.
- N-Version Programming with Coding Agents (arXiv 2606.20158): **48** implementations, **1M** inputs,
  single-version **387.44** failures → three-version majority **130.99**; co-occurring failures traced to
  specification ambiguity.
- Brilliant, Knight & Leveson, *The Consistent Comparison Problem in N-Version Software* (1986 / IEEE TSE
  1989) — that it exists, and that it says an N-version system may fail to reach consensus **when no
  version has failed**.
- Knight & Leveson 1986 — 27 versions, 1M tests, coincident failures above the independence assumption.
- MAST multi-agent failure taxonomy: **Specification Issues 41.77%**.
- Agentic Diagrammatica (arXiv 2603.26990) — the quoted diagnosis about "implicit mathematical
  conventions", and that its remedy is a **trusted backend**, not N versions.
- de Bruijn factor ≈ 4.

⛔⛔ **If you have no web access, say so explicitly and confine yourself to the repository, to internal
consistency, and to what you can establish by reasoning. ⛔ DO NOT verify a citation from memory and do not
report a citation as confirmed on that basis.** A leg that says *"I recall this paper and it is right"* is
producing exactly the defect class this project exists to remove: an assertion with no computation or
source behind it. Saying "unverifiable without web access" is a **correct and useful** answer.

⭐ If you *do* have web access: fetch the source, quote the sentence, and report any place the document
overstates, understates, or attributes a claim to the wrong paper. ⚠ Several citations are already marked
in the document as unread (ⓢ) versus fetched (ⓕ) — check whether that labelling is honest, i.e. whether
anything marked ⓕ is in fact second-hand.

### ③ Does the recommendation follow?

The document's central proposal is to replace **name-based matching + symbolic simplification** with
**value-based comparison**: evaluate both engines' emitted objects at a shared set of **exact** points
(rationals or a finite field, never floating point), hash the evaluation vector into a fingerprint, and
**join on the fingerprint rather than on the tag name**.

⭐ Attack this. Specifically:
- Does it actually dissolve the naming problem, or does it relocate it into the symbol→value map?
- The document claims the project's objects (Lagrangians, system matrices, determinants, polynomial
  dispersion relations, roots) are the "good case" for Schwartz–Zippel because they are rational/algebraic.
  ⭐ **Check that against the actual artifacts** — e.g. `research/pde_ledger_v3/scripts/` and
  `research/pde_ledger_v3/mathematica/`. Are there objects that are *not* rational/algebraic, where
  probing would be unsound or would produce the DLMF false-positive failure mode?
- The proposed fix for "11 DISAGREE rows on nullspace-basis representation" is to compare **subspaces**
  (RREF, or kernel membership plus rank) rather than bases. Is that correct and complete? Is there a case
  where two engines' kernels legitimately differ and the check should fire?
- ⭐⭐ §4 contains a counterweight claiming value-based comparison **does not** fix the project's measured
  "typed conclusion with no computation behind it" defect, because a typed literal fingerprints just as
  well as a computed one. **Is that right?** And is the document consistent — does any earlier section
  oversell numeric probing in a way §4 then has to walk back?

### ④ Rule 16 — prior art is an oracle, never a premise

The project's governing rule: *check our computed result against prior art; never assume a published
result holds for our object; its conditions may not be ours.* ⭐ Find every place the document **adopts** a
method because it is standard, or transfers a published conclusion to this project's object without
establishing that the conditions match. §5 makes a closing disclaimer — check whether the body honours it.

### ⑤ Absence of evidence

§5 reports several NOT-FOUND results. The project rule is that **NOT-FOUND is never evidence of
originality**. ⭐ Check that framing. ⭐ **And do your own search: name anything relevant the document
missed.** Candidate areas it may have underweighted — decide for yourself whether they matter:
property-based testing, golden-master/characterisation testing, interval arithmetic and validated
numerics, the literature on CAS bug-finding, reproducibility infrastructure for computational mathematics,
and anything on comparing intermediate results rather than final ones.

## Do not read

⭐⭐ **Reading order is the blindness mechanism here, and it is load-bearing.** Form your own view of what
the sources establish **before** you read the artifact or the document that commissioned it.

- ⛔ **`/var/projects/toy_physics/docs/method_prior_art_brief.md`** — this is the commissioning directive.
  ⛔ Do not read it until you have completed ① ② and formed your own view. It carries the orchestrator's
  framing of the problem, including its own numbers, and reading it first anchors you to the framing that
  is under test. ⭐ You may read it **last**, solely to judge whether the findings document answered what
  it was asked; report that as a separate, clearly-labelled section.
- ⛔ `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/` — raw build transcripts.
- ⛔ Any file matching `*_review*.md` or `*_review*.txt` under `docs/` or
  `research/pde_ledger_v3/directives/` — prior review outputs, including the other leg's.
- ⛔ This prompt file is not evidence: ⛔ do not treat the claims restated above as confirmed. They are
  restated so you know what to check, ⛔ not so you can quote them back.

## Required method — DOCUMENT branch

⭐ Read the **source of truth first** — for this artifact that means the repository files named in ①, and
the cited papers themselves if you can reach them — form your own view of what they establish and what
they do not, and **only then** open
`/var/projects/toy_physics/docs/method_prior_art_findings.md`.
⛔ Do not read them in the other order.

⭐ **Quote both sides for every finding**: what the source says, and what the document says about it.

⛔⛔ **A PROSE RE-DERIVATION IS WORTH NOTHING.** For every claim you confirm or refute, the evidence is a
**quotation with a locus** — a file and line number, or a sentence from the fetched source with its URL.
⛔ *"I checked and it is right"* will be discarded. Where a claim is checkable by running something (a line
count, a grep, reading a specific line), ⭐ **run it and paste the literal output**.

## Physics filter

⚠ This artifact is a **method** document, so the usual physics filter is adapted, ⛔ not dropped:

> Report a finding only if it could **change what the project computes, what it may claim, or what method
> it adopts**. ⛔ Do not report style, formatting, tone, markdown conventions, or "this section could be
> clearer". ⛔ Do not report that a citation could be more complete unless the incompleteness changes what
> the document concludes.

⭐ **A finding that says "this recommendation is wrong / does not follow / would break X" is worth more
than ten citation-formatting notes.**

⚠ **A leg that returns "nothing survives the filter" is weak evidence on its own.** If that is genuinely
your conclusion, say so plainly — ⛔ but then state explicitly **what you checked, what you could not
check, and what would have had to be true for you to find something.** ⛔ Do not manufacture findings to
fill a quota; the honest answer, with its scope stated, is acceptable.

## Operational constraints

- ⛔ **Do not modify the repository.** This is a read-only review. If you want to run something, run it
  read-only or copy to `/tmp` first.
- ⛔ No Mathematica kernels are needed for this artifact. If you start one anyway, wrap it in
  `timeout 600` and never run more than one at a time (the licence has two seats).
- ⭐ Save any script you write and its literal stdout to named absolute paths, and report those paths.
- ⭐ Report findings ranked most-severe first, each with: the claim, the evidence (quotation + locus), and
  what would have to change.

## Quarantine rule

Not applicable — there is no parallel builder and no sibling implementation of this artifact. Blindness
here comes from **reading order** and from the do-not-read list above.
