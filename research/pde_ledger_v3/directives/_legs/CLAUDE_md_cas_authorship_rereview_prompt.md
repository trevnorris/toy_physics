## ⚠ Operational constraint (resource, this run): STATIC inspection only

This is a document review of governing prose. ⛔ Do NOT execute any `verify_*.py` or anything importing
`S11c_c2_stdout_loader` (they load a ~499 MB `.out`; a prior run was killed doing that). Read files only.

---

# Independent RE-REVIEW (round 2) — CLAUDE.md CAS-instrument authorship clarification

You are ONE of two independent legs. You will not see the other leg's output. This is **round 2**: round 1
(Codex-sol + Grok) blocked commit with the defects listed below, and the orchestrator has folded corrections.
Your job: verify each round-1 defect is now correctly fixed, that the fix introduced **no new** defect or
control-weakening, and that the file is internally consistent. A real defect is a finding, not a failure; report
only what catches a real way the current wording is wrong. ⛔ Do not rubber-stamp.

## Artifact under review
The **current** working-tree diff of `CLAUDE.md` vs `HEAD`:
`git -C /var/projects/toy_physics --no-pager diff CLAUDE.md`
Read the full current `CLAUDE.md` (working tree) and `git -C /var/projects/toy_physics show HEAD:CLAUDE.md` for
preservation. Touched: **E1** (§ "E · Produce evidence"), **G4** (§ "G · Review and advance"), **at-a-glance**
item 2, and the **L-CAS** ledger entry.

## The round-1 defects that were folded (verify each is now fixed)
1. **Grounding went conditional.** Round-1 wording ("this binds the orchestrator too — but for mechanical
   fact-lookups only") could be read as: the command+stdout grounding duty applies only to mechanical claims, so
   a "CAS instrument" claim need not be grounded. FIX INTENT: grounding is **unconditional** (a CAS-backed claim
   still carries its command + literal stdout in `_measurements/`); the routing rule is separate. **Check the
   at-a-glance item 2 and E1 now say grounding is unconditional and do NOT give "it's a CAS instrument" as an
   excuse to skip the measurement.**
2. **Subjective boundary.** Round-1 "mechanical / carrying no physics framing" asked the orchestrator to judge
   its own framing — the exact failure mode; a one-line `subs`/`simplify`/dimensional check would be self-classed
   as a "lookup." FIX INTENT: an **operation-based** definition — a mechanical fact-lookup is only
   existence/verbatim-retrieval/literal-match-count/cardinality-shape of a named stored object, with NO
   algebraic/numeric/symbolic/unit/`subs`/simplify/solve/diff/limit/tolerance/truncation step; anything else is an
   instrument regardless of language, length, or filename. **Check the current E1 boundary is objective enough
   that a one-line `subs`/`simplify`/dimensional/numerical diagnostic is unambiguously an instrument.**
3. **"Grok + I review" contradicted G1.** Round-1 said the orchestrator co-reviews the instrument, but G1 says
   Codex-written → fresh-Claude + Grok, and "my verification is never a leg." Also the author-vs-run conflation:
   round-1 said "I may run only lookups," which (with "Codex writes it") left no one cleanly executing the
   reviewed instrument. FIX INTENT (the orchestrator adjudicated the Grok/Codex split): the instrument is
   **Codex-written and follows G1** (fresh-Claude + Grok) — the orchestrator is **not** a review leg; the
   orchestrator **neither authors nor adjudicates from a privately-run** instrument; who executes is handled
   under the normal G1/build path. **Check E1 + G4 now route the instrument to G1 (not an orchestrator+Grok
   pairing), do not make the orchestrator a leg, and are coherent about author-vs-run (no step where the
   instrument is written but no one may run it, and no step where its author privately judges its own output).**
4. **Proxy-question root.** The actual over-clear was the *question* (`subs(σ_W→0)==0 ⇒ "rep-invariance holds"`
   was a proxy, not the retained-order test). **Check the current wording makes the review confirm the question
   is asked at the retained order and is not a convenient proxy — not merely that the script answers whatever
   question it was given.**

## Additional checks
5. **No new defect / no control weakened.** Every pre-existing E1 obligation intact (print-not-assert;
   emit-then-guard; interpretation in the record; claim carries command+stdout); R11 not loosened. No place a
   review **leg** or **builder** (who MUST author scripts) is swept up by the orchestrator-only prohibition; G1's
   cell allowing orchestrator-written build scripts (reviewed Codex+Grok) not deleted. Internal consistency
   across E1 / G4 / at-a-glance / L-CAS / the two skills (`build`, `review-legs`).
6. **War-story (L-CAS) still accurate** against the real artifacts (static read of the two verify-script
   headers + the physics adjudication) — orchestrator-authored scripts; E over-clear / σ_W→0-projection
   correction recorded. ⛔ Flag any overstatement.

## Output
For each of 1–6: finding + quoted `CLAUDE.md` line(s). Separate CONFIRMED defects from unsettled. End with: is
the change NOW correct to commit, or the exact remaining wording that must change.
