# Independent review — a governing-prose clarification to CLAUDE.md (orchestrator authorship of CAS instruments)

You are ONE of two independent legs reviewing a change to the project's **governing process file**, `CLAUDE.md`.
You will not see the other leg's output. This is orchestrator-written governing prose; an error here — a weakened
control, a new loophole, or an ambiguity — propagates to all future work, so review it adversarially. A real
defect is a finding, not a failure; there is no quota, report only what catches a real way the change is wrong.

## What changed and why
The orchestrator repeatedly mis-read two existing controls — **E1** ("this binds the orchestrator too"; "a claim
carries the command"; commands + stdout go in a `_measurements/` file) and **G4/rule-13** ("verify it myself") —
as *licence to hand-author its own CAS `verify_*.py` scripts*, run them, and use them to overrule review legs.
That makes the orchestrator both the author of the instrument (it picks the framing, method, and truncation) and
the judge of its output — the single-engine bias the two-engine architecture exists to remove. The measured
relapse (2026-09-06): the orchestrator wrote `verify_F.py`/`verify_EG.py`, and **over-cleared** the E/N6
rep-invariance finding (it holds only in the σ_W→0 projection, not at the retained order claimed).

The change makes the intended rule explicit: **the orchestrator NEVER authors a CAS analysis/verification
script; Codex writes it, and Grok + the orchestrator review that script blind-to-output before it runs; only
mechanical fact-lookups (`len`, `grep`, a count) are the orchestrator's to run.** The orchestrator still owns the
*question* (run by Codex first) and the *adjudication* (G4).

## Artifact under review
The working-tree diff of `CLAUDE.md` vs `HEAD`:
`git -C /var/projects/toy_physics --no-pager diff CLAUDE.md`
Read the **full current** `CLAUDE.md` (working tree) and the **prior** version (`git -C /var/projects/toy_physics show HEAD:CLAUDE.md`) so you can judge preservation, not just the hunks. The touched controls are **E1** (§ "E · Produce evidence"), **G4** (§ "G · Review and advance"), the **at-a-glance** item 2, and a new **L-CAS** entry in the evidence ledger.

## What you are handed (read any; nothing you may use is withheld)
- the diff + both versions of `CLAUDE.md`;
- the skills the file defers runbooks to: `.claude/skills/build/SKILL.md` and `.claude/skills/review-legs/SKILL.md`
  (check the change does not now contradict them — e.g. their demand that a *leg/Codex* "show the script and its
  literal stdout" must remain coherent; that demand is on the script's author/reviewer, unchanged here);
- the two orchestrator-written CAS scripts named in the new war-story,
  `research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py` and `…_verify_EG.py`, and the physics
  adjudication `research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md` (verify the L-CAS
  war-story is factually accurate: those scripts are orchestrator-authored, and the adjudication records the E
  over-clear / σ_W→0-projection correction).

## The checks (this is a DOCUMENT review — read the sources, form your own view, quote both sides)
1. **No control weakened.** Does every pre-existing E1 obligation survive intact — a script PRINTS never asserts;
   emit-both-operands-and-residual THEN guard; interpretation belongs to the record; **a claim still carries the
   command that produced it, with commands + literal stdout in `_measurements/`**? The new "mechanical fact-lookup
   only" clause must ADD a routing rule for CAS instruments, ⛔ NOT give anyone (orchestrator or leg) an excuse to
   stop grounding a claim. Flag any reading where "it's a CAS instrument" becomes a way to skip the measurement
   entirely, or where R11 (no control dropped for cost) is loosened.
2. **Loophole actually closed, and no new one opened.** Is the authorship rule unambiguous — orchestrator NEVER
   authors a CAS analysis/diagnostic/verify script (not in `_measurements/`, not "under rule 13"); Codex writes it;
   Grok + orchestrator review it blind-to-output; the orchestrator runs only mechanical fact-lookups? Can it still
   be re-misread the way the old wording was? In particular: **is the mechanical-lookup vs CAS-instrument boundary
   crisp enough to apply under load**, or is there a gray zone (a one-line `sympy` simplify, a `subs`, a
   dimensional check) that will breed the same drift? If the boundary is fuzzy, say exactly where and propose the
   sharper line.
3. **Internal consistency.** Do E1, G4, the at-a-glance item 2, and the L-CAS ledger entry now agree with each
   other and with the rest of CLAUDE.md and the two skills? Any place the file still tells (or implies) the
   orchestrator should write a verification script, now contradicting the new rule? Any place a *leg* or *Codex*
   is wrongly swept up by the "never author a CAS script" prohibition (they MUST author scripts — the prohibition
   is on the ORCHESTRATOR only)?
4. **War-story accuracy.** Is the new L-CAS evidence-ledger entry factually correct against the real artifacts
   (the two verify scripts exist and are orchestrator-authored; the E/N6 over-clear + σ_W→0-projection correction
   is what the adjudication records)? ⛔ Flag any invented or overstated claim.

## Method / filter
Report a finding only if it catches a real way this change weakens a control, opens or leaves a loophole,
contradicts another part of the file/skills, or misstates the war-story. ⛔ Not style/wording preference unless
the wording is genuinely ambiguous in a way that will cause the drift again. Quote the specific `CLAUDE.md`
line(s) (old vs new) for each finding.

## Output
For each of 1–4: your finding + the quoted evidence. Separate CONFIRMED defects from unsettled questions. End
with: is the change correct to commit as-is, or the exact wording that must change (and why).
