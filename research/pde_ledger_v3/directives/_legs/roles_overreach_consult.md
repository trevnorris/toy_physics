# Process question — how do we stop multi-agent role over-reach?

You are being consulted for your opinion on pipeline discipline. Answer concisely (no need to read files).

## The intended roles (a two-engine PDE-derivation pipeline)
- **ORCHESTRATOR (Claude):** writes directives / specs / decision-lists, launches reviews, adjudicates findings. Owns the QUESTION and the judgment. Should NOT author build-code or instruments.
- **BUILDER (codex/astra):** implements the deliverable from its directive, runs it once, reports. Nothing else.
- **REVIEW LEG (codex / grok / fresh-Claude):** independently reviews ONE artifact and returns a verdict. Nothing else.

## The observed failure — every agent reached outside its lane
1. The **builder**, handed a directive that mentioned "reviewed by fresh-Claude + Grok," spawned its OWN grok + claude review subprocesses and ran a review-until-clear on its own output — invalid self-review (whatever writes must not review itself), ~1h45m wasted.
2. A **review leg** (grok) built its own `run_all`/`watch` ablation-orchestration scripts and DEADLOCKED twice, instead of running a bounded review and reporting a verdict.
3. The **orchestrator** (Claude) reacted by AUTHORING load-bearing tooling (a kill-watchdog, monitor loops, diagnostics) and self-validated them with NO independent review — the exact single-engine bias the pipeline exists to remove.

Note: Claude Code called one-shot (`-p`) does NOT charge per token, so this is NOT a cost problem. It is about **role discipline, reliability, and not letting any one engine be builder+reviewer+judge at once.**

## Questions (concrete, minimal rules please)
1. What is the SIMPLEST set of structural constraints that keeps each agent strictly to its assigned task — no self-review, no orchestration-building, no self-cleared tooling?
2. Should the orchestrator's OWN operational scripts (watchdog/monitor/diagnostic) be codex-written + reviewed like everything else, or is there a lighter, honest rule that avoids the orchestrator both writing and clearing them?
3. What is the simplest reliable way to CALL a review leg so it performs a BOUNDED review and returns a verdict, without it building its own orchestration or spawning sub-agents?
