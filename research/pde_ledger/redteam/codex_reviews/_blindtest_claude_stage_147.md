---
test: detection-parity blind review
stage: 147
reviewer: claude-agent (general-purpose, clean context)
blind: true   # agent was forbidden from reading redteam/codex_reviews/
date: 2026-05-28
purpose: >
  Calibration test — can a clean Claude review agent independently detect the
  same defects Codex's review found, without seeing Codex's output? Result used
  to decide the remediation role split (who fixes vs who verifies).
---

# Blind detection-parity test — stage 147

The agent received the same inputs Codex had (script .py/.wl, saved outputs,
directive, paper card) but was explicitly forbidden from reading the Codex review
under `redteam/codex_reviews/`. It produced 5 findings independently.

## Overlap vs Codex (stage_147.md)

| Issue | Codex | Claude blind agent |
|---|---|---|
| chain-rule check = reassociated closed form (X−X) | R1 (py:69) | R1 (py:66–75) — MATCH |
| Wcenter re-expands its own definition | R2 (py:84) | R2 (py:88–99) — MATCH |
| g_*/S_* resubstitution of identical expression | R4 (py:107) | R3 (py:101–113) — MATCH |
| .wl is a transliteration of .py | R5 | R5 — MATCH |
| .wl wCenterConst samples one point x→1/2 | R3 (insufficient) | observed in reasoning, not raised separately |
| root condition g_*=g_minus asserted nowhere | (missed) | R4 (paper_misalignment) — EXTRA, Codex missed |

## Conclusion

Parity — arguably better on this sample. A clean Claude agent is a competent
independent detector of the tautological/transliteration/paper-misalignment
defects in the orchestrator-direct era, not a rubber stamp. Supports using a
clean Claude agent as the independent VERIFIER of Codex-applied remediation fixes
(fixer = Codex, verifier = Claude; fixer ≠ verifier).

(Agent's full 5-finding report is preserved in the conversation transcript.)
