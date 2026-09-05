# CLAUDE.md restructure — review record (rule 2)

The 2026-09-05 CLAUDE.md restructure (17 rules → M/E/G/S + quick-ref + artifact→review-discipline table +
evidence ledger). **Astra designed** it (user-authorized, one-off); **orchestrator authored** the rewrite;
**Codex `gpt-5.6-sol` xhigh + Grok `grok-4.6`** reviewed (whatever-writes-doesn't-review: Astra≠reviewer).
Design + preservation map: `CLAUDE_streamline_proposal_2026-09-05.md`.

## Commands

```
# Design (authorized astra):
codex exec -m gpt-6-astra -c model_reasoning_effort=high --sandbox danger-full-access \
  "$(<directives/_legs/claude_md_streamline_prompt.md)"        # → CLAUDE_streamline_proposal_2026-09-05.md
# Preservation check, each round (OLD=git HEAD:CLAUDE.md snapshot vs NEW=working tree):
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
  "$(<directives/_legs/claude_md_preservation_check_prompt.md)"
grok --prompt-file directives/_legs/claude_md_preservation_check_prompt.md --cwd /var/projects/toy_physics \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain
```

## Round-by-round (both legs each round; raw reports ephemeral in scratchpad, verdicts transcribed here)

| Round | Grok | Codex sol | Real findings → fix |
|---|---|---|---|
| r1 | PRESERVATION-COMPLETE (nits) | **BLOCK** | ⚠ **Both** caught real defects the orchestrator introduced: (a) **fabrication** — 2 invented ledger entries (L-R15/L-R16) asserting historical outcomes not in the original; (b) **R7 weakened** in the table ("normally O"); (c) R17 framing truncated; (d) non-physics row contradiction; (e) gate record added as new authority. **All fixed.** |
| r2 | PRESERVATION-COMPLETE | **BLOCK** | Codex: ledger claimed "verbatim" but L-R12/L-R14/L-R17 were excerpts/edited. **Fixed** (intro re-worded to honest "excerpt where the rest is in the rule"; L-R12 restored verbatim). |
| r3 | PRESERVATION-COMPLETE | **BLOCK** (editorial) | "material issues" undefined (→ tied to R10 scope); R17 not a continuous verbatim quote (→ full account moved verbatim into L-R17, M3 terse, no duplication); "links" overstatement (→ fixed); gate-record "hence" (→ softened). **All fixed.** |
| r4 | — | — | Final confirmation of the r3 fixes; **killed twice by a mobile-navigation interrupt** before completion (not a resource kill — no orphans, 23 GB free). |

## Verdict (rule 10 stopping point; user-approved commit)

Across the 3 **completed** rounds, **both legs confirm all 17 controls are preserved**, the "Repository
infrastructure" annex is **byte-for-byte identical**, the **A1 disambiguation** (decision-list one-pass vs
spec/script/record review-until-clear) is correct and consistent across the quick-ref/table/G2/G4, the
reviewer pairings match old R7 (O→Codex+Grok; Cx→fresh Claude+Grok), and the two commit meanings coexist.
Grok returned **PRESERVATION-COMPLETE** in r2 and r3. Every specific item Codex raised across r1–r3 was
**folded and self-verified**; r3's residual BLOCK was editorial (the "review manufacturing findings" regime,
rule 10). The r4 confirmation could not complete due to repeated mobile interrupts; the user
(AskUserQuestion, this session) approved committing at the rule-10 stopping point. ⇒ **committed as a
presentation-only restructure — no control dropped, weakened, or invented.** Revertable if a later clean
confirmation round finds otherwise.
