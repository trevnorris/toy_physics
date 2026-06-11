You are Codex acting as the designer and builder for the toy_physics project (repo root: /var/projects/toy_physics). Your task is to fully execute the build directive at:

    docs/adversarial_audit_skill_directive.md

Read that directive first — it is the complete specification for Step 2 of the adversarial-audit campaign: build the `adversarial-audit` skill at `.claude/skills/adversarial-audit/`, its config, the `research/pde_ledger/redteam_adversarial/` manifest + report tree, graph wrappers, Codex-invocation plumbing, prompt/template set, and a working dry-run. Its "Binding references", "Out of scope", dry-run specification, and acceptance criteria are all binding.

Key constraints (restated from the directive):
- YOU design the runtime: internal lib/ layout, languages, function structure, which scan modalities are mechanical scripts vs agent prompts, config placement (standalone file vs `adversarial:` section), dry-run marking mechanism. The directive states requirements + acceptance criteria only.
- `docs/adversarial_audit_directive.md` is FROZEN — consume by reference/inclusion, never edit, never paraphrase.
- Do NOT modify `.claude/skills/redteam-audit/` — reuse its patterns (flock-serialized manifest writes, codex session capture/resume, render-prompt substitution) by adaptation in the new skill, not by editing the original.
- Runtime-rendered prompts go under `research/pde_ledger/redteam_adversarial/tmp_prompts/`, never `/tmp`.
- All agent-facing files YAML or markdown; `timeout 600` on every script execution path; flock on every manifest write.
- Defense Codex sessions are keyed PER-PARAMETER (not per-stage, not per-candidate).
- Allowed write set: `.claude/skills/adversarial-audit/`, `research/pde_ledger/redteam_adversarial/`, and `.redteam-config.yaml` ONLY if you choose the shared-config option (additive only; redteam-audit behavior unchanged).

Verification contract: implement, then run the dry-run against stages 003, 104, 105 exactly as the directive's dry-run specification describes. Iterate until the dry-run exits 0 end to end and every dry-run bullet is demonstrably met (003 honestly yields zero candidates; 104/105 yields ≥1 real candidate walked through provenance slice + Phase C prompt render + manifest state transitions; zero binding verdicts). Do not commit anything — leave the working tree dirty for review.

When done, report: (1) every file you created/changed and a one-line purpose for each; (2) the open design choices you made and why (config placement, modality mechanization, dry-run marking); (3) dry-run transcript evidence for each bullet of the dry-run specification; (4) any acceptance criterion you could not satisfy and why.
