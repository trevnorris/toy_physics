You are Codex acting as the fix-applier for the toy_physics project (repo root: /var/projects/toy_physics). Your task is to fully execute the directive at:

    docs/adversarial_audit_paper_integration_directive.md

Read that directive first — it is the complete specification (three items: resolve the paper canonical-source contradiction; add a falsification-stack disclosure section to the reader verification summary; add a definition-only audit-status vocabulary passage next to the claim-status firewall). Its "Path convention", "Out of scope", and acceptance-criteria sections are binding.

Key constraints (restated from the directive):
- All `paper/...` paths are relative to `research/pde_ledger/`.
- Zero edits under `research/pde_ledger/paper/stages/` — descriptions and definitions only.
- Do NOT edit `docs/adversarial_audit_directive.md` (frozen) or this directive.
- No appendix stubs, no per-stage prose patching, no script/notes/graph changes.
- Wording is your call; the directive states requirements and acceptance criteria only.

Verification contract: after applying the edits, build BOTH papers from `research/pde_ledger/paper/`:
    latexmk -pdf pde_ledger.tex
    latexmk -pdf pde_ledger_reader.tex
Iterate on your edits until both builds exit 0. Do not commit anything — leave the working tree dirty for review.

When done, report: (1) every file you changed and a one-line summary of each change; (2) build results (exit codes); (3) any acceptance criterion you could not satisfy and why.
