---
unit_id: 117
batch: IV.3
created_at: 2026-05-27
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 117

See `/var/projects/toy_physics/research/pde_ledger/redteam/tmp_prompts/audit_directive_117.md` for the full directive body. The structure is:

- **F1 — mathematica_transliteration** (Blocked): the `.wl` script is a line-by-line port of the `.py` script. Codex cannot safely invent an independent second-engine derivation; user must specify the alternative path.
- **F2 — tautological_check** (Blocked): D/N tube `kappa_c = 1/3` check at `sympy_audit.py:155-156` and `mathematica_audit.wl:112-113` is circular. Needs forward expression for `kappa0` from the D/N tube BVP (upstream physics); Codex cannot invent it.
- **F3 — tautological_check** (Blocked): `gamma_c = 1/9` check at `sympy_audit.py:157` and `mathematica_audit.wl:114` is pure algebraic substitution. Needs forward expression for `gamma0` from the bare mixed-channel construction; Codex cannot invent it.
- **F4 — tautological_check** (Apply): capstone classification table at `sympy_audit.py:170-182` and `mathematica_audit.wl:126-137` hardcodes the survivability flags. Wire flags to computed booleans drawn from sections 1-5 residuals. Concrete patch shape provided in the full directive.

Apply F4. Skip F1, F2, F3 with `## Blocked` blocks until user provides direction.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT run python or mathematica. Do NOT touch paper.tex, notes/, or any prose documents.
