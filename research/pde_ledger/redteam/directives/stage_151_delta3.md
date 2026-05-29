---
unit_id: 151
batch: IV.6
iteration: 4
created_at: 2026-05-28T23:20:00-06:00
parent_directive: stage_151_delta2.md
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-28T23:22:19-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex delta-directive — unit 151 (iteration 4): add anti-footgun note

delta2 is correct and verified-runnable (exact 5-sample multi-point cross-check). One addition only: the script's docstring explains *what* it is but does not explicitly warn future editors NOT to "improve" it back into a fully-symbolic SymPy proof (which hangs). Add that explicit prohibition so a future human/AI does not repeat the 35-minute hang.

## F1-fix — insert an explicit "do not undo" comment block

Insert the following comment block at the top of `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`, immediately **after** the module docstring (i.e., right after the closing `"""`, before `from __future__ import ...`). Keep it verbatim (ASCII). Change NO logic — this is a comment-only edit.

```python
# === DO NOT "fix" this into a fully-symbolic SymPy proof. =======================
# SymPy CANNOT evaluate  integral_0^1 e^(-Pi_star*x) * {cos,cosh}(...) * x^n dx
# with a SYMBOLIC Pi_star -- it hangs indefinitely (confirmed 2026-05-28; the
# attempt was killed at 35 min and again at 19 min). This script is therefore an
# intentional EXACT, MULTI-POINT cross-check at concrete rational Pi_star values
# (symbolic in r1, r2, A_T, B_T, gprime). The full all-Pi_star symbolic proof is
# carried by the Mathematica engine:
#   mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl
# If you think you can make the SymPy side fully symbolic in Pi_star: you cannot
# (tried definite-integrate, indefinite+bounds, and trig->exp rewrite -- all hang).
# ================================================================================
```

After inserting, RUN `timeout 600 python3 scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py` once to confirm the comment did not break anything and it still exits 0 with all `[Pi=...] M1..M7` lines printing `= 0`. Append an `## Applied: F1-fix` block when done.

## Applied: F1-fix

- files_changed:
  - `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- summary: Inserted the verbatim anti-footgun comment block immediately after the module docstring; no logic was changed.
- deviation: none
