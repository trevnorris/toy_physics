---
unit_id: 150
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 150

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py:37`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl:37`

**Issue:**
`S_q` is defined as `simplify(diff(T_q, x).subs(x, 0))` (line 37 of the `.py`; analogous line 37 of the `.wl`). The later assertion `T_q'(0) - S_q == 0` (line 45 of the `.py`; line 44 of the `.wl`) is therefore `X - simplify(X) == 0`, an algebraic identity that cannot fail. Replace the definition of `S_q` with the hand-derived closed form `S_q = A_q*k - C_q*Pi` obtained by differentiating `T_q(x) = A_q sinh(k x) - C_q cosh(k x) + C_q exp(-Pi x)` at `x=0`. The downstream assertion line is unchanged; only the definition of `S_q` changes.

**Required change:**

(1) In `scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py`, replace line 37.

Before (line 37):
```python
Sq = sp.simplify(sp.diff(Tq, x).subs(x, 0))
```

After (line 37, plus a one-line comment above it for clarity):
```python
# Hand-derived closed form: T_q'(0) = Aq*k - Cq*Pi
# (differentiate Tq = Aq*sinh(k*x) - Cq*cosh(k*x) + Cq*exp(-Pi*x) at x=0).
Sq = Aq*k - Cq*Pi
```

Lines 38–62 are unchanged.

(2) In `mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl`, replace line 37.

Before (line 37):
```mathematica
sQ = FullSimplify[D[tq, x] /. x -> 0, Assumptions -> $Assumptions];
```

After (line 37, plus a one-line comment above it for clarity):
```mathematica
(* Hand-derived closed form: T_q'(0) = aq*k - cq*p
   (differentiate tq = aq*Sinh[k*x] - cq*Cosh[k*x] + cq*Exp[-p*x] at x=0). *)
sQ = aq*k - cq*p;
```

Lines 38–63 are unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 150` and `redteam exec-mathematica 150`. Both scripts must (a) exit 0, (b) still report `T_q'(0)-S_q = 0` (and `PASS: T_q'(0)-S_q` in the Mathematica transcript), (c) display `S_q(Pi) =` showing a compact form involving `Aq` and `Cq` (e.g., `pi*Aq/2 - Cq*Pi` or `aq*k - cq*p`) rather than the previous expanded derivative, and (d) still report `R''(0) - target = 0` (load-bearing check is not affected by this edit).
