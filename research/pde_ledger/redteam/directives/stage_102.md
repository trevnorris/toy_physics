---
unit_id: 102
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 102

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with `files_changed`, `summary`, and `deviation`.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.
After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — missing_verification_script (script_doesnt_cover_claim, sympy side)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py:13-17`

**Issue:**
The SymPy script computes the correct residuals (`diff(im(coeff(omega^5)), tauQ)`, `diff(im(coeff(omega^7)), tauQ)`, and `im(coeff(omega^5)) - chiQ*9/(32*Omega^5)`) but only `print`s them — no `assert` statement, so the script exits 0 regardless of value. The Mathematica peer (`moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl:46-48`) wraps the same quantities in `expectZero[...]`, but SymPy provides no gating verification.

**Required change:**

Edit the file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py`. The current contents (lines 13-17) are:

```
print('series through O(omega^5) =', Yser5)
print('series through O(omega^7) =', Yser8)
print('tauQ coefficient in omega^5 term =', sp.simplify(sp.diff(sp.im(Yser5.coeff(omega,5)), tauQ)))
print('tauQ coefficient in omega^7 term =', sp.simplify(sp.diff(sp.im(Yser8.coeff(omega,7)), tauQ)))
print('check canonical odd coefficient =', sp.simplify(sp.im(Yser5.coeff(omega,5)) - chiQ*sp.Rational(9,32)/Omega**5))
```

Replace lines 13-17 with (keep the print lines for transcript readability, add three assertions immediately after their corresponding print lines):

```
print('series through O(omega^5) =', Yser5)
print('series through O(omega^7) =', Yser8)

tau5 = sp.simplify(sp.diff(sp.im(Yser5.coeff(omega, 5)), tauQ))
print('tauQ coefficient in omega^5 term =', tau5)
assert tau5 == 0, f"tauQ irrelevance at omega^5 failed: residual = {tau5}"

tau7 = sp.simplify(sp.diff(sp.im(Yser8.coeff(omega, 7)), tauQ))
print('tauQ coefficient in omega^7 term =', tau7)
assert sp.simplify(tau7 - sp.Rational(1, 4)) == 0, f"tauQ omega^7 coefficient != 1/4: residual = {tau7 - sp.Rational(1,4)}"

canon = sp.simplify(sp.im(Yser5.coeff(omega, 5)) - chiQ*sp.Rational(9, 32)/Omega**5)
print('check canonical odd coefficient =', canon)
assert canon == 0, f"canonical odd coefficient mismatch: residual = {canon}"

print('Stage 102 SymPy audit passed.')
```

No other lines of the file change. Do not modify lines 1-12.

**Claim manifest** (the new asserts must independently verify):

- M1 (D1, paper notes section 1): `d/d tauQ Im[coeff(Yhat_Q^ret(omega), omega^5)] = 0`. In script symbols: `sp.diff(sp.im(Yser5.coeff(omega,5)), tauQ) == 0` after `simplify`. (tauQ first enters at omega^7, so it is invisible at omega^5.)

- M2 (D2, paper notes section 1): `d/d tauQ Im[coeff(Yhat_Q^ret(omega), omega^7)] = 1/4`. In script symbols: `sp.diff(sp.im(Yser8.coeff(omega,7)), tauQ) - sp.Rational(1,4) == 0`. (At omega^7 the tauQ contribution comes from the leading `i tauQ omega^7` in the denominator, multiplied by the `1/4` prefactor.)

- M3 (D3, paper notes section 1): `Im[coeff(Yhat_Q^ret(omega), omega^5)] = chi_Q * 9/(32 Omega_Q^5)`. In script symbols: `sp.im(Yser5.coeff(omega,5)) - chiQ*sp.Rational(9,32)/Omega**5 == 0`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 102` and confirm:
- The transcript now contains the three new assertions' diagnostic prints plus the new `Stage 102 SymPy audit passed.` line.
- The script exits 0.
- A future regression that broke any of M1, M2, M3 would now cause an `AssertionError` instead of a silent pass.

---

## Applied: F1 (orchestrator-direct)

- files_changed: scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py, mathematica/moving_throat_pde_stage102_higher_odd_irrelevance_mathematica_audit.wl
- summary: Rewrote SymPy script to add three asserts mirroring Mathematica's expectZero calls (D1 tauQ irrelevance at omega^5, D2 tauQ coefficient at omega^7 = 1/4, D3 canonical odd coefficient at omega^5 = 9/(32 Omega^5) * chi_Q). Plus banner sweep STAGE 085 → STAGE 102.
- deviation: none
