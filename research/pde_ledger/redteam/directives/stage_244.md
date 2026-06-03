---
unit_id: 244
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-03T09:48:16-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 244

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Do NOT touch paper.tex or any prose document, EXCEPT the single notes line authorized in `## F0` below (user-approved 2026-06-03).

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F0 — notes correction (USER-APPROVED 2026-06-03)

**Target:** `notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md:366`

**Change (apply verbatim):** `196\sqrt2` → `128\sqrt2` (NOTES-ONLY numerical typo; the published paper card is UNAFFECTED). The SymPy script's `W_bulk`/`W_sess` closed forms reduce to the `128\sqrt2` coefficient (`E0 = 16·η_leak/π²` ⇒ `√2·E0²/2 = 128√2/π^{9/2}`), and the new `.wl` (F2) independently recomputes the leakage/work coefficients → cross-engine corroboration. Edit ONLY notes:366; do NOT change the script or the card. Codex applies; Claude reviews.

## Applied: F0

- files_changed:
  - `notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md`
- summary: Corrected the notes-only bulk-work coefficient typo from `196\sqrt2` to `128\sqrt2`.
- deviation: none

## F1 — insufficient_verification (support/orbit-split self-test trap)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py:130-142`

**Issue:** Section 4's support-versus-orbit split is verified by `sp.diff(expr, R_tr/R_target/eps_eta) == 0`, where `R_tr`, `R_target`, `eps_eta` are fresh symbols (declared L131) that never appear in `S_pull_varrho`, `W_bulk_pull_varrho`, or `W_sess_pull_varrho`. The derivatives are therefore identically zero by construction and the assertions can never fail — a vacuous proof of the paper's central structural claim (eq:app-part08-stage244-support-orbit-split). The split must instead be verified as a free-symbol containment statement.

**Required change:**
Replace the loop body (currently L133-142) so the load-bearing assertion is a free-symbol set test, not a derivative. Keep the orbit-symbol declaration on L131 (it is still used to *name* the excluded set). Use exactly this structure (you may keep the `diff` prints for human readability, but the `assert`s below are the load-bearing checks):

Before (L131-142):
```python
    R_tr, R_target, eps_eta = sp.symbols("R_tr R_target eps_eta", real=True)

    for expr, name in [(S_pull_varrho, "S_leak"), (W_bulk_pull_varrho, "W_bulk"), (W_sess_pull_varrho, "W_sess")]:
        d_Rtr = sp.simplify(sp.diff(expr, R_tr))
        d_Rtarget = sp.simplify(sp.diff(expr, R_target))
        d_epseta = sp.simplify(sp.diff(expr, eps_eta))
        print(f"d/dR_tr    {name} =", d_Rtr)
        print(f"d/dR_target{name} =", d_Rtarget)
        print(f"d/deps_eta {name} =", d_epseta)
        assert d_Rtr == 0
        assert d_Rtarget == 0
        assert d_epseta == 0
```

After:
```python
    R_tr, R_target, eps_eta = sp.symbols("R_tr R_target eps_eta", real=True)
    orbit_syms = {R_tr, R_target, eps_eta}
    support_syms = {Lam, varrho, eta_leak}

    for expr, name in [(S_pull_varrho, "S_leak"), (W_bulk_pull_varrho, "W_bulk"), (W_sess_pull_varrho, "W_sess")]:
        free = expr.free_symbols
        print(f"{name} free symbols      =", sorted(free, key=str))
        print(f"{name} orbit overlap     =", sorted(free & orbit_syms, key=str))
        print(f"{name} support coverage  =", sorted(free & support_syms, key=str))
        # Support-blindness: no orbit variable appears in the compiled observable.
        assert orbit_syms.isdisjoint(free)
        # Positive control so the disjointness is not vacuously true for a constant:
        assert support_syms.issubset(free)
```

Rationale for the positive control: each of the three compiled observables genuinely contains `Lam`, `varrho`, and `eta_leak` (S_leak is degree-1, the works degree-2 in `eta_leak`; all carry `Lam` and `varrho`). The `issubset` therefore passes for the correct expressions and would fail if a future edit collapsed an observable to a constant, guaranteeing the `isdisjoint` test is meaningful. The `isdisjoint` fails the instant any orbit symbol leaks into a compiled formula.

**Self-test confirmation (auditor):**
- `S_pull_varrho` free symbols = `{Lam, eta_leak, lam, mu_w, rho0, varrho}` → contains `{Lam, varrho, eta_leak}` ✓, disjoint from `{R_tr, R_target, eps_eta}` ✓.
- `W_bulk_pull_varrho` / `W_sess_pull_varrho` add `q` → still contain `{Lam, varrho, eta_leak}` ✓, disjoint from orbit set ✓.
- Neither `isdisjoint` nor `issubset` is vacuous: both observables are non-constant and orbit-free, so both asserts pass for the true expressions and would fire on a real regression.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 244` and confirms (a) the new `free symbols` / `support coverage` lines appear, (b) the `assert orbit_syms.isdisjoint(...)` and `assert support_syms.issubset(...)` checks are present and the derivative-based `assert d_* == 0` lines are gone, and (c) the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py`
- summary: Replaced the orbit-derivative support split check with free-symbol containment and support positive-control assertions.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl`

**Issue:** Stage 244 is SymPy-only, but every claim is independently verifiable in Mathematica (Gaussian integrals, algebraic substitutions, parity, free-symbol containment). The dual-engine rule requires a second engine wherever it is *possible*. Create a NEW independent-route `.wl` — native Mathematica primitives, a different decomposition than the `.py` (do NOT port the `.py` line-by-line; e.g. evaluate the Gaussian moments directly with `Integrate[..., {w, -Infinity, Infinity}]` and verify the works via the second-moment integral, ordering the pullback substitutions differently). It must independently verify the claim manifest below.

**Required change:**
Create the file with native Mathematica. Use an `expectZero[label, expr]` helper that calls `FullSimplify` and `Exit[1]` on nonzero, and an `expectTrue[label, bool]` helper that `Exit[1]` on `False`. Declare positivity assumptions matching the SymPy domains: `lam>0`, `E0>0`, `muw>0`, `rho0>0`, `q>0`, `Lam>0`, `eps` and `varrho` real with `0<eps<1`, `etaleak` real (NOT positive — parity needs the sign flip). Do not over-assume `etaleak>0`; the SymPy script's `eta_leak` positive assumption is acceptable for the algebra but for the parity check you must allow `etaleak < 0`, so declare `etaleak \[Element] Reals` and use `Assuming` only for the positive physical scales.

**Claim manifest** (each must have a corresponding `expectZero`/`expectTrue`):

- **M1 (boundary term):** With `W = Exp[-w^2/lam^2]/(lam Sqrt[Pi])`, `phi = 2 w Exp[-w^2/lam^2]/(Sqrt[Pi] lam^3)`, `Ew = -E0 phi`, `jw = muw rho0 Ew`: `Limit[W jw, w -> Infinity] - Limit[W jw, w -> -Infinity] == 0`. → `expectZero`.
- **M2 (leakage compiler, eq:sleak-e0):** `Integrate[D[W,w] jw, {w, -Infinity, Infinity}] - Sqrt[2] E0 muw rho0/(2 Sqrt[Pi] lam^3) == 0`. → `expectZero`.
- **M3 (bulk work, eq:bulk-work):** with `Jw = q jw`, `Integrate[Jw Ew, {w, -Infinity, Infinity}] - Sqrt[2] E0^2 muw q rho0/(2 Sqrt[Pi] lam^3) == 0`. → `expectZero`.
- **M4 (work-leak relation, eq:work-leak-relation):** `Wbulk - q E0 Sleak == 0` where `Wbulk` and `Sleak` are the integrals from M3/M2 (NOT the closed forms — derive both by integration so this is a genuine cross-relation, not a tautology). → `expectZero`.
- **M5 (session scalar, eq:session-work + eq:quadratic-law):** `2 Sqrt[2 Pi] lam Wbulk - 2 E0^2 muw q rho0/lam^2 == 0` AND `2 Sqrt[2 Pi] lam Wbulk - 4 Pi q lam^4 Sleak^2/(muw rho0) == 0`. → two `expectZero`.
- **M6 (pullback, eq:pitr + eq:e0-pullback + eq:work-epsilon):** with `Cmix = 8 Lam (1-eps)/Pi^2`, `Pitr = (4/3) Cmix`: `Pitr - 32 Lam (1-eps)/(3 Pi^2) == 0`; substituting `1-eps -> (3/2) varrho`: `Pitr /. ... - 16 Lam varrho/Pi^2 == 0`; with `E0pull = etaleak Pitr`, substitute into the M2/M3/M5 closed forms and verify the compiled `(1-eps)` forms: `Spull - 16 Sqrt[2] etaleak muw rho0 Lam (1-eps)/(3 Pi^(5/2) lam^3) == 0`, `Wbulkpull - 512 Sqrt[2] etaleak^2 muw q rho0 Lam^2 (1-eps)^2/(9 Pi^(9/2) lam^3) == 0`, `Wsesspull - 2048 etaleak^2 muw q rho0 Lam^2 (1-eps)^2/(9 Pi^4 lam^2) == 0`. → `expectZero` each.
- **M7 (support/orbit split — DO NOT use D[...]):** verify free-symbol containment, mirroring the corrected SymPy F1. For each compiled observable `Spull, Wbulkpull, Wsesspull`: `expectTrue["Spull orbit-free", FreeQ[Spull, Rtr] && FreeQ[Spull, Rtarget] && FreeQ[Spull, epseta]]` and a positive control `expectTrue["Spull support-dependent", (! FreeQ[Spull, Lam]) && (! FreeQ[Spull, varrho]) && (! FreeQ[Spull, etaleak])]`. (Here `Rtr, Rtarget, epseta` are symbols that simply never appear; `FreeQ` returning `True` is the meaningful pass, and the positive control prevents a vacuous "constant is free of everything" pass.) Do NOT use `D[Spull, Rtr]` — that reproduces the exact vacuous trap F1 corrects.
- **M8 (parity, eq:parity):** `(Spull /. etaleak -> -etaleak) + Spull == 0` (odd); `(Wbulkpull /. etaleak -> -etaleak) - Wbulkpull == 0` and `(Wsesspull /. etaleak -> -etaleak) - Wsesspull == 0` (even). → three `expectZero`. This is why `etaleak` must be `\[Element] Reals`, not assumed positive.
- **M9 (recovery slice):** `(Spull /. etaleak -> 0) == 0`, `(Wbulkpull /. etaleak -> 0) == 0`, `(Wsesspull /. etaleak -> 0) == 0`. → three `expectZero`.

**Independence requirement:** Derive `Sleak` and `Wbulk` by *integration* (M2, M3) and then build M4/M5 from those integrated quantities, so the relations are genuine engine-internal cross-checks rather than re-statements of the SymPy closed forms. Do not import or echo the SymPy intermediate-variable names/choreography; structure the file around the claim manifest M1–M9.

**Self-test confirmation (auditor):**
- M7 variable-independence: `Rtr/Rtarget/epseta` are absent from the compiled observables ⇒ `FreeQ[...]` returns `True` (meaningful), and the positive control `!FreeQ[..., Lam/varrho/etaleak]` returns `True` because each observable genuinely contains those ⇒ not vacuous. No `D[...]` used, so the F1 trap is not reproduced.
- M2/M3 parity of integrands: `W'` is odd × `jw` odd = even ⇒ M2 integral nonzero (matches the `√2/(2√π)` moment); `Jw Ew` = odd² = even ⇒ M3 integral nonzero. Both `expectZero` targets are `integral − nonzero_closed_form`, which correctly reduces to 0.
- M8 parity needs `etaleak` to take negative values; declaring it `\[Element] Reals` (not `> 0`) is required, else `etaleak -> -etaleak` under a positivity assumption is inconsistent.
- M6 constants `512√2/9`, `2048/9` and `Pitr=32/(3π²)` are mutually consistent under `1-eps -> (3/2)varrho` and match both the `.py` (L105,107) and the paper card (eq:work-epsilon, eq:pitr).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 244` and confirms the new `.wl` exists at the path above, contains the M1–M9 checks, and exits 0 with all `expectZero`/`expectTrue` passing. The engine cross-check (SymPy vs Mathematica closed forms) then becomes available.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M9 with integral-derived leakage/work checks, pullbacks, free-symbol support/orbit checks, parity, and recovery.
- deviation: none
