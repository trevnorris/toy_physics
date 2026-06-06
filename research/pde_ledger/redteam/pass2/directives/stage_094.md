---
unit_id: 094
batch: IV.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T05:13:26Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 094

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py:59-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:73-85`

**Issue:**
The "static-limit" block tied to paper Check #1 is tautological. `K_g2`/`K_g4` are assigned the literal `0` (`K_g2 = sp.Integer(0)` py:60; `Kg2 = 0` wl:73), so `eps_2 = Omega_Q^2 * K_g2 / K_pole` and `eps_4` are zero *by construction* and the `assert eps_2 == 0` / `expectZero["eps_2 (static limit)", eps2]` checks can never fail regardless of whether the orthogonality theorem holds. They do not consume the angular-overlap integrals proven zero in the block above (py:32-54, wl:52-69), so deleting or breaking that block would not affect these assertions. Additionally `assert c_pole + c_geom == 1` (py:68) / `expectZero["c_pole + c_geom - 1", ...]` (wl:83) only tests the arithmetic identity `1/4 + 3/4 == 1`; it never pins the named deliverable `c_pole = 1/4` (it would pass for any complementary pair). The fix must make the static-limit block CONSUME the computed orthogonality result and assert the two constants individually — without changing any emitted value (still `eps_2=0`, `eps_4=0`, `c_pole=1/4`, `c_geom=3/4`).

**Required change (SymPy):**

In the generic-cross-coefficient loop (py:47-54), accumulate the actually-computed `l=0/l=2` mass and Laplacian overlap integrals into running symbols so they are available below. Concretely, before the loop (after py:46) initialize:
```
Kg2_overlap = sp.Integer(0)
Kg4_overlap = sp.Integer(0)
```
and inside the loop, after `I_lap = ...` (py:50), add:
```
    Kg2_overlap += I_mass   # l=0/l=2 overlap moment (proven 0 mode-by-mode above)
    Kg4_overlap += I_lap    # l=0/(-Delta)l=2 overlap moment (proven 0 above)
```
Then replace the static-limit block (py:59-68) so that `K_g2`/`K_g4` are the DERIVED overlap sums, not bare literals:
```
Omega_Q, K_pole = sp.symbols('Omega_Q K_pole', positive=True)
# K_(g,2), K_(g,4) are the l=0 <-> l=2 overlap moments accumulated from the
# explicit S^2 integrals above; they are 0 BECAUSE those integrals vanished,
# not by assignment.
K_g2 = sp.simplify(Kg2_overlap)
K_g4 = sp.simplify(Kg4_overlap)
assert K_g2 == 0   # carries the orthogonality result, can fail if any overlap != 0
assert K_g4 == 0
eps_2 = sp.simplify(Omega_Q**2 * K_g2 / K_pole)
eps_4 = sp.simplify(Omega_Q**4 * K_g4 / K_pole)
assert eps_2 == 0
assert eps_4 == 0
c_pole = sp.Rational(1, 4)
c_geom = sp.Rational(3, 4)
assert c_pole == sp.Rational(1, 4)   # named deliverable (paper Check #1)
assert c_geom == sp.Rational(3, 4)
assert c_pole + c_geom == 1
```
Leave the print at py:69 unchanged (it already emits the correct values).

**Required change (Mathematica):**

Mirror identically. Before the `Do[...]` loop (after wl:47) initialize accumulators:
```
kg2Overlap = 0;
kg4Overlap = 0;
```
Inside the loop, after `lapCross = dOmega[y00*(-lapS2[y])];` (wl:59), add:
```
  kg2Overlap = kg2Overlap + overlap;
  kg4Overlap = kg4Overlap + lapCross;
```
(`overlap` and `lapCross` are the per-mode integrals already computed at wl:55,59.) Then replace the static-limit block (wl:73-83) so `Kg2`/`Kg4` are the derived sums:
```
(* K_(g,2), K_(g,4) carry the accumulated l=0 <-> l=2 overlap moments from the
   explicit S^2 integrals above; 0 BECAUSE those integrals vanished. *)
Kg2 = FullSimplify[kg2Overlap, Assumptions -> $Assumptions];
Kg4 = FullSimplify[kg4Overlap, Assumptions -> $Assumptions];
expectZero["K_(g,2) overlap moment", Kg2];
expectZero["K_(g,4) overlap moment", Kg4];
OmegaQ = Symbol["OmegaQ"];
Kpole = Symbol["Kpole"];
eps2 = OmegaQ^2 * Kg2 / Kpole;
eps4 = OmegaQ^4 * Kg4 / Kpole;
cPoleStatic = 1/4;
cGeomStatic = 3/4;
expectZero["eps_2 (static limit)", eps2];
expectZero["eps_4 (static limit)", eps4];
expectZero["c_pole - 1/4", cPoleStatic - 1/4];
expectZero["c_geom - 3/4", cGeomStatic - 3/4];
expectZero["c_pole + c_geom - 1", cPoleStatic + cGeomStatic - 1];
```
Leave the final `Print["eps_2 = ", ...]` (wl:84-85) unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 094` and `redteam exec-mathematica 094` and confirm: (a) `K_g2`/`K_g4` (py) and `Kg2`/`Kg4` (wl) are defined via the accumulated `domega`/`dOmega` overlap sums (no bare `sp.Integer(0)` / `= 0` literal assignment to them); (b) new assertions `c_pole == 1/4` / `c_geom == 3/4` (py) and `expectZero["c_pole - 1/4"]` / `expectZero["c_geom - 3/4"]` (wl) appear; (c) both scripts exit 0 and the transcripts still print `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`.

**Self-test (auditor pre-check):** Each per-mode `I_mass`/`overlap` and `I_lap`/`lapCross` was already asserted/expected zero above, so the accumulated sums are 0; substituting them gives `eps_2 = Omega_Q^2 * 0 / K_pole = 0` — the `assert`/`expectZero` still passes, but now ONLY because the integrals vanished (break any overlap and `K_g2 != 0` would propagate to `eps_2 != 0` and fail). The constants `1/4`, `3/4` are exact rationals; `c_pole - 1/4 = 0` and `c_geom - 3/4 = 0` hold. No `sp.diff`/`D` is introduced, so no variable-independence trap. No paper round-trip issue: the emitted values are unchanged and continue to match card:22 (`c_pole=1/4`) and notes:120,122.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl`
- summary: Static-limit contamination checks now derive `K_g2`/`K_g4` and `Kg2`/`Kg4` from accumulated l=0/l=2 overlap integrals and assert the two conservative split constants individually.
- deviation: none
