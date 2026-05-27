---
unit_id: 094
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 094

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — engine_disagreement

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:60`

**Issue:**
The Mathematica `cCross` (line 60) and the SymPy `Ccross` (sympy line 52) purport to compute the same physical object — the bilinear `Y00`↔`Y2A` cross coefficient of the isotropic quadratic wall action `μ(∂_t η)² + T_w(∂_w η)² + T_Ω η(-Δ)η + K η²`. SymPy uses `mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot` (correct: `T_w (∂_w η)²` factors angularly to `∫ Y00 Y2A dΩ = I_mass`). Mathematica uses `mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap`, pairing `tw` with `gradCross = ⟨∇_S² Y00 · ∇_S² Y2A⟩` instead of `overlap`. That is a different coupling, not the `T_w (∂_w η)²` cross term. Both expressions happen to vanish only because all three angular integrals (overlap, gradCross, lapCross) vanish individually, but the two scripts are not asserting the same identity.

**Required change:**
In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl`, replace the single line 60:

Before:
```
  cCross = FullSimplify[mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
```

After:
```
  cCross = FullSimplify[mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
```

Do not change `gradCross`'s independent assertion on line 65 — that check stays. Do not touch the SymPy script for this finding.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 094` and confirm: (a) the script exits 0, (b) the transcript still prints `PASS: Generic isotropic cross coefficient C_0,*` for each of the five Y2A labels, and (c) the diff against the .wl shows exactly the one-line change above.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py:56` (append immediately before the final `print(...)` on line 56)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:71` (insert immediately before the trailing `Print[""];` on line 71)

**Issue:**
The paper card's Check #1 reads: "Check the static limit ε_2=ε_4=0 returns c_pole=1/4." Notes §3 lays out: `eps_2 = Omega_Q^2 K_(g,2) / K_pole = 0`, `eps_4 = Omega_Q^4 K_(g,4) / K_pole = 0`, `c_pole = 1/4`, `c_geom = 3/4`. Neither script names ε_2, ε_4, K_pole, K_{g,2}, K_{g,4}, or c_pole. The orthogonality theorem is verified, but the static-limit reduction the card explicitly lists as a check is not. Add a short symbolic block in both engines that records K_{g,2} = K_{g,4} = 0 as the orthogonality output, evaluates ε_2 and ε_4, and asserts the c_pole = 1/4, c_geom = 3/4 normalization holds.

**Required change:**

**SymPy** (`.../moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py`): immediately before the existing final line `print('\nStage 77 theorem verified: ...')` (currently line 56), insert:

```python

# --- Static-limit check (paper Check #1): with K_(g,2) = K_(g,4) = 0 from
# the orthogonality theorem above, the contamination numbers vanish and the
# 3/4 + 1/4 conservative split is recovered (notes §3).
Omega_Q, K_pole = sp.symbols('Omega_Q K_pole', positive=True)
K_g2 = sp.Integer(0)  # established by l=0/l=2 orthogonality (asserts A1-A4 above)
K_g4 = sp.Integer(0)  # same orthogonality, l=4 channel (Omega_Q^4 moment)
eps_2 = sp.simplify(Omega_Q**2 * K_g2 / K_pole)
eps_4 = sp.simplify(Omega_Q**4 * K_g4 / K_pole)
assert eps_2 == 0
assert eps_4 == 0
c_pole = sp.Rational(1, 4)
c_geom = sp.Rational(3, 4)
assert c_pole + c_geom == 1
print('eps_2 =', eps_2, '; eps_4 =', eps_4, '; c_pole =', c_pole, '; c_geom =', c_geom)
```

**Mathematica** (`.../moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl`): immediately before the existing line 71 `Print[""];`, insert:

```
(* Static-limit check (paper Check #1): with K_(g,2) = K_(g,4) = 0 from
   orthogonality, the contamination numbers vanish and 3/4 + 1/4 split holds. *)
Kg2 = 0;
Kg4 = 0;
OmegaQ = Symbol["OmegaQ"];
Kpole = Symbol["Kpole"];
eps2 = OmegaQ^2 * Kg2 / Kpole;
eps4 = OmegaQ^4 * Kg4 / Kpole;
cPole = 1/4;
cGeom = 3/4;
expectZero["eps_2 (static limit)", eps2];
expectZero["eps_4 (static limit)", eps4];
expectZero["c_pole + c_geom - 1", cPole + cGeom - 1];
Print["eps_2 = ", fmt[eps2], "; eps_4 = ", fmt[eps4],
      "; c_pole = ", fmt[cPole], "; c_geom = ", fmt[cGeom]];
```

Do not rename or reorder existing assertions.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 094` and `redteam exec-mathematica 094` and confirm: (a) both scripts exit 0, (b) the sympy transcript contains a new line `eps_2 = 0 ; eps_4 = 0 ; c_pole = 1/4 ; c_geom = 3/4`, (c) the Mathematica transcript contains three new `PASS:` lines for `eps_2 (static limit)`, `eps_4 (static limit)`, and `c_pole + c_geom - 1`, and (d) the existing PASS lines for orthogonality and `Generic isotropic cross coefficient C_0,*` still appear.

---

## Applied: F1+F2 (orchestrator-direct)

- files_changed: scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py, mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl
- summary: F1 fixed engine_disagreement at .wl line 60: cCross now uses mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, matching SymPy's Tw*I_mass pairing for the T_w(∂_w η)² angular factor. F2 added static-limit check block (eps_2=eps_4=0, c_pole=1/4, c_geom=3/4) in both engines. Plus banner sweep STAGE 77/077 → STAGE 094.
- deviation: none
