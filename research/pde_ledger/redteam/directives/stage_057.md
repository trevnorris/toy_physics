---
unit_id: 057
batch: III.2
created_at: 2026-05-26T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-26T03:04:39-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false  # Q2=(a) applied — see Applied: F1 block below
---

# Codex directive — unit 057

This unit has one `paper_misalignment` finding (F1) that requires user resolution before Codex may proceed on it. Findings F2, F3, F4 are script-side only and may be applied normally.

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit `paper/`, `notes/`, or scripts to "fix" F1 unless the user has explicitly chosen a direction in a follow-up directive.

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/`, `notes/`, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

---

## F1 — paper_misalignment

**Subtype:** `script_missing_paper_claim`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_057.tex:23` quote: "It is monotone in the constructive Peclet direction, with zero-bias value"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:140-148` quote: "Increasing in transport bias `Pe` — From Stage 39, `dOmega_Pe/dPe > 0` on the constructive branch `Pe >= 0`, so `partial_Pe zeta_0^(Pe+R) > 0`."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage057_physical_parameter_map.md:303` quote: "this map is monotone increasing in `Pe` and monotone decreasing in both `eta` and `kappa`"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:62-73` — only `partial_kappa` and `partial_y` are computed; no `sp.diff(zeta_phys, Pe)`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:56-67` — only `D[zetaPhys, kappa]` and `D[zetaPhys, y]`; no `D[zetaPhys, pe]`.

## Resolve before fix_loop

The paper card and notes claim `partial_Pe zeta_0^(Pe+R) > 0` on the constructive branch (`Pe >= 0`); the notes mark this as carry-forward from Stage 056 (Stage 39 in old numbering), but the stage-057 scripts neither check it nor reference its upstream source. Which direction should the resolution take?

Possible directions (the user picks one):
- (a) Stage-057 scripts must check Pe-monotonicity locally → in both scripts, after the existing `partial_kappa`/`partial_y` blocks, add a numerical sweep (e.g. `Pe in {0.1, 0.5, 1, 2, 5, 10}`, `kappa = 1`, `y = pi/4`) confirming `D[zeta_phys, Pe] > 0`, and a comment explaining the carry-forward from Stage 056. Re-run sympy+mathematica.
- (b) Pe-monotonicity stays as Stage 056 carry-forward → follow-up directive authorizes a paper/notes edit adding a single line ("Pe-monotonicity established at Stage 056; not re-verified here"). No stage-057 script change. The Codex pass does NOT include this edit; the user must issue a separate directive after choosing direction (b).
- (c) Both → user adds the carry-forward note AND a confirming numerical spot-check at this stage.

The orchestrator will not invoke Codex on F1 until the user has chosen a direction.

## Applied: F1

- direction: (a) — local Pe-monotonicity numerical sweep added at Stage 057.
- files_changed:
  - `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:75-83` — added `dPe` symbolic derivative + sweep over `Pe in {1/10, 1/2, 1, 2, 5, 10}` at `kappa=1, y=pi/4` with `if val <= 0: raise AssertionError` guard; new PASS line: `partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS`.
  - `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:69-80` — analogous `Module` block using `D[zetaPhys, pe]` + same sweep; new PASS line: `partial_Pe zeta > 0 on constructive branch (numerical sweep)`.
- summary: Stage-057 scripts now anchor the paper's Pe-monotonicity claim locally; carry-forward comment cites Stage 056 §4 covariance identity.
- deviation: none.
- post-edit verification: both scripts re-ran exit 0; new PASS lines appeared in transcripts.

Codex must NOT re-apply F1.

---

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:70` (kappaMax definition)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:83` (kappaReq definition)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:82,96-101` (yReqSq usage / y_req defining equation)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:107-112` (y_req defining equation)

**Issue:**

(a) The Mathematica `kappa_max identity` check (L78) compares `kappaMax` to its own literal definition (L70). Tautological.
(b) The Mathematica `kappa_req identity` check (L92-95) compares `kappaReq` to its own literal definition (L83). Tautological.
(c) The Mathematica/SymPy `y_req defining equation` check substitutes `y^2 -> y_req_sq` into the defining equation, but `y_req_sq` was defined as exactly the expression that makes the substitution collapse to `zeta_req = zeta_req`. The cancellation is structural (`kappa + y_req_sq = (Omega_Pe^2/zeta_req)(kappa+pi^2/4)` by construction, so the ratio is `zeta_req` regardless of whether `y_req_sq` is the correct solution). This is the same tautology v1's directive thought it had fixed — it merely swapped a self-subtraction tautology for a substitution-cancellation tautology.

The real fix is to derive `kappaMax`, `kappaReq`, and `y_req_sq` via `Solve`/`sp.solve` instead of writing them as literal closed forms, then compare those derived expressions to the closed forms.

**Required change:**

**(a) Mathematica `kappa_max` — derive via `Solve`.** At line 70 of `mathematica/.../stage057_..._mathematica_audit.wl`, the current block is:

```mathematica
kappaMax = FullSimplify[Pi^4/(4 (4 zetaReq - Pi^2)), Assumptions -> zetaReq > Pi^2/4];
```

Replace with:

```mathematica
kappaMaxSol = Solve[zetaReq == (Pi^2/4) (kappa + Pi^2/4)/kappa, kappa, Reals];
kappaMax = FullSimplify[kappa /. First[kappaMaxSol], Assumptions -> zetaReq > Pi^2/4];
```

The existing `kappa_max identity` check at L78 (which asserts `kappaMax - Pi^4/(4 (4 zetaReq - Pi^2)) == 0`) remains unchanged — it now becomes a real solve-vs-closed-form comparison.

**(b) Mathematica `kappa_req` — derive via `Solve`.** At line 83, the current block is:

```mathematica
kappaReq = FullSimplify[(omegaPe^2 Pi^2/4 - zetaReq y^2)/(zetaReq - omegaPe^2), Assumptions -> $Assumptions];
```

Replace with:

```mathematica
kappaReqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + y^2), kappa, Reals];
kappaReq = FullSimplify[kappa /. First[kappaReqSol], Assumptions -> $Assumptions];
```

The existing `kappa_req identity` check at L92-95 remains unchanged — it now becomes a real solve-vs-closed-form comparison.

**(c-sympy) SymPy `y_req` — derive via `sp.solve`.** At lines 107-112 of `scripts/.../stage057_..._sympy_audit.py`, the current block is:

```python
expect_zero(
    "y_req defining equation",
    zeta_req - sp.simplify(
        Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y_req_sq)
    ),
)
```

Replace with:

```python
y_req_sq_solved = sp.solve(
    sp.Eq(zeta_req, Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y**2)),
    y**2,
)[0]
expect_zero(
    "y_req identity",
    y_req_sq - y_req_sq_solved,
)
```

Leave the definition of `y_req_sq` at lines 92-93 unchanged. Leave the `print("y_req^2 = ...")` line unchanged.

**(c-mathematica) Mathematica `y_req` — derive via `Solve`.** At lines 96-101, the current block is:

```mathematica
expectZero[
  "y_req defining equation",
  zetaReq - FullSimplify[
    omegaPe^2 (kappa + Pi^2/4)/(kappa + yReqSq),
    Assumptions -> $Assumptions
  ]
];
```

Replace with:

```mathematica
yReqSqSol = Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + ysq), ysq, Reals];
yReqSqSolved = FullSimplify[ysq /. First[yReqSqSol], Assumptions -> $Assumptions];
expectZero["y_req identity", yReqSq - yReqSqSolved];
```

Leave the definition of `yReqSq` at line 82 unchanged. Leave the `Print["y_req^2 = ", ...]` at line 86 unchanged.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 057` and `redteam exec-mathematica 057` and confirm:
1. Both scripts exit 0.
2. SymPy transcript prints `y_req identity = 0` (replacing the earlier `y_req defining equation = 0`).
3. Mathematica transcript prints `PASS: y_req identity` (replacing `PASS: y_req defining equation`).
4. Mathematica transcript still prints `PASS: kappa_max identity`, `PASS: kappa_req identity` — but the underlying source now derives `kappaMax`, `kappaReq` via `Solve` (verified by grep for `Solve[zetaReq == (Pi^2/4)`, `Solve[zetaReq == omegaPe^2`, `Solve[zetaReq == omegaPe^2 (kappa + Pi^2/4)/(kappa + ysq)`).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
  - `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
- summary: Replaced tautological closed-form threshold checks with Solve/sp.solve-derived comparisons for kappa_max, kappa_req, and y_req_sq.
- deviation: none

---

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:35` (symbol assumptions) and after L73 (new sign check)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:35-36` ($Assumptions) and after L67 (new sign check)

**Issue:**

Notes deliverable (4c) claims `partial_kappa zeta < 0` because `0 < y < pi/2`. The scripts assert the algebraic form of `partial_kappa zeta` but never its sign, and neither declares `y < pi/2`. A literal sign sweep across `y in (0, pi/2)` (with sample `Pe`, `kappa`) is the minimum acceptable verification.

**Required change:**

**(a) SymPy.** No change to L35; add after the existing block at L70-73 (so it becomes the new L74 onward, BEFORE the existing `# Constructive-branch closure ceiling.` comment at L75):

```python
# Sign check on partial_kappa over 0 < y < pi/2 (from y tan y = eta, eta finite).
# Notes deliverable (4c) requires partial_kappa zeta < 0 on the constructive branch.
for y_val in (sp.pi / 8, sp.pi / 6, sp.pi / 4, sp.pi / 3, sp.Rational(7, 16) * sp.pi):
    val = sp.simplify(dkappa.subs({Pe: 1, kappa: 1, y: y_val}))
    if val >= 0:
        raise AssertionError(f"partial_kappa zeta sign failed at y={y_val}: {val}")
print("partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS")
```

**(b) Mathematica.** At line 35-36, replace
```mathematica
$Assumptions =
  Element[{pe, kappa, y, zetaReq, x}, Reals] && pe > 0 && kappa > 0 && y > 0 && zetaReq > 0 && x > 0;
```
with
```mathematica
$Assumptions =
  Element[{pe, kappa, y, zetaReq, x}, Reals] && pe > 0 && kappa > 0 && y > 0 && y < Pi/2 && zetaReq > 0 && x > 0;
```

Add after the existing block at L64-67 (before line 69 `zetaMax = ...`):

```mathematica
(* Sign check on partial_kappa over 0 < y < Pi/2 (from y tan y = eta, eta finite).
   Notes deliverable (4c) requires partial_kappa zeta < 0 on the constructive branch. *)
Module[{yvals, signOk},
  yvals = {Pi/8, Pi/6, Pi/4, Pi/3, 7 Pi/16};
  signOk = AllTrue[yvals, TrueQ[N[(D[zetaPhys, kappa] /. {pe -> 1, kappa -> 1, y -> #})] < 0] &];
  If[signOk,
    pass["partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)"],
    fail["partial_kappa zeta sign sweep"]
  ]
];
```

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 057` and `redteam exec-mathematica 057` and confirm:
1. Both scripts exit 0.
2. SymPy transcript prints `partial_kappa zeta < 0 on 0 < y < pi/2 (numerical sweep): PASS`.
3. Mathematica transcript prints `PASS: partial_kappa zeta < 0 on 0 < y < Pi/2 (numerical sweep)`.
4. Mathematica `$Assumptions` line now includes `y < Pi/2`.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- summary: Added the requested partial_kappa negative sign sweeps and tightened the Mathematica y-domain assumption to y < Pi/2.
- deviation: none

---

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:37` (insert before the existing `omegaPe` line)

**Issue:**

The Mathematica `aKX`/`xSub`/`aKKappa` substitution chain (L42-44) mirrors the SymPy `A_K_x`/`x_sub`/`A_K_kappa` chain. The Mathematica script should derive `A_K` independently from the physical operator quantities described in the notes (section 2), not by porting SymPy's `x`-substitution algebra.

**Required change:**

Insert the following block into the Mathematica script immediately AFTER the existing `omegaPe = FullSimplify[...]` block at L38-41 (i.e., as the new block starting around L42) and BEFORE the existing `aKX = FullSimplify[...]` line:

```mathematica
(* Independent derivation of A_K from the physical support operator.
   Notes section 2: K_W^(eff) = (T_X/L^2)(kappa + Pi^2/4),
                    K_phi,0^(eff) = (T_X/L^2)(kappa + y^2),
   so A_K = K_W^(eff)/K_phi,0^(eff) = (kappa + Pi^2/4)/(kappa + y^2). *)
Module[{KW, Kphi0, aKPhys},
  KW = KX + Pi^2 TX/(4 LL^2);
  Kphi0 = KX + TX y^2/LL^2;
  aKPhys = FullSimplify[KW/Kphi0,
    Assumptions -> KX > 0 && TX > 0 && LL > 0 && y > 0 && y < Pi/2];
  aKKappaFromPhys = FullSimplify[aKPhys /. KX -> kappa TX/LL^2,
    Assumptions -> kappa > 0 && TX > 0 && LL > 0 && y > 0 && y < Pi/2];
  expectZero["A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)",
    aKKappaFromPhys - (kappa + Pi^2/4)/(kappa + y^2)]
];
```

This adds an independent physical-operator derivation in the Mathematica script without touching the existing `aKX`/`xSub`/`aKKappa` chain (which now functions as a consistency check rather than the primary derivation). Leave L42-44 unchanged.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 057` and confirm:
1. Mathematica exits 0.
2. Transcript contains `PASS: A_K(physical) - (kappa+Pi^2/4)/(kappa+y^2)`.
3. Source contains the symbols `KX`, `TX`, `LL` (which do not appear in the SymPy script).

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- summary: Added an independent Mathematica derivation of A_K from the physical support operator quantities KX, TX, and LL before the existing substitution chain.
- deviation: none
