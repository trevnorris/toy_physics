---
unit_id: 013
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 013

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For the `paper_misalignment` finding (F1), do nothing. The orchestrator is holding for user resolution. Do not edit `paper/stages/stage_013.tex`, `notes/`, or the scripts to "fix" F1 unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding requires. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch `paper/stages/stage_013.tex`, `notes/`, or any prose documents. The red-team only modifies scripts.

---

## F1 — paper_misalignment

**Subtype:** `paper_missing_script_claim`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_013.tex:43-45` quote:
  > "Stage~013 exports the mouth-Taylor primitive map \eqref{eq:stage007-z0n0}--\eqref{eq:stage007-z4}."
  (The card's body declares only `z_0, n_0, z_2, z_4` and the narrative `Xi_load = n_0/N_0 + z_0/D_0`. No mention of `deltaP_2, deltaP_4, n_2, n_4`, or a mechanism-sieve test.)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:48` quote:
  > "Mouth-local Taylor map for projected source and flux corrections."
  (broad enough to plausibly cover the script's extra content, but does not enumerate `deltaP_2`, `deltaP_4`, or a sieve test by name)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:106-119` defines `deltaP2`, `deltaP4` and tests `assert_zero("d(delta P2)/dGprime", sp.diff(deltaP2_der, Gx) + 2*P/(D0*Delta**2))` at line 122.
- Same script lines 126-151 defines and tests the `qd_matrix`, `sh_matrix` mechanism sieve (determinants nonzero + linear solves trivial).
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:80-98` defines `n2, n4` and asserts their literal polynomial forms (the Mathematica script at lines 87-106 independently derives and verifies them).

## Resolve before fix_loop

The stage 013 paper card lists only four equation-numbered deliverables (`z_0, n_0, z_2, z_4`) plus the narrative `Xi_load` formula, but the scripts also verify `deltaP_2, deltaP_4, n_2, n_4`, the substantive coefficient `d(deltaP_2)/dG_w' = -2P/(D_0 Delta^2)`, and a 2x2 mechanism sieve on `(K_1, H_even)`. Which is authoritative?

Possible directions (the user picks one):
- (a) Paper card is incomplete -> add `deltaP_2, deltaP_4, n_2, n_4` and the mechanism-sieve claim to `paper/stages/stage_013.tex` as explicit deliverables (with brief formulas or named cross-references). No script change. Follow-up directive will instruct Codex to edit the paper card.
- (b) Script is overreaching -> trim the SymPy and Mathematica scripts to verify only `z_0, n_0, z_2, z_4`, `Xi_load`, and the one-sided Taylor lemma. Move `deltaP_2`, `deltaP_4`, the `n_2, n_4` checks, and the mechanism sieve to whichever downstream stage owns them (likely stage 014, where `K_1` and `H_even` are explicitly defined). Follow-up directive will instruct Codex to perform the moves.
- (c) Mixed -- keep `n_2, n_4` and `deltaP_2, deltaP_4` here (they are natural extensions of the same chain rule) but move only the K_1/H_even sieve to stage 014. Resolve by editing both paper and downstream scripts in a coordinated follow-up.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

---

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py` (insert a new block between the existing line 102 `H_even = ...` and line 104 `assert_zero("dXi/dPprime", ...)`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` (insert a new block between the existing line 124 `HEven = ...` closing `]` and line 126 `deltaP2 = ...`)

**Issue:**
The paper card declares `Xi_load = n_0/N_0 + z_0/D_0` (eq stage013-xi) as a deliverable. The scripts define `Xi` by writing its expansion in terms of `p1, P, d1, Delta, q1, Q, D0` directly, which silently assumes `N_0 = P^2/Delta^2`. The only test on `Xi` is `dXi/dPprime - 2/P == 0`, which exercises one partial derivative of one summand and would still PASS if the `z_0/D_0` half of `Xi` were dropped entirely. The paper's closed-form equality is not directly asserted.

**Required change:**

### F2 step 1 — SymPy

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`, insert the following block between line 102 (the existing `H_even = sp.simplify(...)` line) and line 104 (the existing `assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)` line). Note: `Xi` is already defined on line 100, so it is in scope here.

```python
    # Paper round-trip: verify Xi matches the paper's closed-form Xi_load = n0/N0 + z0/D0,
    # with the natural identification N0 = P^2/Delta^2 (so that n_0/N_0 = d/ds log(P^2/Delta^2)).
    z0_form = (Delta * q1 - Q * d1) / Delta**2
    n0_form = 2 * P * (Delta * p1 - P * d1) / Delta**3
    N0_form = P**2 / Delta**2
    Xi_paper = sp.simplify(((n0_form / N0_form + z0_form / D0).subs(subs_der)) / mu1)
    assert_zero("Xi matches paper closed form n0/N0 + z0/D0", Xi - Xi_paper)
```

### F2 step 2 — Mathematica

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`, insert the following block between line 124 (the closing `]` of the `HEven = FullSimplify[...]` definition, followed by `;`) and line 126 (the existing `deltaP2 = FullSimplify[...` line). `Xi` is already defined on lines 113-118.

```mathematica
(* Paper round-trip: verify Xi matches Xi_load = n0/N0 + z0/D0 with N0 = P^2/Delta^2. *)
z0Form = (Delta q1 - Q d1)/Delta^2;
n0Form = 2 P (Delta p1 - P d1)/Delta^3;
N0Form = P^2/Delta^2;
XiPaper =
  FullSimplify[
    ((n0Form/N0Form + z0Form/D0) /. subsDer)/mu1,
    Assumptions -> $Assumptions
  ];
assertZero["M5 Xi matches paper closed form n0/N0 + z0/D0", Xi - XiPaper];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 013` and `redteam exec-mathematica 013` and confirm (a) both scripts exit 0, (b) the Mathematica transcript adds a new line `OK M5 Xi matches paper closed form n0/N0 + z0/D0 residual = 0`. Self-test: with `subs_der` mapping `q1 -> mu1*Qx, d1 -> mu1*Dx, p1 -> mu1*Px`, both `Xi` (script-side) and `Xi_paper` (this directive) collapse to the same symbolic expression in `{P, Q, Delta, D0, Px, Dx, Qx}`, so the residual is 0 by construction; the new check fails only if a future edit changes one definition without the other.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`
- summary: Added direct Xi_load round-trip checks against `n0/N0 + z0/D0` with `N0 = P^2/Delta^2` in both audit scripts.
- deviation: Inserted after the current `Xi` definitions because the directive's named `H_even`/`deltaP2` insertion anchors are not present in this checkout.

---

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py` (insert the new block immediately after the F2 step 1 block, still before the existing line 104 `assert_zero("dXi/dPprime", ...)`).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` (insert the new block immediately after the F2 step 2 block, still before the existing line 126 `deltaP2 = ...`).

**Issue:**
`K_1 = -(z_2 + z_0/9)` and `H_even = -z_4 + (2/3) z_2 - z_0/27` use literal coefficients `1/9, 2/3, -1/27` that match stage 014's paper definitions, but the only assertions involving K_1 and H_even — the sieve determinants and trivial-solve checks — are insensitive to the specific coefficient values. A typo turning `1/9` into `1/91` would not fail any existing assertion.

**Required change:**

### F3 step 1 — SymPy

After the F2 step 1 block, insert:

```python
    # Coefficient sanity (catches typos in the 1/9, 2/3, -1/27 weights).
    # With (Qx, Dx, Sx, Hx, Px, Gx) = (1, 0, 0, 0, 0, 0) only q1 = mu1 is nonzero, so:
    #   z0 -> 1/Delta, z2 -> 0, z4 -> -1/Delta^2 (from the -Delta^2*q1 term of z4).
    # Hence K_1 = -(z2 + z0/9) = -(0 + 1/(9*Delta)) = -1/(9*Delta),
    #       H_even = -z4 + (2/3)*z2 - z0/27 = 1/Delta^2 + 0 - 1/(27*Delta).
    K1_only_q = (K1 * mu1).subs({Qx: 1, Dx: 0, Sx: 0, Hx: 0, Px: 0, Gx: 0})
    H_only_q = (H_even * mu1).subs({Qx: 1, Dx: 0, Sx: 0, Hx: 0, Px: 0, Gx: 0})
    assert_zero("K1 weight 1/9 anchor (q'=1 probe)", K1_only_q - sp.Rational(-1, 9) / Delta)
    assert_zero(
        "H_even weights -1, -1/27 anchor (q'=1 probe)",
        H_only_q - (1 / Delta**2 - sp.Rational(1, 27) / Delta),
    )
    # With (Qx, Dx, Sx, Hx, Px, Gx) = (0, 0, 1, 0, 0, 0) only s1 = mu1 is nonzero, so:
    #   z0 -> 0, z2 -> Q/Delta^2,
    #   z4 -> (-Delta^2*Hport*1 + 2*Delta*Q*S2*1)/Delta^4 = -Hport/Delta^2 + 2*Q*S2/Delta^3.
    # Hence K_1 = -(Q/Delta^2 + 0) = -Q/Delta^2,
    #       H_even = -(-Hport/Delta^2 + 2*Q*S2/Delta^3) + (2/3)*(Q/Delta^2) - 0
    #              = Hport/Delta^2 - 2*Q*S2/Delta^3 + (2/3)*Q/Delta^2.
    K1_only_s = (K1 * mu1).subs({Qx: 0, Dx: 0, Sx: 1, Hx: 0, Px: 0, Gx: 0})
    H_only_s = (H_even * mu1).subs({Qx: 0, Dx: 0, Sx: 1, Hx: 0, Px: 0, Gx: 0})
    assert_zero("K1 z2-weight anchor (S2'=1 probe)", K1_only_s - (-Q / Delta**2))
    assert_zero(
        "H_even weight 2/3 anchor (S2'=1 probe)",
        H_only_s
        - (Hport / Delta**2 - 2 * Q * S2 / Delta**3 + sp.Rational(2, 3) * Q / Delta**2),
    )
```

### F3 step 2 — Mathematica

After the F2 step 2 block, insert:

```mathematica
(* Coefficient sanity (catches typos in the 1/9, 2/3, -1/27 weights). *)
(* q' = 1, all others 0 -> K1 = -1/(9 Delta), H_even = 1/Delta^2 - 1/(27 Delta). *)
K1OnlyQ =
  FullSimplify[
    (K1 mu1) /. {Qx -> 1, Dx -> 0, Sx -> 0, Hx -> 0, Px -> 0, Gx -> 0},
    Assumptions -> $Assumptions
  ];
HOnlyQ =
  FullSimplify[
    (HEven mu1) /. {Qx -> 1, Dx -> 0, Sx -> 0, Hx -> 0, Px -> 0, Gx -> 0},
    Assumptions -> $Assumptions
  ];
assertZero["M6 K1 weight 1/9 anchor (q'=1 probe)", K1OnlyQ - (-1/(9 Delta))];
assertZero[
  "M6 H_even weights -1, -1/27 anchor (q'=1 probe)",
  HOnlyQ - (1/Delta^2 - 1/(27 Delta))
];

(* S2' = 1, all others 0 ->
   K1 = -Q/Delta^2,
   H_even = Hport/Delta^2 - 2 Q S2/Delta^3 + (2/3) Q/Delta^2. *)
K1OnlyS =
  FullSimplify[
    (K1 mu1) /. {Qx -> 0, Dx -> 0, Sx -> 1, Hx -> 0, Px -> 0, Gx -> 0},
    Assumptions -> $Assumptions
  ];
HOnlyS =
  FullSimplify[
    (HEven mu1) /. {Qx -> 0, Dx -> 0, Sx -> 1, Hx -> 0, Px -> 0, Gx -> 0},
    Assumptions -> $Assumptions
  ];
assertZero["M6 K1 z2-weight anchor (S2'=1 probe)", K1OnlyS - (-Q/Delta^2)];
assertZero[
  "M6 H_even weight 2/3 anchor (S2'=1 probe)",
  HOnlyS - (Hport/Delta^2 - 2 Q S2/Delta^3 + (2/3) Q/Delta^2)
];
```

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 013` and `redteam exec-mathematica 013` and confirms (a) both scripts exit 0, (b) the Mathematica transcript adds four new `OK ... residual = 0` lines for the four new `assertZero` calls. Self-test: the expected RHS values were derived by hand from the script's existing `z_0, z_2, z_4` literals with the named probe substitutions. If the coefficient `1/9` were perturbed to (say) `1/8`, the `K1 weight 1/9 anchor (q'=1 probe)` residual would become `(-1/8 + 1/9)/Delta = -1/(72*Delta)` and `assert_zero` would raise.

## Blocked: F3

- reason: The current stage 013 SymPy and Mathematica scripts no longer define `K1`, `H_even`, `deltaP2`, or `deltaP4`, so the requested coefficient-anchor blocks have no valid insertion point or in-scope symbols.
- question: Should the K1/H_even coefficient-anchor checks be omitted for the trimmed stage 013 scripts, or should a follow-up directive restore/move those definitions from the downstream stage that owns them?
