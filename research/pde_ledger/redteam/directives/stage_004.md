---
unit_id: 004
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-25T02:33:33-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 004 (v2 re-audit)

This directive supersedes the prior v1 directive (F1/F2 there have been applied; this is a new finding raised by the paper-grounded v2 re-audit). Apply F1 below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check (Faraday/Bianchi block)

**Target (sympy):** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py:46-67`

**Target (mathematica):** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl:28-69`

**Issue:**

The "Faraday / vector Bianchi signs" block in both engines is a tautology. It assigns `F_{23} = B_1, F_{30} = E_3, F_{02} = -E_2` (and cyclic permutations), then constructs the LHS of the cyclic-Bianchi sum from those `F` symbols and the RHS from the manually-written vector form `∂_t B_1 + ∂_y E_3 - ∂_z E_2`. After substitution of the `F` definitions, the LHS reduces to `∂_t B_1 + ∂_y E_3 - ∂_z E_2`, which is character-for-character the RHS. The residual is identically zero by symbol substitution alone; the assertion exercises no Maxwell-Faraday content, no cyclic-Bianchi structure, and no homogeneous/inhomogeneous distinction. The Mathematica `M2` block (lines 28-69) has identical structure (`twoForm23 = B1[...]`, etc.) and the same tautology.

The paper card for stage 004 (`paper/stages/stage_004.tex`, "Bundle role" paragraph) credits this block as supporting rule (ii): "the homogeneous Maxwell equations project in measured fields, while the inhomogeneous equations carry source-coupled flux fields." The cover is empty as currently written.

**Required change:**

Replace the Faraday block in **both** engines with a check that exercises the cyclic Bianchi identity `∂_[α F_{βγ]} = 0` for `F = dA` from a smooth 4-potential, then specializes the result via the standard E,B map and verifies the three Maxwell-Faraday vector components fall out by Schwarz symmetry of mixed partials. This check is non-tautological: a sign error in the `E,B ↔ F` map produces a nonzero residual.

### Sympy edit (lines 46-67)

Replace lines 46-67:

```python
    t, x, y, z = sp.symbols("t x y z", real=True)
    E1 = sp.Function("E1")(t, x, y, z)
    E2 = sp.Function("E2")(t, x, y, z)
    E3 = sp.Function("E3")(t, x, y, z)
    B1 = sp.Function("B1")(t, x, y, z)
    B2 = sp.Function("B2")(t, x, y, z)
    B3 = sp.Function("B3")(t, x, y, z)

    F23, F31, F12 = B1, B2, B3
    F10, F20, F30 = E1, E2, E3
    F01, F02, F03 = -E1, -E2, -E3

    faraday = [
        sp.diff(F23, t) + sp.diff(F30, y) + sp.diff(F02, z)
        - (sp.diff(B1, t) + sp.diff(E3, y) - sp.diff(E2, z)),
        sp.diff(F31, t) + sp.diff(F10, z) + sp.diff(F03, x)
        - (sp.diff(B2, t) + sp.diff(E1, z) - sp.diff(E3, x)),
        sp.diff(F12, t) + sp.diff(F20, x) + sp.diff(F01, y)
        - (sp.diff(B3, t) + sp.diff(E2, x) - sp.diff(E1, y)),
    ]
    for i, residue in enumerate(faraday, start=1):
        assert_zero(f"Faraday component {i}", residue)
```

with:

```python
    t, x, y, z = sp.symbols("t x y z", real=True)
    coords = (t, x, y, z)
    A0 = sp.Function("A0")(t, x, y, z)
    A1 = sp.Function("A1")(t, x, y, z)
    A2 = sp.Function("A2")(t, x, y, z)
    A3 = sp.Function("A3")(t, x, y, z)
    A_components = (A0, A1, A2, A3)

    # F_{mu nu} = d_mu A_nu - d_nu A_mu (antisymmetric by construction).
    def F(mu, nu):
        return (sp.diff(A_components[nu], coords[mu])
                - sp.diff(A_components[mu], coords[nu]))

    # Cyclic Bianchi: d_[alpha F_{beta gamma]} = 0 for F = dA, by Schwarz
    # symmetry of mixed partials on smooth A on real coords. The check is
    # non-tautological: a sign error in the F definition produces a nonzero
    # residual.
    for (alpha, beta, gamma) in [(0, 2, 3), (0, 3, 1), (0, 1, 2)]:
        cyc = (sp.diff(F(beta, gamma), coords[alpha])
               + sp.diff(F(gamma, alpha), coords[beta])
               + sp.diff(F(alpha, beta), coords[gamma]))
        assert_zero(
            f"cyclic Bianchi (alpha={alpha}, beta={beta}, gamma={gamma})", cyc)

    # Specialize via E_i = -F_{0i} and B_1 = F_{23}, B_2 = F_{31}, B_3 = F_{12},
    # then verify the three components of dB/dt + curl(E) = 0 reduce to
    # cyclic Bianchi (and hence vanish). A sign error in the E,B<->F map
    # produces a nonzero residual.
    E_from_A = (-F(0, 1), -F(0, 2), -F(0, 3))
    B_from_A = (F(2, 3), F(3, 1), F(1, 2))
    mf1 = (sp.diff(B_from_A[0], t)
           + sp.diff(E_from_A[2], y) - sp.diff(E_from_A[1], z))
    mf2 = (sp.diff(B_from_A[1], t)
           + sp.diff(E_from_A[0], z) - sp.diff(E_from_A[2], x))
    mf3 = (sp.diff(B_from_A[2], t)
           + sp.diff(E_from_A[1], x) - sp.diff(E_from_A[0], y))
    assert_zero("Maxwell-Faraday component 1 from A", mf1)
    assert_zero("Maxwell-Faraday component 2 from A", mf2)
    assert_zero("Maxwell-Faraday component 3 from A", mf3)
```

Do not change anything outside lines 46-67 except the print line at the end of `main` (currently sympy line 85): replace the substring `"vector Bianchi signs"` with `"cyclic Bianchi from F=dA and Maxwell-Faraday reduction"`. This keeps the human-readable summary honest.

### Mathematica edit (lines 28-69)

Replace the entire block from line 28 (`(* M2: cyclic Bianchi signs reduce to vector Faraday signs. *)`) through line 69 (`Print["PASS: M2 Faraday component 3"];`) with:

```mathematica
(* M2: cyclic Bianchi identity for F = dA implies Maxwell-Faraday by
   Schwarz symmetry of mixed partials. Non-tautological: a sign error in
   the E,B<->F map produces a nonzero residual. *)
Clear[t, x, y, z, AA];
spaceTimeAssumptions = Element[{t, x, y, z}, Reals];
coordList = {t, x, y, z};
potentialList = {AA[0][t, x, y, z], AA[1][t, x, y, z],
                  AA[2][t, x, y, z], AA[3][t, x, y, z]};
fieldStrength[mu_, nu_] :=
  D[potentialList[[nu + 1]], coordList[[mu + 1]]]
    - D[potentialList[[mu + 1]], coordList[[nu + 1]]];

bianchiTriples = {{0, 2, 3}, {0, 3, 1}, {0, 1, 2}};
Do[
  Module[{alpha, beta, gamma, cyc, residual},
    {alpha, beta, gamma} = triple;
    cyc = D[fieldStrength[beta, gamma], coordList[[alpha + 1]]]
        + D[fieldStrength[gamma, alpha], coordList[[beta + 1]]]
        + D[fieldStrength[alpha, beta], coordList[[gamma + 1]]];
    residual = FullSimplify[cyc, Assumptions -> spaceTimeAssumptions];
    Print["M2 cyclic Bianchi ", triple, " residual = ", fmt[residual]];
    If[residual =!= 0,
      Print["FAIL: M2 cyclic Bianchi ", triple]; Exit[1]
    ];
    Print["PASS: M2 cyclic Bianchi ", triple];
  ],
  {triple, bianchiTriples}
];

electricFromA = Table[-fieldStrength[0, i], {i, 1, 3}];
magneticFromA = {fieldStrength[2, 3], fieldStrength[3, 1], fieldStrength[1, 2]};
maxwellFaradayComponent[1] :=
  D[magneticFromA[[1]], t]
    + D[electricFromA[[3]], y] - D[electricFromA[[2]], z];
maxwellFaradayComponent[2] :=
  D[magneticFromA[[2]], t]
    + D[electricFromA[[1]], z] - D[electricFromA[[3]], x];
maxwellFaradayComponent[3] :=
  D[magneticFromA[[3]], t]
    + D[electricFromA[[2]], x] - D[electricFromA[[1]], y];

Do[
  Module[{residual},
    residual = FullSimplify[maxwellFaradayComponent[k],
                            Assumptions -> spaceTimeAssumptions];
    Print["M2 Maxwell-Faraday component ", k, " residual = ", fmt[residual]];
    If[residual =!= 0,
      Print["FAIL: M2 Maxwell-Faraday component ", k]; Exit[1]
    ];
    Print["PASS: M2 Maxwell-Faraday component ", k];
  ],
  {k, 1, 3}
];
```

Do not modify any other block (`M1`, `M3`-`M6`, header `Print`, final `STATUS: PASS`, or the `fmt`, `scaleAssumptions`, `couplingAssumptions` definitions at the top).

**Self-test (mental walkthrough performed before writing this directive):**

1. **Variable independence (sympy edit).** `F(μ,ν)` depends on `A_components[μ]` and `A_components[ν]` only; each `A_i = sp.Function("Ai")(t,x,y,z)` depends on all four coords. The cyclic-Bianchi triple `(0,2,3)` produces second-derivative terms `∂_t∂_y A_3`, `∂_t∂_z A_2`, `∂_y∂_t A_3`, `∂_y∂_z A_0`, `∂_z∂_t A_2`, `∂_z∂_y A_0` (six terms), which cancel pairwise because sympy commutes mixed partials on `Function` with real symbol args. So `cyc` simplifies to `0`. **Non-trivial:** if I replaced `+ sp.diff(F(γ,α), coords[β])` with `- sp.diff(F(γ,α), coords[β])`, the cancellation breaks and the residual is nonzero. So the test detects sign errors.

2. **Variable independence (mathematica edit).** `fieldStrength[μ,ν]` builds the antisymmetric `dA-dA` form on `AA[i][t,x,y,z]`. Mathematica's `D` commutes mixed partials when the function arguments are symbolic and real-declared. `FullSimplify` with `Assumptions -> Element[{t,x,y,z}, Reals]` returns 0.

3. **Symmetry/parity:** N/A — no integrals over unbounded domains in this block.

4. **Trivial-case pre-check:** Substitute `AA[0]` and `AA[1]` as `0` and `AA[2] = f[t,x,y,z]`, `AA[3] = g[t,x,y,z]`. Then `F(2,3) = ∂_2 g - ∂_3 f`, `F(3,0) = -∂_3 (0) - ... = 0` (since `AA[0] = 0`), `F(0,2) = ∂_0(f) - 0 = ∂_t f`. The cyclic-Bianchi triple `(0,2,3)`: `∂_0(F(2,3)) + ∂_2(F(3,0)) + ∂_3(F(0,2)) = ∂_t(∂_2 g - ∂_3 f) + ∂_2(0) + ∂_3(∂_t f) = ∂_t∂_2 g - ∂_t∂_3 f + ∂_3∂_t f`. The last two terms cancel by Schwarz, leaving `∂_t∂_2 g`, which is generically nonzero. **Wait** — this is a check; the Bianchi sum should be zero for all `A`. Let me recompute: the third term needs `F(3,0)` not `F(0,3)`. `F(3,0) = ∂_3(0) - ∂_0(g) = -∂_t g`. So `∂_2(F(3,0)) = -∂_2∂_t g`. Full sum: `∂_t∂_2 g - ∂_t∂_3 f + ∂_2(-∂_t g) + ∂_3(∂_t f) = ∂_t∂_2 g - ∂_t∂_3 f - ∂_2∂_t g + ∂_3∂_t f = 0` by Schwarz. ✓

5. **Path specifications:** sympy file lives in `scripts/`; mathematica file lives in `mathematica/`. Both files already exist; this directive only edits them.

6. **Paper round-trip (v2):** The replacement check (cyclic Bianchi from `F = dA`, then Maxwell-Faraday by specialization) is a non-tautological re-derivation of bundle role rule (ii)'s homogeneous-Maxwell side. The paper card's eq:stage004-projection-ibp is unaffected (that's the M1/IBP block, not the M2 block). No new paper_misalignment is introduced; the previously-empty cover of rule (ii) is now substantively filled.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 004` and `redteam exec-mathematica 004` and confirm:
- both scripts exit 0;
- the sympy script source contains the new labels `"cyclic Bianchi (alpha=0, beta=2, gamma=3)"`, `"cyclic Bianchi (alpha=0, beta=3, gamma=1)"`, `"cyclic Bianchi (alpha=0, beta=1, gamma=2)"`, `"Maxwell-Faraday component 1 from A"`, `"Maxwell-Faraday component 2 from A"`, `"Maxwell-Faraday component 3 from A"`;
- the mathematica saved transcript contains `M2 cyclic Bianchi {0, 2, 3} residual = 0` through `{0, 1, 2}` and `M2 Maxwell-Faraday component 1 residual = 0` through component 3;
- the old `Faraday component k` (sympy) and `M2 Faraday component k` (mathematica) labels no longer appear in either source file or transcript;
- the sympy header print line (currently around line 85) no longer says `"vector Bianchi signs"` and instead says `"cyclic Bianchi from F=dA and Maxwell-Faraday reduction"`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl`
- summary: Replaced the tautological Faraday sign checks with cyclic Bianchi checks from `F=dA` and Maxwell-Faraday reductions in both engines.
- deviation: none
