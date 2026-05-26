# Codex — batch I.2 remediation (two math-decision items)

The verifier wave on batch I.2 surfaced two items that need your math authority. Trevor explicitly wants you (not Claude) deciding the math.

## Item 1 — Stage 021 F1 rework

**Status:** the Q6=a apply you did in the prior session ADDED a composed assertion, but the verifier flagged it as tautological. Concretely:

In `scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py` (and mirrored in mathematica), you added something like:
```python
expect_zero(
    "delta D_2^(odd) composed",
    Dcorr.subs(Gamma_port, a**5/(27*c_s**5)) - (-I * N0 * a**5/(27*c_s**5) * omega**5)
)
```

The problem: `Dcorr = -I * Gamma_port * omega^5 * N0`, so after `Gamma_port → a^5/(27 c_s^5)`, the LHS is exactly `-I * N0 * a^5/(27 c_s^5) * omega^5`. The assertion is `X - X = 0`, which holds trivially without exercising N0's closed form OR Γ_5^port's closed form.

**What the paper claims (eq:app-stage021-wall-odd):**
> δ D_2^{odd}(ω) = -i N_2(0) a^5/(27 c_s^5) ω^5 + O(ω^7)

where N_2(0) is the transfer factor at zero frequency, computed in Section III as the closed form for the reduced one-port:
```
N(0) = (Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² - R²)²    (from sympy line ~210 region)
```

**What you need to do:**

Replace the tautological assertion with one whose RHS uses the **explicit closed form** for both N(0) and Γ_5^port — not the symbol `N0`. The check should look like:

```python
expect_zero(
    "delta D_2^(odd) composed from Section III N(0) closed form and Section IV Gamma5 = a^5/(27 c_s^5)",
    Dcorr.subs(Gamma_port, a**5/(27*c_s**5))
    - (-I * ((OA**2 * gW + R * gA)**2 / (OA**2 * OW**2 - R**2)**2) * a**5/(27*c_s**5) * omega**5)
)
```

This non-tautologically requires N0 (Section III) to actually equal the rational expression, AND Γ_5^port (Section IV) to equal `a^5/(27 c_s^5)`. A future regression in either piece would break this.

Mirror in the Mathematica script using the local symbol names (`radius`, `cS`, `oA`, `oW`, `r`, `gA`, `gW`).

**Authorization:** edit `scripts/moving_throat_pde_stage021_*.py` and `mathematica/moving_throat_pde_stage021_*.wl`. Run both scripts to confirm exit 0. Iterate until clean.

After applying, append `## Applied: F1 (iter2)` to `redteam/directives/stage_021.md` (or note the rework in the existing `## Applied: F1` block if you'd rather just supersede the prior text).

---

## Item 2 — Stage 023 F2 sign decision

**Status:** the F2 directive prescribed adding a Schur-complement derivation block in both engines for the reduced mass-spring matrix. You correctly blocked F2 in iter1 because the directive's prescription has a sign conflict:

> The exact requested Schur block uses off-diagonal `+Rmix`/`+rMix`, whose adjugate gives a `-2*g_U*g_W*Rmix` cross term, but the existing `Q_expr`/`qExpr` rational numerator uses `+2*g_U*g_W*Rmix`; applying the block as written would introduce a failing assertion.

**What you need to do:**

Read the actual math source so the sign is determined by physics, not by the directive's possibly-typo'd prescription. Specifically:

1. Read `paper/stages/stage_023.tex` — what does the paper's §2 (or whichever section) say about the sign of the off-diagonal `R_mix` term in the U-W coupling Lagrangian? Cite the equation label.
2. Read `notes/stages/moving_throat_pde_stage023_*.md` if a derivation note exists — that's likely where the Schur-complement derivation was first written down with the convention.
3. Read the existing `Q_expr` / `qExpr` in `scripts/moving_throat_pde_stage023_*.py` and `mathematica/moving_throat_pde_stage023_*.wl` to confirm which sign convention the existing code uses.
4. Decide which sign is correct based on the paper/notes (the source of truth for the math), then either:
   - **(option α)** Apply the Schur block with the corrected sign so the new derivation matches existing `Q_expr` — if the paper supports `+Rmix` and the existing `Q_expr` is right.
   - **(option β)** Fix the existing `Q_expr` cross term sign AND apply the Schur block as originally prescribed — if the paper supports `-Rmix` and the existing `Q_expr` is wrong.

Whichever direction, you're the math authority — derive the sign from the Lagrangian and apply consistently.

**Authorization for Item 2:**
- READ: `paper/stages/stage_023.tex`, `notes/stages/moving_throat_pde_stage023_*.md`, `scripts/moving_throat_pde_stage023_*.py`, `mathematica/moving_throat_pde_stage023_*.wl`. Also read prior stages 021/022 if relevant for the coupling convention.
- EDIT: `scripts/moving_throat_pde_stage023_*.py`, `mathematica/moving_throat_pde_stage023_*.wl`. Plus append a `## Applied: F2` block under `redteam/directives/stage_023.md` documenting which option you chose and why.
- DO NOT edit: paper/tex (if the paper's sign needs changing, mark F2 still blocked with a one-line note for Trevor).
- Run both scripts to confirm exit 0 after editing. Iterate until clean.

---

## Both items

Mathematica single-seat: only one `math -script` at a time. Run sympy in parallel if you like, but mathematica serially.

When done, both stage 021 and stage 023 scripts (sympy + mathematica) should exit 0, and the two `## Applied:` blocks should appear in their respective directives.
