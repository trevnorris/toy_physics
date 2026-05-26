You are applying the user-approved resolutions for batch III.2 v2 user-gate items. The user has explicitly authorized both directions (Q1=a, Q2=a) from the recommendations you wrote earlier in `redteam/resolutions/batch_III2_paper_alignment.md`.

# Per-question scope and authorization

## Q1 — Stage 050 paper-card edit: advertise the enhancement ceiling `S_n^(max)`

**Authorized to edit:**
- `paper/stages/stage_050.tex` — add a fifth boxed equation for the enhancement ceiling AFTER the existing four boxed equations and BEFORE `\stagefield{Output}`; update the Output line to mention the ceiling.

**NOT authorized to edit:**
- Scripts (already verify the ceiling correctly — sympy:88-107, wl:82-95)
- Notes (already own the ceiling in Section 5)
- Any other paper stage

**Specific edits:**

1. **Add a new boxed equation block** after the last boxed equation in the paper card body (just before the `\stagefield{Output}{...}` line at paper:44). The current boxed equations cover (in order): the lowest-twin saturation, the lowest-twin success condition, the higher-harmonic exclusion identity, and `x_max` (the `\stagefield{Output}` reference target). Insert the new block in a form matching the existing boxed-equation style:

   ```latex
   For \(n\ge 1\), the higher-harmonic enhancement is bracketed below by its support-ceiling value
   \begin{equation}
   \label{eq:app-stage050-Sn-max}
   \boxed{
   S_n^{\rm twin}(x;\varepsilon) < S_n^{\rm max}(\varepsilon)
   := 1 + \frac{1-\varepsilon}{(2n+1)^2-\varepsilon}}.
   \end{equation}
   ```

   Place this block such that it appears immediately before `\stagefield{Output}` (so it's the last math content of the body). Match the indentation/whitespace of the existing boxed equations.

2. **Update `\stagefield{Output}`** at paper:44. Current text:
   ```
   \stagefield{Output}{Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success} and higher-harmonic exclusion/softness thresholds.}
   ```
   Change to:
   ```
   \stagefield{Output}{Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success}, higher-harmonic exclusion/softness thresholds, and the enhancement ceiling \eqref{eq:app-stage050-Sn-max}.}
   ```

3. **No script change.** The SymPy script's `S_n_max` definition at sympy:90 and the Mathematica's `sNMax` at wl:83 already match the paper-side definition exactly.

4. **Verify** the new `\label{eq:app-stage050-Sn-max}` is unique by grepping for it across `paper/`. The label should not appear anywhere else.

**Destination cross-check (MANDATORY before editing):** confirm `paper/stages/stage_050.tex:11-44` does not already contain a `S_n^(max)` boxed equation under a different label name. Grep for `S_n.*max`, `S_n^.*max`, `Sn.*max`, `enhancement.*ceiling` in stage_050.tex. If found, this is a sign of partial prior work that needs investigation — record as `## Blocked: Q1` and halt instead of double-adding.

---

## Q2 — Stage 057 script edits: add local Pe-monotonicity check

**Authorized to edit:**
- `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py` — add a numerical Pe-monotonicity sweep after the existing `partial_kappa`/`partial_y` block (around line 73, before the constructive-branch closure ceiling block at line 75).
- `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl` — add the analogous block after the existing `D[zetaPhys, kappa]`/`D[zetaPhys, y]` block (around line 67, before the existing `zetaMax = ...` line at L69).

**NOT authorized to edit:**
- Stage 056 scripts (the user explicitly chose direction (a) and declined to also backfill Stage 056; that's deferred to a later batch)
- Paper TeX or notes files for Stage 057 (the paper card already states monotonicity; we're adding script-side verification, not paper edits)
- Any other stage's files

**Specific edits:**

**(a) SymPy:** add after the existing `expect_zero` block at L70-73, before the constructive-branch closure ceiling comment at L75:

```python
# Pe-monotonicity sweep — carry-forward from Stage 056 (notes §4: dOmega_Pe/dPe > 0 on
# the constructive branch via Cov_Pe(chi_0, s) > 0). Stage 056's scripts verify the
# covariance identity (sympy:79) but not the sign, so we anchor the sign locally here.
dPe = sp.simplify(sp.diff(zeta_phys, Pe))
for Pe_val in (sp.Rational(1, 10), sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(5), sp.Integer(10)):
    val = sp.simplify(dPe.subs({Pe: Pe_val, kappa: 1, y: sp.pi / 4}))
    if val <= 0:
        raise AssertionError(f"partial_Pe zeta sign failed at Pe={Pe_val}: {val}")
print("partial_Pe zeta > 0 on constructive branch (numerical sweep): PASS")
```

Use `sp.Rational(1, 10)`, `sp.Rational(1, 2)`, etc. rather than Python floats — the script already uses `sp.simplify` and rational arithmetic; mixing floats triggers `Limit::indet` style warnings.

**(b) Mathematica:** add after the existing `expectZero` block at L64-67, before the `zetaMax = ...` line at L69:

```mathematica
(* Pe-monotonicity sweep — carry-forward from Stage 056 (notes §4: dOmega_Pe/dPe > 0 on
   the constructive branch via Cov_Pe(chi_0, s) > 0). Stage 056's scripts verify the
   covariance identity (wl:65) but not the sign, so we anchor the sign locally here. *)
Module[{pevals, signOk},
  pevals = {1/10, 1/2, 1, 2, 5, 10};
  signOk = AllTrue[pevals,
    TrueQ[N[(D[zetaPhys, pe] /. {pe -> #, kappa -> 1, y -> Pi/4})] > 0] &];
  If[signOk,
    pass["partial_Pe zeta > 0 on constructive branch (numerical sweep)"],
    fail["partial_Pe zeta sign sweep"]
  ]
];
```

**Note on Mathematica:** the Mathematica script already declares `pe > 0` in `$Assumptions`, so the substitutions are sound. Use rationals (`1/10`, `1/2`) for consistency.

**Post-edit:** run both scripts:
- `python3 scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
- `math -script mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`

Confirm both exit 0 and that the new PASS line appears in each transcript.

**Mathematica single-seat:** orchestrator has NO other `math -script` running. You may run the math script as part of your apply session.

**Destination cross-check (MANDATORY):** before editing, grep both stage 057 scripts for `Pe` derivative usage to confirm no other `dPe`/`D[zetaPhys, pe]` block already exists in the script. If found, that's a sign of partial prior work — record `## Blocked: Q2` and halt.

---

# Output format

Append to `redteam/resolutions/batch_III2_paper_alignment.md` a new `## Apply log` section listing for each Q:

```
### Q<N> applied
- direction: <a|b|c|skip>
- files modified: <list with paths and line ranges>
- destination_verified: <yes — file:line | n/a>
- post-edit checks: <script exit codes, or "n/a — paper-only">
- notes: <anything surprising>
```

# Critical rules

1. **Per-question scope is strict.** Do NOT edit files outside the per-Q authorization list.
2. **Mathematica single-seat.** One `math -script` invocation at a time across the system.
3. **No JSON output.** YAML or markdown only.
4. **No fake commentary scripts.** Read and reason directly.
5. **For Q1 — duplicate-label cross-check (MANDATORY)**: grep `eq:app-stage050-Sn-max` across `paper/` before adding the new label. Record as `destination_verified: yes — no duplicate found`.
6. **For Q2 — pre-existing Pe-derivative cross-check (MANDATORY)**: grep both stage 057 scripts for any existing `Pe` derivative or sign check before adding the new block.
7. **Report blocked items explicitly** — don't silently skip.

# Working directory

`/var/projects/toy_physics/research/pde_ledger`
