---
sweep_date: 2026-05-22
files_scanned: 196
files_flagged: 3
total_suspect_blocks: 4
risk_breakdown:
  high: 0
  medium: 0
  low: 4
---

# lRed-class continuation defect sweep (P4-32)

## Summary

Scanned all 196 `.wl` files under `mathematica/` (182 top-level audit scripts + 14 in `numerical/`) for the multi-line assignment continuation defect described in PAPER_CLEANUP_TRACKER item P4-32. The primary grep `^\s+[+\-]\s` returned 52 hits across 11 files; after manual context inspection, every hit is either (a) inside enclosing brackets `(...)`, `[...]`, or `{...}` that protect newlines from being treated as statement terminators, or (b) already patched (stage 003 with explicit re-add block, stage 021 with enclosing parens). The only residual cases worth flagging are four multi-line RHS assignments living inside `Module[...]` / `Function[...]` bodies (stage001, stage185_187_orbit_stress); these are conservatively low-risk because they sit inside outer brackets where Wolfram treats newlines as whitespace, but they match the user's explicit "flag Module/Block/With bodies with the same pattern" guidance. No top-level instance of the unprotected `name =\n  term\n  + term\n  - term;` pattern survives. Recommendation: P4-32 is closed for the audit corpus; the two known prior cases remain the only true defects.

## Reference defect

Defect signature: an assignment whose RHS spans multiple lines where (a) the assignment line ends with a complete term (no trailing binary operator), and (b) the next non-blank line begins with a unary `+` or `-`, and (c) the multi-line RHS is NOT enclosed in `(...)`, `[...]`, or `{...}` brackets. At top-level Wolfram scripts, the newline terminates the assignment after the first term, and the `+ ...` / `- ...` continuation lines become orphan top-level expressions that do nothing. Downstream `Coefficient[...]` and `expectZero[...]` checks against the truncated assignment can still pass, making the truncation silent.

Known prior cases (both already patched):
- `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl` lines 54-59 (patched in commit `10d4027` via an additive `lRed = lRed + ( ... )` re-add block at lines 62-67).
- `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl` lines 69-74 (patched by wrapping the multi-line RHS in `lRed = ( ... );`).

## Flagged files

### `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`

**Risk:** low

**Blocks:**
- lines 132-137: `lap` inside `lapEig[l_] := Module[{Ylm, lap}, ...]`. The body assigns
  ```
  lap = (1/Sin[theta]) D[Sin[theta] D[Ylm, theta], theta]
        + D[Ylm, {phi, 2}]/Sin[theta]^2;
  ```
  Line 134 ends with `theta]` (complete operand), line 135 begins with `+ D[Ylm, {phi, 2}]/Sin[theta]^2`. The Module body sits inside `Module[...]` brackets, so newlines are whitespace and the parser keeps reading; however the local form of the RHS matches the defect pattern verbatim. Low risk because of the enclosing `Module[`.

**Recommended fix:** wrap the RHS in parens — `lap = ( (1/Sin[theta]) D[Sin[theta] D[Ylm, theta], theta] + D[Ylm, {phi, 2}]/Sin[theta]^2 );` — to make the protection local and immune to refactors that move the body out of Module.

### `/var/projects/toy_physics/research/pde_ledger/mathematica/numerical/stage185_187_orbit_stress.wl`

**Risk:** low

**Blocks:**
- lines 116-121: `dmu` inside `fullTangentFromFree[ref_Association, free_] := Module[{...}, ...]`. The assignment is
  ```
  dmu =
    deta
    + 2 kW
    - 2 lam
    - ref["E"] (2 gamma + 2 lam - kU - kW)
    - ref["F"] ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);
  ```
  Line 117 holds the complete first term `deta`; lines 118-121 are pure `+`/`-` unary continuations. Protected by the enclosing `Module[`.
- lines 321-326: `dmuExpected` inside `stage187Check[...] := Module[..., Scan[Function[{item}, ...], ...]]`. The shape is identical to the block above:
  ```
  dmuExpected =
    detaExpected
    + 2 kW
    - 2 lam
    - ref["E"] (2 gamma + 2 lam - kU - kW)
    - ref["F"] ((1 + ref["delta"])/(1 + ref["chi"])) (gamma + ceta - kU);
  ```
  Protected by the enclosing `Function[...]`/`Module[...]` brackets.

**Recommended fix:** wrap each RHS in parens — `dmu = ( deta + 2 kW - 2 lam - ... );` — to make the protection independent of the outer Module/Function brackets. Both `dmu` and `dmuExpected` are numerically cross-checked against `muVal` via `nearQ`, so a silent truncation here would be caught at runtime, which further reduces risk; the wrap is a belt-and-suspenders hardening.

## Clean files

Files containing multi-line `+`/`-` continuation lines that were inspected and confirmed safe (all because the continuation lives inside `[...]`, `(...)`, or `{...}` brackets, OR the assignment uses trailing-binary-operator continuation, OR the assignment has been previously patched):
- `moving_throat_pde_stage003_bdg_mathematica_audit.wl` — patched (additive re-add + paren-wrapped Module bodies at lines 62-74, 86-90).
- `moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl` — continuations inside `FullSimplify[...]` brackets (lines 67-68, 85-86).
- `moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl` — `polyInScaledU =` (line 35) uses trailing-operator continuation (line 36 ends with `+`).
- `moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl` — patched with `lRed = (...)` paren wrap (lines 69-74).
- `moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl` — continuation inside `FullSimplify[...]` (line 173).
- `moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl` — continuation inside `FullSimplify[...]` (line 100).
- `moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl` — `deltaPerpExpected = Expand[ ... ]` (lines 124-132) — wrapped in `Expand[...]` brackets.
- `moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl` — `deltaPerp = Expand[ ... ]` (lines 39-47) — wrapped in `Expand[...]` brackets.
- `moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl` — `bCombFormula = Expand[ ... ]` (lines 54-57) — wrapped in `Expand[...]` brackets.
- `moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl` — continuation inside `FullSimplify[...]` (lines 87-92).
- `moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl` — six multi-line blocks (lines 148-157, 169-185, 196-205, 211-215, 216-221), all inside `normalizeExpr[...]` brackets or `(...)` parens within `expectZero[...]`.

All remaining 185 files contain no multi-line `+`/`-` unary continuation hits at all; they use single-line assignments, lists/associations `{...}`/`<|...|>`, or function-call brackets exclusively.

## Methodology notes

Primary first-pass grep:
```
grep -rEn '^[[:space:]]+[+\-][[:space:]]' mathematica/ --include='*.wl'
```
Returned 52 hits across 11 files. Each hit was inspected with `Read` on a ~5-25-line window around the candidate line to determine:
1. What the immediately enclosing bracket context is (top-level vs `(...)` vs `[...]` vs `{...}` vs `Module[...]`/`Block[...]`/`With[...]`/`Function[...]`).
2. Whether the assignment line itself ends with `=` (binary, forces continuation) or with a complete operand (defect-prone).
3. Whether the preceding non-blank line ends with a binary operator `+`/`-`/`*`/`/`/`,` (forces continuation) or with a complete operand.

Secondary grep for unary `*` and `/` continuations:
```
grep -rEn '^[[:space:]]+[*/]' mathematica/ --include='*.wl'
```
Returned only a comment terminator `*)` and one `/. theta -> ...` replacement-rule continuation inside `[...]`. Both safe.

Tertiary grep for assignments where the LHS line ends with `=`:
```
grep -rnP '[a-zA-Z][a-zA-Z0-9]*[[:space:]]*=[[:space:]]*[\n]' mathematica/ --include='*.wl'
```
Returned ~25 hits, all of which then either (a) had a trailing-binary-operator-style continuation, (b) had the first RHS term begin with a bracket-opening function call (`Element[...]`, `FullSimplify[...]`, `Expand[...]`, etc.), or (c) were already covered by the primary grep.

Edge cases skipped with rationale:
- Lines `{a,\n b,\n c}` or `<| key -> ..., \n key -> ... |>` — protected by `{...}`/`<|...|>` brackets, the parser never sees a statement boundary there.
- Block comments `(* ... *)` — handled by `grep`'s literal-line semantics; one match (`stage025:83 *)`) was a comment terminator and ignored.
- The `:=` definition heads (e.g. `f[x_] := body`) — the `:=` itself is a binary operator forcing continuation; the body, when a single `Module[...]` / `Block[...]` / `With[...]` call, is bracket-protected.
- Multi-line `expectZero["msg",\n  expr - (...)]` — inside `expectZero[...]` brackets, safe.

The sweep is exhaustive against the defect's grep signature (`^\s+[+\-]\s`); a false negative would require a defect form that does NOT start a continuation line with a unary `+`/`-` (e.g., a juxtaposition like `lRed = expr1\n  expr2` where `expr2` is interpreted as a separate top-level statement). I did not find any such pattern in the corpus during the cross-checks above, but note that grep-based detection cannot fully exclude that variant.
