# Build review — S11c-b background-order MULTIGRADE INSTRUMENT (a SCRIPT)

You are independently reviewing a **Codex-written script**. Do not trust its author's self-report; derive
from first principles, run it, and **ablate its load-bearing checks**. Report a findings list; if it is
sound, say so and say exactly what you ablated and what the literal output was.

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_background_multigrade.py`

## What this instrument must do (the correctness claim to test)
For twenty aligned operand pairs the committed adjudication layer routed as genuine cross-engine differences
(twelve `FLAG`: ADVECTIVE×4, KINETIC×4, ADMISSIBILITY BODY_FORCE/THETA×4; eight `RESIDUAL_BULK`:
COUPLING_KERNEL×8), it must compute the **exact `(eta_bg, sigma_W)` background-bookkeeper multigrade** of the
SymPy operand (A), the Wolfram operand (B), and their residual (A−B), treating every other symbol as an atom.
The grade coefficient at `(a,b)` is the exact `∂ᵃ_{eta_bg} ∂ᵇ_{sigma_W} leaf |_{eta_bg=sigma_W=0} / (a! b!)`
over a fixed inclusive window `a,b ∈ [0,N]` (`N=4`), plus an **exact remainder** `R = leaf − Σ_window` for the
content beyond the window (the coupling operands are rational in `eta_bg` via `1/(1+eta_bg·w1_profile)`, so
the series is unbounded and the remainder must be exact, not dropped). It is a **measurement**: it must state
no conclusion about which engine is right, and it must re-derive **no physics** — it consumes the layer's own
aligned operands.

The correctness that matters downstream: (i) the grades are the TRUE grades of the actual aligned operands;
(ii) the residual graded is the SAME object the layer routed as `A_minus_B` (for KINETIC, via the layer's
`_kinetic_pairs`/`_kinetic_residual` adapter, with A/B/residual keyed by the SAME semantic labels); (iii) the
guards (`RECONSTRUCTION`, `WINDOW_CLEAN`, `GRADE_DIFFERENCE`, `REMAINDER_DIFFERENCE`) are genuine
self-consistency checks, not tautologies that pass for a wrong extraction.

## What you are handed
- The script above.
- The layer it must reuse: `.../scripts/S11c_b_adjudicated_comparison.py`
  (`transform`, `_bridge_d`, `_family_cases`, `_arithmetic_residual`, `_kinetic_pairs`, `_kinetic_residual`,
  `WL_TO_PY_RENAME`, `PROFILE_GRADE_SUBS`); comparator `.../scripts/S11c_b_cross_engine_comparator.py`
  (`load_py`, `load_wl`, `DEFAULT_PY`, `DEFAULT_WL`, `serialise`); engine
  `.../scripts/S11c_b_brane_operator_sympy_audit.py` (`eta_bg`, `sigma_W`, both `real=True`).
- The build directive (its intended contract): `.../directives/S11c_b_multigrade_instrument_build_directive.md`.
- The full committed transcripts it parses: `.../scripts/out/S11c_b_brane_operator_sympy_audit.out`,
  `.../mathematica/out/S11c_b_brane_operator_mathematica_audit.out`.
- A prior full adjudication run's stdout (ground-truth `operand_A`/`operand_B`/`A_minus_B` srepr for the 20
  cases): `/home/trevnorris/.s11_build/S11c_b_adjud_fullres_run.out`.

## Required method — DERIVE, RUN, ABLATE (a prose re-derivation is discarded)
1. **Independent re-derivation (mandatory, with a script + its literal stdout saved to named absolute paths).**
   Write your OWN small script that, for at least the ADVECTIVE and ADMISSIBILITY-THETA cases (whose operands
   are compact) and at least one COUPLING_KERNEL case (rational), pulls the aligned operand independently
   (you may reuse the layer's `transform`/`_bridge_d` — that is the operand source, not the thing under test)
   and extracts the `(eta_bg,sigma_W)` coefficients your own way (e.g. `sp.Poly` on the polynomial part,
   `sp.series`/limits on the rational part, or direct `.diff().subs()`), then compare YOUR coefficients to the
   instrument's emitted `MULTIGRADE_*`. Report the literal diff. ⛔ Do not verify the instrument against its
   own output or its own identities — derive the coefficients by a route the script does not use.
   Save your script and its stdout to named paths and report them; without them your derivation is discarded.
2. **Ablate every load-bearing check and report the LITERAL before/after output.** At minimum:
   - Corrupt the coefficient formula (e.g. drop the `1/(a!b!)`, or use `a→a+1`) on a copy and confirm the
     `RECONSTRUCTION`/`GRADE_DIFFERENCE` guards go NONZERO. If they stay zero, the guards are decorative.
   - **FORM ablation (mandatory — this is the only thing that catches the worst defect).** Change the
     STRUCTURE of a load-bearing object on a copy: swap the roles of `eta_bg` and `sigma_W` in the extraction,
     OR collapse both bookkeepers to one symbol, OR feed the wrong operand (B where A is expected) into one
     multigrade — then re-run and report the LITERAL diff in the guards and in a couple of grade cells. A
     COEFFICIENT rescale tests arithmetic; only a STRUCTURE change tests that the instrument measures what it
     claims. If the output is byte-identical under a structural change, that check is not implemented.
   - Ablate the remainder path: on a copy, set `R = leaf` with no coefficients extracted, and confirm the
     acceptance/`WINDOW_CLEAN` guard REJECTS it (a builder must not be able to hide all content in R).
   - For KINETIC: confirm A/B/residual multigrades are keyed by the SAME `_kinetic_pairs` labels — corrupt the
     alignment (key A by raw tuple index) and confirm `GRADE_DIFFERENCE` breaks. If mismatched keys silently
     become 0, the kinetic guard is vacuous.
3. **Which line computed this?** For every emitted object, point to the line that computed it; flag any tag
   that is a typed literal with no computation behind it. Report any `assert`/guard that PRECEDES the value it
   guards (a perturbation strong enough to flip it would crash before emission, hiding the defect).
4. **Reuse / no-re-derivation.** Confirm the instrument obtains operands from the layer (Bridge A + Bridge D
   via `transform`+`_bridge_d`, `bridge_d` applied once) and reconstructs no energy/operator/kernel of its own.
   Confirm the residual it grades equals the layer's routed `A_minus_B` (spot-check against the run output).
5. **Leak check (rule 5).** Confirm the script encodes no expected grade population, no "which engine is
   right", no "first/higher order" verdict, no coefficient target — the grades must be discovered by
   computation. Confirm the guards reference no expected physics value (they are decomposition self-checks).

## Physics filter
Report a finding only if it catches a way the instrument would MIS-MEASURE the multigrade, manufacture false
agreement (or false disagreement), silently drop content, or misalign the kinetic leaves — i.e. a way the
downstream background-order adjudication would be corrupted. Do not report style. Do not offer the physics
answer (which engine is right) — that is not this instrument's job and not yours here.

## Ablation sandbox and ops
- ⛔ Copy the script to your own scratch directory and ablate the COPY. ⛔ Never modify the working tree.
- Pure SymPy; **no Wolfram/Mathematica kernel** — no licence seat, no `timeout` needed.
- ⚠ The instrument RUN parses two large committed transcripts and does exact high-order differentiation on
  rational operands; a full run takes several minutes. **Run it in the FOREGROUND and WAIT for it. Do NOT
  background it, do NOT poll it in a monitor loop, do NOT yield before it finishes** — just block on the
  foreground command and read its output when it returns. Each ablation is another foreground run.
- If `N=4` exact differentiation on the 21KB rational coupling operands proves pathologically slow, say so and
  report the observed per-case time — that is itself a finding worth reporting (the observed polynomial
  degrees are small, so a large `N` may be unnecessary cost), but first confirm whether the run completes.

## Output
A numbered findings list — each: the defect, the literal ablation output that exposes it, the failure it
causes downstream, and the minimal fix — followed by a one-line verdict and the absolute paths of every
derivation/ablation script and stdout you saved.
