# Independent build review — S11c-b WL admissibility repair (Codex-written)

You are one of two independent review legs on a **Codex-written** change to a blind Mathematica engine.
Derive independently; a prose re-derivation is worth nothing — save every derivation/ablation **script and
its literal stdout** to named absolute paths under `/tmp` and report them, or the claim is discarded.

Work in `/var/projects/toy_physics`; paths are relative to it.

## Artifact under review
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (working tree, **uncommitted**).
The change is 3 lines in `constructFullFieldBackgroundEnergy` (≈L528–577): the mixed thickness-gradient
invariant `Dot[gradTheta, gradLocalEw]` was changed to `Dot[gradTheta, gradFullEw]` with
`gradFullEw = anchoredWidth^(-1) gradient[fullWidth]` (was `gradLocalEw = gradient[localEw]`), where
`fullWidth = anchoredWidth + WZero ewVariation` and `localEw = (WZero/anchoredWidth) ewVariation`.
`git diff research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` shows it.

## What it claims to establish
The background-order admissibility operand `ADMISSIBILITY_OPERATOR_OPERAND` (built by
`backgroundBalanceFromModel`, ≈L1328–1351, from `constructFullFieldBackgroundEnergy`) is the §3d
background-order (ε⁰) first variation of the §3a full-field energy functional. §3d requires the thickness
gradient content be that of the **full** thickness field `∇(W_bg+δW)`. The change lifts the mixed
`∇θ·∇e_W` invariant's thickness-gradient factor to the full field, completing the treatment already given to
the pure-thickness invariant (invariant [[7]] on `fullWidth`).

## Sources of truth (read these; do NOT read the SymPy engine or any run output for an answer)
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`: §1c (L93–149), §2a/§2b (L172–235), §3a
  (L242–270), §3d (L325–356), §1d (L150–171).
- You are NOT given the sibling (SymPy) engine's θ body-force value, and must not derive toward it. Do not
  open `scripts/S11c_b_brane_operator_sympy_audit.py` or any committed run to lift an answer.

## Required checks

1. **FORM ablation (mandatory — the only test that has ever caught the worst defect).** Copy the `.wl` to
   `/tmp` and ablate the COPY (never the working tree). Change the **structure** of the repaired object —
   e.g. revert `gradFullEw` to the perturbation-only `gradient[localEw]`, or replace `fullWidth` inside it
   by `anchoredWidth` (background only) or by `WZero ewVariation` (perturbation only) — re-run, and report
   the **literal** change in the θ component of `ADMISSIBILITY_OPERATOR_OPERAND`. Confirm the θ body force
   **genuinely depends** on the full-field lift (a byte-identical θ output under the ablation would mean the
   term is not doing the work it claims). A COEFFICIENT rescale tests only arithmetic; the FORM change tests
   the physics. Report which line/emit computed the θ component.

2. **Independent derivation of the correct lift — including its normalization.** Derive from §3a/§3d the
   §3d-correct full-field form of the mixed `∇θ·∇e_W` invariant's thickness-gradient factor, and its
   background-order θ Euler covector, in SymPy (script + stdout). The implementation uses
   `gradFullEw = anchoredWidth^(-1) gradient[fullWidth]` (i.e. `∇(W_bg+δW)/W_bg`). Determine whether that
   normalization is the §3a/§3d-correct one, or whether the correct full-field representative carries a
   different normalization (e.g. `/W_0`) — these differ off the uniform limit by `W_0/W_bg` and change the
   θ body force's `η` dependence. State the mandated normalization and whether the implementation matches it.
   Confirm the uniform limit (`W_bg→W_0`, `∇W_bg→0`) still reduces to the S11b `∇θ·∇e_W` invariant.

3. **N15 untouched (double-count check).** Confirm the change does not lift the `N15` spurion invariants
   (`newInvariantExpressions`); those already carry a background jet, and lifting their `∇e_W` would inject a
   second background jet (second order in `σ_W`/`η`). Verify `newInvariantExpressions` is byte-identical.

4. **Scope byte-identity.** Confirm the following are byte-identical vs `HEAD`: both `kineticEw = muW
   WZero^2 …` (L838/L923), `constructEnergyData`, the §3b slab operator + origins, the §3c coupling kernel,
   `backgroundBalanceFromModel`. Confirm every emitted tag other than the admissibility operand (and its
   direct residual/embeddings) is unchanged. `git diff` is the reference.

5. **Blindness.** The engine still imports nothing and re-derives from the specs; no sibling import/read was
   added.

## Mathematica run discipline (both legs get this identically)
⛔ Wrap EVERY kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it and move on; NEVER
raise the timeout. ⛔ Run only ONE kernel at a time (the licence has TWO seats and another job may hold one).
⛔ Copy the artifact to `/tmp` and ablate the COPY; NEVER modify the working tree. Use
`--sandbox danger-full-access` semantics as needed. Save every ablation script AND its literal stdout to
named absolute `/tmp` paths and report them. Emit runs write only to `/tmp`, never `mathematica/out/`.

## Physics filter
Report a finding only if it catches a way the physics or method could be wrong — the θ body force not
depending on the lift, a wrong normalization, N15 double-counted, a scope regression, or blindness broken.

## Output
Lead with a one-line verdict: `SOUND` or `DEFECT(S) FOUND`. Then per check: the finding with the ablation/
derivation script + stdout paths. If a check is clean, say so. Do not modify the working tree or commit.
