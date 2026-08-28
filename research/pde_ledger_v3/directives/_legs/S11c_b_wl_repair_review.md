# Independent physics review — S11c-b Wolfram engine, REPAIR round 1 (a SCRIPT, Mathematica)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
(repo root `/var/projects/toy_physics`). A Codex-repaired blind Wolfram engine. Its first build had three
defects (W1–W3) found by FORM ablation; those were just repaired against the corrected spec. **Verify the
repairs actually hold, and that nothing ablation-clean regressed.** No expected numeric answer (the spec
withholds every value).

## What you are handed
- The script above.
- Physics authority: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  (corrected §3a/§3c/§3d/§5a), and the specs it inherits (`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`).
- Prior finding record (what W1–W3 were):
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_b_wl_build_review.md`.
- ⛔ NOT the sibling SymPy engine's output — the cross-engine check is a separate downstream comparator.

## Required method — SCRIPT; derive independently and ABLATE (prose is discarded)
Write your OWN derivation for the load-bearing physics BEFORE studying the artifact; save it + literal stdout
to named /tmp paths and report them. Copy the `.wl` to /tmp and ablate the COPY; re-run only the affected
task; report the LITERAL diff.

**Verify each repair by FORM ablation:**
1. **W1 (BLOCKING) — admissibility is the background-order (ε⁰) balance, NOT vacuous.** Confirm
   `WL_S11CB_ADMISSIBILITY_OPERATOR_OPERAND` is the first variation of the FULL-field background energy at 𝔅⁰
   with the profile's gradients (`∇(W_bg+δW)`) retained — genuine data dependence on `W_bg`/its jets (including
   the generated second spatial derivative). ⛔ Not identically zero; ⛔ not `D[perturbation-energy(scaled
   fields), backgroundOrder]`. Ablate a background jet → the operand must MOVE.
   ⚠ **Also probe the uniform limit** (zero all background jets `η,σ_W`): a uniform unloaded slab must source
   ZERO background force — report any jet-INDEPENDENT survivor (e.g. a force ∝ moduli with no `W_bg` jet). The
   correct construction keeps the scalar perturbation fields (`θ,e_W`) as perturbations (∝ background order,
   vanishing at first variation) and full-fields only the gradient content; a scalar over-promotion
   (`θ→1+θ`) would leave a spurious jet-independent survivor.
2. **W2 — coupling kernel extracted via the §3c weak variational restriction** (solenoidal/irrotational
   trial+test), off-diagonal block → §1d decoupled zero in the uniform limit. Ablate: break the operator →
   kernel must break; zero the jets → the off-diagonal block must collapse (no diagonal-thickness or
   gradient-independent survivor in the emitted block).
3. **W3 — adjointness is the pairing-based operator-block residual**, ⛔ not the scalar-Hessian Clairaut
   (`D[D[E,·],·]−D[D[E,·],·]`). If adjoint by construction, it emits the two blocks / "no independent second
   route", not a structural zero.
4. **No regression:** the N15 energy-basis spurion construction, the operator, and the operator-extraction
   data dependence must still be ablation-clean.

Also report general defects: a guard that crashes instead of emits; a hand-typed payload; a tautological
residual; a VERDICT/PASS/FAIL or native/evaluated boolean residual; which line computes each §4 object.

## Physics filter
Report a finding only where the physics could be WRONG (still-vacuous or spurious-force admissibility, a
kernel not reducing to the uniform zero, a still-tautological control, a regression). If clean, state which
objects you ablated and the literal diffs.

## Operational constraints (IN THIS PROMPT — Mathematica)
- ⛔ Wrap EVERY kernel run in `timeout 600`; ⛔ never raise it. ⛔ Run at most ONE kernel at a time (two-seat
  licence; the sibling leg may hold one). ⛔ Kill any kernel past 6 GB RSS by PID and report.
- ⛔ Copy the `.wl` to /tmp and ablate the COPY. ⛔ Never modify the working tree or write under `mathematica/out/`.
- Save every ablation script + its literal stdout to named absolute paths and report them.

## Output
Findings (script quote + line, ablation + literal diff, why wrong, one-line repair), or the explicit
ablated-and-clean list with the repair-verification diffs. Read-only on the working tree.
