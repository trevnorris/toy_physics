# Independent physics review — PY constraint-fold build (S11c-b slab operator, spec-pin B)

You are one of two independent build-review legs. The artifact is a Codex-written change to a SymPy engine. Derive
independently, ablate load-bearing objects on a COPY, and report literal stdout. A prose re-derivation is discarded:
write your own script and save both the script and its literal stdout to named absolute paths, and report those
paths. Report a finding only if it catches a way the physics could be WRONG (not "would be wrong on a different
input").

## Artifact

The working-tree change (vs `git HEAD`) to:
- `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (the SymPy engine)
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (§3b amendment paragraph)

Read the diff: `git -C /var/projects/toy_physics diff HEAD -- research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`

## What the build was supposed to do (the pinned physics — implement, not re-adjudicate)

The S11c-b slab U-momentum and thickness rows are the CONSTRAINT-REDUCED equations under S11b's binding material
virtual-displacement rule: the material virtual constraint `δ_vΣ_mat = 0` (the LIVE-`W_bg`, non-uniform object,
imported as the S11c-a `virtual_constraint` substrate — NOT the uniform three-term formula on
`S11c_b_SHARED_PHYSICS.md:143`) eliminates virtual `θ`, so the same held-fixed `μ_θ` feeds an in-plane reaction and
a thickness reaction. The `θ`-row is the imported sourced mass-evolution operand (`evolution_mass_balance`), not
`μ_θ=0`. `μ_θ = (δU/δθ)|held-fixed` stays a SEPARATE constitutive operand. Depths were cascaded (strong→3,
coupling→4, control→5). The build directive is
`research/pde_ledger_v3/directives/S11c_b_py_constraint_fold_directive.md` and the pin record is
`research/pde_ledger_v3/directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`. You are NOT asked
whether pin (B) is correct — that is settled.

## What you are handed

The engine, the S11c-a substrate engine `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
(source of the imported `virtual_constraint` and `evolution_mass_balance` operands), the S11b governing text
`research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:337,341-345,356-357,426`, and the spec
`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §1c/§3b. Copy the engine to `/tmp` and ablate the COPY;
never modify the working tree.

## Required probes — each is an ABLATION with literal stdout, not a code-read

⛔ A FORM ablation is MANDATORY (change the STRUCTURE of a load-bearing object, re-run, report the literal diff). A
coefficient rescale tests arithmetic; only a form change tests physics. For each, state WHICH LINE computes the
object.

1. **The reaction is COMPUTED from the imported constraint, not hand-typed.** This is the primary risk (a hand-typed
   CAS `Add` is still hand-typed). FORM-ablate the imported `virtual_constraint` SOURCE the fold consumes — e.g.
   replace the non-uniform imported constraint with the uniform three-term relation, or corrupt one of its
   coefficients on the COPY — and show the constrained `U`/`e_W` rows MOVE. If the reduced rows are byte-identical
   under a changed constraint, the reaction was frozen/typed (the defect this build exists to remove). Also: does the
   reduced `U`-row carry the non-uniform coefficient content the live constraint implies, or only a bare gradient of
   `μ_θ` (which would mean the uniform constraint was frozen in)?
2. **The two independent routes are genuinely independent.** The engine builds the reduced rows two ways
   (constraint-elimination and an independent Lagrange multiplier) and prints their residual. One-sided corruption:
   break ONE route only on the COPY and show the residual goes nonzero and only THAT route's object moves. If
   corrupting one route moves both, or the residual is zero by construction, it is operand theatre, not a check.
3. **#88 raw reference intact.** Confirm `operator_from_density` and `committed_strong_rows` still return the RAW
   held-fixed Euler–Lagrange rows (no constraint reaction, `THETA_BALANCE` still the held-fixed `μ_θ`-EL), so #88's
   standing measurement is not silently redefined. The fold must live in a SEPARATE layer.
4. **Depth cascade retains, not truncates.** Ablate a depth cap on the COPY (e.g. lower the coupling/control depth)
   and show a background jet the reaction generates is LOST — i.e. the raised depths are load-bearing, and nothing in
   the reduced `U`-row or the kernel is frozen below the order the reaction carries (rule 17).
5. **θ-row = imported mass-evolution, no double-count.** Confirm the emitted `θ` strong row is the imported
   `evolution_mass_balance` `(branch, DELTA_W, representative)` case and that `ADVECTIVE_MASS_OPERAND` is NOT added
   into it (that would double-count background advection). Confirm `μ_θ` (`MU_THETA_OPERATOR`/`MU_THETA_FACE_BINDING`)
   is unchanged and still held-fixed.
6. **Provenance + controls.** Is the reaction attributed under `BULK_ENERGY` (not a new top-level origin)? Does the
   per-term `energy_origins` loop avoid repeating the global mass row per energy term? Do `corrupt_material_constraint`
   and the §5b form ablations perturb the constraint/evolution SOURCE (not just zero an already-built operand)? Is the
   §5c uniform-limit comparison like-with-like (`U`/thickness vs the S11b `inplane_eom`/`thickness_eom`), and does
   `uniform_transverse_dispersion` stay free of an injected longitudinal force?
7. **No acceptance-value leak / no fix-until-match**, and every emitted payload is a CAS object, not a typed
   conclusion; report any `assert` that precedes the value it guards.

## Operational

Every kernel-free SymPy ablation is fine to run directly; there is no Mathematica seat here. Run a smoke build on a
single `(branch, representative)` case (e.g. `S11CB_PRIMARIES_ONLY=1`); ⛔ do NOT run the full package loop (it OOMs
~3.5h). Save every ablation script AND its literal stdout to named absolute paths and report them.

## Output

Per probe: CLEAR / FINDING (with the failure it catches), each backed by your script path + literal stdout and the
computing line. If everything is CLEAR, say what your ablations would have caught had it been wrong. Do not edit the
working tree.
