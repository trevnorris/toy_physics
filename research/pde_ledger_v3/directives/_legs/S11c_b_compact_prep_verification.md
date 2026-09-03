# Codex verification — S11c-b session state (compact-prep "are we in the clear")

Verify that the session's committed work and the durable records are CORRECT and internally consistent, against the
actual code. This is a sanity pass before a context compaction. Report any error, overstatement, or loose end. A
claim you cannot confirm against the code/diff is a FINDING.

## What was done this session (all committed + pushed to origin+GIN; HEAD `5f9cfe00`)
1. **PY constraint-fold (pin B)** — commit `82f53828`. The slab U/thickness rows are the constraint-reduced
   equations; θ-row = imported `evolution_mass_balance`; `μ_θ` a separate held-fixed operand; depth cascade 3/4/5.
2. **#90 §3c coupling content** — commit `7677aa18`. PY's coupling was bulk-only; now it computes the face
   generalized-force rows INTO the §3b operator (reversible tilted-face + irreversible response, Λ symbolic) and
   weak-restricts the full operator, plus a §0 clarity pin.

## Verify these claims against the code (spot-check; cite lines)
Read `directives/_measurements/S11c_b_{py_constraint_fold,90_coupling_content}_directive.md` (the two build-review
records), `directives/S11c_b_90_coupling_content_directive.md`, the top two sections of `STATUS.md`, and the diff
`git diff 82f53828~1 7677aa18 -- research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`. Confirm or refute:

1. **Constraint-fold:** the slab `U_BODY_BALANCE`/`E_W_BALANCE` `EXPANDED` rows are the constraint-reduced rows
   (raw held-fixed EL + the reaction from the IMPORTED non-uniform `virtual_constraint`, not the uniform three-term
   formula); `THETA_BALANCE` `EXPANDED` = the imported `evolution_mass_balance`; `MU_THETA` stays separate; the
   two-route residual is honestly relabeled `CONSTRAINT_FOLD_TRANSCRIPTION_RESIDUAL` (tautological, not independence).
2. **#90 coupling:** the emitted `S11CB_COUPLING_KERNEL` is the weak restriction of the FULL operator (bulk + face
   rows computed into the operator via the consumed virtual work), NOT a weak restriction of the raw
   `FACE_FLUX_BOUNDARY_OPERANDS` bundle (the §3c-forbidden parallel route); `build_kernel` does not read that bundle;
   live `μ_θ` is bound (`bind_mu_theta_operand`); reserved `mu_theta_L/M` do not survive in the blocks; `Λ_I(ω)` is
   symbolic (no `Z`/DtN/impedance/`.solve`); the skipped `closure_shape_deriv` (Λ_A/Λ_V) is now folded; `ζ_c` is not
   added; the dead `bulk_kernel_from_density`/`paired_kernel_from_density` are uncalled.
3. **#88 preserved:** `operator_from_density` and `committed_strong_rows` are byte-identical to their pre-fold form
   (the #88 raw-EL reference is intact) across BOTH commits.
4. **The two #90 flags are correctly characterized as CROSS-ENGINE / step-record items, NOT build defects:**
   (a) the closure-fold sign/magnitude (`mass_balance − Σ closure_shape_deriv`, unit coefficient) is only
   PY↔WL-verifiable; (b) the §5c uniform-limit residual is nonzero and Λ-bearing — §5c is a smoke test with no
   supplied value, so this is not a build gate; whether the response legitimately survives at uniform vs violates
   §1d decoupling is a cross-engine + step-record adjudication. Confirm the build legs did NOT gate on these.
5. **Terminology / plan:** the records correctly state that closing S11c-b produces the **S11c-b** step card (not an
   S11c-wide roll-up) and that S11c-c (curved-bulk closure/DtN) + S11c-d/e still remain (`directives/S11c_decisions.md`);
   the cross-engine integration pass is doable on this box (single-case WL operator ~0.9 GB, not ≥64 GB-gated), and
   ≥64 GB is off the critical path except the belt-and-suspenders full in-band control run.
6. **In the clear:** the engine parses; no acceptance-value leak was introduced; the records' commands match what
   they claim; nothing is overstated (e.g. no claim of cross-engine agreement that was not computed).

## Output
For each of 1–6: CONFIRMED (with the line/commit you checked) or FINDING (what is wrong + the fix). If a record
overstates or a claim is unverifiable, say so. Do not edit any file.
