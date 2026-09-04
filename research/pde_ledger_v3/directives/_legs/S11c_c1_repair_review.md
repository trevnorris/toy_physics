# Independent review — the S11c-c1 engine REPAIR (Codex-written; round 2)

## Context
The first c1 engine build was 2-leg-reviewed: the core physics (two-momentum DtN kernel, operator-inverse face
response, 44-row own-rows delta) was confirmed SOUND, but five self-checks could not fail (were toothless). A
repair directive (`directives/S11c_c1_repair_directive.md`) folded those five findings. Your job: verify each of
the five controls now GENUINELY BITES (a real error moves its residual), the new constructions are physically
correct, and the repair introduced NO regression to the confirmed-sound core.

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py` (repaired; runs to
exit 0, regenerates `S11c_c1_exports.py`, 63 tags). The baseline before repair is git `65afa1cd`; the diff added
`outgoing_farfield_poynting`, `hanzawa_first_kernel`, `reduce_on_shell`, `intended_own_bind_closure`,
`simplified_object_difference` and changed `port_matrix`/`port_hermitian`/`build_model`/`build_delta`. Confirmed
already: `dtn_first_kernel`, `closed_coefficients`, `IMPORT_KEYS`, and the 44 exported rows are byte-unchanged.

## What you are handed
- The repaired engine and the baseline (`git show 65afa1cd:…`, `git diff 65afa1cd -- …`).
- The PHYSICS AUTHORITY: `directives/S11c_c1_SHARED_PHYSICS.md` (§3a, §3b obj. 1–3, §5a, §5d) — derive from this.
- The prior review record `directives/_measurements/S11c_c1_build_review.md` (the 5 findings), the repair
  directive, `scripts/ledger_fold.py`, the frozen base `scripts/S11c_b_exports.py`.

## Required method — this is a SCRIPT
Derive independently; write your own derivation/ablation scripts and save them + literal stdout to named /tmp
paths (a prose "it now bites" is discarded). Ablate load-bearing checks on a /tmp COPY, never the working tree.
The engine runs in ~9 min; capture stdout to an absolute path and grep the tags you need.

## The five repairs — for EACH, prove it now bites (the check must FAIL on a real error)
1. **R1 energy route (§3b obj. 3).** The bulk operand must now be the outgoing far-field Poynting flux from φ at
   `|w|→∞` (`outgoing_farfield_poynting`), NOT `δp` at the face. **Ablate:** (a) corrupt the `q_out` branch sign
   in the bulk/φ route only → the residual `ENERGY_RESIDUAL` must go NONZERO (a branch error is now caught); (b)
   corrupt the `t_s` sign only → only the FACE operand moves. Confirm the two operands do NOT share a common
   `z_energy` factor (the old defect). Confirm the far-field flux is the physically correct outgoing acoustic
   intensity of the radiated wave — derive it yourself and compare.
2. **R2 rep-invariance §5a (`hanzawa_first_kernel`).** EULERIAN and HANZAWA must now be STRUCTURALLY DIFFERENT
   constructions of the same DtN (a real coordinate change / layer potential with its own outgoing kernel), not
   the same rational function. **Ablate:** one-sided-corrupt the Eulerian route only → `REP_INVARIANCE_RESIDUAL`
   must go NONZERO. Confirm `EULERIAN − HANZAWA` genuinely simplifies to 0 for the CORRECT kernel (the routes
   agree) but the residual is not a structural `A−A` (they are two constructions). Verify the Hanzawa route
   preserves the §1b radiation branch (⛔ not the secular global scaling `w′=[w−ζ_c]/[W_bg+δW]`).
3. **R3 on-shell reduction (`reduce_on_shell`).** `DTN_RIGID_SHIFT_RESIDUAL` must now be emitted ON-shell (the
   rigid-shift cancellation actually demonstrated, not an off-shell expression), and the §5d `ZERO_JET` operand
   must no longer carry the spurious off-shell term. **Check:** grep the real run for `DTN_RIGID_SHIFT_RESIDUAL`
   and `ZERO_JET`; confirm the on-shell reduction fires (is representation-robust — `sp.factor` no longer breaks
   it). Confirm it still cancels for the correct kernel and would be nonzero for a wrong one.
4. **R4 Hermitian/port fed Z_0+Z_1.** `DTN_HERMITIAN_PART` and the port Hermitian form must now include the flat
   `Z_0`, so at `(η,σ_W)=0` they reduce to `H_a[Z_0]` (the leading bulk radiation) and the leading closed-port
   map — NOT vanish. **Check:** evaluate at `(η,σ_W)=0`; confirm the Hermitian part is the nonzero `Re(ρ_m ω/q)`
   propagating-subspace object, not 0. Confirm the sign-definiteness restriction + `NOT_ESTABLISHED_…` token on
   the zeroth-order nullspace, and no sign stated in prose.
5. **R5 minimum-mode (`intended_own_bind_closure`).** `assert_delta_is_minimal` must now check the bind-closure
   of the five INTENDED export roots (derived from the roots, not read back from the delta). **Ablate:** inject a
   stray inherited row into the delta on a /tmp copy → the check must RAISE. Confirm `own_closure` is NOT the
   delta's own key-set.

## No-regression checks
- The core is byte-unchanged, but confirm the REPAIRED run still emits the sound objects: the two-momentum DtN
  kernel with both legs live (form ablation `q_input=q_output` still MOVES it), the operator-inverse response,
  `Λ_X`-only-in-traction, opaque `μ_θ`.
- The delta is still the 44-row own-rows delta; `check_consumer` + the lookup smoke-test still bite; the 5
  model-level rows are byte-identical to `65afa1cd` (only the self-digest changed).
- ⛔ No NEW toothless check, no `A−A`, no hand-typed CAS payload with no data dependency, no conclusion emitted
  as prose, no tag name encoding a sign/parity/regime.

## Physics filter & sandbox
Report a finding only if a control still cannot fail, a new construction is physically wrong, the core regressed,
or a new defect was introduced. Copy anything you execute-ablate to /tmp; never modify the working tree. Pure
SymPy, no Mathematica.
