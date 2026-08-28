# Independent review — S11c-a step record (interpretation layer)

You are reviewing a **step record** — the prose interpretation layer for S11c-a. It makes physics claims
about a dual-CAS-engine derivation that will be read by later sub-steps (S11c-b..e). A prior draft of this
record **overclaimed** ("the two engines AGREE; every residual representational") on the back of a flawed
classifier that hid two real findings; this rewrite is supposed to correct that. Your job is to catch any
claim that overstates, misstates, or is not backed by the artifacts.

## Artifact under review
`research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md`

## Sources of truth — read these and form your own view BEFORE reading the record
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`.
- The two engines: `scripts/S11c_a_interface_geometry_sympy_audit.py`,
  `mathematica/S11c_a_interface_geometry_mathematica_audit.wl`.
- The final comparator run accounting: `~/.s11_build/comparator_final_cccb4f9e.out` — grep `RUN_ACCOUNTING`
  and every `^ACCOUNTING` line. This is the authoritative cross-engine state.
- Commit history: `git log --oneline -12` (the five fix shas the record cites: c36beac4, 6fae82b8,
  49b5c525, 8c1a5ed1, cccb4f9e).
- Measurements: `directives/_measurements/S11c_a_faceshift_nonjoin.md` (the FACE_SHIFT fix + the
  material-transport adjudication), `directives/_measurements/S11c_a_control_battery_result.md`,
  `directives/_measurements/S11c_a_T7_adjudication_verdicts.md`.

## What to check (quote the record and the artifact for each finding)
1. **Overclaim.** The record's central claim is "no genuine T7 physics disagreement survives … the 11
   families still showing unpaired cases are all pre-existing keying/schema asymmetries in the CONTROL and
   bookkeeping families." Verify against `RUN_ACCOUNTING` and the per-family `ACCOUNTING` lines: is
   `families_with_unpaired=11` right, are those 11 exactly the control/bookkeeping families the record
   names, and is it true none is a physics-bearing family with a genuine nonzero residual? If any
   physics-bearing family (geometry, projection, evolution, traction, FACE_SHIFT, uniform-limit) has
   unpaired cases or a nonzero residual the record calls clean, that is a finding.
2. **FACE_SHIFT resolution.** Does `FACE_SHIFT` show `join=160, py_only=0, wl_only=0, axis_set_mismatch=0`
   in the final run? Does the record's account of the fifth fix (free-premise `rhoBulkBackground` → grounded
   `rho4Profile`, all three sites) match the engine diff at `cccb4f9e`?
3. **The five fixes.** Are all five shas real and do they do what the record says (git show --stat each)?
   Is the substantive/free-premise split honest (which actually change a computed physics object vs compute
   to a forced value)?
4. **Material-transport adjudication.** Does the record faithfully represent the adjudication verdict (T-e
   correctly omits the transport; it lives in δh_s + T-f/T-g/T-h; not a shared blind spot)? Is the
   "EXACT_TRACE_SOURCE shows ρ(x) — changes no computed object" claim consistent with the measurements? Does
   it overstate the adjudication's certainty?
5. **Control battery.** The record retracts "every object bites" and says the ∂W_bg CONTROL_FORM ablation
   bites the geometric objects but not those riding on abstract symbols, and flags a post-fix change for
   FACE_SHIFT. Is that consistent with `S11c_a_control_battery_result.md` and the accounting? Any remaining
   overclaim about the controls?
6. **μ_θ / δρ / CONORMAL representational identities.** Are these still accurate after the #3/#4 current
   fixes and #5 density grounding? In particular, is the μ_θ "zero jet-suffixed / zero Derivative" probe
   still true in the current transcript, and does the δρ claim (PY primitive symbol = full §3a perturbation,
   EVOLUTION 8/8 to 0) hold?
7. **Any unbacked claim** — a number, a commit, a file path, or a physics statement the artifacts do not
   support (rule 2). Flag citations to measurement files that do not exist or do not say what the record
   says.

## Method
- Read the source of truth first; form your own view of the cross-engine state from the accounting before
  reading the record's summary of it.
- Where you assert a computed fact (a residual, a join count, a probe result), show the command and its
  literal output. A prose "I checked and it's fine" is discarded.
- Report a finding only if it changes what the record should say (an overclaim, a wrong number/commit, a
  physics misstatement, an unbacked claim). "Could be clearer" is not a finding.

## Return
A short list of record defects (or "accurate as written"), each with the quoted record text, the artifact
it contradicts (with command/output or file:line), and the correction.
