# S11c-b row-residual instrument — build review (2 legs, Codex-written → fresh Claude agent + Grok)

Artifact: `scripts/S11c_b_row_residual.py` (Codex build, uncommitted). Legs (raw):
- fresh Claude agent (35-min baseline, all 4 FORM ablations, 216k tok): verdict **instrument SOUND**; one defect
  (nondeterministic certificate witness sign, verdict-preserving). Evidence under `/tmp/s11cb_review/`.
- Grok (`~/.s11_build/S11c_b_row_residual_buildreview_grok.txt`): 3 findings (F1 coupling residual, F2 mass
  activation, F3 certificate bridge convention). Evidence under `/tmp/s11c_b_row_residual_review/`.

## Both legs CONFIRM sound (independent ablation): truncation wired both directions (bound ≤0/≤2 move the
right rows; crosscheck flips false); assembly-accounting fires on a dropped operand (guard, not decorative);
ε order from family metadata (c=1 wave/coupling, c=0 admissibility) not symbol-counting; no `assert` before
any residual emit; no residual-value-conditional emission; both coupling cross-sector blocks + relabelled
adjoints emitted (not "±forward"); admissibility 52 componentwise (U×3/θ/e_W + PLUS/MINUS tractions, 4 cases).

## Leg DISAGREEMENT resolved (rule 1/6/13) — F2
Grok F2 said the mass/advective residual is BLIND to divergence activation (byte-identical under disabling the
instrument's `_expand_held_divergences`; mechanism = face subtraction `wl − face`). Fresh Claude's Ablation #2
(disabling the LAYER's `_held_divergence`) showed ALL 4 mass rows MOVE. The two ablated DIFFERENT functions.
RESOLUTION (rule-13, my reads): the mass row's `face_attributed = origins[FLUX] = projectedFaceFlux`
(`…row_residual.py:648-651`), and `projectedFaceFlux = Total[face["MEASURE"]·face["RELATIVE_FLUX"]]`
(`…_mathematica_audit.wl:782`) — a GENUINE boundary flux through the tilted faces, correctly excluded from the
bulk residual. WL advective (`divergence[σ_E u_t]`, active, L694) stays IN the residual. So Grok's ablation
moved only the correctly-EXCLUDED face-flux HeldDiv (residual correctly invariant); Claude's moved the bulk
content (residual correctly moves). ⇒ **Grok F2 = FALSE ALARM; the mass residual is the right object** (bulk
advective+accumulation compared, genuine face flux excluded). No instrument fix. (Validate at adjudication
that PY's mass row lacking ACCUMULATION — PY books ∂_tρ in the θ/density eq — is a representational split.)

## Confirmed real defects → one fix round
- **F3 (Grok; confirmed by code read, UNCONTESTED by Claude).** The coupling certificate is classified on the
  POST-`_bridge_d` residual with `apply_bridge_d=False` (`…row_residual.py:304, 406-409`), but the reviewed
  layer's weak route classifies the PRE-bridge anchored residual `_arithmetic_residual(left_pre_bridge_d,
  right_pre_bridge_d)` with `apply_bridge_d=True` (`…adjudicated_comparison.py:985-992`) — a different Euler
  object (profile-grade substitution can change Euler content before classification). The instrument must
  match the REVIEWED convention (established by the layer's v2 build legs). FIX: for coupling, form the
  anchored residual from the PRE-bridge operands with PY truncated to the requested truncation, classify with
  `apply_bridge_d=True`, and emit route + `anchored_remainder` + `euler_signature` + `V_ANCHORED` as the layer
  does.
- **F1 (Grok=defect, Claude=acceptable-by-design) — resolved by the F3 fix.** Both legs agree the FACT: the
  coupling `ROW_RESIDUAL = _exact_scalar_residual(py_trunc, wl_complete)` (the pre-quotient exact diff), with
  the modulo-divergence verdict living in the certificate. Claude: acceptable (certificate is a genuine
  second-route object). Grok: should be the quotient representative. The F3 fix makes `ROW_RESIDUAL` the
  layer's `anchored_remainder` (R − div V) = the quotient representative, satisfying both.
- **Claude Finding 1 (nondeterministic certificate witness sign, verdict-preserving).** For a byte-identical
  input residual, the emitted `euler_signature` flips overall sign between processes (cause: `sp.cancel` in the
  layer's `_normalise_exact` L617-631 is hash-ordering-sensitive). ROUTE (RESIDUAL_BULK vs REPRESENTATIONAL) is
  stable, so no verdict impact — but the committed `.out` is not byte-reproducible and a hand-check of the
  witness sign can spuriously mismatch (it also produced Grok's "phantom" 6/20 sign-only cert diffs). FIX:
  canonicalize the emitted witness sign deterministically in the instrument (e.g. fix the sign of a chosen
  leading monomial) OR pin `PYTHONHASHSEED` for the committed run. Prefer the in-code canonicalization
  (self-contained, robust). ⛔ Do NOT modify the committed/reviewed layer.

## Physics observation for the step record (Claude; NOT an instrument fault, VERIFY after the F3 fix)
Under the current (pre-fix) convention all 20 coupling certificates classify RESIDUAL_BULK (incl. the 4
ADJOINTNESS_RESIDUAL), i.e. genuine bulk not total divergence. ⚠ This used the WRONG convention (F3) — RE-READ
after the fix under the layer's pre-bridge/apply_bridge_d=True convention before trusting it for the coupling
adjudication.

## Decision
One focused Codex fix (coupling weak-route convention + deterministic witness), scoped 2 build legs on the
CHANGED coupling handling + determinism (the rest was already legged sound, unchanged), then commit → run →
adjudicate. First build-leg round found real defects → fix once (not a rule-15 repeat-fail).

## FIX-DIRECTIVE decision-list review (2 legs on the fix directive; both REJECTED v1/v2 draft, CONVERGED on the construction)
Legs: Codex `~/.s11_build/S11c_b_fix_dirreview_codex.txt` (7 defects), Grok `~/.s11_build/S11c_b_fix_dirreview_grok.txt`
(constructive). Both SOUND: (a) Euler must run PRE-bridge (independent jets); (b) truncation only acts POST-bridge;
(c) my "classify pre-bridge, then bridge+truncate the Euler SIGNATURE" is INCOHERENT — `_bridge_d` does NOT commute
with `formal_divergence`/`formal_dx` (layer fixture L1304-1308; Grok direct check: `div(W_bg u_T)` is
REPRESENTATIONAL pre-bridge but post-bridge/apply_bridge_d=False gives false RESIDUAL_BULK Euler `−W_0 η ∂w1 + σ_W ∂w1`
— exactly the current instrument bug). Also: (Codex D2/D3/D5) the Euler signature is a field-keyed obstruction,
early-break INCOMPLETE (L664), and absent on RESIDUAL_BULK (L735-739 short-circuit, no V/remainder) & on homotopy-throw
DIVERGENCE_INCOMPLETE — it is NOT a residual representative; (Codex D6/§382) ROW_RESIDUAL must stay the A−B density;
(Codex D4/Grok) full-order route ≠ in-scope verdict EITHER direction (`R=div(W_bg u_T)+(W_bg−W_0)²u_T`: route=RESIDUAL_BULK
full, but in-scope content is the divergence IMAGE (nonzero ησ_W) while genuine bulk is η² truncated away → a naive
`truncate(bridge(R))` FALSELY reads in-scope bulk); (Codex D7/Grok Q5) sign-flip determinism unsafe on symbolic coeffs
& must not touch ROW_RESIDUAL/remainder — use deterministic SERIALIZATION on the witness copy only.
SOUND (both): no rule-5 leak; strong/mass/admissibility correctly kept out of the quotient; layer immutability.

⇒ PINNED CONSTRUCTION (both legs' convergent correction; I rule-13-cross-checked both reports):
  R_pre := A._arithmetic_residual(left_pre_bridge_d, right_pre_bridge_d)          # pre-bridge A−B (layer L985-992 operands)
  FULL_PREBRIDGE_ROUTE, euler_signature := A.classify_total_divergence(R_pre, apply_bridge_d=True)   # reviewed full-order diagnostic
  ROW_RESIDUAL := requested_truncation(A._bridge_d(R_pre))                         # in-scope A−B DENSITY (§382; Codex D6)
  IN_SCOPE_WEAK_REMAINDER := requested_truncation(A._bridge_d(_normalise_exact(R_pre − A.formal_divergence(A._homotopy_vector(R_pre, REG)))))
      # the modulo-total-in-plane-divergence representative, computed IN THE INSTRUMENT even when the layer short-circuits
      # (RESIDUAL_BULK) — Grok's fix; on `_homotopy_vector` throw emit a NO_CLEAN_QUOTIENT flag + the raw ROW_RESIDUAL.
  determinism := deterministic serialization of the emitted euler_signature/witness COPY only (leading-monomial sign /
      pinned signsimp), AFTER route+residual computed; never applied to ROW_RESIDUAL / remainder.
  Verdict reading (off-instrument, adjudication): ROW_RESIDUAL=0 ⇒ engines agree in-scope; ROW_RESIDUAL≠0 ⇒ then
  IN_SCOPE_WEAK_REMAINDER=0 ⇒ representational (in-scope difference is a total in-plane divergence);
  IN_SCOPE_WEAK_REMAINDER≠0 ⇒ genuine in-scope bulk difference.
Strong slab rows / mass / admissibility: UNCHANGED (both legs sound). This is the legs' verified construction
transcribed; the 2 BUILD legs on the fix gate the IMPLEMENTATION (no third decision-list round — not iterating to green).

## FIX BUILD legs (2 legs on the coupling fix; both CLEAN — no defects)
- Grok `~/.s11_build/S11c_b_fix_buildreview_grok.txt`: no findings. Ablations bite (pre→post bridge Euler moves;
  homotopy peel: synthetic div(W_bg u_T)+(W_bg−W_0)²u_T peels to 0 with-peel, nonzero without; 20/20 real
  ROW_RESIDUAL≠IN_SCOPE_WEAK_REMAINDER; 2 runs byte-identical incl. all 20 witnesses). Layer git diff = 0 bytes.
- fresh Claude agent (162k tok): no defects. Independent non-commutation check reproduced the defect; FIX witness
  SHA present in delivered output & ABL SHA absent (proves pre-bridge convention); homotopy on pre-bridge,
  truncation on bridged (L849-893); canonicalization scoped to euler_signature only (raw residual/remainder
  byte-identical across seeds); NON_COUPLING byte-identical via independent stripper (sha match); NO_CLEAN_QUOTIENT
  reachable (demonstrated a throw) but false on all 20 production cases; no asserts, no hardcoded verdicts; layer
  unchanged. Caveat: full run ~2h so used faithful harnesses + the builder's authenticated preserved full outputs.
⇒ instrument REVIEWED-SOUND (coupling fix correct, load-bearing, deterministic, no regression, layer immutable).
COMPUTED (both legs, for adjudication): FULL_PREBRIDGE_ROUTE=RESIDUAL_BULK ×20 (genuine full-order bulk under the
correct convention); ROW_RESIDUAL≠IN_SCOPE_WEAK_REMAINDER ×20; NO_CLEAN_QUOTIENT=false ×20. The IN-SCOPE verdict
(is IN_SCOPE_WEAK_REMAINDER zero?) is READ during adjudication, not by the legs (rule 5).
