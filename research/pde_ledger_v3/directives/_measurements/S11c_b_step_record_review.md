# S11c-b step record — document-review record (2 legs, convergent → folded)

Artifact `steps/S11c_b_variable_coefficient_operator.md`. Orchestrator-written prose ⇒ 2 document legs (Codex + Grok,
rule 7). Prompt `directives/_legs/S11c_b_step_record_review.md`; logs `~/.s11_build/S11c_b_step_record/{codex,grok}.log`.

**Both legs CONFIRMED the central close/defer split is HONEST** — the record does not claim the full residual ran; it
correctly carries the deferred items (full `row_residual`, P2a/P2b builds, #88 re-adjudication, the 2 hardenings, the
two sign conventions, #90's two flags), scopes the coarse consistency as coarse, and does not revive the superseded
step-2/4 family verdicts. All cited SHAs resolve.

**~16 correctable accuracy/provenance/completeness defects folded (v1→v2, one pass):**
1. **θ-row** (both): the current θ-row is `evolution_mass_balance − Σ closure_shape_deriv` (#90), NOT plain
   mass-evolution — that was the state AT `82f53828`; the plain description hid the #90 closure-fold flag object.
2. **`9f40c18e`** (both): it CLEARED #89, did not PRODUCE the engine (engine = checkpoint `f655ea65`; `git show
   --name-status 9f40c18e` = status/review files only). Relabeled.
3. **frozen-26 mechanism** (both): the common 26 was undercomplete by TWO mechanisms — PY froze the Hessian; WL
   hand-coded 8 with no quotient. Distinguished.
4. **#88 scope** (Codex): limited to the PY-side witness, LAB_HELD stored-energy rows only. Added.
5. **"operator rows"/"equations of motion"** not "Euler–Lagrange rows" (Codex): U/e_W are constraint-reduced.
6. **`u_L`** not `U` in the coupling sector `{θ,e_W,u_L}` (Codex).
7. **`Λ_A𝒜_s+Λ_VV_s / Λ_X𝒜_s`** not the rejected `A_T·Λ` shorthand (A_T collides with PY's solenoidal test
   potential); S11c-a exported the full T-a…T-i substrate, not "geometric boundary operators only" (Codex).
8. **Deferred recipe** must also name the fresh PY `S11CB_PRIMARIES_ONLY` `.out` + the out-of-band-verified #89 PY /
   #89b WL controls (Codex).
9. **#89b legs DISAGREED** (Codex): Grok caught the blocker, the other cleared it (insufficient presence test),
   orchestrator resolved computationally; only the REPAIRED artifact got 2 CLEAR re-legs. Narrative corrected.
10. **Rule-2 OOM** (Codex): only the production-alone OOM is a measurement; the concurrent double-kill is a
    diagnosis (STATUS), not a formal measurement. Reworded to match.
11. **P1-WL composition exception** (Grok): 2 fresh Claude agents (Grok's leg died on an xAI 500); no separate
    build-review `_measurements/` record — stated explicitly, not glossed as a standard 2-leg CLEAR.
12. **"spans verified"** (Grok): WL #89a = span-union identity; PY #89 = independent rank/count + form ablation.
    Distinguished in the ESTABLISHED bullet.
13. **Coarse consistency** (Grok): evidence is a SCRATCH file (not git / not `_measurements/`), single-case — named
    as such, so it does not fail the record's own rule-2 promise.
14. **Stale in-tree `.out`** (Grok): the committed WL `.out` (`d4adbd99`) is pre-#89b frozen, PY is pre-fold; a
    reader of them won't see the established objects. Stated in the close note.
15. **Exports** (Grok): `S11c_b_exports.py` regenerated (faithful — digests match the committed folded+#90 engine/spec)
    and later COMMITTED `af560257` (per user direction, satisfying N1's per-sub-step obligation; coupling
    per-engine-verified, cross-validation deferred) — the step record was noted here as "held" pre-commit;
    N1's per-sub-step exports obligation deferred with the residual. Stated.
16. **Deferral was a CHOICE** (Grok): a lighter core-only residual (~8 GB, fits) was OFFERED and NOT taken — the
    ≥64 GB deferral is the user's choice, not forced by physics. Stated.

⇒ v2 folded once. The close is honest per-engine, with the cross-engine residual explicitly owed.
