# S11c-b #84 — CAS verification of the bulk-core basis claim (rule 13; "no rederivation trusted until in CAS")

Codex (round 2) claimed the ~118-term bulk-core coupling residual is a GENUINE basis-content disagreement
(not a coefficient-renaming artifact), because both engines' §3a first-jet independence test "froze the
spurion during IBP", undercounting the §1d-mandated invariants. I verified this myself in CAS.

## What the reviewed comparator already localized (grounding)
`scripts/S11c_b_handcoded_comparison.py:205-213` (committed, leg-reviewed) states: WL's
`gamma{Width,Modulus}DivGrad{Theta,Ew}` multiply PY's OMITTED candidates 08/11, NOT PY's selected quotient
representatives 07/10; "A scale substitution or quotient identity would be physics-bearing, so all of these
names remain visible in the measured residuals." ⇒ the bulk-core residual = WL keeping reps {08,11} vs PY
keeping {07,10}; whether these are equivalent is exactly the §1d "physics-or-representational" question (#85).

## The mechanism, code-verified (PY)
`scripts/S11c_b_brane_operator_sympy_audit.py`: `enumerate_new_candidates(g_vector)` builds each first-jet
invariant as `spurion_i · (DOF bilinear)` with `spurion = g_vector` (the background first jet ∂W_bg).
`basis_euler_signatures` builds the Euler–Lagrange total-derivative `derivative_maps` ONLY from `basis_fields`
(the DOF fields u,θ,e_W + their jets). The background spurion `g_vector` is ABSENT from `derivative_maps`, so
`basis_dx` treats ∂W_bg as CONSTANT — it never generates the ∂²W_bg (Hessian) term. That is the UNIFORM
total-divergence quotient, which §1d says does NOT lift to variable coefficients.

## The probe (`/tmp/.../scratchpad/basis_rank_probe.py`; reconstructs PY's exact enumeration + EL machinery verbatim)
```
# candidates enumerated: 15   (FIRST_JET_CONTRACTION_01..15, per profile source)
(A) FROZEN spurion (PY's code):  rank=8   selected(1-based)=[1,4,5,6,7,9,10,13]  omitted=[2,3,8,11,12,14,15]
(B) HESSIAN-RETAINED (§1d):      rank=15  selected=all
=> PY selects 8 = [1,4,5,6,7,9,10,13] — MATCHES the emitted selection exactly (reconstruction is faithful).
```
- (A) reproduces PY's ACTUAL selected indices exactly ⇒ my reconstruction of PY's quotient is faithful, and
  PY provably uses the FROZEN-spurion (uniform) quotient.
- ⚠ (B) rank=15 is an OVER-COUNT of my probe, NOT the corrected §1d count: my probe treats the three spurion
  components `g_i` as independent symbols, so it misses the genuine null-Lagrangians `∂_i(W_bg·h_i)` with
  `∂_i h_i = 0` (which are g-linear divergences). The CORRECT variable-coefficient count is between 8 and 15
  (Codex estimated the two engines' selection-union at 10). Pinning it exactly = the #85 computation (needs a
  quotient expressed via W_bg explicitly, or g enforced as ∇W_bg).

## VERIFIED conclusion (direction solid; magnitude open)
1. PY's §3a first-jet independence test FREEZES the background spurion (uniform quotient) — code-verified and
   the frozen rank reproduces PY's exact 8-selection. §1d says this quotient is INVALID at variable coefficients
   (IBP of a variable coefficient generates first-jet/Hessian terms that are PHYSICS). ⇒ PY's §3a basis is
   UNDERCOMPLETE (drops invariants genuinely independent at variable coefficients).
2. ⇒ the ~118-term bulk-core coupling residual is a GENUINE §1d physics difference (PY reps {07,10} vs WL reps
   {08,11}, inequivalent at variable coefficients), NOT a benign coefficient relabeling. Codex's direction
   CONFIRMED; the exact "undercomplete by N" is not yet pinned (my probe over-counts).
3. Both engines kept DIFFERENT frozen 8-subsets ⇒ both are affected. WL's Mathematica basis reduction should be
   checked for the same frozen-spurion quotient (very likely, given it kept a different frozen representative).
4. ⛔ SCOPE ESCALATION: this is the deferred §1d energy-basis quotient (#85), now root-caused, and it is
   FOUNDATIONAL — the §3a basis feeds the §3b operator and EVERY family (slab operator, admissibility, coupling).
   The coupling bulk core cannot be cleanly closed without fixing the §3a basis quotient in BOTH engines.
   ⚠ Whether it reopens the already-adjudicated families (admissibility/kinetic/advective — which may have
   agreed only because both engines froze the SAME way) is an open risk to assess.

## Relation to the §3c CONTENT verdict (independent)
The INCLUDE/INCLUDE verdict (WL carries the correct face + response coupling; PY under-extracts BOTH) is a
SEPARATE axis from this bulk-core/basis issue: face/response = `A_T`/`Λ` content PY's §3c extraction drops;
bulk core = the frozen-spurion §3a representative difference. Both are real; they need different repairs.
