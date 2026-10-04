NEEDS REVISION FOR THIS SUPPORT AND CENTRE NONUNIFORM WORK METHOD

The plan is sound in structure, and I found no physics that has to be dropped. It needs one substantive addition about the background centre force, plus a clarification about the face-traction pairing.

**Substantive changes**

1. **The background centre force has no registry entry.** Section 2a lists only body U/θ/E_W and per-face tractions.
   - `admissibility_operator` (`S11c_b_brane_operator_sympy_audit.py:4057-4078`) and `admissibility_support` (`:4081-4098`) have no centre/ZETA_C component. The three registry entries therefore say nothing about a background centre force R0_c.
   - The plan pairs R0 with rates of f^(1) and f^(2), and 2b pairs the centre row with ζ_c,t. But 2a never says what R0_c is.
   - Add an explicit record for it. It is derived from the original traction normal components via the face geometry, or it is `UNSUPPLIED`. It must not be inferred as zero from the registry's silence.
   - R0_c·ζ_c,t at order ε¹ and ε² must stay a named term. This matters because the centre sector is never solved.

2. **Face traction is not independent of the body E_W force.**
   - In the original operator, face s has traction `(0,0,0, s*e_force/W0)` (`:4062-4074`). Its fourth component is the E_W body force divided by W0. In the support bundle, the matching `t_hold_*_0_1..4` symbols (four components, `:342-353`) are free and independent.
   - Operator-minus-support compared componentwise would therefore count E_W twice if the components are summed into work. The plan's "do not count twice" rule is correct but generic.
   - State this specific join. Also say which virtual rate the fourth component pairs with. In c2 that is the height velocity `(v − n·tangential)/n[3]`, with a tilted-normal area factor (`c2:957-963`).

**What checked out against source**

- **Legacy residual.** `object_difference` returns `sp.Equivalent(...)` whenever either operand is Boolean (`native-legacy-residual-definition.txt:12-13`). The residual view contains Equivalent/Not nodes (2 matching lines). Preserving it as provenance only and forming a new arithmetic difference is correct.
- **Support symbols.** `f_hold_u_a_0`, `f_hold_theta_0`, `f_hold_e_W_0` and `t_hold_*` are declared `"PREMISE"` symbols with no law attached. Treating them as `UNBOUND_SUPPORT`, never zero, is faithful.
- **Centre drive.** `CENTER_FACE_GENERALIZED_ROW` is built from physical DELTA_W with virtual ZETA_C (`:2154-2178`) and exported in `FACE_GENERALIZED_FORCE_ROWS` (`:3111-3114`).
- **Five-field restriction.** c2 reads only `U` and `E_W` from the face rows (`c2:964-965`) and `expanded_rows` takes only U/THETA/E_W. The centre row is silently dropped.
- **Selection is not a constraint.** `dof_fields` returns zero centre fields for DELTA_W, and `REPRESENTATION='DELTA_W'` is set at `c2:47`. The plan is right that this is a selection, not a closure.
- **Centre velocity.** The plan's outward form V_s = s·ζ_c,t + W0·eW_t/2 matches `dof_fields`: the lab velocity is `face*W0*eW_t/2`, and ζ_c enters unsigned.
- **Refusal rule.** "NO_SUPPLIED_CENTRE_CLOSURE, stop if the drive is nonzero or unknown" is a legitimate finite check. It invents no inertia, clamp or six-field inverse.
- **Order bookkeeping.** The ε versus defect-λ distinction (ε² field f^(2) versus λ² field f2) is clear.
- **Q_port.** It is kept formal and unevaluated. Reflection counts as survival, and the O(λ) weighted-remainder and G2 duties are explicit.
- **Background versus work.** Stationarity, incremental support work and mean power are kept distinct. The plan correctly avoids claiming that a static force times a periodic velocity is dissipative.
- **Missing inputs.** Honest outcomes such as `BACKGROUND_SUPPORT_UNRESOLVED` and `UNSUPPLIED` are allowed results. I found no real missing material or agency choice yet. The two changes above are plan edits, not user decisions.

**Optional wording**
- Say that "V_s candidate" is verified against the saved ZETA_C `face_velocity` table, not assumed.
- Say that the weak-form centre test should use the physical ZETA_C virtual-work table as a separate operand. The plan already says this; it could be restated in section 3.

**Reading coverage and uncertainty**
- Read in full: `plan.txt`, `guide.txt`, `review-prompt.md`, and the support, centre-field, five-row, face-selection, legacy-residual and publication views.
- Read in part:
  - `packet-index.json`, first 80 lines only.
  - The residual view, first 6 of 12 lines. The rest was truncated, so the Equivalent/Not presence rests on a count search.
  - The b source, only lines ~325-360, 2150-2250, 3090-3130 and 4020-4160.
  - The c2 source, only lines 955-985 plus line-count searches.
- Not read: the a source beyond searches, the exports, all saved numerical/current/receiving JSON, the `face_velocity`/traction/closure ZETA_C tables, hash pins, and the 8,000+ other packet files.
- I could not confirm the 4th-component pairing, the ZETA_C table contents or the a/b source pins, and the saved-operand joins are untested. Nothing was executed.