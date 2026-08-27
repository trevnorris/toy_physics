# Independent review — S11c-a PY projection-current FIX DIRECTIVE (decision list, pre-build)

## Your task
Review an orchestrator-written fix directive **before** a builder acts on it. It instructs a change to one
computed object in the SymPy interface-geometry engine: the dynamic-window projection. A defect in this
decision list becomes a defect in a physics-bearing engine object, so find what is wrong, ambiguous, or
missing. Derive your own view from the spec and the code; a review that finds nothing is weak evidence.

## Read (source of truth first)
- The directive under review:
  `research/pde_ledger_v3/directives/S11c_a_py_projection_current_fix_directive.md`
- The spec it must honor: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — **§1b** (the current
  law `j=ρ_4D v_bulk`, conservation `∂_tρ_4D+∇₄·j=0`, "integrating this conservation law against Ω and
  integrating by parts in w"), **§3c** (the shifted-trace law and the sentence scoping zero-background to
  *traced* bulk fields at the face), **§4 T-f** (the `S11CA_PROJECTION_*` objects), **§6/§7** (script
  discipline; parallel tag grammar).
- The engine: `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` — read the current
  declarations (~lines 176–181), `projection_terms` (~1114–1171), and the trace construction that uses the
  normal jets (~line 651). Confirm for yourself whether `dw_delta_j_bulk` enters the projection.

## The questions that matter most
1. **Is the defect characterization accurate?** Does `projection_terms` in fact use bare `j_bulk[i]`
   constants and in-plane jets only, with the perturbation current's normal jets `dw_delta_j_bulk` absent
   from the projection (while present in the trace at line 651)? Quote the lines. If the code does something
   other than the directive claims, say so.
2. **Is the §1b obligation correct and does §3c genuinely NOT override it here?** Verify from the spec that
   the projected conservation law carries `∂_w δj_w` and that §3c's zero-background/trace language is scoped
   to traced face fields, not the bulk current under the `w`-integral. If the spec is actually ambiguous on
   this point, say so plainly (a one-engine fix is a spec question first).
3. **Is the correction well-scoped — names the object without dictating a wrong recipe?** The directive says
   to make the projection integrate the full `∇₄·δj` using the engine's existing `affine_bulk_perturbation` +
   `dw_delta_j_bulk` mechanism. Is that the right mechanism, or would applying it there double-count,
   mis-place a jet, or conflict with the IBP-in-`w` the projection already performs (the `window_normal`
   term)? Flag any way a literal application would be wrong.
4. **Do the guards protect the right things?** The directive forbids touching the trace construction, the
   (zero) background current, the window, other T-objects, and the `PROJECTION_*` tag set. Is anything
   load-bearing left unguarded, or any guard over-broad enough to block the actual fix?
5. **Does the directive leak an expected result the builder could iterate toward (rule 5)?** It must not hand
   a target residual, a "reduce to S11b" or cross-engine-match criterion, or the identity of which normal-jet
   terms survive. Flag any sentence that does.
6. **Anything ambiguous or missing** that would breed a build defect.

## Output
A ranked list of concrete defects (most severe first), each citing the directive line/section and the
spec/code evidence. If you judge an axis sound, say so and state what you checked. Point at artifacts and
cite them; a prose-only claim is discarded.
