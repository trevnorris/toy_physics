# Directive review — S11c-a cross-engine comparator build directive

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_comparator_build_directive.md`

This is an orchestrator-written build directive for a cross-engine comparator that will difference the two
independent S11c-a engines' tag streams (SymPy `PY_S11CA_*` ↔ Wolfram `WL_S11CA_*`) and report, per shared
object, AGREE / DISAGREE / UNDECIDED / UNCOMPARED. It reuses the frozen "T7 contract" comparator
`research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py` and adds S11c-a-specific reconciliation
"bridges." Your job is to decide whether this directive is CORRECT, CLEAR, and — above all — SAFE, i.e.
whether it can be handed to a builder without any bridge silently masking a real cross-engine disagreement.

## What you are handed
- The directive (above).
- The frozen comparator precedent: `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`.
- The two REAL tag streams, to verify the directive's schema claims against reality:
  WL `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`;
  PY — regenerate if needed by running `python3 research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
  (⚠ slow; it rewrites `S11c_a_exports.py` in Dummy-index counters only — `git checkout` it after), or inspect
  the sympy source's emit calls. You do NOT need to run it if you can read the shapes from the source + the WL out.

## The decisive question — can ANY bridge turn a real VALUE difference into AGREE? (CLAUDE.md rule 6)
The whole point of two engines is that they can disagree; a comparator that reconciles too much destroys the
test. For EACH bridge in the directive's §"S11c-a delta" (retarget S11B_→S11CA_; outer case-key alignment;
inner field-name alias; held-integral head + bound var; nested-order derivatives; Dummy canonicalization;
window-function head rename), decide with a concrete argument or counterexample:
1. Is it purely syntactic / a fixed table, and **value-preserving** — does it rename only a symbol, head, key,
   or bound variable, and never move/drop/coerce a payload VALUE or map two structurally different objects
   together? Give a case where a defective (over-reaching) version WOULD mask a difference, and confirm the
   directive forbids it.
2. **Case-key alignment (delta #2):** is dimension-classification by disjoint token vocabularies actually
   unambiguous against the REAL keys (check FACE_PLUS/FACE_MINUS↔1/-1, DOF_ strip, branch, density)? Could a
   token be misclassified, or two distinct cases be force-paired? Does an unmatched case correctly become a
   surfaced KEY DISAGREE rather than being dropped or positionally paired?
3. **Field alias (delta #3):** is the two-entry table (VALUE→EXPRESSION, MULTIGRADE→MULTIGRADE_EPSILON_ETA_SIGMAW)
   correct and injective against the real field vocabularies? Are WL-only `EXACT_SOURCE` and PY-only
   `OPERAND_A/B` correctly left UNaliased and surfaced (not force-mapped, not silently dropped), while the
   matching fields still residual? ⚠ Does aliasing a KEY ever risk hiding that the two engines put a
   different STAGE of the object under `VALUE` vs `EXPRESSION` (e.g. graded result vs exact source)? If so, is
   that risk surfaced?
4. **Window function (delta #7):** the two engines parameterize the window differently — PY `O_window(G_+,G_-)`
   (2-arg) vs WL `windowFunction[single arg]`. The directive permits renaming the HEAD but ⛔ forbids bridging
   the ARITY/parameterization, so a genuine window difference surfaces as a projection DISAGREE. Is that the
   right call, and is the prohibition tight enough that no builder could "helpfully" introduce an argument map
   that reconciles them toward AGREE?
5. **Held integral / Dummy / nested derivative:** are these canonicalizations touching only head/bound-var/
   index (safe), with the integrand/order/name still residualled so a real difference DISAGREEs?

## Also check
- **Frozen contract preserved:** are the four verdicts + precedence (DISAGREE > UNCOMPARED > UNDECIDED >
  AGREE), native-boolean rejection at the leaf, UNDECIDED reserved for an authoritative status token, the
  no-terminal-verdict rule, and the exit-code policy (0 ran / 1 operational / 2 malformed) carried over intact
  and not weakened by any delta?
- **Fixtures (rule 5):** does the acceptance section require a SYNTHETIC fixture for EACH bridge that a
  defective (masking) comparator FAILS — in particular a "same aliased field, DIFFERING value → DISAGREE"
  fixture, an "unmatched case → KEY DISAGREE" fixture, and a "2-arg vs 1-arg window → structural DISAGREE"
  fixture — plus the carried-over repoint ablation? Is any acceptance criterion value-tuned to a real payload
  (forbidden)?
- **Rule-5 leak:** does the directive state any real computed VALUE or any expected real cross-engine result
  (which objects actually agree)? It must not.
- **Anything missing** that the builder needs to reconcile the real shapes, or any instruction that would let a
  real disagreement be masked or a real agreement be manufactured.

## Output
For each bridge: SAFE / UNSAFE (with the masking counterexample) / CLEAR-UP-NEEDED. Then: is the frozen
contract intact? Are the fixtures adequate and value-free? Any rule-5 leak? End with: safe to hand to the
builder as-is, or the exact edits needed. Quote directive line numbers and, where you checked a schema claim,
the `.out`/source evidence.
