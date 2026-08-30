# Independent decision-list review — S11c-b STEP 3 WL-repair directive (v2, one item)

You are one of two independent review legs on an **orchestrator-written** directive. **No builder runs until
this review clears.** A v1 of this directive had two items; two independent decision-list derivations withdrew
the second (thickness kinetic: WL's `μ_W W_0²` is spec-correct because §1a defines `e_W ≡ δW/W_0` with
constant `W_0`) and re-diagnosed the first. This v2 has **one** item. Your job is to verify v2 is correct and
buildable — do not rubber-stamp; a wrong or over-specified directive here corrupts a blind engine.

Work in `/var/projects/toy_physics`; paths are relative to it.

## Artifacts
- Directive under review: `research/pde_ledger_v3/directives/S11c_b_wl_admissibility_repair_directive.md` (v2).
- Its measurements: `research/pde_ledger_v3/directives/_measurements/S11c_b_wl_admissibility_repair_directive.md`.
- Shared physics (source of truth): `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  (§1c L93–149, §2a/§2b L172–235, §3a L242–270, §3b L272–288, §3d L325–356, §1d L150–171).
- The WL engine the directive modifies (you MAY read it — it is the target):
  `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`.

You are NOT handed the sibling (SymPy) engine's values. Derive from the shared physics; do not lift an answer
from `scripts/S11c_b_brane_operator_sympy_audit.py` or any run output. **Save any derivation as a
SymPy/Python script + its literal stdout to named `/tmp` paths and report them** — a prose derivation is
discarded. Prefer Python (avoid the two-seat Mathematica licence); if you must use Mathematica, `timeout 600`,
one kernel.

## What to verify

1. **The named object is the §3d-correct one (independently).** v2 says the background-order admissibility
   operand requires the §3d full thickness-field gradient `∇(W_bg+δW)` in **every** §3a invariant carrying a
   thickness gradient, and that WL applies this to the pure-thickness invariant (`fullWidth`, L554) but not to
   the mixed `∇θ·∇e_W` invariant (`gradLocalEw`, L555/L543). Derive from §3a/§3d whether that is the correct
   object and the correct locus of the under-retention — i.e., that lifting the mixed invariant's thickness
   factor to `∇(W_bg+δW)` yields the §3d-mandated background-order variation, and that the perturbation-only
   `gradLocalEw` does not. Confirm or refute with CAS.

2. **Completeness of the requirement.** v2 names the requirement at the level of "every invariant carrying a
   thickness gradient." Enumerate the §3a invariants (uniform list + the `N15` new invariants) and confirm
   which carry a thickness gradient, so the uniform requirement covers them all and none is silently missed
   (e.g., are any `N15` spurion invariants also affected?). Report any invariant the requirement should cover
   that the directive's framing might leave ambiguous.

3. **The scope claim (the key new assertion).** v2 asserts the fix belongs to the background-order
   admissibility functional (`constructFullFieldBackgroundEnergy`) **only**, and that the §3b wave operator
   (`constructEnergyData`) and §3c coupling kernel must stay **byte-identical** because the wave operator's
   mixed invariant is legitimately perturbation-only (the operator is first order in the perturbations).
   Verify this: is the full-field thickness lift genuinely inappropriate for the wave operator, or would a
   faithful repair have to touch it too? A wrong scope guard here is a regression. Argue from §3b/§3d.

4. **Leak (rule 5) and over-specification (rule 3).** Confirm v2 states no value, coefficient, sign, order, or
   residual for the repaired object, and does not say whether it is currently zero or nonzero. Confirm it
   names the object (the §3d full-field requirement) rather than prescribing a construction recipe (e.g., it
   must not dictate a `ρ_4D` density multiplier — §1c/§3a do not name one). Flag any residual leak or
   over-specification. Re-check the measurements file's v2 leak-gate claim yourself.

5. **Blindness & buildability.** The engine must keep importing nothing and re-deriving independently. The WL
   line cites (L540/L541/L543/L554/L555; L1334–1340; L838/L923) must be accurate. A Mathematica engineer must
   be able to execute the item from the directive alone.

## Physics filter
Report a finding only if it catches a way the physics or the method could go wrong.

## Output
Numbered list. Lead with a one-line gate verdict: `CLEARS` or `DOES NOT CLEAR`. Then, per check, your finding
with the derivation path where you ran CAS. If the directive is clean on a check, say so. Do not edit any file.
