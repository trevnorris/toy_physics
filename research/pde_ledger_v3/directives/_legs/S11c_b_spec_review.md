# Independent physics review — S11c-b SHARED PHYSICS spec

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`

This is a **specification** (an obligation-to-compute), not a script and not a results record. It is the
physics authority that TWO blind engines will read and independently execute: a SymPy engine (imports the
S11b `LEDGER` + the S11c-a exports) and a Wolfram engine (imports nothing; re-derives every object from this
spec and its two cited sibling specs). Because BOTH engines read this one document, an error HERE makes both
engines agree on the same wrong thing and the cross-engine residual cannot catch it. **Finding such a
shared-blind-spot error is the single most important job of this review.**

## What this sub-step is
S11c-b builds the **variable-coefficient brane (slab) operator** and its **off-diagonal transverse→{θ,e_W,u_L}
coupling kernel** — the first S11c physics object — on top of S11c-a's tilted-face shape-derivative substrate.
On a uniform background the transverse (photon) sector is decoupled from the thickness sector; a non-uniform
in-plane background (varying thickness `W_bg(y)` and modulus `μ_R,bg(y)`) is what can turn the coupling on.

## What you are handed (read these; nothing else is withheld — there is no "expected answer," this is a spec)
- The artifact above.
- The physics sources it inherits, to check every inherited equation and every cited line anchor against:
  - `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (the sibling
    S11c-a spec — §§1–3 are imported by S11c-b by reference; §4 T-a..T-i are the substrate).
  - `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` (the S11b spec —
    the source of the energy `U`, the symmetry group, the DOF, and the interface laws).
  - `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` (the S11b result
    — the uniform transverse decoupling and the non-unique energy-basis quotient).
  - `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md` (the S11c decision list,
    N-series — the requirements this spec must honour: N1,N3,N4,N5,N6,N8,N11,N12,N13,N14,N15).

## Required method (this is a DOCUMENT review)
Read the four sources FIRST and form your own view of (a) what physics S11b/S11c-a actually established and
(b) what a correct variable-coefficient coupling-kernel spec must ask for. ONLY THEN read the artifact and
compare. Quote both sides (artifact text + source line) for every finding. You do not need to run a full CAS
derivation, but where you claim an inherited equation is mis-stated, or a requested object is ill-posed, or a
cited line anchor is wrong, verify it against the named source file and quote the line. A prose assertion
with no source quote is discarded.

## Check, and report a finding only where the physics could be wrong or the spec would mislead both engines
1. **Shared-blind-spot / fidelity of inherited objects.** Is the stored energy `U` (§1c), the kinetic energy,
   the symmetry group (§1c), the DOF set (§1a), the interface laws (§1c), and the background ansatz (§2a)
   stated correctly and consistently with the cited S11b/S11c-a source lines? A wrong inherited equation is
   the worst possible defect. Verify the cited line anchors (`S11b_SHARED_PHYSICS.md:…`, `S11c_a_…`,
   `S11b_interface_coupling_law.md:…`) actually say what the spec attributes to them.
2. **The transverse/longitudinal split (§1a).** S11b never defined `u_T`/`u_L` — it carried the 3-vector `u`
   and used `∇×u` / `∇·u`. Is defining `u = u_L + u_T` (in-plane Helmholtz, `∇×u_L=0`, `∇·u_T=0`) as a
   supplied definition correct and sufficient to make the "off-diagonal transverse→{θ,e_W,u_L} block"
   well-defined? Does it misrepresent the S11b sector structure in any way?
3. **rule 5 — no leaked values / expected results.** The spec must state WHAT to compute, never what anything
   equals, is expected, or was measured. Flag ANY place it leaks: the coupling being nonzero / its sign /
   coefficient / grade; a specific new invariant or its coefficient; the S11b transverse-stiffness identity
   `μ_⊥ = μ_R + μ_S/2` (this must NOT appear); the admissibility residual being zero; the S11c-b basis count.
   Conversely, flag any place it OVER-withholds a genuinely SUPPLIED input (S11b under-specification —
   e.g. a missing governing equation — has cost this ledger more than contamination).
4. **rule 3 — name the object, not the recipe.** Does §3a (the N15 variable-coefficient basis) / §3b (the
   divergence-form operator) / §3c (the off-diagonal kernel) prescribe a derivation PATH (over-specify), or
   correctly name the object and let the engine compute it? Is the N15 "background first jets as spurion
   data" construction rule the right altitude — enough for a blind engine to build a well-defined basis,
   without dictating the answer? (Compare S11b task B0-energy at `S11b_SHARED_PHYSICS.md:272–288`.)
5. **The (ε,η,σ_W) power counting (N12/N5).** §2a requires every object multigraded and fixes the truncation
   (first order in ε, first shape order in η and σ_W) without stating any grade. Is that correct, or does it
   under-constrain so an engine could discard the leading coupling as "second order," or over-constrain by
   implying a grade? N12 states the coupling is O(εη) and leakage O(ε²η²) — is withholding the specific grade
   as a computed result the right call here, given σ_W (the first-jet bookkeeper) is separate from η?
6. **The admissibility computation (§2b, N12).** Is `S11CB_ADMISSIBILITY_{OPERATOR,SUPPORT,RESIDUAL}` a
   well-posed, can-fail test? Is the ε→0 stationary reduction of the §3 operator the right definition of the
   operator operand? Does the declared support bundle `𝒮_hold⁰ = {f_hold⁰, t_hold,s⁰}` carry enough to form
   the support operand? Is the spec right to forbid the "insert `W_bg−W_0` into the uniform equations"
   implementation (S11c-a §2d)?
7. **The controls (§5, N4/N6).** Is the representation-invariance route (Eulerian vs material-coordinate
   flatten, map back, compare) + the one-sided corruption the genuine control N6 demands? Is the uniform-limit
   correctly demoted to a smoke test (N6: `(η,σ_W)→0` is NOT an accepted corruption; a "coupling ∝ ∇W_bg ⇒
   vanishes" statement is the vacuous uniform limit renamed)? Is the §5b form ablation a FORM change (not a
   coefficient rescale)?
8. **Names / F9 (N14/G3).** Do the new objects take fresh names and avoid reusing the constant keys
   `W_0`,`mu_R`,`rho_br`,`e_W`,`v_0`? Is `v_bulk_normal_0` kept distinct from `v_0`?
9. **Scope boundary (N5/N10/N11).** Is anything that belongs to S11c-c (bulk DtN/impedance), S11c-d
   (spectrum/scattering — no global `ω(k)`), or S11c-e (falsification) wrongly pulled into S11c-b, or any
   S11c-b-owned object wrongly deferred? Is the carry-in `v_bulk_normal_0` standing-limit handled per N11a?
10. **The inheritance-by-reference architecture.** The spec inherits S11c-a §§1–3 by reference rather than
    restating them, and says the blind Wolfram engine re-derives the S11c-a substrate from the S11c-a spec.
    Is that clean — does the Wolfram engine have everything it needs from the two specs to re-derive without
    importing the sibling engine's output? Or does any needed supplied input go missing under the by-reference
    inheritance?

## Physics filter
Report a finding only if it catches a way the physics could be wrong, the engines could be misled into
agreeing on a wrong thing, or a rule (2,3,5,7 / an N-item) is violated. Do not report style. A leg that finds
nothing is weak evidence — if you conclude the spec is clean, say specifically WHICH shared-blind-spot classes
you checked and ruled out (inherited-energy fidelity, symmetry-group fidelity, the u-split, the admissibility
definition, the anchor citations).

## Output
A short list of findings, each: the artifact quote + the source line it contradicts (or the rule/N-item it
violates) + why it could make the physics wrong or mislead both engines + a one-line suggested repair. If
clean, the explicit ruled-out list above. Do not modify the working tree; this is read-only.
