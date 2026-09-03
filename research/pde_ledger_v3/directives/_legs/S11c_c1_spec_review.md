# Independent physics review — S11c-c1 SHARED PHYSICS spec (decision leg)

You are one of two independent reviewers of a **physics specification**. It will be read by two blind
computer-algebra engines (a SymPy engine and a from-scratch Wolfram engine) that each build the S11c-c1 physics
from it. An error in this spec is the worst kind: both engines would faithfully compute the same wrong thing and
their agreement would hide it. Find every such error **before** any engine is built.

⛔ This is a DOCUMENT review, not a script review. A prose claim that you "re-derived X and it works" is worth
nothing unless you show the computation. Where a question is settled by a short symbolic derivation, **write a
small script (SymPy), save it and its literal stdout to a named absolute path, and cite the path** — otherwise
your derivation claim is discarded. `sympy` is available. Save scripts under
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/` or `~/.s11_build/`.

## The artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md`

## What S11c-c1 is, and the split you must respect
S11c (staged family) computes the non-uniform (curved-brane) transverse↔thickness coupling. S11c-a produced the
tilted-face interface shape derivatives (T-a..T-i); S11c-b produced the variable-coefficient slab operator and
its off-diagonal coupling kernel, leaving the bulk pressure `δp_s` and the flat-face response kernels `Λ_I(ω)`
**symbolic** in every face/flux slot (it did NOT solve the bulk).

The old "S11c-c" was split into **two** build units (this is a ratified refinement — check it is respected):
- **S11c-c1 (THIS spec):** the **curved bulk closure** — solve the perturbed curved two-face outgoing bulk
  acoustic problem for the nonlocal DtN/impedance operator, compose with the interfacial mass balance + face
  closure (S11b's B0c on curved faces) for the permeable curved face response `(δp_s, J_s, t_s)(V_s, μ_θ)`, its
  degenerate/Fredholm loci, and its dissipation audit. c1 **exports** the closed face response + DtN operator.
- **S11c-c2 (NOT this spec):** fold the closed response into S11c-b's slab operator, re-extract the closed
  off-diagonal kernel from the closed full operator, produce the coupled nonlocal self-energy operator. c1 must
  NOT do this fold or re-derive any slab row.

## The sources of truth — read these and form YOUR OWN view of what S11c-c1 must require, BEFORE opening the spec
Build your own list of the objects, setup, controls, and hazards a correct curved-bulk-closure spec must contain;
THEN read the artifact and compare. (Reading order is method, not a blindness control.)
1. `research/pde_ledger_v3/directives/S11c_decisions.md` — the ratified N-series. N1 (staged family / chain
   wiring / own comparator+exports), N5 (⛔ no global dispersion `ω(k)`; emit the operator/kernel, not its
   spectrum; profile-conditioning is S11c-d), N8 (T7 + reconciliation schema), **N11a** (`v_bulk_normal_0`
   inherited as a standing rest-frame limit; ⛔ the convective bulk operator is NOT an S11c task), N12
   (multigrade), N13/N7 (confinement/falsification are S11c-e), N14 (fresh names).
2. `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` §1, §1b, §2, tasks B0a/B0b/B0c, and the energy
   check `:463-474`, the locus protocol `§8`, and the two traps `:760-816` (propagating `Re Z` is bulk
   RADIATION not interfacial transfer; branch re-selection turns a sink into a source). This is the FLAT bulk
   closure c1 generalizes: `Z ≡ δp_face/(bulk OUTWARD normal velocity)`, the radiation condition, the branch
   object `q_out` and its sound-cone branch points, the three regimes (bulk normal wavenumber² +/−/0=grazing),
   per-face inertial loading, the permeable face response, the degenerate loci.
3. `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — the tilted-face substrate c1 consumes: T-a
   normal (orientation law `s(n̂_s·ŵ)>0`, NOT `sign(∇₄F)`), T-a′ conormal, T-a″ measure, T-b velocity, T-c
   flux, T-c′ kinematic balance, T-d traction, T-e shifted trace, T-i closure shape-derivative (which EXCLUDES
   the bulk solve — the seam c1 fills).
4. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` §1c/§3b and
   `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md` — the slab operator with symbolic
   `δp_s`/`Λ` slots (c2 will close), the θ-row = `evolution_mass_balance − Σ closure_shape_deriv`, and the
   Λ-channel placement (`Λ_A𝒜+Λ_V V` in the flux closure `J_s`; `Λ_X𝒜` in the traction `t_s`, NOT the flux and
   NOT T-i). S11c-b is per-engine-verified with its cross-engine residual DEFERRED — c1 imports `μ_θ` as such.

## What to check — report a finding only where it catches a way the physics could be wrong
Derive your own answer first, then check the spec. In particular:

1. **The DtN as a two-momentum OPERATOR (rule 17 / N5).** c1 §3a insists the curved DtN is an operator
   `Z_s(ω;k,k′)` (flat Fourier symbol + first-shape-order two-momentum kernel / physical-space `W_bg(x)`
   product), and forbids collapsing it to a single-`k` multiplier `Z(k;∇W_bg(x))`. Independently derive the
   first-shape-order half-space DtN correction for a wavy boundary `w=s W_bg(x)/2` and confirm (or refute) that
   it depends on BOTH `q_out(k)` and `q_out(k′)` and is not Fourier-diagonal. Is the §3a operator framing
   correct and sufficient to prevent the freeze, and is it consistent with N5 (operator/kernel, not spectrum)?
   Is a rigid-shift (`k=k′`) cancellation the right sanity check?

2. **Face parity (rule 5 / N3).** c1 §2a states the lab-`w` graph slopes are opposite but refuses to infer a
   parity of any outward face object, because the OUTWARD normal's in-plane tilt is `−½∇W_bg` for both faces
   (even). Independently compute the outward normal `n̂_s` from `s(n̂_s·ŵ)>0` and confirm its in-plane tilt
   parity. Is c1's refusal-to-state correct, or does it (still, anywhere) leak or misstate the `(δW,ζ_c)` mixing?

3. **The permeable closure and the Λ_X placement.** c1 §1d puts `Λ_A/Λ_V` in the flux closure `J_s` and `Λ_X`
   in the traction `t_s` only, and §3b makes the bulk elimination an OPERATOR inverse `[I+(Λ_A/ρ_m²)Z]^{-1}`.
   Independently solve the flat B0c closure and confirm the Λ-channel placement and the operator-inverse
   structure. Is c1 right that c1 computes `t_s` (with `Λ_X`) but the routing of `J_s`/`t_s` into slab rows is
   c2's, not c1's?

4. **The degenerate loci vs the Fredholm condition (N5).** c1 §3b/§6 emits a FORMAL noninvertibility (Fredholm)
   condition for the operator, and restricts the algebraic CAS locus protocol to the flat Fourier-diagonal
   symbol and finite-dimensional coefficient degeneracies, deferring the profile-conditioned resolvent/singular
   set to S11c-d. Is this the right boundary — does forcing the operator's singular set through the algebraic
   protocol require an illegal generic spectral solve (N5)? Is the Fredholm condition the right object to emit
   in c1?

5. **The dissipation audit — Hermitian part + independent flux route.** c1 §3a/§3b define the reactive/radiative
   split as the operator HERMITIAN part under the true-area boundary pairing `⟨f,g⟩=Σ_s∫a_s f_s* g_s`, and add
   an INDEPENDENT energy route (slab pressure work vs independently-computed outgoing bulk flux, real ω /
   propagating / impermeable / Λ_X⁰=0). Independently check: is entrywise `Re` wrong for a mode-mixing operator?
   Is the two-route audit (operator-Hermitian + independent flux) sufficient to catch a shared sign error in the
   curved `t_s`, given S11c-b's face sign is still cross-engine-unvalidated?

6. **The radiation-preserving second route (§5a route 2).** c1 forbids the global scaling
   `w′=[w−ζ_c]/[W_bg+δW]` over the whole half-space (secular at infinity → wrong radiation branch) and requires
   a cutoff/Hanzawa flattening (= face map near the boundary, identity at infinity, transformed radiation
   condition stated) OR a boundary-integral/layer-potential route. Independently check the secularity claim
   (transform the half-space Laplacian under the global scaling and look for coefficients growing with the
   normal coordinate / a secular first variation of an outgoing wave). Is the cutoff/Hanzawa alternative a
   genuinely independent, radiation-preserving second construction of the SAME operator?

7. **The one-sided controls and their teeth.** §5a runs BOTH mutations (tilt slope-flip MANDATORY on the DtN;
   advection `Σ_E` on the response — structurally absent from `Z` since `ρ_m` is constant). §5b form-ablation.
   §5c uniform-limit regression (cannot validate first-shape-order curvature). §5d ZERO-JET regression
   (`σ_W→0`, `η` retained — catches a finite-gap-cavity `O(η)` error invisible to §5b/§5c). §5e branch/momentum
   liveness (sign-flip one leg; momentum-freeze the output leg). Do these controls actually bite the
   first-shape-order curvature, the two-momentum structure, the branch selection, and the `η`-order cavity
   error? Is any control tautological (`A−A≡0`) or a schema mismatch? Is any needed control missing?

8. **The rest-frame limit and the non-uniform grazing domain (N11a).** c1 §0/§2b keep the convective operator
   out of scope (legitimate at first shape order — the drain-projection correction is `O(σ_W²)`), but state the
   validity domain is non-uniform near grazing: grazing = the strict `v_bulk_normal_0=0` result, and away from
   grazing each result is conditional on `|q v₀/ω|≪1` AND `|ω v₀|/(c_s0²|q|)≪1` plus the subsonic condition.
   Independently check: does `|q v₀/ω|≪1` become vacuous as `q→0` while the correction it bounds does not? Is
   the boundary-layer exclusion + explicit limit order the right fix, and is excluding the convective operator
   itself legitimate at the requested order?

9. **The split seam (c1↔c2).** Does c1 export exactly what c2 needs (the closed face response `(δp_s,J_s,
   t_s)(V_s,μ_θ)`, the DtN operator, the flat symbol) and no more? Is anything that belongs to c1 (the bulk
   solve, the closure) wrongly deferred to c2, or anything that belongs to c2 (the fold into the slab operator,
   the re-extraction of the closed coupling kernel) wrongly pulled into c1? Is the c1/c2 boundary reviewable —
   can c1 be built and checked WITHOUT c2, and does c2's later isolation from S11c-b's deferred residual survive
   this split?

10. **N-series + chain output.** N1/N8 (own `scripts/S11c_c1_exports.py`, `BUILD_INPUT_DIGESTS`, D3 round-trip,
    `_RELATIONALS`, own T7 comparator — §7), N5 (operator not spectrum), N11a, N12 (multigrade; first order in
    ε, first shape order in η and σ_W — is this the right order to capture the leading curvature?), N13/N7
    (deferred), N14 (fresh `S11CC1_*` names, no imported-key reuse). Any leaked value (the two inherited S11b
    caveats — propagating `Re Z` = radiation, per-face inertial loading against outward acceleration — are they
    legitimate setup caveats or a leak of the curved answer?).

## Physics filter
Report a finding only if it catches a way the physics could be wrong, a leaked answer, a frozen varying quantity,
a missing/incorrectly-consumed input, a tautological/toothless control, a broken split seam, or a genuine
tractability blocker. Do NOT report style/wording unless a concrete misread produces wrong physics. For each
finding: quote the spec line, quote the source it contradicts (or your saved derivation path), and state the
concrete wrong-physics consequence and a proposed correction. A clean axis stated explicitly is a finding too.

## Output
A ranked list of findings, most-severe first, each with the spec location, the contradicted source or saved
derivation path, the wrong-physics consequence, and a concrete proposed correction. End with an overall verdict:
is c1 safe to build two blind engines against as-is, safe after the listed folds, or does it need a re-author?
