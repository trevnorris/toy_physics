NEEDS REVISION FOR THIS NONUNIFORM WORK AND LEADING SURVIVAL METHOD

This is a Claude-only method assessment. It is not worker or result clearance. I ran nothing, restored no opaque objects, and evaluated no leakage coefficient.

## Scope and actual coverage

**Read in full:**
- `guide.txt`, `plan.txt` and `review-prompt.md`.
- The full nonuniform `SOURCE_ENERGY` (`views/unbound-energy-source.json:5`).
- `original/source/physical-input.json` and `extended-binding-context.json` lines 1–80.
- `uniform-material-constraint.json` and `native-acoustic-ports.json` lines 1–60.
- Frozen-engine `UniformSlabCurrent` and `SlabEnergyBalance` (`frozen-engine.py:1685–2046`), and `ClosedCurrentPairing.construct` (`:2367–2549`).

**Read in part:**
- `LEFT-original.json`, at selected keys only: mass row, constraint remainders, balance residuals, surface density, chemical driver and memory kernels.
- Four saved local rows: U1 grade-00, U2 grade-10, THETA grade-10 and E_W grade-10.
- The all-local-cells header, with 200 cells carrying σ-grade 1.
- The first 40 text lines of `plus-UA-affine.json`.

**Not read:**
- `native-slab-work.json`, the minus-face records and the other three affine views.
- `flux/`, `matching/`, `regularity/` and `source-projection/` (about 8,166 files).

I performed no row-by-row Euler comparison, so the actual rows-versus-energy agreement is undecided. Everything below rests on these partial reads.

## Physics blockers

**1. The λ² identity cannot close in the retained model as proposed.**
- The engine takes the retained Taylor projection of the very quantities that make up the balance (`frozen-engine.py:1729–1734, 1967, 1978`). `ENERGY_BALANCE_RESIDUAL=0` therefore holds only through grades 00/10/01/11.
- `SOURCE_ENERGY` also looks already truncated at first order in η. For example, the factor `(1−2·η·w)` multiplies `γ15·σ·w'·eW·eW_d1`, which looks like a truncated `1/(1+ηw)²`.
- On the diagonal path the omitted grades weigh η²:ησ:σ² = 1 : 1/10 : 1/100. The retained model keeps the 1/10 term and drops the weight-1 term, so the retained λ² is not a controlled approximation of the physical one.
- The saved zero residuals at `LEFT-original.json:560, 592, 718, 730` come from the homogeneous background. The engine sets η=σ=0 when `end is None` (`frozen-engine.py:1749–1750`), and the constraint coefficients at `:534–560` contain no η. They say nothing about η-dependent constraint elimination or end currents.
- The plan's observable, `1 − (J_t+J_r)/(a†G0a)`, needs the right-end current matrix at order η². Only G0 and K1=G1+… are saved. Using the truncated `G0+ηG1` to define `J_t` puts a G2-type artifact into the deficit at the same order as L2.
- The packet also has no grade-20 or grade-02 rows or currents, so a finite check cannot test the physical-model λ² balance.

*Correction:*
- Declare that the target is the retained-model deficit, or the exact-in-η deficit of one stated functional. The currents, rows and energy must all belong to that one model.
- State explicitly that the right-end current in `J_t` is the conserved current of that model, including its η² part.
- Replace the plan's general instruction to show which omitted terms could enter with one narrow check: the contribution of the omitted 20/02 terms contracted with the transverse baseline, where θ0=eW0=0.
- If that contraction is not purely a total x-derivative or end current, stop and report the obstruction.

**2. The kinetic and thickness-inertia action, and the chemical normalization, are unsupplied for the nonuniform slab.**
- `SOURCE_ENERGY` has no time derivatives, so it carries no inertia.
- The only kinetic action is the homogeneous one in `frozen-engine.py:2423–2431`. It has constant `ρ_br` and `μ_W`, with thickness taken from the face-displacement difference.
- The saved local rows do carry a profile inertia. U2 grade-10 has `−9(1+tanh)/2` on `u3`, i.e. `−ω²·w`, which matches `ρ_br_bg=ρ_br(1+ηw)` at `extended-binding-context.json:77`.
- The chemical normalization is constant `ρ_br` in the engine (`LEFT-original.json:894–900`). The profile binding would need `1/ρ_br_bg`, whose 1/(1+ηw) expansion is truncated again.
- The plan's "derive the kinetic contribution from the native rows/action" risks circularity. If the inertia form is read off the rows, the energy-versus-rows identity is tautological.

*Correction:*
- State the nonuniform kinetic functional as an explicit hypothesis: `½ρ_br_bg|u_t|²`, the thickness inertia with its profile dependence, and the density-gradient term at grade σ.
- State which chemical normalization is used, then test it against all five saved rows with the W0 conjugacy factors.
- Any mismatch stays a named residual. It is not absorbed into LAB_HELD work.

**3. The saved rows are not purely energy-derived, so the residual bookkeeping must separate face response from stored-energy terms.**
- THETA grade-10 has complex coefficients proportional to `(1+tanh)`, such as `(300−1000i)/1090000`. These look like the `1/(1−iωτ_A)` memory kernel at ω=3, with x-dependent weight.
- The affine maps also carry the eliminated bulk response, with `q/(q+(30+9i)/109)`.
- *Correction:* subtract face memory and acoustic response using the same affine maps before forming the stored-energy residual, and show that the remainder is zero. The plan gestures at this but does not say which complex terms belong to the face memory. This does not change D_mem at λ² if χ0=0, but without the subtraction the residual cannot be classified.

**4. The incident wave appears to sit exactly at the exterior sound-cone edge, and first-order regularity does not give o(λ²).**
- The receiving shift is `l−√595/10` (`plus-UA-affine.json:73`), which is p. The shear modulus is μ_R+μ_S=3/2 (`SOURCE_ENERGY`), ω=3 and ρ_br=1. That gives |k|²=6, so kn²=5.95=p².
- I inferred this from displayed numbers and did not verify the dispersion relation. If it holds, the incident transverse wave has q(kn)=0.
- The right-end momentum then shifts by O(η), since the grade-10 inertia is `−ω²w`. The coalescence scale is |l−p|~λ, where q~λ^{1/2}.
- A fixed-order bound on P̂1 and L1 slab decay do not control `∫q|P̂(λ)−λP̂1|²dl = o(λ²)` near l=p.
- The plan's "if the expansion exists" hedge is correct, but its headline `+o(λ²)` is unsupported.

*Correction:*
- State that only the formal second-order functional of the first-order fields is derivable.
- Make existence and uniformity of the remainder a separate stated assumption.
- Record the ω=3, c=√6/2 binding as a model assumption. With the bare input (ω=1, c_s0=10) the exterior is evanescent at all l, because ω²/c²=0.01<|k_t|²=0.05. Radiative P_rad would then be identically zero.

**5. Do not require far-bulk and side-flux equivalence for the definition.**
- The slab finite-domain balance already gives the exterior power as the face-plane port work. Parseval on that plane gives the plan's `∫q|P̂|²/(4πρω)` over |l|<p. This needs only the outgoing choice of q and `q^{1/2}P̂1∈L²`.
- The plan's requirement to "prove … control side surfaces before claiming equivalence to total far-bulk flux" (plan lines 191–193) is stronger than the method needs. It also leaves the depth-versus-x limit-order problem open.
- *Correction:* define leakage as the two-face plane power. Treat far-bulk escape as an optional later statement.

**6. L1 decay needs to be matched to the maps that use it.**
- L1 plus boundedness of `E3χ3` gives L² by Plancherel, which is enough for D_mem if the face map only uses the derivative orders for which decay was shown.
- The face maps use jets up to order 2. Show decay for exactly those orders, or state the one that fails.
- For P_rad only boundedness on |l|<p is needed, so the plan's grazing argument is adequate there.

## What is sound and reusable

- The pointwise identity `(P+Xχ)V+μj = Pv+χj+XχV` follows algebraically from `j=ρ_m(v−V)` and `χ=μ−P/ρ_m`. It still has to be shown to be the slab-bound work, as the plan says.
- The memory relation `χj = j²/Λ + d_t[τj²/(2Λ)]` is right, and `Re A_mem=Λ/(1+ω²τ²)` matches the stored kernel.
- The Fourier normalizations check out: D_mem coefficient `Re A/(4π)`, and `v̂=qP̂/(ρω)` from `P=−ρφ_t`.
- The face lift `W0·eW/2` and velocity `−iωW0·eW/2` match the plus-face affine map.
- The structural point that all sinks are face-quadratic with zero baseline face fields avoids u2·baseline work. This is correct if the plan's own derivations hold.

## Missing physical input: none proven yet

Nothing here needs the user yet. The decisions in blockers 1, 2 and 4 are method choices or finite checks against the saved rows. They become a user question only in two cases:
- A residual between the saved rows and the stored-energy-plus-inertia-plus-face-response form survives, which would mean an unspecified external agency.
- The user wants the leakage defined for the original ω=1, c_s0=10 input instead of the substituted ω=3, c=√6/2.

## Optional wording, separate from the blockers

- Plan line 122: `(1/2)Re(A* conjugate(B))` double-conjugates. It should read `½Re(Ā·B)`.
- Controls (plan lines 215–218): add a mutation of `ρ_br` versus `ρ_br_bg` in the chemical normalization.