**CLEAR FOR THIS BOUNDED SOURCE-COMPOSITION METHOD**

I read the packet in the source/JSON form only. I ran no code and imported nothing. No peer report was supplied. The method holds up against the original c2 routing and the saved raw operands. I found no scientific, domain or claim blocker. Items 1–5 below are required guards and operands for the instrument; they are not blockers.

## What I checked against the raw evidence and c2 routing

- **Normal jet and face sign:** `build_face` builds the jet as `d/dN[exp(i·f·qo·(N−ref))·F_ref]` at the reference. That gives `i·f·q_out(l)·F_ref`, matching the method.
  - `retained-response-census.json` has `normalJet` `+i·qo` for the plus face and `−i·qo` for the minus face.
  - The plus and minus jet kernels differ by exactly that sign.
- **Reference pressure routing:** `plus-final-native-slot-routing.json` and `minus-final-native-slot-routing.json` give value coefficient 1, no constant operand, and `P_ref = P_phys − (±η·ĥ)·N`.
  - So the reference pressure is already trace-converted, as the method says.
  - The height is `f·η·ĥ`, with an opposite sign on the minus face.
  - The pointwise product with `N` appears in the saved kernels as a convolution through the middle variable `m`.
- **Fourier normalization:** `kernel_apply` uses `1/(2π)^3` on the source integrals and the profile hat is normalized (profile forward power −3, source forward power 0).
  - Reducing the two edge deltas gives plain `dk dl` with a normalized `1/(2π)` 1D transform.
  - Constants give `b·δ`, and the output multiplier gives `ĉ(r−l)`. This matches §5.
- **Source and consumer coefficients:**
  - Source side: the face velocity is `½W0·e_W_t·ε` with the identical constructor on both faces, and the plus and minus source grade-split hashes are identical.
  - Consumer side: the pressure coefficients are grade (0,0) and the normal-slot coefficients are grade (1,0), both from the saved consumer rows.
  - The consumer coefficients are single exact terms: `c^p` is constant and `c^n ∝ ε·η·w1_profile` (undifferentiated). The enumeration in §4 will therefore be mostly explicit zeros.
- **Pressure census:** THETA has one occurrence per slot and E_W has two per slot. U0, U1 and U2 have none.
  - `FACE_GENERALIZED_FORCE_ROWS.U` is the literal `(0,0,0)`.
  - The method already frames this as a census statement about those rows and not as decoupling.
- **Grade triples:** there are 1+3+3+9 = 16.
  - The direct (1,1) response enters only with source (0,0) and consumer (0,0).
  - The normal slot has no (0,0) consumer coefficient, so direct (1,1) enters only through the pressure slots.
  - The mixed iteration carries both height/slope assignments in the saved joins.
- **Speed independence:** `effective-speed-reuse-domain.json` shows no `cs` in the restored source or consumer symbols. `cs` enters only through `q`.
- **Controls:** each is applicable at its addressed route, with caveats.
  - Lower normal-jet sign reversal acts through `c^n(1,0)×F00×S00`. It drives the face sum to zero, which is a real movement and not an inert control.
  - `q(l)→q(r)` needs `r≠l`, which the nonconstant `ĉ` supplies.
  - `p→k` is inert at source grade (0,0) because the coefficient is constant. It needs a source grade of at least (1,0) with spatial jets.
  - Doubling the direct addend works through the pressure slots only.

## Required guards and operands (tooling and evidence, not blockers)

1. **The old selected increment is out of the rectangle.**
   - `native-pressure-consumer-omission.json` gives the THETA and E_W `selectedIncrement` with `η²·σ_W·w1_profile` and the whole kernel. That is consumer (1,0) × normal jet (1,1), which is grade (2,1) and is dropped by `retained_shape`.
   - The method's preface and §6 call the old selected check "(1,1)*(0,0)*(0,0)", but this normal-slot piece does not fit that label.
   - Treat it as excluded data. Do not match the new direct (1,1) block to it "by argument identity".
   - Match only a pressure-slot `c^p00` route. The same applies to the direct-addend control.
2. **The `whole_*` tags have no arguments.**
   - `whole_height_slope_convolution`, `whole_iterated_density_integral` and `certified_closed_direct_whole_convolution` are bare symbols in `retained-response-census.json`.
   - The direct certificate uses different variable names (`grazing_output`, `grazing_transfer`, `k`, ω+iδ).
   - Before any address can be written, each tag needs an explicit `(l,k[,m])` signature and a variable map to `reference_l/k/m`.
3. **The source grades are lumped.** `plus-source-grade-split.json` and `plus-chemical-grade-split.json` hold only `full`, `zeroGrade` and one lumped `higherGrades`.
   - The (1,0), (0,1) and (1,1) source pieces must be derived exactly from the saved `full` expression, with exact reconstruction.
   - The jet and derivative-order decomposition by `(A,α)` is also not saved and must be derived.
   - The saved zero-grade chemical source includes second-order spatial jets.
4. **No `Inputs` object means a reimplementation.**
   - A reimplemented single-pass `xreplace`, `at_source`, `retained_shape` and the `L^n` profile-jet scaling are needed.
   - `xreplace` is single-pass, so a leftover `W_bg` inside a replacement would be silently graded as (0,0). The fail-closed check for leftover background atoms must be real.
   - Validate by exact reproduction of the saved `full` and grade-zero expressions and the constructor hashes.
5. **Lower-face evidence is provenance only.**
   - The census limitations say lower-face source and consumer entries are unchanged provenance.
   - `units-and-scope.json` lists lower-face correction and cancellations as deferred, and `selected-fourier-contraction.json` has `lowerFaceCorrectionConstructed: false`.
   - Support for both-face equality here is the identical constructor text and the identical hashes, so state it as native-inherited.

## Unresolved obligations to carry forward

- **Branch points inside the certified window.** The `q(l)` branch points sit at `|l| = √(9/cs² − 1/20)`, which is 1.48 to 2.99 for cs in [1,2]. That is inside [−3,3]. Compositions integrate across them.
- **Normal multiplier growth.** `q(l)` grows like `|l|` at large `l`, so the normal-jet kernel does not decay by itself.
- **Nondecaying profile.** `w` is nondecaying, so its transform carries δ plus PV tails. Products of PV tails in `l` are defined only distributionally.
- **Real-frequency limit.** The saved response objects are for ω+iδ with δ>0. The real ω=3 prescription is explicitly uncomputed. The direct certificate covers the kernel only, on compact `k,l`.
- **Test space and transverse fixing.** The test-space hypothesis and the fixed transverse edge momenta (not compactly supported trials) are unresolved.

The method already marks all of these UNRESOLVED and makes no composed-grazing claim.

## What a successful instrument would establish

- A complete, hash-addressed, exact-algebra both-face ordered inventory with reconstruction against the actual rows.
- All 16 triples per face, slot and channel, including explicit zeros, with ε counted once.
- The routing of derivatives, normals and momenta through the algebraic multipliers.
- The applicable controls.

It would not establish any value, convergence, or a grazing limit.

## Still open before a finite near-unity defect pilot

- A domain, test-space and quadrature statement for the global momenta.
- A branch and PV treatment.
- The δ→0 handling.
- Control of the finite inverse and the cutoffs.
- The deferred lower-face cancellations.
- The matrix route.

## Runtime evidence needed

- Pinned census hashes.
- Exact reproduction of saved `full` and grade-zero operands.
- Exact derivation of the three missing source grades and the jet table.
- An exact affine-residual-zero check for each of the 12 pressure children.
- An enumeration of all 16 triples per face, slot and channel, with reasons for each zero.
- Tag signatures and variable maps.
- Algebraic Fourier routing at distinct rational `p,k,l,r` away from branch points.
- The controls run only on routes that are nonzero at the addressed face, slot and channel.
- A record of the old (2,1) increment as excluded.