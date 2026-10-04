**CLEAR FOR THIS FIRST-ORDER RECEIVING-BLOCK AND FACE BUILD**

I found no blocker. This is a read-only source review. I ran nothing, so nothing here is runtime proof. The verdict covers the finite prerequisite only, and T1, R1, G1, K1, reflected flux, the radiation tail and all-real regularity remain required later.

**What I checked**
- **Local00 and C00 R00 S00 correspondence:**
  - Local00 is rebuilt from the 400 saved cells as `coeff·(il)^xOrder`, with duplicate and `cancelled==0` guards.
  - Source00 jets use `(−iω)^t` and `(il)^n`, a linearity check, and the saved reconstruction operands.
  - All 20 pressure pieces (5 rows × 4 slots) use the saved consumer and the whole flat factor with `q` independent, and are compared against `S·native·Sᵀ` for the full 5×5 off-wave matrix.
  - The field names in the assembly (`row`, `column`, `face`, `slot`, `consumers`, `flatNormalFactor`) and the flat-argument files (`plus`/`minus` × `pressure`/`normal`) match what the worker reads.
- **Coupling and thresholds:** both coupling blocks are compared against zero entry by entry, with the full transformed matrix emitted first and no Schur path. The threshold stage cancels denominators, saves them, and requires a finite nonzero determinant at both `l=±p`, `q=0`. The geometric columns satisfy `k·t = k·b = 0`, so a transverse–scalar decoupling is plausible. Pressure acts only on the x row, and `b` has x-component `H2`, so a nonzero b↔k coupling is possible on shell. That would be a legitimate stop, and it is preserved because the matrix is emitted before the check.
- **RIGHT end:**
  - The worker uses `nativeSource` with only the original map, so the old finite-1% `origin` pairs are not applied.
  - It replaces `eta_bg` with `lam` and `sigma_W` with `lam/10`. `nativeSource` has no unevaluated `Limit` atoms.
  - The LEFT binding carries `eta_bg` only in its `origin` metadata, so comparing LEFT with `Pweak` has no λ dependence to leak in.
  - The sign chain is correct: `−P1(p)U = localRight = 9U`, then `δ = −(DB·P1·U)₀₀ / D_T'(p) = 9/(3p) = 3/p`, then all five rows of `P1U + δ(P'U + P·B')`.
  - `B'` is carried through `B = E·C`.
- **Native faces:**
  - The face records carry orientation ±1, outward velocity, lift `W0/2`, pressure, mass flux, affinity and load.
  - Native pressure has the form `ρω/(q+ωΛ_A/(ρ(1−iωτ_A)))`, which matches the worker's `R` and `beta` check.
  - The load is `lift·(P+X·affinity)`, so pressure stays in it at `X=0`.
  - The direct `mu1` tag carries the saved envelope with the 1/10 factor and argument `l−p`. Pressure enters once through `R·(V+Aμ/ρ)`.
  - Orientation drives the face name, the lab-normal join `sign·i·q·R`, the velocity context, and the `liftEW` and V1 zero operands and return.
- **Pins and containment:**
  - The gate, argv, manifest, helper, method, build-review and authority pins are exact, and the output check precedes `mkdir`.
  - Helpers are loaded by name from the unchanged source, and `sympy` is imported only after containment.
  - Posthash covers every declared input, and the launcher uses no deadline and no retry.
  - The `VERDICT` string matches the requested wording.

**Advisories (not blockers)**
1. **Reversed-lab-normal control:** it compares `s·iqP` with `−s·iqP` on locally built expressions, so it is nearly tautological. It never passes through the saved normal factor, so it cannot catch a join error. A mutant built from the saved normal factor would be stronger. The `normal-flat-arguments` join and the actual-orientation check are the real protection.
2. **Direct-chemical control:** it needs the envelope at `T=0` to be nonzero. If it is zero, the control fails honestly rather than passing falsely.
3. **End-sign identity at free `q`:** `zeros('end-RHS-minus-operator-sign')` and the five-row end check require `−P1(p)U` to have zero θ and e_W rows identically in `q`. `P1`'s rows 3 and 4 hold q-dependent pressure and chemistry terms. I could not confirm that they cancel for both `U` columns without computing. If they don't, the run fails and preserves the operands. That would be a physics finding about the saved local force versus the end symbol, not a worker defect.
4. **Orientation and naming assumption:** the worker assumes `ORIENTATION=+1` means `plus`, and checks it only through the normal-factor join. A wrong assumption fails loudly.
5. **First failure aborts the stage:** a nonzero coupling or an unknown threshold would stop the end-sign and face-map joins too. That matches the required stop, but those later joins would then be unevidenced.
6. **Fourier convention:** the forward-transform sign and 2π convention are inherited from the saved S01 returns. The worker adds no independent check, and the tag is never evaluated.
7. **Memory and journal size:** `sympy.cancel` on the full 5×5 with native rationals, plus large journal emits, may be heavy within 4 GiB. A memory stop would be an honest failure.

**Coverage**
- I read the three named files in full, plus the launcher, manifest head, guide and tooling results.
- I spot-checked the evidence JSON for field names and structure, not every one of the 912 files. That included the RIGHT binding, the pressure-assembly structure, the native face and chemical subobjects, and the return and operand keys.
- I did not read every native leaf symbol, so unmapped leaves would raise `ValueError` at runtime rather than pass.
- The 14 tests are tooling-only. I did not hash-verify the pins.