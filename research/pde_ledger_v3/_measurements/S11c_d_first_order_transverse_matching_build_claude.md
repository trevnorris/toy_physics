**CLEAR FOR THIS FINITE MATCHED TRANSVERSE AMPLITUDE BUILD**

This verdict covers the transverse amplitude extraction in saved incident coordinates only. It does not cover flux normalization, leakage, all-real coupled regularity or other-channel scattering. I used Read/Glob/Grep only. I found no blocker in the worker, build notes or method plan. The manifest, gate, launcher, helper source, guard and supervisor were not in the packet, so I could not inspect them.

**What I checked by hand against the saved operands**
- **Forces and step identity:**
  - Both `etaForce` columns equal `A·U` with `A=9w−6m=2.5+4.5T+2T²`. The right limit is `9U` and the left limit is 0.
  - `physicalPathForce = eta + sigma/10` matches in the constant and T² coefficients that I spot-checked.
  - All five rows of the sigma force vanish at `T=±1`. I checked several entries, including rows 3 and 4 and the unexpanded row-3 product. So the single approached-end factor division is valid, including the scalar rows the projector discards.
- **Pressure:** The pressure assembly has 40 literal-zero consumer entries for the displacement rows. That is consistent with the `require`, and nothing pressure-related enters the projection.
- **Half-line moments:**
  - The substitution `dx=L dT/(1−T²)` and the division by `(1±T)` are correct.
  - The remaining simple fraction integrates to `−r·log(1−T)` on the left and `+r·log(1+T)` on the right.
  - Each primitive is differentiated back to its integrand, and the endpoint `log 2` terms are finite.
- **Green function and Abel step:**
  - The jump `−(3/2)[G′]=1` holds for `G=i e^{ip|x−y|}/(3p)`, and `δ=3/p=6√595/119`.
  - The Abel integral `∫₀^∞ e^{2ipy}=−1/(2ip)` gives `−δ/(2p)` for both the reflection and the forward finite term.
  - The linearized exact `t` and `r` agree: `r=(p−p_R)/(p+p_R)` and `t=2p/(p+p_R)`, both with coefficient `−δ/(2p)`.
  - The right-end secular term is `iδxU`.
- **Sign of `δ·D_B′(p)U`:**
  - Expanding `D_B(p+Q)·9U/(iQ)` gives the regular constant `−9i·D_B′U`.
  - Multiplying by `i/(3p)` gives `+δ·D_B′U`.
  - The output-side derivative gives `δB′`, so the total is `δΠ′U=δ(B′+B·D_B′U)`, as coded.
  - The code's `dprojected` split and the final join `δB′+U·T1` are the same algebra.
- **Saved `D_B′U`:**
  - From the saved `DBprimeAtIncident` and `C`, `D_B′(p)U = [[−p/6, 0],[−p²/3, 0]]`.
  - I confirmed the (0,0), (0,1) and (1,0) entries from the saved numbers.
  - The (1,0) entry is `≈−1.983`, nonzero and real.
  - So the derivative-omission control is applicable.
- **Imaginary and off-diagonal sigma terms:** They enter in full through `D_B(p)·moment` on all five rows. No Hermiticity is assumed.
- **Reflection:**
  - It is `−δ·D_B(−p)U/(2p)` plus `i/(3p)·∫e^{2ipx}D_B(−p)F_loc`, with `Q=−2p`.
  - The phase `e^{ipy}·e^{ipy}` matches, and the `D_B(−p)` evaluation follows from integrating the multiplier by parts.
  - The integral stays as a tagged, unevaluated integral and is absolutely convergent because of the `2/L` tails.
- **Coordinate-tail argument:**
  - The only coordinate poles are `l=±iH` with `H²=1/20`, and the kernels are `e^{−H|x|}/(2H)`.
  - The tail rates are about 0.22 from the kernel and 0.2 from `2/L`, so both decay.
  - Both off-wave blocks are inherited as zero and the dual rows annihilate `k`, `θ` and `e_W`, so no acoustic `q` enters.
  - The worker honestly sets `fullCoupledFieldAsymptoticsProved: False`.

**Inherited inputs and runtime**
- Source replay is absent. Only saved operands, sympy arithmetic and the helper census are used.
- The `D_T=(3/2)(l²−p²)` operand, the end-forcing operand and the incident-chart operand are compared as structural equalities with the saved operands. A mismatch fails closed.
- Argv and output are pinned by exact command and `outputDirectory` equality.
- The gate checks worker, manifest, launcher, source, review, method, authority, guard and supervisor hashes.
- The resource tuple is fixed with `durationLimits: None`. Post-hashes cover every copied saved file, the manifest and the gate.

**Coverage limits and advisories (none blocking)**
1. **The controls are weak.**
   - The `D_B′` omission control is a value-sensitivity check on `T1[1,0]`. It is not run through the `right-field-finite-matching` join.
   - The reversed scalar-step control compares `s` with `−s`, so it passes for any nonzero `s`.
   - The real sign checks are the exact `t` and `r` zero tests, my derivation above, and the join.
   - The join is largely an identity once `B(p)=U` and `D_B(p)U=I`. It validates the chart algebra and the inherited `B′` and `D_B′`, not the extraction logic.
   - Do not report these controls as independent evidence of the sign.
2. **Helper pin gap.** The worker does not assert that `manifest['helperSource']` is in `sourcePins` or hashed. When the manifest and gate are assembled, confirm it is pinned, along with the method bytes matching this `method-plan.txt`. I could not verify either.
3. **File count mismatch.** The packet holds 751 files, but the notes describe 732 complete JSON inputs. That is probably the `-raw` files, but I cannot reconcile it without the manifest. Confirm that the manifest `savedFiles` covers every consumed alias.
4. **Possible fail-closed stops.** These depend on `J.zero` and `J.nonzero` simplifying sqrt(595), `I` and `log 2` expressions, which I could not inspect. Examples are the `rem` zero checks, the primitive derivative checks and the `is_finite` requirement. A failure there would be preserved, not a wrong result.
5. **Pending.** These items remain, as the notes state:
   - `G0`, `G_ref`, `G1` and `K1`.
   - Origin and gauge checks.
   - The face maps.
   - Evaluation of the reflection integrals.
   - All-real coupled regularity.
   - Other-channel radiation, including the longitudinal source from the `1/Q` step.

Source-only tests are not native scientific proof. `T1` in saved coordinates must not be read as attenuation or leakage.