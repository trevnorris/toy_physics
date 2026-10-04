**CLEAR FOR THIS BOUNDED KERNEL WAVE/MEASURE UNIT-TRANSPORT BUILD**

This is a source-level verdict only. It is not execution or result acceptance, and no evaluator, action, current, loss, inverse or sweep is cleared. I found no concrete blocker. The advice below is optional.

## What I read

**Read completely:**
- `build.md`
- `evidence-guide.md`
- `worker.py`
- `library.py`
- `source-contracts.json`
- `review-prompt.md`
- `tooling-tests.py`
- `tooling-test-record.json`
- `runtime-source/S11c_d_defect_packet_source_units_lib.py`
- `saved/native/binding.json`
- `saved/native/lower-geometry.json`
- `saved/native/fourier-contract.json`
- `saved/native/profile-scale.json`
- `saved/pressure/whole-tags.json`
- `saved/fourier/transform-convention.json`

**Read in part:**
- `inputs.json`: lines 1–70 and 12660–12723, plus greps of the alias counts.
- `saved/native/effective-registry.json`: lines 1–3246 of 4273, plus greps for `w1_profile`, `sigma_W`, `eta_bg` and `W_0`.
- `method.md`: lines 60–148.
- Original source:
  - `preflight.py` lines 195–285
  - `inner_lib.py` lines 40–170
  - `inner.py` lines 105–175
  - `reference_grazing.py` lines 160–180 and 285–345
  - `raw_increment.py` lines 425–460
- `native_units_lib.py`: only the grep hits for `unit_tuple`, `add`, `scale` and `Units`.
- `first-shape-native-transport.json`: lines 1–120.
- `fourier/family-0.json`.

**Sampled only, not read in full:**
- Address 10023, 10024 and 8619 input/return files.
- The count of `consumerRequiredUnit` values over all 544 address returns.
- The counts of live, zero, flat and `measure: dl` addresses.

**Not read:**
- The 17 factor-proof sets, the `preflight/numeric-factor-*` files other than a head, and the `inner/*` and `reference/*` proof files.
- `raw-direct/*`, `whole-origin/*` and the Fourier library source.
- The shared guard, supervisor, launcher and execution authority.
- The remaining registry lines, which hold the gamma entries.
- `packet-action-method.md` and the tooling log.

I did not execute anything. All dimensional checks below are my own hand walks of the pinned source text.

## Checked against original source and data

- **mu, a, beta:** The registry gives omega (0,-1,0), rho_m (-4,0,1), Lambda_A_0 (-5,1,1) and tau_A (0,1,0).
  - mu comes to (-4,-1,1), a to (3,1,-1), and beta to (-1,0,0), so `q+beta` is homogeneous.
  - omega·tau_A is dimensionless, so `1-i·omega·tau` is homogeneous.
  - The worker also separates input omega=1 from held omega=3.
- **Generic `kernel_components`:** J, reflected, height and quadratic each come to (-1,-1,1) as separate terms.
  - They keep distinct qi, qo, qh and qs, with no root product.
  - Whole J and D then carry (-2,-1,1) after dt.
  - W·L is forced to be length², so the unit check would catch a wrong W/L identification.
  - The saved Jwhole and Dwhole definitions carry the same factor structure.
- **First shape:** I walked `first-shape-native-transport.json` and got (0,-1,1). Removing the two edge deltas gives (-2,-1,1), which equals the 1D template unit.
  - The native map `hat_transfer → 2ĥ/W` and `jet_hat → 2ĵ` is consistent with hhat=L² and jhat=L.
- **H contact and PV:**
  - Hsub is L², the H density is L³, and the contact is L².
  - The W/4 coefficient alone is L, and its delta adds one more L to give hhat.
  - The `C` in the reference source is a local expression, not the registry `C`.
  - The contact route and the paired-PV route both total (-2,-1,1).
- **Templates:** All five are homogeneous.
  - Flat is (-3,-1,1) and the rest are (-2,-1,1).
  - The normal slot adds one q(l).
  - Consumer units split 272 at (-1,1,0) for the pressure slot and 272 at (0,1,0) for the normal slot, as `method.md` requires.
  - Jet orders tie to `sourceRequiredUnit` in the sampled cases: for 10024, (2,-1,0) plus the jet (-1,0,0) gives (1,-1,0).
- **Addresses:** There are 102 live and 442 explicit zeros. 288 are flat, with 288 `measure: dl` and 288 non-null `flatSupport`; the counts are consistent. Totals for flat, normal and height come out at (-2,-1,1).
- **Worker and library logic:** I traced the `require` bool usage, journal-name uniqueness and the AST walker's handling of the pinned expressions. I found no runtime defect.

## Optional advice (none of it blocks)

1. **Controls can't fail.**
   - The missing-dt controls and the flat-delta control compare `unit` against `unit±L`, which always differs.
   - Only the `a`-unitless control runs the actual kernel AST.
   - It would be stronger to feed each mutant into the real requirement (for example, whole J must be (-2,-1,1)) and show a refusal.
   - There is also no missing-normal-factor control. Dimensions can't catch a wrong normal depth, as `build.md` already says.
2. **Pin the W/L literal origins.**
   - `reference_grazing.py:165` has `W=sp.S.One; L=sp.Integer(10)`, and line 169 requires them to equal the physical `W_0` and `L_W`.
   - `raw_increment.py:443` has `L=values['L_W'];W=values['W_0']`.
   - Neither is in `source-contracts.json`.
   - The call sites that pass literals (`c.one`, `c.mpf(10)`, `sp.Integer(1)`, `sp.Integer(10)`) are only pinned and emitted.
   - The numeric A lambdas at `inner.py:138` and `preflight.py:237` are also unpinned.
   - Pinning these, with equality requires, would make the "original source assignments" join literal. Today it rests on the accepted algebra, which is a declared dependency.
3. **Bind the caller environment from the call ASTs.**
   - The worker hard-codes qi=q(k), qo=q(l), qh=q(k+t), qs=q(l−t) and the A1/A2 dimensionless assignment.
   - It parses neither pinned call nor `inner_lib.py:116`, which defines those depths.
4. **Assert the first-shape join.** `worker.py:173–175` computes the reduced unit (-2,-1,1) but never `require`s it equal to the template unit.
5. **Use the address's own measure and flat-support fields.**
   - Measure is chosen by the `NATIVE_FLAT` label.
   - Non-flat addresses are not required to have `flatSupport is None` or `measure == 'dl dk'`.
   - The current data are consistent.
6. **Check the normal sign.** The sign is emitted from the face label (`worker.py:245`) rather than compared with `a['normalMultiplier']` or the template. The inherited preflight zero proofs do bind it.
7. **Fix the `normalMap` key.** `worker.py:270` uses `a.get('completeNormalMap')`, but the key lives at `a['fullFactorProof']['completeNormalMap']`. Every summand-decision record will carry `normalMap: null`. The map is still preserved in the summand-input and template records.
8. **Confirm the single middle measure from the whole-definition tags.** The "plain dt, one middle integral" assumption is hard-coded. The tags' `boundVariable` entries (t, td, s) would confirm it.

## Dependencies this verdict leaves open

These are dependencies, not findings:
- The inherited gamma registry inference.
- The accepted numeric-to-symbolic algebra joins, including the A lambda, W/L and `qm→qh`.
- The Fourier/distribution conventions.
- The prior accepted source-transport required units.

Unit consistency does not distinguish q(l) from another momentum of the same unit. That rests on the original argument and normal proofs, which I did not re-read in full.