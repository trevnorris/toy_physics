I'll start by reading the build guide, method, and worker in the packet, then check the launcher, helpers, and saved operands against the stated method.The worker is in view. Next I’ll read the launcher, tests, manifest, and the native source files the method says the build must join.The manifest pins a large inherited-pressure surface. I’ll read the helpers, tests, and the native wave and binding sources next, focusing on joins that could pass without the actual operands.The review is still open. I am reading the remaining joins inside the packet: the pressure-wave definition, local denominator bases, helper runtime, and the addressed-assembly metadata.The remaining joins are the guard limits, the non-trivial pressure factor, local memory denominators, and the addressed-assembly metadata. I will read those packet files before choosing a verdict.NEEDS REVISION

The memory-domain join stops on actual local children. In `worker.py` the post-derivation check requires every negative-power base that contains `tau_A`, `tau_V`, or `tau_X` to have free symbols inside `{omega, rho_br, that tau}` and to equal `scale*I*rho_br*(1-I*omega*tau)`.

`THETA_BALANCE` local child 0 (`saved-inputs/inventory/THETA_BALANCE-full-native-partition.json`, `childIndex` 0, no pressure hit) has the base

`W_0*omega*rho_br*tau_A + I*W_0*rho_br`.

That is `W_0*I*rho_br*(1-I*omega*tau_A)`. The time sign matches, with constant scale `W_0`, but `W_0` is outside the allowed set, so `native local memory denominator form` raises before the identity is recorded.

`THETA_BALANCE` local child 9 has the base

`L_W**2*W_0**4*rho_br**2*(omega*tau_A+I)**2`,

which is `-L_W**2*W_0**4*rho_br**2*(1-I*omega*tau_A)**2`. `L_W` and `W_0` fail the same set check. The degree-one template is also unequal to that square. After the saved numeric binding these bases are finite constants, so the later quotient-domain rule would accept them. The raw symbol filter rejects the real negative-time factors first. `E_W_BALANCE` child 0, whose base is only `omega*rho_br*tau_X + I*rho_br`, would pass. The failure is specific to prefactors and powers that the saved local rows actually contain.

## Nonblocking

- `exact_nonzero_number` is loaded in `main` and never called. The live certificate is `constant()`.
- `effectiveSpeedOnlyInPressure` is written true in the binding emit before the raw speed census. A failed run can preserve that early flag. The success path still refuses a raw speed name.
- Final flags such as `nativeEpsilonOnce` and `wholeDirectOnce` are literals placed after the requires that establish them.
- Whole-tag handling checks operand names and counts. Representative factor operands carry the tag functions themselves: factor 10 is `Hwhole` plus `Jwhole_minus`, factor 4 is `Dwhole_plus`, factor 11 is `Dwhole_minus`, and factor 12 is the minus normal `-I*reference_qo` times the flat response.
- The supervisor still requires the historical endpoint progress file and a verified pool manifest. That dependency fails closed and is outside this packet.
- The gate file is absent, with inputs status `PREPARED_NO_SCIENCE_OR_GATE`. `verify_gate` still requires `READY_FOR_ONE_FULL_WEAK_INSTRUMENT` and both literal clear verdicts.

## Coverage and limits

The worker does select the local children before the pressure slots, reconstructs epsilon once and one wave, keeps eta and sigma independent, binds extra gammas from this physical input, and checks pressure child hashes, affine slot sums, and inherited cancel proofs by input plus saved zero. Factor 12 joins the normal argument `-I*reference_qo` to `-I*common_outgoing_q(composition_l)` and the product to the saved mapped factor. Containment in the helper runs before `import sympy`. The launcher arms the hook before the guard. The shared guard uses pool `s11c-near-unity`, 4 GiB, 32 tasks, zero swap, infinite `RuntimeMaxUSec`, and no retry. Local `T=±1` endpoints are labeled as local bookkeeping. No integral, inverse, current, or loss is computed.

The six ordered-address arrays and `weak-address-coverage.json` were not read. Their omission is recorded in `saved-input-map.json`. Tests were not run. They read those arrays from the original measurement paths. No scientific object was restored.