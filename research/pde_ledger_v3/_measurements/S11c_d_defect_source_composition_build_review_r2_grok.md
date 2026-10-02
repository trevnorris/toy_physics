I'll review the bounded source-composition build from the packet inputs and the saved source they cite, then return one verdict with concrete blockers.The packet pins a source-only review. I'll read the worker, the evidence guide, and the cited local sources before judging the composition.The worker is long enough that the joins live in helpers and saved operands. I'll read those next and check the typed factors against the method.The adapter and grade joins depend on the native assignments and the exact nonzero test. I'll read those next, then the pinned operands they cite.The remaining load-bearing operands are still unread, so I am checking those source and JSON records before choosing a verdict.The route file stops before `kernel_bridge`, and the adapter namespace may not define every name the extracted assignments use. I am checking those joins against the helper and the remaining saved records.The remaining source checks are the full row-density join, the E_W omission tail, unit-join keys, source-jet names, and the census symbols the worker looks up. I’ll read those heads directly and then give the verdict.CLEAR FOR THIS BOUNDED SOURCE-COMPOSITION BUILD

The packet’s worker and the typed-factor method specify the bounded both-face inventory: native c2 pressure response plus the separate direct whole correction, at real omega 3, rest bulk `LAB_HELD`/`RHO4_CONSTANT`, saved edges 1/5 and 1/10, and effective bulk speed in [1, 2]. There is no scientific blocker in that scope. This is a source reading. It does not accept a run, a restored result, or an execution gate. `inputs.json` is still `PROPOSED_BUILD_NOT_EXECUTION_READY`, and the return forces `scientificAcceptance` false.

## Typed objects

Three objects stay separate in `input/worker.py` around lines 350–386 and 461–467.

- `Rprod = qi*qo/((qi+beta)*(qo+beta))` with `beta = omega/(10-I*omega)`. The bare closure factor, the isolated reference-factor pair, and the linear coefficient of `closed[0,2]` in the bare direct symbol are joined to this.
- `E = 1/((qi+beta)*(qo+beta))`, so `Rprod = qi*qo*E`. `E` is only that denominator factor.
- `Dwhole` is the saved complete closed density, with its own numerator, transfer, and sinh factors. It is not `Rprod` and not `E`. The inventory component is the tag `Dwhole_face`; it is not multiplied by `Rprod`, `E`, or `qi*qo`.

`kernel_bridge` in `input/source/native-c2.py` lines 409–411 puts a literal integer 0 at `z_three[0,2]`. The worker never calls `kernel_bridge`, `kernel_apply`, `reference_pressure_kernels`, `build_face`, or `upper_triangular_solve`. The formal adapter executes only the four saved assignments `extension`, `jet_transfer`, `normal_jet`, and `reference_pressure`, with diagonal and second fixed at 0 and source fixed at 1. Its probe `DeltaP = eta*sigma*Rprod*D` is a routing placeholder. The retained direct address is `INHERITED_DIRECT_WHOLE_OFF_DIAGONAL`.

Both new trace matrices have unit diagonal and `[0,2] = 0`. An increment confined to `[0,2]` therefore has no new `[1,2]` contribution. Plus height is `+eta_bg*reference_height_hat(-k+l)` and plus normal is `+I*qo`. Minus height is `-eta_bg*reference_height_hat(-k+l)` and minus normal is `-I*qo`. The unprojected reference is `P - height*I*sign*qo*P`, support `{(1,1),(2,1)}`. The retained grade is `P`; the excluded height correction is the `(2,1)` term.

H’s bound symbol is `reference_left_height_transfer`, checked against `variableChange[1] = l - hs`. J is bound by `composition_middle_transfer`, and D by `composition_direct_transfer`. Those three variables stay distinct.

## Row, address, and density joins

The native row substitution puts independent response tokens times source-jet coefficients into the current face’s pressure and normal slots, then extracts exact quotient grades. The address sum is `consumerOriginal * responsePlaceholder * sourceOriginal * sourceAtom` over the same grade. The normal factor `I*sign*q(l)` is stored on the token index and is already absent from the affine row coefficients, so it is not applied twice. Componentwise grades in G make `b_i+c_i≤1` automatic. U rows are identically zero. The separate noncommuting product is an enumeration check, not that row identity.

Saved full-density rights match the worker’s current operand, `density(at frequency 3) * (plus pressure-00 + minus pressure-00) * source-00`:

- THETA pressure-00 is `-I*epsilon_shape*(3-10*I)/109` on each face. Twice that, times the common source constant term, reproduces the printed THETA right, including the jet polynomial `-15*I*e_W + 6*I*e_W_d1d1 + …` and the pole `30/109+9*I/109`.
- E_W pressure-00 is `epsilon_shape/2` on each face. Their sum times the same source constant term reproduces the printed E_W right, whose jet polynomial is the negative of the THETA jet polynomial.
- Both canonical returns have `cancelled = Integer(0)`. Those zeros are inherited. The isolated reference-factor pairs are a different object and are joined only to `Rprod`.

Historical THETA and E_W omissions keep `removedSlot = delta_p_plus`, `mixedPerSource = 0`, and `ablatedRowPerD = 0`. Both `selectedIncrement` expressions are grade `(2,1)` in `(eta_bg, sigma_W)` and are not added to the new direct sum. Binding `Lambda_X_0 = 0` and `W_0 = 1` leaves the E_W pressure and normal coefficients that the grade splits already record.

Unit joins exist for all four THETA slots and all four E_W slots, and each saved total equals its expected total. With `Lambda_A_0 = 1/100`, `omega = 3`, `rho_m = 1/10`, and `tau_A = 1/10`, the THETA unit coefficients bind to the grade-split values, including `±I*epsilon_shape*eta_bg*w1_profile*(6-20*I)/436` on the normal slots. The normalized source dimension is the inherited `[1,-1,0]`.

The source constant term contains spatial jets, including `e_W_d1d1` and `u_1_d1`. The `(1,0)` quotient inherits a nonzero profile multiple of those jets through the `eta_bg*w1_profile` denominator, so a d1-addressed control candidate exists. Epsilon is absent from the saved source. Nonzero consumer grades are exactly one explicit factor of `epsilon_shape`; that count is homogeneity, not a nonzero-field proof. Controls then require an exact finite signed component through `exact_nonzero_number`, which emits the value before the guard and rejects an unknown zero, finite, or nonzero flag.

Flat support replaces the whole diagonal pole by `q(l)`. The `q(l)→q(r)` control scales the already-bound normal jet by `q(r)/q(l)`, moving the prefactor and leaving the output pole at `q(l)`. Direct omit and double move the numeric consumer–source coefficient by 0 and 2 while the `Dwhole` signature stays in the context record, with external resolvent multiplier 1. The positive outgoing sheet is `sqrt(9/cs^2-1/20-p^2)` when `p^2-9/cs^2 < -1/20`, and `I*sqrt` of the opposite radicand when greater. Zero-profile source, consumer, and nonflat response-tag reductions are checked. The nonzero constant-height full-response reduction is explicitly `NOT_COMPUTED`. That is enough for this inventory and is not a claim that the operator is globally applicable at constant height.

The guard path writes containment before SymPy, requires `RuntimeMaxUSec=infinity` and `Restart=no`, uses a 4 GiB job inside a 16 GiB pool with a 4 GiB host reserve, zero swap, one CPU, one thread, and 32 tasks, and does not impose a scientific deadline. The launcher arms the hook before `verify()` and the guarded command.

## Tooling and wording

These do not change a computed quantity. No wording revision is requested.

- `input/runtime-source/launcher.py` `verify()` uses `assert`, which `python -O` removes. The scientific gate in `worker.py` `verify_gate` uses `require`.
- The launcher’s 30-second `select` is a startup handshake.
- `input/evidence-guide.md` still says no scientific worker has been written. The implementation guide and `worker.py` are the instrument.

## Coverage still outside this inventory

No integral value, global composed momentum, grazing limit of the composed operator, cutoff, finite matrix, current, loss, defect sweep, or production replacement is claimed. `globalComposition`, `globalTestSpace`, and `composedGrazingLimit` are emitted `UNRESOLVED`. Internal response certificates remain `k,l` in `[-3,3]`; nonconstant multipliers can leave that interval. Effective bulk speed stays the depth-sheet parameter: the bound source and pressure consumers read here do not contain it, and a profiled coefficient that retained any other symbol would fail the `{x}` closure. There is no separate emitted name list of every speed-like symbol. U rows skip the full-density product because their pressure consumers are zero; their raw kernel operands are still joined. Consumer slot totals are inherited unit declarations joined to the bound coefficients, not a new dimensional derivation.

## Runtime evidence still required

A later guarded run must produce literal `checks.json`, journal inputs and raw residuals written before each zero or nonzero guard, posthash equality, and containment samples for the 4 GiB cap, zero swap, 32 pids, one CPU, and address-space limit. It must also produce the 320 coverage rows, the address sums, the control certificates, and the source-jet applicability records. `verify_gate` cannot open until a later review record carries both literal clear verdicts, the current worker and method hashes, and one authorized no-deadline run. That record does not exist in this packet.