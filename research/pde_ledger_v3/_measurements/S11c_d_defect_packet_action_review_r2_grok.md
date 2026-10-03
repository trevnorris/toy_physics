I'll review the packet-action method from the named inputs and the sources they point to, using only those files.The method and evidence boundary are in. Next I’ll read the packet index and the source files it names, without running any computation.The index lists the exact sources. I’ll read the review prompt and the load-bearing definitions next, then check the address and cell selections against them.CLEAR FOR THIS BOUNDED PACKET-ACTION METHOD

No substantive blocker. This review clears only the bounded method for two complex-bilinear Schwartz pairings of the corrected retained operator, input `e_W` into row `THETA_BALANCE`, at the one saved left speed. It does not clear a worker, a build, a quadrature result, or the old finite benchmark.

The calculation stays inside one fixed pairing unit. It keeps both native faces, all selected addresses and grades, and the local cells. It is not a transverse mode, a current transfer, a scattering solve, or a loss run. Earlier source, response, local, and weak-end records are restored as saved objects.

## Scope, selection, epsilon, grades, derivative order

The speed and branch point are the saved left grazing operands in `ends/left-match.json`: `weak_end_cs = sqrt(6)/2`, `weak_end_p = sqrt(595)/10`, `weak_end_q = 0`. With real frequency 3, `omega^2/cs^2 = 6`, so `q(p) = sqrt(6 - 1/20 - p^2)` vanishes at that `p`. The development-file parameter `omega = 1` is the old input value; `local/context.json` and both weak methods override it with held real frequency 3. `c_s0 = 10` stays out of this pairing. Tangents `1/5` and `1/10`, `W_0 = 1`, and `L_W = 10` match `physical-input.json`.

`selected/pressure-addresses.json` records `originalCount` 2652, `selectedCount` 544, and `formalCount` 102, under `row THETA_BALANCE` and `jet.channel == e_W`. `selected/local-cells.json` records `originalCount` 400 and `selectedCount` 16, under `row == THETA_BALANCE`, `field == e_W`, all `xOrder` and grade. That is the full local rectangle: `x` orders `0..3` and grades `(0,0)`, `(1,0)`, `(0,1)`, `(1,1)`. Sampled headers include `(xOrder 0, grade (0,0))`, `(0, (1,1))`, and `(1, (0,1))`.

Epsilon is the coefficient of `epsilon_shape` once. Formal address 7956 has `epsilonCount` 1. Exact-zero addresses have `epsilonCount` 0, with the saved meaning that an exact zero carries no power. Independent output grades stay `(0,0)`, `(1,0)`, `(0,1)`, `(1,1)`. In `local/context.json` the grade order is `eta_bg`, then `sigma_W`. The optional display `B00 + (1/100) B10 + (1/1000) B01 + (1/100000) B11` is those origin values from `ends/left-source-binding.json` key `origin`: `eta_bg = 1/100`, `sigma_W = 1/1000`, and their product. `sigma_W` is not a separate parameter in `physical-input.json`. The display is after the four grade tables and does not add `eta^2` or `sigma^2`.

Derivative order matches the saved jets. Address 70020 area, jet `e_W_t`, has `waveMultiplier = -3*I`. Address 90018, jet `e_W_d2d2` with spatial orders `(0,2,0)`, has `waveMultiplier = -1/25 = (i/5)^2`. The source field is stored separately, including when it is zero, so the source is `b(x) D_j u(x)`.

## Flat support, contact, PV, and typed wholes

Flat address 7956, `plus` / `pressure` / `NATIVE_FLAT`, has `deltaSupport` `k=l`, `l=p`, `r=l`, `flatSupport` `k=l`, and output depth `q(l)`. On that support the pole and normal depth collapse together. Normal flat address 8892 carries `normalMultiplier = I*q(l)`; the last minus flat address 10576 carries `-I*q(l)`. That is `i*f*q(l)` with face signs `f = +1` and `f = -1`.

The height action in the new method is the same contact-plus-integral prefactor as pressure-method equation (7):

`(W/4) f_k(0) + W/(2i) ∫ A(Q) [f_k(Q) - f_k(0) chi(|Q|<=1)] / Q dQ`.

The half-line form with `[f(Q)-f(-Q)]` is the even/odd reduction of that formula because `A(Q) = L Q / [4 sinh(pi L Q / 2)]` is even and `A(0) = 1/(2pi)`. For `L = 10` this is the saved profile `5 t / (2 sinh(5 pi t))`. The method requires that equality be derived before use, and it forbids sampling the quotient at `Q = 0`.

Saved `H` contact `5*(-10k+10l)/(16 sinh(5 pi (l-k)))` equals `(W/4) j(l-k)` with `j(Q) = L A(Q)/2`. The candidate profiles `h(x) = W(1+tanh(x/L))/4` and `j(x) = sech^2(x/L)/4` are the `1/(2pi)` representatives of those kernels: `h' = (W/L) j`, and `hat{j}` reproduces `L A/2`. They are a join check against saved `H`, not a replacement. `Jwhole` has middle momentum `k+t` and reflected momentum null. `Dwhole` is a separate complete integral, entered once, and is not passed back through the native middle integral. Bare `Rprod`, factored external denominator, and whole direct density stay distinct in `pressure/typed-direct.json`.

## Branch geometry and the direct sum

Outgoing `q` is positive real for `|p| < kappa`, positive imaginary for `|p| > kappa`, and zero at `±kappa`. Inner splits are `t = ±kappa - k`, `t = l ± kappa`, and the removable points `t = 0`, `t = l-k`. Coincident endpoints merge with both labels kept. No hole is left, and a `0/0` depth ratio is not filled with zero.

Equating a middle root endpoint with a reflected root endpoint gives `k+l = (s1-s2) kappa` for `s1,s2 in {±1}`, hence `k+l = 0, ±2 kappa`. The line `l = k` is `t = 0 = l-k`. Root/profile hits reduce to `k = ±kappa` or `l = ±kappa`. The stated `l` lines are `-K`, `K`, `±kappa`, `-k`, `-k±2 kappa`, and `k`, plus the carrier and profile resolution lines and every in-box intersection, including the box edges. Parallel families correctly contribute no intersection. Route A’s two-ended square substitution does not sample those endpoints.

The four nongrazing nodes hit those internal lines and no external grazing line:

- `(kappa/3, -kappa/3)` is `k+l = 0`
- `(kappa/3, 5 kappa/3)` is `k+l = 2 kappa`
- `(-kappa/3, -5 kappa/3)` is `k+l = -2 kappa`
- `(kappa/3, kappa/3)` is `l = k`

Shifts `±kappa/64` and `±kappa/128` leave every line in that set. The method does not treat those offsets as a derivative or an extrapolated limit.

Saved direct `Bc` is a sum, not a pinch product. In `pressure/typed-direct.json` and pressure-method section 4 the three terms are

`k(2l-t)/(qo+qs) + qi^2/qh + k(t+2k) qi / [qh (qh+qi)]`,

with `qh = q(k+t)` and `qs = q(l-t)`. There is no factor `1/(qh qs)`. `J` contains only `qm`. The inherited continuity argument is a sum of simple-root envelopes and does not claim differentiability at a collision.

## Fourier envelope, tails, routes, and controls

The strip `|Im x| <= 5` lies inside the `tanh(x/10)` poles at `|Im x| = 5 pi`, and `|tanh| <= 1` there because `5/10 < pi/4`. Shifting by `-i*5*sign(nu)` produces the Gaussian factor `exp(25/(2 s^2) - 5|nu|)`. Absolute-coefficient and absolute-moment majorants, with the real center shift included, are a concrete tail recipe. The build still has to evaluate `C_(field,jet)` and the finite `K`, `T` tails. Inherited polynomial constants remain existence bounds. The method requires a new truncated-tail inequality below `1e-11` per unweighted grade and component, separately for the `(k,l)` square and the `(k,Q)` height domain, and it forbids treating cutoffs 4 or 6, or a larger box, as that bound.

Fourier routes do not share nodes, weights, or arrays. The stop is an absolute difference below `1e-12`, with route B’s internal target at most `1e-13`. The method states that contour damping does not prove relative accuracy at large `|nu|`, and a miss stops the action. Complete routes are independent on the same line set. Declared settings are Gauss orders 24/48, Fourier orders 96/160, route B at 50 decimal digits, and the enlargement `K → K+2`, `T → T+2`. The acceptance width `tau = 1e-9 + 1e-7 |refined value|` is empirical. The analytic tail is separate. Their sum is not a flux or leakage bound, and a component may not cancel out of the comparison.

Controls are addressed and predeclared: drop the surviving `H` contact, replace reflected `q(l-t)` by `q(k+t)`, replace normal `q(l)` by `q(k)`, and move a nonconstant coefficient inside a nonzero derivative. Off-diagonal normal height is a real target for the depth swap: `plus-(1,0)-NATIVE_HEIGHT` has `flatSupport` null, `qi = q(k)`, `qo = q(l)`, and normal factor `I*q(l)`. Local cells include a nonconstant `tanh` polynomial on a positive `x` order, so the Leibniz control has a metadata target. Silence in both packet settings stays unestablished. A nonzero movement must exceed ten times the combined empirical envelopes. Formal movement is not a numerical success. The transverse lift in `ends/left-match.json` is zero on the `theta` and `e_W` rows, so this block is outside what that end comparison tested.

## Blockers and optional notes

Substantive blockers: none. No correction is required.

Optional, not a method change: the inner spot checks named “zero transfer”, “both external grazing approaches”, and “nonzero reflected momentum” do not have coordinates. The lines themselves are already in the arrangement, and the four collision pairs plus their two-sided shifts are fixed. The phrase “16 ordered source/response/consumer grade triples” means the 16 local `(xOrder, grade)` cells. Pressure grade triples remain the per-address source, response, and consumer grades across all 544 entries.

The exact `h`/`j` adapter joins, numerical tail constants, and execution of the controls are concrete-worker obligations under this prescription. They do not require a larger method change. The symbolic tag note “UNRESOLVED outside `k,l` in `[-3,3]`” is not a cutoff: the later global envelopes are all-real existence bounds, and this method still requires a derived finite-`K`/`T` tail.

## Files inspected

Read in full: `input/method.md`, `input/evidence-guide.md`, `input/review-prompt.md`, `input/packet-index.json`, `input/physical-input.json`, `input/local/assembly.json`, `input/pressure/typed-direct.json`, `input/pressure/whole-definitions.json`, `input/pressure/whole-tags.json`, `input/pressure/Fourier-order.json`, `input/pressure/source-jets.json`, `input/pressure/global-parameter-domain.json`, `input/pressure/profile-envelope.json`, `input/pressure/whole-envelopes.json`, `input/pressure/normal-growth.json`, `input/pressure/shift-root-bound.json`, `input/pressure/H-bound.json`, `input/pressure/J-numerator-envelope.json`, `input/pressure/D-height-envelope.json`, `input/pressure/D-reflected-envelope.json`, `input/pressure/plus-height-PV.json`, `input/pressure/minus-height-PV.json`, `input/ends/left-match.json`, `input/ends/correspondence-summary.json`, `input/background/pressure-weak-method.md`, `input/background/full-weak-method.md`.

Read in full but only the load-bearing keys were used from the long symbolic matrices: `input/ends/left-source-binding.json` (through `origin` and `halfHeights`). `input/local/context.json` was read through the frequency override, `sigma_W = W_0*eta_bg/L_W`, and `independentGrades = [eta_bg, sigma_W]`; the later repeated physical-parameter block was not re-read line by line.

Partial: `input/selected/pressure-addresses.json` (header counts; addresses 7956, 7957, 8892, 10576; `NATIVE_HEIGHT` maps `plus-(1,0)` and `minus-(1,0)`; jet multipliers `-3*I` and `-1/25`). Statuses seen: `FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED`, `EXACT_ZERO_SOURCE_JET`, `EXACT_ZERO_CONSUMER`. The 106/336 zero split was not re-counted entry by entry. `input/selected/local-cells.json` (header and sampled cell headers and polynomials, not all 16 headers). `input/pressure/fields.json` and `input/pressure/coefficient-certificates.json` (opening constant and nonconstant `tanh` records). `input/background/numerical-route-inventory.json` (opening legacy-route index only).

Not read line by line: `input/source/composition.py`, `input/source/pressure-weak.py`, `input/source/full-weak.py`, `input/source/weak-ends.py`, `input/source/legacy-integration.py`, `input/source/legacy-finite.py`. The method restores the saved JSON definitions and does not replay those instruments.