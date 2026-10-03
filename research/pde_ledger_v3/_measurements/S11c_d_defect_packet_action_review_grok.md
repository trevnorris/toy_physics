I'll start by reading the method packet and evidence guide, then inspect only the sources they name.The packet map is in. Next I’ll inventory the submitted sources and read them in the order the evidence guide names.I’ll read the review prompt, physical input, and background methods first, then the smaller certificates.Next I’ll read the remaining pressure certificates, ends files, and the smaller source inventories.I’ll inspect the large selected inventories and the LEFT matching operands in bounded reads.I’ll sample later address records, remaining local cells, and the source-interface formulas.I’ll locate cell headers and distinct pressure-address components by sampling further through those files.I’ll read the source-interface files for the kernel formulas and then sample remaining address types.**CLEAR FOR THIS BOUNDED PACKET-ACTION METHOD**

This packet specifies two Schwartz pairings of the restored retained operator, \(e_W\to\) `THETA_BALANCE`, at the saved LEFT matching speed. It is a numerical-method proposal. It does not clear a worker, a build, a scattering/loss claim, or the old \(c_s=10\) finite benchmark.

No substantive blocker. Adapter joins, explicit \(C_{\mathrm{field,jet}}\) numbers, \(K,T\), and addressed control applicability remain concrete worker obligations, with stop/limitation branches already written.

---

## Scope and units

The object is the \(\varepsilon\) coefficient of \(B_g(v,u)\) for \(g\in\{00,10,01,11\}\), with only \(u_{e_W}\) and \(v_{\mathrm{THETA}}\) nonzero. Packets

\[
u_{e_W}(x)=e^{-(x-x_u)^2/(2s^2)}e^{ip_0 x},\qquad
v_{\mathrm{THETA}}(x)=e^{-(x-x_v)^2/(2s^2)}e^{-ip_0 x}
\]

at \(s=8\), \(x_u=-L/4\), \(x_v=L/4\), \(p_0=\kappa\) and \(p_0=0\), are complex bilinear (no conjugation). Results stay in the \(e_W\)–THETA pairing unit.

LEFT operands in `ends/left-match.json` are \(c_s=\sqrt6/2\), \(\kappa=\sqrt{595}/10\), \(q=0\). Then \(\omega^2/c_s^2=6\) and

\[
q(p)=\sqrt{6-1/20-p^2},\qquad \kappa=\sqrt{6-1/20}.
\]

`physical-input.json` still carries development \(c_{s0}=10\) and \(\omega=1\); `local/context.json` records `frequencyOverride` \(1\to3\) and `effectiveSpeedOnlyInPressure: true`. The RIGHT matching speed is unused. Jet rule

\[
D_j=(-3i)^{n_t}(i/5)^{n_2}(i/10)^{n_3}\partial_x^{n_1}
\]

matches `background/full-weak-method.md` and address `waveMultiplier` samples (`I*composition_p` on \(e_W_{d1}\); \(-3I\) on \(e_W_t\); \(-1/25=(i/5)^2\) on \(e_W_{d2d2}\)).

Source is \(b(x)D_j u(x)\). Epsilon is taken once; grade weighting is last. The optional \(B_{00}+(1/100)B_{10}+\cdots\) display is not the object.

---

## 544 addresses and 16 local cells

`selected/pressure-addresses.json`: `originalCount` 2652, `selectedCount` 544, `formalCount` 102, predicate `row THETA_BALANCE, jet.channel == 'e_W'`. That matches \(102+106+336=544\).

Inherited grade-route contract (`source/pressure-weak.py`): 16 triples \((c,r,s)\) with \(c+r+s\in G\), 2 faces, 2 slots, components

| \(r\) | components | # triples |
|---|---|---|
| \((0,0)\) | `NATIVE_FLAT` | 9 |
| \((1,0)\) | `NATIVE_HEIGHT` | 3 |
| \((0,1)\) | `NATIVE_SLOPE` | 3 |
| \((1,1)\) | `NATIVE_MIXED_ITERATION`, `INHERITED_DIRECT_WHOLE_OFF_DIAGONAL` | 1 |

gives \(2\times2\times17=68\) cells per jet. Eight \(e_W\) jets \(\times68=544\).

Sampled addresses:

- `addressId` 7956: plus, pressure, `NATIVE_FLAT`, \(e_W\), `FORMAL`, \(\varepsilon\)-count 1, `flatSupport` \(k=l\), depths \(q(l)\), transfer \(l-p\).
- 8275: plus, pressure, `NATIVE_HEIGHT`, \(e_W_t\), transfer \(k-p\), depths \(q(k)\) and \(q(l)\), measure \(dl\,dk\), `EXACT_ZERO_SOURCE_JET`.
- 8818: plus, normal, `NATIVE_SLOPE`, \(e_W_{d2d2}\), map `plus-(0,1)-NATIVE_SLOPE`, \(q(k)\) and \(q(l)\) both present.
- 8892: plus, normal, `NATIVE_FLAT`, nonconstant \(\tanh^5\) source, `normalMultiplier` \(I\,q(l)\).
- \(\sim70032\): `minus-(0,0)-NATIVE_FLAT`.
- Last record: `EXACT_ZERO_CONSUMER`.

Zeros are kept. `globalComposition: UNRESOLVED` is the frozen inventory tag; global bounds are restored from the later certificates.

`selected/local-cells.json`: `selectedCount` 16, predicate `row == 'THETA_BALANCE' and field == 'e_W'; all xOrder and grade`. Headers include \((x,g)=(0,00),(0,10),(0,01),(0,11),(1,00),(1,10),(2,01),(3,11)\). Several cells are identically zero (`epsilonPower` 0); those stay.

---

## Flat, contact/PV, H, J, direct, depths, types

Fourier convention (`method.md` §2, `pressure/Fourier-order.json`):

\[
\hat f(k)=\frac1{2\pi}\int e^{-ikx}f(x)\,dx,\quad
X=\widehat{b D_j u}(k),\quad
Y=2\pi\widehat{c v}(-l),
\]

so \(F=\delta(l-k)\) returns \(\int c\,v\,b D_j u\,dx\). Coefficient \(\times\) packet is transformed; \(\hat b\) is not used as a free distribution.

On flat support, pole, jet, and normal factor reduce at \(k=l\). Off-diagonal HEIGHT/SLOPE keep source \(p\), input \(q(k)\), output \(q(l)\).

Height action (method §3), with \(f_k(Q)=B(k+Q,k)Y(k+Q)\):

\[
\mathrm{height}(k)=\frac W4 f_k(0)+\frac W{2i}\int_0^\infty A(Q)\frac{f_k(Q)-f_k(-Q)}{Q}\,dQ.
\]

Saved chi form (`background/pressure-weak-method.md` (7), `pressure/plus-height-PV.json`):

\[
\frac W4 f_k(0)+\frac W{2i}\int A(Q)\frac{f_k(Q)-f_k(0)\chi_{|Q|\le1}(Q)}{Q}\,dQ.
\]

\(A\) even and \(A\chi/Q\) odd imply the two agree as principal values. The worker must derive that before use. \(Q=z^2\) at the origin; branch coincidences stay Hölder-\(1/2\), not a fabricated \(f'(0)\). Normal height is bounded directly (`qTimesTestAssumedC1: false`). Plus/minus normal signs are \(+I q(l)\) and \(-I q(l)\) (`pressure/plus-height-PV.json`, `minus-height-PV.json`).

H uses the saved contact \(j(Q)/4\) and subtracted integrand \(A(t)[j(Q-t)-j(Q)]/(2i t)\) (`source/pressure-weak.py` joins). The product check against \(h(x)=W(1+\tanh(x/L))/4\) and \(j(x)=\mathrm{sech}^2(x/L)/4\) is a new adapter join; mismatch stops; it does not replace saved H.

J is native iteration once (middle \(k+t\) only). `Dwhole` is the closed direct density once (middle \(k+t\), reflected \(l-t\)); it is not sent through another middle integral. `pressure/typed-direct.json`: `Rprod=qi*qo*E`, `barePlaceholderIsWholeTag: false`, `multiplyWholeTagByResolvents: false`. `pressure/whole-definitions.json`: H, Jwhole, Dwhole joined, each once.

---

## Real outgoing quadrature

\(q\) is positive real for \(|p|<\kappa\), positive imaginary for \(|p|>\kappa\), continuously \(0\) at \(\pm\kappa\). Sources stay at real \(\omega=3\); no finite-\(\delta\) extrapolation.

Middle splits: \(t=\pm\kappa-k\), \(t=l\pm\kappa\), removable \(t=0\) and \(t=l-k\). Square maps at endpoints; coincident nodes merge with ancestry; near-coincident nodes stay separate unless coordinates alias. No omitted hole. Removable \(A\) by sinhc. A \(0/0\) depth is not set to \(0\); it needs a saved closed limit or an unsampled endpoint with an integrable transformed limit (`typed-direct.json` `endpointsNotAssignedPointwise: true`).

Outer \(k,l\) split at \(\pm\kappa\). Height \(Q\) also splits at \(Q=|\pm\kappa-k|\). Collisions \(k=\pm l\) and external/inner endpoint coincidences are recorded. Inner H/J/GD tests at zero transfer, grazing approaches, reflected momentum, and \(k=\pm l\) use closed/two-sided limits and do not clear the packet action.

---

## Fourier envelopes and tails

Strip \(|\mathrm{Im}\,x|\le5\): nearest \(\tanh(x/10)\) poles at \(|\mathrm{Im}\,x|=5\pi\); \(|\mathrm{Im}(x/10)|\le1/2<\pi/4\) gives \(|\tanh|\le1\). Shift \(x\mapsto x-i5\,\mathrm{sign}(\nu)\) produces \(\exp(25/(2s^2)-5|\nu|)\). Bound \(P(\tanh)\) by the absolute-coefficient sum; bound Hermite/Gaussian polynomials on the shifted line by absolute coefficients and translated Gaussian moments, including the \(1/(2\pi)\) and \(2\pi\) factors. Those inequalities are proposed certificates, not executed ones.

Outer tails: those exponential envelopes times inherited polynomial growth \(P^2\) (pressure) / \(P^3\) (normal) (`pressure/whole-envelopes.json`). Height uses the Hölder/PV bound and \(A\)’s exponential \(Q\) tail, not an \(L^1\) bound on the undistributed kernel. Middle: \(T\ge K+4\) beyond shifted branch endpoints (\(\kappa\approx2.44\), so \(K+\kappa<K+4\)); profile \(A\le11 e^{-|t|}\) and numerator envelopes in `J-numerator-envelope.json`, `D-height-envelope.json`, `D-reflected-envelope.json`. Old \(|H|\le100\), \(K_J\), \(K_D\) are existence bounds. The finite-\(K,T\) truncation inequality must be derived and saved before quadrature. If the envelopes cannot meet \(10^{-11}\) per unweighted grade/component under containment, that is a reported limitation, not a silent return to cutoffs \(4/6\).

---

## Two routes, tolerances, error meaning

Route A: contour-shifted Gauss–Hermite, orders 96/160. Route B: independent real-axis adaptive Fourier at 50 decimal digits, own subdivision on the original variables, own branch/PV subtraction, no shared nodes, weights, transform arrays, or response values. Inner Gauss 24/48. Panel cuts at the carrier and \(\pm n/s\), \(\pm n/L\), kept distinct from branch points. Adaptive inner target \(10^{-11}\) per addressed component. Failed library estimates stop.

Agreement, refinement, and \(K\to K+2\), \(T\to T+2\) use

\[
\tau=10^{-9}+10^{-7}\,|\text{refined reference}|
\]

in the pairing unit, with raw absolute differences reported. Analytic tails (conditional on the assessed inequalities) plus empirical route/refinement/enlargement differences form a declared conservative envelope. That envelope is not a flux or leakage bound. An unresolved component is not hidden in a cancelled total. No retry campaign.

Assembly consistency of grouped vs addressed sums is labeled as such.

---

## Numerical controls

Same selected pairings, metadata chosen before values, changed-minus-baseline:

1. drop surviving H contact (`NATIVE_MIXED_ITERATION`; prior formal selector used constant \(e_W\) source/consumer);
2. replace \(q(l-t)\) by \(q(k+t)\) on `INHERITED_DIRECT_WHOLE_OFF_DIAGONAL`;
3. replace \(q(l)\) by \(q(k)\) on an applicable normal slot (prior formal selector: `NATIVE_SLOPE`, nonconstant consumer);
4. Leibniz: move a nonconstant source coefficient inside a positive \(x\)-derivative (pressure and/or local \(x\)-order \(\ge1\); local \((1,00)\) and \((2,01)\) exist).

Silence in both packet settings is reported as unestablished. A nonzero numerical movement must exceed ten times the sum of empirical envelopes. Invalid-domain mutants are domain refusals. Prior end-equation agreement (`ends/correspondence-summary.json`) does not test lift-annihilated pressure entries; this \(e_W\to\) THETA pairing does (LEFT lift rows for \(\theta\) and \(e_W\) are zero). No current transfer.

---

## Optional improvements (not blockers)

These are worker-document polish, not a method change:

- Name the local pairing rule in one sentence (direct \(x\)-space \(\int v\,a\,\partial_x^n u\) on the same two routes, or \(Y X\) with \(F=1\)).
- Write one displayed inequality for \(C_{\mathrm{field,jet}}\) including the \(1/(2\pi)\) and \(2\pi\) factors.
- Name the exact-vs-float coalescing receipt (the method already requires ancestry and the union of intervals).

Exact H-product joins, numerical \(C_{\mathrm{field,jet}}\), \(K,T\), and control applicability are already runtime obligations with stop/limitation branches.

---

## Files inspected

Complete or essentially complete: `input/method.md`, `input/evidence-guide.md`, `input/review-prompt.md`, `input/packet-index.json`, `input/physical-input.json`, `input/local/assembly.json`, `input/ends/correspondence-summary.json`, `input/ends/left-match.json`, `input/background/pressure-weak-method.md`, `input/background/full-weak-method.md`, `input/background/numerical-route-inventory.json`, all `input/pressure/*.json` listed in the index except the two large inventories below.

Partial (structure, counts, headers, and targeted addresses/cells; not every record):

- `input/selected/pressure-addresses.json` (114780 lines). Sampled FLAT/HEIGHT/SLOPE, plus/minus, pressure/normal, FORMAL and both zero statuses. MIXED/DIRECT not line-sampled; they are required by the 544 count and the inherited grade-route contract.
- `input/selected/local-cells.json` (~5058 lines). Eight of sixteen \((x,g)\) headers plus several polynomials/zeros; not every cell body.
- `input/local/context.json`: start and end (frequency override, `effectiveSpeedOnlyInPressure`).
- `input/source/composition.py`, `pressure-weak.py`, `full-weak.py`, `weak-ends.py`, `legacy-integration.py`, `legacy-finite.py`: defining formulas, grade/component/depth contracts, controls, and legacy \(4/6\) Gaussian route. Not every helper line.

Not read: anything outside this packet, including peer reports (none are included). No execution, CAS, or peer communication.