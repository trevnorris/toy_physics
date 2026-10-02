I'll start from the frozen packet index and evidence guide, then inspect the named source and method files only.The index and guide are in place. Next I’ll read the method and the saved factorization, closure, and physical-input files they point to.The factorization is on the page. I’ll now check the Fourier weights, both-face closures, and the beta/sheet records against the proposed rearrangement.Closure and iteration residuals match so far. I’ll now inspect the lower-face mirror, uniform matches, grade-zero rows, and the Fourier construction in the worker.I’ll pull the exact grazing matches, the worker’s Fourier/A(t) construction, and the native outgoing/beta law next.The saved both-face direct mixed term has a source-bound locally integrable closed limit at the two selected matches. Identity (1), the native `Omega=3+i delta` continuation, and the moving-endpoint `L1` argument check against the frozen operands.

## Verdicts

**(A) Source/evidence applicability: SUPPORTED WITH STATED LIMITS**

**(B) Proposed method: CLEAR FOR BOUNDED CLOSED-GRAZING IMPLEMENTATION**

This is a method clearance for a later join/limit instrument. It is not an implementation clearance, an executed convolution, or a physical-loss result.

---

## 1. Ordered density, factors, faces, rearrangement (1)

Rearrangement (1) follows from the saved operands.

`raw/physical-factorization-return.json` gives

\[
B=-i\omega\rho_m\left(\frac{k(-H+2Q+2k)}{q_i q_o(q_o+q_s)}+\frac{k(H+2k)}{q_h q_o(q_h+q_i)}+\frac{q_i}{q_h q_o}\right)
\]

on the regular domain, with exact residual `0` after the two dispersion identities in `raw/factorization-before-reduction.json`. With \(t=H\), \(l=k+Q\), this is the method’s \(B\). Both-face contact at \(t=0\) is `0` (`raw/upper-complete-contact-return.json`, `raw/lower-complete-contact-return.json`).

The native two-leg factor on both faces is

\[
R=\frac{q_i q_o}{(q_i+\beta)(q_o+\beta)},\qquad \beta=\frac{30+9i}{109}
\]

(`raw/plus-closure-operands.json`, `raw/minus-closure-operands.json`). Multiplying \(B\) by \(R\) produces (1) termwise:

- \(k(2l-t)/(q_i q_o(q_s+q_o))\cdot q_i q_o=k(2l-t)/(q_s+q_o)\)
- \(k(t+2k)/(q_h q_o(q_h+q_i))\cdot q_i q_o=k(t+2k)\,q_i/(q_h(q_h+q_i))\)
- \(q_i/(q_h q_o)\cdot q_i q_o=q_i^2/q_h\)

`raw/plus-external-depth-cancellation-return.json` and the minus copy already give residual `0` for \(R\cdot C\) against the same numerator over \((q_i+\beta)(q_o+\beta)\). A later instrument still has to emit the residual of (1) against each saved `*-closed-raw-increment.json` physical density.

The ordered ordinary density matches `raw/raw-ordered-before-cancel.json`:

\[
\frac{W}{2}\,B\,A(t)\,\frac{L A(Q-t)}{2}\,\frac{1}{i}=\frac{WL}{4i}A(t)A(Q-t)B
\]

with \(W=1\), \(L=10\), \(A(t)=Lt/(4\sinh(\pi L t/2))\), and \(A(0)=1/(2\pi)\). That \(A\) is the saved transform join in `source/raw-worker.py` (the `saved-transform-A` identity). The uncancelled contact/PV form is the worker’s `uncancelledHeight` string; it stays alongside the factored density. The reduced 1-D factor \(1/(2\pi)\) is `raw/edge-delta-reduction.json`, consistent with \(A(0)\).

Face conventions match the saved lower derivation: the scalar mixed coefficient is equal on the two faces (`raw/lower-upper-outward-mirror-return.json` residual `0`), while normal jets are \(+iq_o\) and \(-iq_o\) (`raw/plus-trace-domain.json`, `raw/minus-trace-domain.json`). Zero-grade reference trace is `1`. First-shape iteration is left once: `raw/plus-iteration-unchanged-return.json` and the minus copy give `closed[0,2]-closed[0,2]|_{D=0}=D\cdot R`.

The permeability/memory law is the saved native definition

\[
I+\frac{\Lambda_{A0}}{\rho_m^2(1-i\omega\tau_A)}Z,
\]

with live \(\Lambda_{A0}\), \(\tau_A\), and \(\omega\), evaluated at the construction frequency 3 to \(a=(100+30i)/109\) and \(\beta=(30+9i)/109\) on both faces.

---

## 2. Continuation, \(\beta\), roots, external cancellation

The proposed \(\Omega=3+i\delta\) law is the same native outgoing prescription named as uncomputed in `raw/physical-sheet.json` (`limiting omega+i0+ prescription remains uncomputed`). Depths use the first-quadrant square root of \(\Omega^2/c_s^2-1/20-p^2\). As \(\delta\to 0^+\) this recovers the saved piecewise sheet: positive real on a positive radicand, positive imaginary on a negative radicand (`raw/physical-depth-bindings.json`). No conjugated current or parameter swap is introduced. Independent \(\eta\) and \(\sigma_W\) stay independent.

Continuing the same \(a(\Omega)\) and \(\beta(\Omega)=a(\Omega)\rho_m\Omega\), with \(D_\delta=(1+\tau_A\delta)^2+9\tau_A^2\),

\[
\operatorname{Re}\beta=\frac{\Lambda_{A0}}{\rho_m}\frac{3}{D_\delta}>0,\qquad
\operatorname{Im}\beta=\frac{\Lambda_{A0}}{\rho_m}\frac{\delta(1+\tau_A\delta)+9\tau_A}{D_\delta}>0.
\]

The identity \(\operatorname{Re}(\Omega\overline{(1-i\Omega\tau_A)})=3\) is exact, so \(\operatorname{Re}\beta\) is independent of the \(\tau_A\delta\) cross terms. For \(0<\delta\le 1/10\) and \(\tau_A=1/10\), \(D_\delta\le 11101/10000\), hence \(\operatorname{Re}\beta\ge 3000/11101=\beta_{\min}\). All four depths and \(\beta\) lie in the closed first quadrant, \(\beta\) strictly interior. For \(u,v\) in that quadrant, \(|u+v|\ge\max(|u|,|v|)\), so \(|q_i+\beta|\ge\beta_{\min}\) even when \(q_i\to 0\). External cancellation remains multiplication by \(q_i q_o/((q_i+\beta)(q_o+\beta))\) with live \(\beta(\Omega)\).

---

## 3. Compact domain, envelope (2), collisions, \(\kappa=0\)

The declared envelope \(c_s\in[1,2]\), \(k,l\in[-3,3]\) contains both saved matches in `uniform-context.json`:

| match | \(c_s\) | \(k\) | \(\kappa^2=9/c_s^2-1/20\) |
|---|---|---|---|
| LEFT, index 1 | \(\sqrt{6}/2=\sqrt{3/2}\) | \(\pm\sqrt{595}/10\) | \(595/100\) |
| RIGHT, index 8 | \(5\sqrt{606}/101=\sqrt{150/101}\) | \(\pm\sqrt{601}/10\) | \(601/100\) |

Forward/backward pairings \(l=\pm k\) stay in the rectangle (\(|Q|\le 6\)). Nearby radiating and evanescent rows in the same file also sit inside the envelope.

On this set, \(|\Omega|\le 4\) and \(|q_i|,|q_o|\le 5\) hold (actual maxima are about \(3.01\) and about \(3\)). Combined with \(|q_i/(q_h+q_i)|\le 1\) and \(1/|q_s+q_o|\le 1/|q_s|\), this gives envelope (2) exactly as written. \(\kappa_\delta=\sqrt{(9-\delta^2)/c_s^2-1/20}\) is minimized at \(c_s=2\), \(\delta=1/10\), so \(\kappa_\delta\ge\sqrt{879}/20=\kappa_{\min}>0\). The comparison

\[
|q_\delta(p)|\ge\sqrt{|\kappa_\delta^2-p^2|},\qquad |\kappa_\delta^2-p^2|\ge\kappa_{\min}\operatorname{dist}(p,\{\pm\kappa_\delta\})
\]

is valid, because \(\max(|p-\kappa|,|p+\kappa|)\ge\kappa\ge\kappa_{\min}\). Inverse depths are therefore bounded by a *sum* of two inverse square roots. The elementary estimate \(\int_E |t-a|^{-1/2}\,dt\le 2\sqrt{2m}\) is uniform in the moving endpoint \(a\).

Collision, including \(q_h=q_s=0\): the majorant remains a sum, not a product. At simultaneous external grazing \(l=k=\kappa\) (or \(l=-k\)), the point \(t=0\) is a common simple zero of \(q_h\) and \(q_s\), so \(|B_c|\sim |t|^{-1/2}\), which is integrable. The same holds at \(t=\pm 2\kappa\).

The excluded \(\kappa=0\) degeneracy is the merger of \(\pm\kappa\), where \(1/|q|\sim 1/|p|\) fails local \(L^1\). It is the right excluded case. It does not occur on this envelope (\(\kappa\ge\sqrt{11/5}\) at \(\delta=0\)) and is irrelevant to both selected matches (\(|k|\approx 2.44\)).

Tails: every comparison endpoint lies in \([-6,6]\). Outside a large compact interval, \(A(t)A(Q-t)\) decays exponentially while (2) stays \(O(1)\). An explicit \(T\) and constant belong in the instrument; existence of a uniform integrable tail does not require a numerical integral.

---

## 4. \(L^1\) limit, radiating/evanescent paths, contact/PV

Uniform absolute continuity on compact \(t\)-intervals, uniform tightness of the tails, and pointwise convergence off the finite limiting endpoint set give \(L^1(\mathbb{R})\) convergence of \(G\) along every path in the envelope, including \(\delta\to 0\) and \(q_i,q_o\to 0\) together. That is a moving-singularity uniform-integrability argument. It yields a unique \(L^1\) class, hence a unique finite pairing against any bounded test function. It does not yield a pointwise-regular kernel, differentiability in \(c_s\), or an operator-norm bound.

Almost-everywhere limits (3) and (4) are the correct pointwise formulae off the endpoint set. At a finite set of \(t\) values the density is left unassigned. Radiating and evanescent approaches are both covered: the comparison uses \(|\kappa_\delta^2-p^2|\), which is the same majorant on either side of the cut. The joint \(\delta\)/external limit is included because the majorant is uniform in \(\delta\) and in \((c_s,k,l)\), while at finite \(\delta\) one has \(\operatorname{Im}(q^2)=6\delta/c_s^2>0\), so the four real depths never vanish.

Contact/PV: at every \(\delta>0\) the four-leg identities still give \(C(t=0)=0\) (they are algebraic in the depths). The contact piece \((W/2)(\delta(t)/2)\,j_{\hat{}}(Q)\,C(0)\) is therefore zero before the limit. Uniform integrability of the ordinary density \(G_\delta\) prevents a new Dirac mass from forming: \(\int_E |t|^{-1/2}\) vanishes with \(m(E)\), uniformly in the endpoints. A surviving \(\delta(t)\) contribution is incompatible with that bound. A pointwise evaluation at \(t=0\), or separately cancelled summands, is not a substitute for this \(L^1\) statement.

The separate jet \(if q_o G\) inherits the same \(L^1\) conclusion: \(q_o\) is independent of \(t\) and bounded on the envelope. At \(q_o\to 0\) the jet is even smaller.

---

## 5. Reference/jet and grade-zero reuse; execution guards

Scope is preserved if the instrument restores saved operands and does not rebuild them.

Already saved and reusable:

- both-face \(R\), reference, and jet (`raw/*-closed-before-cancel.json`, `raw/*-trace-domain.json`)
- iteration-once residuals
- grade-zero sources `sourcePlus`/`sourceMinus` in the retained increments, with finite `denominatorAtZero = 90900+303000i` (`raw/plus-source-zero-grade-operands.json`)
- U-body pressure census empty, mixed coefficient `0` (`raw/U0-executed-native-census.json`, `raw/U0-retained-increment.json`)
- THETA and E_W affine pressure degree 1, both faces present (`raw/THETA_BALANCE-executed-native-census.json`, `raw/E_W_BALANCE-executed-native-census.json`)
- consumer symbols independent of \(c_s\) (`raw/effective-speed-reuse-domain.json`)

Required execution guards, not already proved as closed-grazing certificates:

1. Exact residual of (1) versus each saved `physical: raw*R` density in `raw/plus-closed-raw-increment.json` and `raw/minus-closed-raw-increment.json`.
2. Actual THETA/E_W slot coefficients of \(D_\pm\) after the saved replacements; confirm those multipliers are independent of \(t\) (they may be polynomials in external \(k,l\)). Unknown or \(t\)-dependent joins stop.
3. Contact identity emitted at continued \(\Omega\), not only at real frequency 3.
4. Domain lemmas: \(\beta_{\min}\), \(\kappa_{\min}\), \(|q_i|,|q_o|\le 5\), endpoint enclosure in \([-6,6]\), and an explicit tail.
5. Construction frequency **3** as in `raw/physical-sheet.json` and the closure numbers \(3/10\), \((100+30i)/109\). `physical-input.json` still lists `"omega": "1"`; that listed value is not the saved kernel frequency.

These guards do not invent an incoming transverse excitation or a finite-deficit theorem. `source/finite.py` remains a reason to keep this term separate from an untruncated finite inverse. `finish/applicability.json` still records `exactMatch: UNRESOLVED`; that is the gap this method is meant to close.

---

## 6. Instrument evidence and controls

The proposed acceptance list is adequate for this bounded object: per-face joins of (1), live \(\beta\) sign/lower bound, four-leg root/domain certificates, simple-root moving-endpoint bound with explicit tail, qi-only / qo-only / simultaneous \(\pm k\) almost-everywhere limits, contact/PV \(L^1\) argument, and separate pressure / jet / grade-zero consumer joins. Operands and residuals before guards; a sampled approach to grazing cannot replace the bound; no convolution evaluation and no finite matrix in this job.

Responsive controls, on the *closed* expressions with the same source bindings:

- drop the external factor \(R\) (or freeze \(\beta=0\))
- flip a root sheet \(q\to -q\)
- flip the lower jet sign \(if q_o\)

Those are the same defect classes already shown movable in `finish/reuse-wrong-sheet-response.json`, `finish/reuse-wrong-lab-face-response.json`, and `finish/lower-reference-jet-sign-return.json`. Those old outputs remain provenance, not independent correctness.

**Precise corrections (instrument obligations, not method blockers):** join (1) to the saved closed densities, not only to \(B\cdot R\) in isolation; attach the ordered-density jet factor \(L A/2\) from `raw/raw-ordered-before-cancel.json` (the worker’s `saved-transform-jet` identity is \(LA\) before that geometric \(1/2\)); pin frequency 3 from the saved closure/sheet, not from `physical-input.json`’s listed `omega`.

**Optional presentation:** record that \(|\Omega|\le 4\) and \(|q_i|,|q_o|\le 5\) are loose; state \(\beta_{\min}\) as \(\min\operatorname{Re}\beta\) over \(0<\delta\le 1/10\); note that `plus-saved-source-zero-grade-input.json` displays identical copies, while the actual grade-zero join is `source0` against saved `source00` in the worker.

---

## Strongest claim after these source-bound checks

After exact joins of (1) to the saved both-face closed densities, the direct mixed kernel of grade \((1,1)\) with grade-zero source/consumer, at LAB_HELD/RHO4_CONSTANT, frequency 3, edge momenta \((1/5,1/10)\), live \(a(\Omega)\) and \(\beta(\Omega)\), admits a unique \(L^1(\mathbb{R})\) limiting ordinary density on the envelope \(c_s\in[1,2]\), \(k,l\in[-3,3]\). That limit exists along radiating, evanescent, and joint \(\delta\to 0\) / external-grazing paths, including the two selected matches \(c_s=\sqrt{3/2}\), \(k=\pm\sqrt{595}/10\) and \(c_s=\sqrt{150/101}\), \(k=\pm\sqrt{601}/10\), and their \(l=\pm k\) pairings. The almost-everywhere density is \(G=(WL/(4i))A(t)A(Q-t)B_c\) with \(B_c\) as in (1); no contact Dirac mass survives; convolution against bounded test functions is finite and path-independent. The same \(L^1\) class passes to the native jets \(\pm i q_o G\) and to THETA/E_W after the saved grade-zero slot joins.

**Exclusions:** this is not an executed integral, a finite-matrix result, a loss or leakage smoothness theorem, a full Green function, a first-shape re-iteration, a pure \(\eta^2\) or \(\sigma_W^2\) block, a primitive-speed calibration, a defect sweep, a pointwise-regular kernel, differentiability in \(c_s\), a global operator-norm bound, an incoming transverse excitation, a finite-deficit protection, the \(\kappa=0\) bulk-edge degeneracy, or an impermeable \(\beta=0\) limit. `finish/` sensitivity certificates and author success flags are not physical acceptance.
