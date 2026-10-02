I'll assess the one-dimensional both-face mixed-increment method from the staged sources, starting with the five construction sections and the packet index, then checking only the exact source and evidence those sections use.The method and index are in hand. Next I’ll read the review request and the saved coefficient, geometry, transform, and closure operands those five rules actually join.The ordered contact, sheet, and face rules are supported on the domains below. The next increment is clear for a bounded nongrazing implementation. The lower mixed coefficient, the closed grazing limit, and any integral value remain objects that implementation has to construct and test.

- ORDERED CONTACT/SHEET/FACE RULES: SUPPORTED WITH STATED DOMAINS
- NEXT INCREMENT METHOD: CLEAR FOR BOUNDED IMPLEMENTATION

No substantive recipe error turned up. The checks below are hand algebra and source joins against the saved operands. No producer was run, and no zero guard is treated as acceptance.

## 1. Ordered transfer and normalization

The saved upper pressure coefficient in `bare/boundary-coefficient-return.json` is the ordered object the method uses:

\[
C_+(k;H,S)=i\omega\rho_m\frac{N}{q_h q_i q_o},\qquad N=-H q_i^2+k q_h q_i+k q_h q_o-k q_h q_s-k q_i^2,
\]

with \(q_h=q(k+H)\), \(q_s=q(k+S)\), \(q_o=q(k+H+S)\). The convolution pairs \(H=t\) with \(h_\mathrm{hat}(t)\) and \(S=Q-t\) with \(j_\mathrm{hat}(Q-t)\) once. \(N\) is not symmetric in \(H\) and \(S\), so averaging \(C\) alone changes the integrand. The half-weight \([F(t)+F(Q-t)]/2\) on the whole integrand, hats included, is the substitution \(t\to Q-t\) and preserves the integral.

Native `shape_source` is \(h=\eta W_0 w_{1,\mathrm{hat}}/2\) and tilt \(=\sigma_W\,\mathrm{jet}/2\). The profile binding is the normalized forward transform \(\mathrm{e}^{-i\Delta k\cdot y}/(2\pi)^3\). A profile that depends only on the profile coordinate therefore contributes one edge delta from each edge integral. Two such hats, integrated over the two middle edge momenta, leave the output-minus-input edge deltas. The profile-direction Jacobian of \(m=k+t\) is \(+1\).

`kernel_apply` is a different convention: the source is transformed without a \((2\pi)^3\) in the forward integral, and one \((2\pi)^{-3}\) sits on the inverse. That factor stays with the source. It is not part of \(D_f\). In the normalized-hat convention the product of two profiles is the plain convolution

\[
D_f(k_\mathrm{out},k)=\int C_f(k;t,Q-t)\,h_\mathrm{hat}(t)\,j_\mathrm{hat}(Q-t)\,dt.
\]

Hand dimension count of that reduced kernel is \([-2,-1,1]\): \(C\) carries \([-4,-1,1]\), the one-dimensional height hat \([2,0,0]\), the jet hat \([1,0,0]\), and \(dt\) contributes \([-1,0,0]\). That is the saved reduced kernel dimension, obtained from the native kernel dimension \([0,-1,1]\) by removing two edge deltas of dimension \([1,0,0]\) each.

At \(k=0\) the factored raw integrand collapses to the saved selected integrand. With \(B_+=-i\omega\rho_m q_i/(q_h q_o)\), \(j_\mathrm{hat}=L A/2\), and \(1/i=-i\),

\[
\frac{W_0}{2}B_+\frac{A(t)}{i}j_\mathrm{hat}(Q-t)=-\frac{W_0\omega\rho_m q_i L}{4 q_o}\frac{A(t)A(Q-t)}{q(t)},
\]

which is the saved prefactor in `bare/profile-action-return.json`. No missing \(2\pi\), \(1/2\), or extra transfer is required on this join. The inherited whole \(D\) stays a comparison symbol: `kernel_apply(..., second=...)` already integrates the middle momentum, so that symbol must not be placed there.

## 2. Complete contact and its domain

The contact and principal-value split follows the height transform recorded in `bare/profile-transform-return.json`:

\[
w_\mathrm{hat}(t)=\frac{\delta(t)}{2}+\mathrm{PV}\!\left[\frac{A(t)}{it}\right],\quad A(t)=\frac{L t}{4\sinh(\pi L t/2)},\quad A(0)=\frac{1}{2\pi}.
\]

Hence \(h_\mathrm{hat}=(W_0/2)w_\mathrm{hat}\) has contact weight \(W_0/4\), and

\[
\mathrm{contact}_+=(W_0/4)\,C_+(k;0,Q)\,j_\mathrm{hat}(Q).
\]

The jet transform \(L A(u)\) is an ordinary function, so only the height delta produces a contact.

On one sheet, \(H=0\) forces \(q_h=q_i\) and \(q_s=q_o\). The whole numerator then cancels:

\[
N\big|_{H=0}=k q_i^2+k q_i q_o-k q_i q_o-k q_i^2=0.
\]

The cancellation uses the \(k q_h q_o\) term against \(-k q_h q_s\). It is not a property of any single summand. For nonzero \(q_i q_o\), the upper contact is zero. This is ordinary pointwise algebra while those external depths stay off zero.

The square identities \(q_h^2-q_i^2=-H(H+2k)\) and \(q_o^2-q_s^2=-H(2k+2Q-H)\) hold whenever \(q(p)^2=\kappa^2-p^2\) and \(S=Q-H\). Dividing by \(q_h+q_i\) and \(q_s+q_o\) reproduces \(N\) exactly, so

\[
C_+=H B_+
\]

with the stated \(B_+\), provided \(q_h,q_i,q_o,q_h+q_i,q_s+q_o\) are nonzero. On the positive sheet the two sums vanish only when a root vanishes, so those conditions sit inside the same exclusion. Where the external depths are regular, \(B_+\) is finite at \(t=0\) and \(t\,\mathrm{PV}(1/t)=1\) gives the ordinary density \((W_0/2) B_+ A(t) j_\mathrm{hat}(Q-t)/i\).

That identity is not a distributional product at a grazing collision. If \(q_i\), \(q_o\), or a factor sum hits zero, the unfactored contact and principal-value pair stay in force until a closed limit exists. The lower numerator is not this \(N\) with a sign typed in; it is still to be derived.

## 3. Sheet and both faces

For \(\omega>0\) and \(c_s>0\),

\[
\kappa^2=\omega^2/c_s^2-b^2,\qquad b^2=(1/5)^2+(1/10)^2,
\]

and \(q(p)\) is the positive square root when \(\kappa^2-p^2>0\) and \(+i\) times the positive square root when \(\kappa^2-p^2<0\). The same rule fixes \(q(k)\), \(q(k+t)\), \(q(k+Q-t)\), and \(q(k+Q)\). A zero radicand is the limit \(\omega\to\omega+i0^+\) at fixed \(c_s\). Solver order, and the symbol `MIDDLE_Q`, are not a substitute for that sheet. At the saved point \(c_s=10\), \(k=0\), this reproduces \(\kappa^2=1/25\), \(q_i=1/5\), \(q_o=\sqrt{3}/10\).

Native Eulerian geometry has \(h_0=f W_\mathrm{bg}/2\) and unit normal \((-f\nabla h_\mathrm{lab},f)/\sqrt{1+|\nabla h_\mathrm{lab}|^2}\). The reference continuation in `reference_pressure_kernels` is \(\exp(i f q(w-f W_0/2))\) with that same \(q\). Its \(w\)-derivative at the reference face is \(i f q\). Sending \(q\to -q\) as well would flip the lower root twice.

The linear joins cancel one face and leave the native face-independent first-shape kernel: lab displacement \(f h\) times derivative factor \(i f q\) is \(i q h\), and lab slope \(f s\) times normal factor \(-f\) is \(-s\). `dtn_first_kernel` has no extra face sign, and the default `shape_source` tilt is unflipped, so those linear pieces are the native check for both faces. The quadratic area factor starts at slope squared and does not create or remove the mixed grade. The minus pressure and normal-jet slots are the actual \(f=-1\) slots. In `continued/consumer-unit-joins.json` the value-slot coefficients agree across faces and the normal-jet coefficients flip, which is the same \(i f q\) pattern. An outward mirror is only a supplemental comparison. The lower mixed numerator itself is not in the packet; the rules above are what the implementation has to expand and then match to the native linear coefficients.

## 4. Closure, insertion, and consumers

Native `z_three[0,2]` is the constant 0, and the trace worker checks that before tagging a direct slot. The first-shape transfers stay on the \((0,1)\) and \((1,2)\) entries. One insertion at \([0,2]\), multiplied by \(\eta\sigma_W\), adds the direct grade without a second copy of the iteration.

For an upper-triangular factor \(I+a_f Z\), the direct \((0,2)\) entry of \((I+a_f Z)^{-1}Z\) is \(z_{02}/(D_\mathrm{out} D_\mathrm{in})\). With \(z_0(p)=\rho_m\omega/q(p)\) and \(D(p)=1+a_f z_0(p)\),

\[
R_f(p)=\frac{q(p)}{q(p)+a_f\rho_m\omega}=\frac{1}{D(p)}.
\]

So the bare direct kernel is closed by \(R_f(k_\mathrm{out})\,C_f\,R_f(k)\). Here \(a_f\) is the coefficient of \(Z\) in the inherited face resolvent, read for that face. The saved numerical factor \((60+91i)/(105+241i+\sqrt{3}(30+250i))\) is the value of this closure at one upper point. It is not the function of \(q\). With the supplied \(\Lambda_{V0}=0\) the velocity-channel factor is 1, but it still comes from the actual source expression, together with density, memory, and \(\epsilon\).

Each \(R_f\) supplies one power of \(q\). Forming \(R C R\) before a limit therefore removes the explicit \(1/q_i\) and \(1/q_o\) while the resolvent denominators stay nonzero. That cancellation is required before any pointwise external-grazing use, and it does not give an integrable collision of the branch point with \(t=0\).

The reference trace inverse and the multiplier \(i f q_o\) stay in the map. At the selected upper point the saved reference and physical factors happen to agree; the minus jet sign shows that agreement is not the general both-face map.

`retained_shape` keeps grades \((0,0)\), \((1,0)\), \((0,1)\), \((1,1)\). `finite.py` assembles the matrix and LU-solves the full block, with an SVD cross-check. A zero direct transverse column does not make a missing scalar entry irrelevant to that inverse. The conditional grade identity is not a substitute, and no separate formal-identity job is needed. Constant-end modes, currents, and trace maps can be reused only after the new both-face zero-jet joins exist. The selected U-row pressure census is zero on the three expanded U rows; theta and \(E_W\) still carry pressure and jet slots. The zero-grade jet coefficients vanish at the retained mixed grade, while the jet slots themselves remain.

## 5. Where the raw increment stops

For \(\kappa^2>0\) the middle branch points of the four depths lie at

\[
t=-k\pm\kappa,\qquad t=k+Q\pm\kappa.
\]

The second pair is a branch of \(q_s\) inside the complete coefficient. For \(\kappa^2\le 0\) the same formula is not a real square root; the domain and any zero endpoint have to be recorded as such. As positive \(k\) approaches input grazing, the colliding root is \(t_b=\kappa-k=q_i^2/(\kappa+k)\), of scale \(|q_i|^2/(2|k|)\). The partition uses that exact root. Negative profile momentum uses the root that actually approaches the contact.

The symbolic closed limit on both sides of the outgoing sheet, the shrinking-layer Jacobian, and an integrable local bound are acceptance conditions for exact-match use. They are not results supplied by this method. A nongrazing raw increment, with \(q_i\) and \(q_o\) nonzero and the branch points off \(t=0\), is the bounded output. Exact-match finite matrices stay blocked until every new contact, principal-value, and branch piece has an explicit integration prescription. If that prescription needs a larger endpoint or exterior method, the bounded task stops there.

## Next action

Build the task-local raw increment for this tanh profile only: both faces, live profile momenta, fixed edge \((1/5,1/10)\), \(\omega=3\), \(c_s>0\), rest bulk, `LAB_HELD`/`RHO4_CONSTANT`. Derive the lower coefficient from the face geometry above, emit the upper factorization residual and both-face linear joins to native c1, keep the contact and principal-value forms beside the factored density, reduce the ordered convolution once, and attach the actual \(R_f\), source factor, and reference map. Restrict the emitted domain to nongrazing external depths. Do not evaluate the middle integral, call the finite solve, replay a producer, or treat exact-match applicability as obtained.