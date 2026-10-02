I'll read the method, evidence guide, and packet index first, then inspect only the source and JSON named in this packet.The proposal is a bounded inventory, so I'll check it against the saved native routes, pressure census, and Fourier operands rather than treating the write-up as a result.The Fourier convention and face routing are the load-bearing parts, so I’m checking those against the saved operands and the original c2 route.The review is still open on the pressure census, the plus/minus source-input difference, and whether the flat-diagonal momentum binding is a real blocker. I will inspect those records directly and then give one verdict.The plus and minus source inputs differ only in the unevaluated face labels. I am checking the pressure-occurrence census and the native flat-kernel momentum binding next.NEEDS REVISION

The ordered rectangle, both-face source join, pressure census, and the refusal to claim a grazing limit or a cutoff are sound. One routing rule is not. As written, the normal-slot flat address does not match the saved jet.

## Blocker

**Location.** `method.md` §2 and §4–§5, the single rule \(N_{f,h}=i f q(l)\,F_{f,h}\) and the single route \(\int dl\,dk\,\hat c(r-l)\,F(l,k)\,\hat s(k-p)\,d_\alpha(p)\). Saved operands: `reference/retained-response-census.json` (`flat`, `jetKernels`, `taggedTotalMixed`) and `native/function-routes.json` (`kernel_apply`, `build_face`). Class label: `continuation/grazing-limit-objects.json`, flat = diagonal delta.

**What the saved records actually are.** The pressure flat coefficient is a rational function of `reference_qi` only:

\[
\frac{\omega}{10\bigl(\mathrm{reference\_qi}+\omega/(1000(-i\omega/1000+1/100))\bigr)}.
\]

The saved plus jet flat is the same function of `reference_qo`, with prefactor \(i\,\mathrm{reference\_qo}\). The minus jet flat carries the extra minus. Height, slope, and mixed iteration do match algebraic multiplication by \(i f\,\mathrm{reference\_qo}\) with their existing poles left alone. Flat is the exception: \(i q(l)\) times the saved flat text leaves the pole at `reference_qi`, while the saved jet pole is `reference_qo`.

Native `kernel_apply` does not put that diagonal piece in the off-diagonal double integral. `p0` multiplies the diagonal coefficient by the c1 factor `DiracDelta(k-k')` and integrates output momentum and the source point only. `build_face` forms the jet by differentiating \(\exp(i f q_o(N-N_{\mathrm{ref}}))\) against the reference diagonal, which attaches \(i f q_o\) and does not by itself rewrite an independent `reference_qi` pole.

**Consequence.** An instrument that follows the written formula on the saved flat text emits a normal-slot flat address with prefactor \(q(l)\) and pole `reference_qi`, or, if it drops the delta and integrates the rational function in \(dl\,dk\), a different operator. The §6 control that replaces \(q(l)\) by \(q(r)\) moves only the explicit multiplier, so a frozen `reference_qi` pole still passes. Replacing \(q(l)\) by \(q(r)\) is a valid control for height, slope, and mixed. It is not a flat-diagonal control.

**Minimum correction.** For the flat piece only, record the saved rational coefficient times the native diagonal delta, bind that coefficient’s momentum to the diagonal point, and write both the prefactor and the pole as \(q(l)\), with \(f=+1\) on the plus face and \(f=-1\) on the minus face, so the address reproduces `jetKernels.flat`. Keep `reference_qi` and `reference_qo` distinct in height, slope, and mixed iteration. Take the direct whole convolution from the separately tagged direct operand once. `taggedTotalMixed` and `jetKernels.*.mixed` already contain that direct addend, so they are sum-checks, not extra summands. This is a distributional identification already named by the flat class. It does not require evaluating a tanh integral or inserting a cutoff.

## What the saved operands already support

These do not need a formula change.

- The rectangle \(G=\{(0,0),(1,0),(0,1),(1,1)\}\) has 9, 3, 3, and 1 triples. Direct \(F_{(1,1)}\) sits only with source \((0,0)\) and consumer \((0,0)\). `directMultiplicity` is 1, and the whole convolution is not passed through the second slot’s middle integral.
- Plus and minus source amplitudes are the same: velocity \(e_{W,t}/2\), chemical coefficient \(10(1+3i/10)/(109(1+\eta w))\). The source-input hash split is only the unevaluated labels `s11cc1_mu_theta_lab_held_{plus|minus}` and `s11cc1_V_lab_held_{plus|minus}`. Join the evaluated amplitudes; keep those raw labels distinct until that join. There is no extra outward-normal sign on the velocity.
- One source division by `epsilon_shape`. At \(\eta=\sigma=0\) the live denominator \(1+\eta w\) is 1. Consumer rows keep a single `epsilon_shape`. U0/U1/U2 have 517 children and pressure counts 0. THETA has 380 children and one degree-1 child in each of the four slots. E_W has 989 children and two degree-1 children in each slot. \(\Lambda_{X,0}=0\) drops the memory child only after both children are stored. Normal-slot consumer signs are row data: THETA \(d_w\) plus is negative, E_W non-memory \(d_w\) plus is \(+\epsilon\eta w/4\).
- Consumer \((1,0)\) is a real nonzero control on the THETA and E_W normal slots. Consumer \((0,1)\) and \((1,1)\) are absent from the pressure children; record that absence. U rows are not a responsive control. Source cross grades remain inside the chemical polynomial and the factor \(1/(1+\eta w)\). The old kin\(=0\), kout\(=1/10\) contraction, and the upper-face `delta_p_plus` omission, do not supply those grades.
- The reduced transform \(\hat b(t)=(1/(2\pi))\int e^{-itx}b(x)\,dx\) matches the edge contract: constants become \(b\,\delta(t)\), two edge deltas are already applied, and the remaining measures are plain one-dimensional \(dk\,dl\). Momenta \(p,k,l,r\) stay distinct, the source derivative multiplies the wave at \(p\), and the middle momentum stays off that list. Products of profiles stay whole transforms.
- Trace matrices agree on both faces because the minus height sign and the minus jet \(-i q_o\) cancel. Reference pressure is already the solved reference value. Effective \(c_s\in[1,2]\) is a response-depth family. The census pole \(\omega/(1000(-i\omega/1000+1/100))\) and the direct pole \(\omega/(100(-i\omega/100+1/10))\) are the same rational function of \(\omega\). No \(c_s\) symbol is present in the source or consumer coefficients. Sections 1 and 5 do not claim a composed grazing limit, a global action, or a projection onto \([-3,3]\).

## What a corrected bounded instrument establishes

A complete unevaluated address list: every face, slot, channel, and all 16 grade triples, including explicit zeros; native child hashes; one epsilon on the consumer and none on the source; rectangular source and consumer coefficients in the quotient by discarded \(\eta^2\) and \(\sigma^2\); the flat diagonal bound to \(q(l)\); direct multiplicity one; and each global-momentum obligation marked unresolved. That is an algebraic and distributional inventory. It is not a value of the composed operator.

## Still open before a finite near-unity defect pilot

The response certificate covers internal \(k,l\in[-3,3]\) on a fixed \(C_c^1((-3,3))\) output test. Nonconstant \(w\) and \(m\) move momentum outside that interval. Grazing integrals are unevaluated, endpoints are not pointwise, and the continuation proof does not give \(q(l)\phi\in C^1\). No finite matrix, quadrature, cutoff, or mode solve is licensed by this inventory. \(\sigma=\eta W/L\) stays out until after grade selection. The rectangular Taylor of \(1/(1+\eta w)\) and the jet-by-jet source split are still to be derived from the saved polynomials; the published `higherGrades` blob is not that split.

## Runtime evidence the instrument has to emit

Under the guarded runner, with no call to `Inputs`, `build_face`, `build_case`, or the old response producers: exact reconstruction of each retained pressure row from independent placeholders; denominator values at \(\eta=\sigma=0\); pressure-name counts against the stored censuses; a symbol-level identity of each flat normal address with `jetKernels.flat`; separated-piece direct count one; control movements for source and consumer cross-grade omission, lower-face jet sign, \(p\to k\), and \(q(l)\to q(r)\) on height, slope, or mixed, each on a nonzero row before any face sum; and an explicit unresolved mark wherever a profile factor leaves \([-3,3]\). Source text and this review are not that evidence.

## Tooling, not physics

`physical-input.json` still holds historical \(\omega=1\) and \(c_{s0}=10\); the composition uses \(\omega=3\) and \(c_s\in[1,2]\) only through the saved response records. Some consumer files predate a formatter failure; their hashes are intact and they are not a finished top-level return. `native/consumer-census.json` is pinned constructor text. The inherited velocity coefficient `1` is the c1 prefactor, distinct from the amplitude \(e_{W,t}/2\). No wording pass is required beyond the flat-address correction above.