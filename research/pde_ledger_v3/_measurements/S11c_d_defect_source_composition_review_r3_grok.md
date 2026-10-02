I'll review the bounded source-composition method from the private packet only, starting with the method, evidence guide, and packet index.The method is a bounded inventory proposal. I'll open the native routing text and the saved closure, slot, and density operands it depends on.The native `build_face` and `kernel_apply` text is in hand. Next I’ll check `kernel_bridge`, the flat-jet convention, and the saved closure operands those formulas join to.The review is still open: the closure coefficient and the closed density use different external factors, and the retained-increment rows, pressure-name census, and lower-face final slot have not been opened yet. I will read those operands before choosing a verdict.The resolvent identification is still the leading issue. I am finishing the pressure-name census and the section 4a operands before choosing the verdict.NEEDS REVISION

The ordered inventory is usable once one factor identification in section 4a is split. The saved operands already support a no-replay formal address inventory. They do not support treating the closed density as \(R(q_o)R(q_i)\).

## Blocker

**Location.** `method.md` section 4a, the definition \(R(q)=q/[q+\beta]\), the sentence that `direct/closed-density.json` already includes those two external resolvents, \(\Delta P=\eta\sigma\,R(q_o)\,D\,R(q_i)\), and the requirement that the routed real-3 reference and normal factors join the closed-density maps with a literal zero residual. `evidence-guide.md` repeats that the closed \(D_{\mathrm{whole}}\) already contains its resolvents.

**What the operands show.** At real \(\omega=3\), \(\beta=30/109+9I/109\).

- The augmented closure factor in `raw/plus-closure-operands.json` and `raw/minus-closure-operands.json` is \(q_i q_o/((q_i+\beta)(q_o+\beta))\). That is \(R(q_i)R(q_o)\). `direct/plus-reference-factor-return.json` prints a cancelled residual of \(0\) against that product. Both `raw/*-closed-before-cancel.json` files store the same unreduced reference factor, with jets \(+I q_o\) and \(-I q_o\) times that factor.
- `direct/closed-density.json` puts \(1/(q_i+\beta)\) and \(1/(q_o+\beta)\) outside a numerator and the two sinh factors. The prefactor has no \(q\)-resolvents. So the density’s external factors are \(R(q_i)R(q_o)/(q_i q_o)\). The same density is copied in `reference/retained-response-census.json` and `reference/restored-direct-argument-join.json`.

**Consequence.** A literal factor residual against the saved density cannot be zero. Following the written stop rule rejects a correct density, or an implementer multiplies \(D_{\mathrm{whole}}\) by \(q_i q_o\), or identifies the bare routing placeholder \(D\) with \(D_{\mathrm{whole}}\). Either move changes the direct address. The full numerator identity between the density and \(R\cdot\mathrm{raw}\cdot R\) is not established by these files and must stay unresolved.

**Minimum correction.** State three separate facts.

1. The linear coefficient of the independent bare symbol, and the saved closed-before-cancel reference factor, equal \(R(q_i)R(q_o)\), with jets \(\pm I q_o\) times that factor.
2. \(D_{\mathrm{whole}}\) is the saved closed density, used once. Its external factors are \(1/(q_i+\beta)\,1/(q_o+\beta)\). It is not rewritten as \(R\cdot\mathrm{raw}\cdot R\) and is not multiplied by \(R\) again. Do not require \(\mathrm{density\ factor}-R R=0\).
3. The formal placeholder \(D\) is not \(D_{\mathrm{whole}}\). \(\Delta P=\eta\sigma\,R(q_o)\,D\,R(q_i)\) is only the closure image used for the `build_face` routing check and the grade \((2,1)\) drop. That factor joins the closure and reference factors. The density join is the saved argument map plus the factor relation above.

## Checks that match the saved operands

Native `kernel_bridge` writes literal \(0\) at `z_three[0,2]` (`source/native-c2.py` line 409). `second` is the later \([0,2]\) entry of the triangular image, then the \(\eta\sigma\) rectangle. The direct correction is only the augmented `*_raw_direct` entry. Section 4a correctly treats the new route as a formal adapter through `reference_matrix[0,1]`, with `jet_diagonal=jet_second=0`, no middle integral, and no call of `kernel_apply`, `reference_pressure_kernels`, or `build_face`.

Both new traces have `valueCoefficient` \(1\). The plus trace is \(N\cdot\eta w/2+P\) with normal jet \(+I q_o\). The minus trace is \(-N\cdot\eta w/2+P\) with normal jet \(-I q_o\). `trace_three[0,2]` is \(0\). With height one power of \(\eta\) and \(\Delta P\) of grade \(\eta\sigma\), the height times \(i f q_o\Delta P\) is grade \((2,1)\), and the retained pressure slot stays \(\Delta P\). Prior final-slot files do not discharge this check: composition was not performed there. The double sign flip explains identical off-diagonal trace structure. The face sign has to be taken from `normalJet` and the extension, not from a signless off-diagonal entry.

The consumer census query is `substring delta_p` on every `Symbol`/`Function`, with child hashes for every top-level child. Opened counts: each `U_BODY_BALANCE` component \(0,1,2\) has \(517\) children and all four pressure counts \(0\); `THETA_BALANCE` has \(380\) children and one of each of the four names; `E_W_BALANCE` has \(989\) children and two of each. The eight `E_W` texts are the \(\Lambda_{X0}\) children plus the simple \(W_0\) and \(W_0^2\) children. `physical-input.json` has \(\Lambda_{X0}=0\) and \(W_0=1\), which is why the bound `E_W` slots are \(\epsilon/2\) and \(\pm\epsilon\eta w/4\). The `THETA` census coefficient at \(\omega=3\), \(\rho_m=1/10\), \(\Lambda_{A0}=1/100\), \(\tau_A=1/10\) prints as \(-I\epsilon(3-10I)/109\) and \(\mp I\epsilon\eta w(6-20I)/436\), matching the grade splits. Historical \(\omega=1\) is not the composed frequency. Epsilon stays once in those row coefficients. The direct \((1,1)\) block pairs with pressure consumer \((0,0)\) and source \((0,0)\). Normal consumers begin at grade \((1,0)\), so normal times direct is \((2,1)\). The omission record’s unprojected \(\eta^2\sigma\) remainder stays excluded, and its removed slot is only `delta_p_plus`.

Flat pressure uses input depth \(q_i\). Both flat jets use output depth \(q_o\) in the prefactor and the pole. The written \(1000\) and \(100\) beta forms are the same \(\Omega/(10-I\Omega)\). \(\delta(l-k)\) may move the whole flat coefficient onto \(q(l)\). It does not identify off-diagonal depths. `kernel_apply` `p0` is already the reduced route. The edge contract reduces two edge deltas with one \(1/(2\pi)\) factor and a plain middle measure, which is the method’s one-dimensional \(\hat b\).

`taggedTotalMixed` and `jetKernels[f].mixed` already contain the direct tag. Using them only as reconstruction targets is the right rule. Row `mixedCoefficient` entries in the retained increments already multiply raw kernels by a reference-factor structure; they are not extra summands on \(D_{\mathrm{whole}}\).

The control point \(c_s=\sqrt{10/7}\), \(p=k=3/2\), \(l=2\), \(r=30/13\) gives \(q(k)=2\), \(q(l)=3/2\), \(q(r)=25/26\) from \(q(n)^2=9/c_s^2-n^2-1/20\), with the edge contribution \((1/5)^2+(1/10)^2=1/20\). Those are proposed supports. Nonzero channel applicability is still an instrument check. Global composition and the grazing limit stay unresolved, and no cutoff is inserted. That is in scope for an inventory.

## Not blockers

`shape_coefficients` drops products once a combined eta or sigma power exceeds \(1\). Section 3’s later quotient sentence controls: grades come from the saved full rational, with numerator, denominator, zero grade, and the unsplit higher remainder kept. The instrument must not call that truncating function as the extractor.

`fullLowerNormal` is the same vector in both face-domain files, ending in \(-1\). Face signs are carried by the height and jet fields, which do flip. The census does not store constructor text for non-hit children, so a `d_w_` name that does not contain `delta_p` is not visible from the occurrence arrays alone. Section 2 already requires that broader scan against source text before the four-name filter. Both face-domain files storing one shared lower normal is consistent with the jet evidence above.

## What a corrected bounded instrument would establish

It would establish an ordered both-face address inventory on \(G=\{(0,0),(1,0),(0,1),(1,1)\}\): native c2 flat, height, slope, and mixed iteration, plus one inherited direct whole convolution that is absent from native `z_three[0,2]`. It would record all \(16\) grade triples, explicit zeros, epsilon count one, the normalized one-dimensional Fourier route, four momenta where the address allows them, source derivatives at incident \(p\), and \(i f q(l)\) immediately before the response. Formal controls would show tagged-coefficient movement on delta support, including removal and doubling of the direct tag, without evaluating a convolution. Compact response momenta \([-3,3]\) would remain an explicit limitation.

## Still unresolved before a finite near-unity defect pilot

Evaluated source, consumer, and response convolutions; domain, test-space, and quadrature for momenta outside \([-3,3]\); a grazing limit of the composed operator; any reuse of cutoffs \(4\) or \(6\); cross-grade runtime joins beyond the old selected direct grade-zero witness; and the unfinished raw-versus-density numerator identity. None of those is required to clear an honest inventory, and none is cleared by it.

## Runtime evidence still required

The new joins themselves: the corrected three-way factor relation and its printed residuals; formal adapter bindings and original `build_face` assignment text for both faces; the \((2,1)\) height term retained as excluded data; quotient grade extraction of both source faces, with the lower face joined from its own operands; the broad `delta_p` or `d_w_` name scan before four-slot filtering, with any extra name stopping the run; material bindings \(W_0=1\), \(\Lambda_{X0}=0\), and real \(\omega=3\), with eta and sigma left independent; and control coefficients on the stated delta supports, reported as formal tag coefficients when the transform or whole kernel is unevaluated.