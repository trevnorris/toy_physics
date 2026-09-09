# Codex-sol (gpt-5.6-sol xhigh) review — S11c-c2 N6 reconcile DISPOSITION (round 1)

⚠ The raw codex transcript (~1.4 MB) echoed the full contents of the source files it read; that noise is trimmed
for tree hygiene (regenerable by re-running the leg on `..._disposition_review_prompt.md`). The review content
below is verbatim.

---

NOT-SOUND: three claim-changing corrections are required. The numerical tallies themselves are faithful, and the load-bearing EL calculation is correct.

### Must-fix findings

1. §1 incorrectly sweeps genuine `FROZEN_RELATIONS` operand agreement into `(0)−(0)`

Artifact: "`N6RC_FROZEN_RELATIONS` 8/0" … "These matched zeros are `(0)−(0)`." — disposition:33, disposition:37
Source: `N6RC_FROZEN_RELATIONS 8 0 692` — tally:24

Those eight matches are not all zero operands. SymPy defines `adv = b.dot(b.u, g)` and
`h = b.dot(b.u, b.grad_W) / b.W_bg if alpha == 'LAB_HELD' else 0` — reconcile_sympy.py:48; WL defines the same
nontrivial operands — mathematica audit:839. Consequently, four of the eight matched leaves are genuine equality
of nonzero operands: `h_alpha` in the two `LAB_HELD` cases and `a_rho` in the two `RHOBR_CONSTANT` cases. The
other four are `(0)−(0)`.

Likewise, the 400 structural-support zeros are agreement between two `NONZERO_WITNESSED` observations: if either
side reports `NO_NONZERO_FOUND`, the comparator returns `UndecidedResidual`, not zero. — comparator:537

What must change: restrict the `(0)−(0)` statement to the covariance/control residual objects, carrier bridge, and
zero-valued advection-absence leaves. State separately that `FROZEN_RELATIONS` contains four genuine nonzero
premise-operand agreements and four zero/zero matches; preserve support as support agreement.

2. §2 omits the three schema-unmatched reconcile-source families

Artifact's exhaustive-looking list says "`R_N6` itself, the channels, and the guards are ENTIRELY UNMATCHED" and
lists `*_OPERAND`, `*_CHANNEL`, `CROSS_CHANNEL`, `DIMENSIONS`, and guards. — disposition:55
But the tally also has `N6RC_SOURCE_BRIDGE_RESIDUAL 0 0 640`, `N6RC_SOURCE_EULERIAN 0 0 640`,
`N6RC_SOURCE_MATERIAL 0 0 640` — tally:27. The framing explicitly requires these reconcile-engine sources to
remain schema-unmatched and UNDECIDED. — question:48

What must change: add `N6RC_SOURCE_EULERIAN`, `N6RC_SOURCE_MATERIAL`, and `N6RC_SOURCE_BRIDGE_RESIDUAL` to §2's
unmatched list and carry them explicitly as schema-UNDECIDED.

3. §3 overstates upstream replay as the unique remedy and permanent claim boundary

Artifact: "The correct instrument applies … UPSTREAM of EL" and such a replay would "NOT retroactively reconcile
the already-emitted `.out` streams." — disposition:92
The astra source is narrower: "Equivalently, a bridge could carry coefficient jets and their chain rules…" —
astra consult:44 [line 48]; and a replay does not retroactively reconcile the original streams "unless their
relation to that replay is also established." — astra consult:48 [line 50]

What must change: say that the proposed post-EL coefficient table is insufficient and that the clean proposed
replacement is upstream of EL. Do not call it the only correct instrument: a derivative-aware
coefficient-field/chain-rule bridge is an acknowledged alternative. Say that replay alone does not reconcile the
old outputs unless their relation to the replay is separately established.

### Independent EL calculation
Both engines form source amplitudes after Euler–Lagrange differentiation: WL `el` (.wl:144), SymPy `variation`
(diagnostic:345). In 1-D, `L = a θ'e' + b eW'θ'`, `T(a)=kR, T(b)=−kR/W, R=W_0/W`.
`T(EL_θ L) = −kR e'' + (kR/W)(e'W' + eW'')`. But `T(L) = kR θ'e' − (kR/W) eW'θ' = k θ'(Re)'`, so
`EL_θ(TL) = −k(Re)''`. Using `R' = −(R/W)W'`, `R'' = (2R/W²)(W')² − (R/W)W''`:
`EL_θ(TL) − T(EL_θ L) = (kW_0/W²)W'e' − (2kW_0/W³)e(W')²`.
The engines encode `W = W_0(1+η w1)`, `W' = σ_W q` before grade extraction (WL:127, SymPy:245). So
`(kW_0/W²)W'e' = (k σ_W q e'/W_0)(1 − 2η w1) + O(η²σ_W)`: retained `σ_W^1` and `ησ_W` coefficients; the second
term is `O(σ_W²)`. Reason (a), including survival at `σ_W^1`, is correct.

All §1/§2/§8 numerical counts match the tally literally. The artifact also correctly avoids claiming agreement for
the 40/76/18 operands, preserves the density-table/source distinction, honors the "not known to be just
thickness" guard, carries the debt and all three premise caveats plus earlier carries, and leaves Reading B and
c1 standing.

**NOT-SOUND**
