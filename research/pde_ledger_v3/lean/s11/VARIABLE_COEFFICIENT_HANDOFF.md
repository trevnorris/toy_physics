# Handoff to the S11c calculation session: variable coefficients and interfaces

This addresses **question 2: which constant-coefficient bulk equivalences
survive spatial profiles and interfaces?** The bounded VC1–VC4 Lean proofs,
mutation controls and compact native checks pass. Claude and Grok independently
cleared the approved fixed packet with no required mathematical corrections.
**The bounded VC1–VC4 contract is complete. Its S11c application limits below
remain obligations for the calculation side.**

The result is precise: the same local D3 and D4 actions acquire the
coefficient-gradient terms below. A density that was a divergence at constant
coefficient also acquires a gradient/current correction when weighted by a
profile. At an interface, the actual integration-by-parts identity retains the
jump in normal momentum flux. Constant-coefficient bulk equality does not by
itself justify deleting these profile or boundary terms.

This formalizes familiar product-rule and variational mathematics for the
specified families. It is not a new physical discovery, a general transmission
theorem, or a verification of the complete S11c reduced operator. No S11c
calculation source or pinned export was edited or regenerated.

## Conventions and scope

Use real fields u and prescribed real coefficients, held fixed under variation.
Write spatial coordinates x1,...,xD, with

\[
G_{ij}=\partial_i u_j.
\]

Thus the derivative is the **row** index. Formulas here use one-based spatial
indices; Lean field components are zero-based and the spacetime point contains
time at index 0 followed by the spatial coordinates. The local theorems assume
smooth fields and profiles, with no positivity, transverse ansatz or background
field-equation premise. Lean permits time-dependent profiles as well: the time
momentum is identically zero because these particular actions use only spatial
derivatives. The compact native tests use spatial profiles.

Momenta are actual jet derivatives of the supplied action. The convention is

\[
p_{ij}=\frac{\partial L}{\partial G_{ij}},\qquad
\mathrm{EL}_j=-\sum_i\partial_i p_{ij}.
\]

The proofs establish local EL and pointwise first-variation identities. They
do not add a general integrated variable-coefficient stationarity theorem.

## D3: the full product-rule correction

For the same classified D3 family,

\[
L=-\frac12\left[a(x)(\operatorname{div}u)^2
                 +b(x)\operatorname{tr}(G^2)
                 +c(x)\sum_{ij}G_{ij}^2\right],
\]

the actual momentum and local EL are

\[
p_{ij}=-a\delta_{ij}\operatorname{div}u
       -b\,\partial_j u_i-c\,\partial_i u_j,
\]

\[
\begin{aligned}
\mathrm{EL}_j={}&(a+b)\partial_j\operatorname{div}u+c\Delta u_j\\
&+(\partial_j a)\operatorname{div}u
 +\sum_i(\partial_i b)\partial_j u_i
 +\sum_i(\partial_i c)\partial_i u_j.
\end{aligned}
\]

The b-gradient contracts with the **transposed** field derivative
\(\partial_j u_i\), unlike the c-gradient term. This distinction is tested
explicitly: b=x1, a=c=0 and u=(x2,0,0) give EL=(0,1,0). Replacing that b-term
by \(\sum_i(\partial_i b)\partial_i u_j\) incorrectly gives zero.

The formerly bulk-null family a=-b,c=0 therefore retains

\[
\mathrm{EL}_j=(\partial_j a)\operatorname{div}u
              -\sum_i(\partial_i a)\partial_j u_i.
\]

It is generally nonzero: a=x1 and u=(0,x2,0) give EL=(1,0,0).
Constants recover the preceding D3 bulk contract. We have not classified every
variable profile/background for which this residual happens to vanish.

For the reviewed current

\[
J_i=\sum_j\left[u_i\partial_j u_j-u_j\partial_j u_i\right],
\qquad \operatorname{div}J=(\operatorname{div}u)^2-\operatorname{tr}(G^2),
\]

the weighted null-family action satisfies the exact local identity

\[
L_{a,-a,0}=-\frac12\operatorname{div}(aJ)
           +\frac12\nabla a\cdot J.
\]

The second term must remain if a varies. The first term still matters at a
boundary or interface.

## D4: the odd density with a profile

In the one-based notation used here, let

\[
F=G-G^T,\qquad P=F_{12}F_{34}-F_{13}F_{24}+F_{14}F_{23}.
\]

This is the same P written as F01 F23-F02 F13+F03 F12 in the Lean sources.
The native computed odd basis has P_D=P; the fully summed contraction
\(\sum_{ijkl}\epsilon_{ijkl}G_{ij}G_{kl}\) is 2P. Keep this normalization
when mapping the coefficient beta. The epsilon identity is prose/hand algebra
in this increment and recorded in the earlier native D4 check; Lean proves
the derivative-defined dual matrix and G:M=2P.

Define \(M_{ij}=\partial P/\partial G_{ij}\). For \(L=-\beta(x)P/2\),

\[
p_{ij}=-\frac{\beta}{2}M_{ij},\qquad
\mathrm{EL}_j=\frac12\sum_i(\partial_i\beta)M_{ij}.
\]

The proof uses the reviewed smooth-field identity
\(\sum_i\partial_i M_{ij}=0\). Thus constant beta gives zero bulk EL, but
variable beta generally does not. For beta=x1 and u=(0,0,0,x3), the actual
response is EL=(0,1/2,0,0).

The reviewed current \(K_i=\frac12\sum_j u_jM_{ij}\) satisfies div K=P,
and the corresponding weighted action identity is

\[
-\frac\beta2 P=-\frac12\operatorname{div}(\beta K)
                +\frac12\nabla\beta\cdot K.
\]

These conclusions concern the unique classified D4 odd family. They do not
provide a full D4 even-family bulk census or a general null-Lagrangian result.

## Interface term: orientation and traction

For scalar normal-coordinate flux representatives p_minus and p_plus and a
common test h, Lean proves the actual two-interval identity

\[
\begin{aligned}
\int_a^c p_-h' +\int_c^b p_+h'
={}&p_+(b)h(b)-p_-(a)h(a)\\
&+[p_-(c)-p_+(c)]h(c)
 -\int_a^c p_-'h-\int_c^b p_+'h.
\end{aligned}
\]

Each flux is a separate real-line function with actual derivatives on its
closed interval, and the derivative data are interval integrable. These are
sufficient extension hypotheses, not a weak one-sided-trace theorem. The two
flux values at c need not agree. The theorem uses oriented interval integrals;
the adjacent physical-interval interpretation uses a<=c<=b. If h(a)=h(b)=0,
the outer terms vanish; compact support is not needed for that corollary.

For a normal n pointing from the minus side toward the plus side, define

\[
\pi_j=\sum_i n_i p_{ij}.
\]

The interface pairing is \((\pi_- -\pi_+)\cdot h\). A separate Lean theorem
proves that it vanishes for **every finite-dimensional test value** iff
\(\pi_- =\pi_+\). For these negative actions, physical stiffness traction is
\(t=-\pi\), so the same continuity condition is \(t_-=t_+\), while the
interface pairing written in t has the sign \((t_+-t_-)\cdot h\).

The sign has a nonzero integral check: p_minus=2, p_plus=5 and h=1-x^2 on
[-1,1], with c=0, give -3. Constant-coefficient bulk-null terms can still
have nonzero momentum traction, which is also explicitly tested.

The scalar slice identity and finite-dimensional pairing criterion are the
proved results. To infer a multidimensional transmission condition from
stationarity, an application must separately justify the tangential integration,
trace regularity and admissible variations, as well as its bulk equations and
absence of an independent surface action/source. Field continuity, solution
existence and realization of arbitrary solution traces are not proved here.
For discontinuous profiles, use separate side data under these hypotheses;
the smooth-profile theorem does not license products of delta distributions
with ambiguous discontinuous traces.

## What to use in the S11c calculation work

1. Identify the actual local action family, coefficient normalization and
   coordinate scaling before importing these component formulas. This packet
   does not identify the nonlocal terms or every sector of the closed operator.
2. Insert profiles into the action or its momentum **before** taking coordinate
   divergence. Substituting profiles into an already simplified constant-
   coefficient EL expression misses the displayed product-rule terms.
3. When using a divergence equivalence, retain both the coefficient-gradient
   correction and the boundary current. Derive interface momentum from the
   original action, not solely from the constant-coefficient bulk operator.
4. Supply side regularity, a common test trace and any surface source/action
   explicitly when using the jump identity. A proof of full multidimensional
   stationarity or weak transmission remains a separate application obligation.

These are concrete identities and checks to apply, not a finding that the
current S11c code has omitted them. No production rerun is requested by this
handoff. The new native instrument derives variable-profile momentum divergence
itself; it does not attribute that extension to the old constant-coefficient
native EL helper.

Evidence: forty recorded builds (34 unchanged local dependencies, five new
modules and the audit root), 31 standard-axiom audits, twelve mathematical
rejections, sixteen positives, seven exact native identities and ten native
controls. Run2 reused all forty successful run1 builds under full source,
dependency and input/output object guards and reran all controls. Lean proves
the stated nonzero components of the D3/D4 witnesses; the native check also
evaluates the full vectors above. No instrument failure counted as a rejection.
Historical proof/evidence and all five T1 objects matched at authorization.

Source entry points, all under research/pde_ledger_v3/lean/s11/:

- S11VariableCoefficients/Common.lean
- S11VariableCoefficients/D3.lean
- S11VariableCoefficients/D4.lean
- S11VariableCoefficients/Interface.lean
- S11VariableCoefficients/Controls.lean
- VARIABLE_COEFFICIENT_COVERAGE.md
- VARIABLE_COEFFICIENT_FIDELITY.md
- VARIABLE_COEFFICIENT_VERIFICATION.txt

The compact instrument is
research/pde_ledger_v3/_measurements/S11_lean_variable_source_check.py.
The approved 61-file packet has aggregate SHA256
8152df86a31ecaabe89713a3fb178c7c469885f905b9e486207bbc695c02cf45.
This requested handoff is outside that frozen packet. Both independent reviews
are CLEAR; their optional observations are resolved or explicitly deferred in
VARIABLE_COEFFICIENT_FIDELITY_REVIEW.md. No proof or instrument changed.
The user authorized committing VC1–VC4 and proceeding to question 3 after
these clearances.
The condition and scope are recorded in
research/pde_ledger_v3/_measurements/S11_lean_variable_continuation_authorization.json.
