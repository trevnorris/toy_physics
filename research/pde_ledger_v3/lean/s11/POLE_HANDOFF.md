# Handoff to the S11c calculation session: nonlinear poles and projectors

This addresses **question 3: which hypotheses justify the residue and contour
construction, particularly at multiple or defective poles?** Local NP1–NP4
verification passes. Following explicit user authorization, Claude and Grok
independently returned CLEAR on the fixed packet with no required mathematical
corrections. **The bounded NP1–NP4 contract is complete. Its physical S11c
application hypotheses remain obligations; this is not a physical pole certificate.**

The result supports your already adopted `nonlinearPoleV2` correction:
`research/pde_ledger_v3/directives/S11c_d_NONLINEAR_POLE_CONTRACT.md`.
It supplies proofs of a bounded algebraic/contour core and explicit examples
that distinguish residues, multiplicity counts, projections and physical
response. It does not discover a new physical pole or supersede your repair.
No S11c source, pinned export or native diagnostic was edited or regenerated.

## 1. The valid local projection formula is conditional

Keep field, source and modal spaces distinct. With

\[
A:X\to Y,\quad V:K\to X,\quad W:Y\to K,\quad D=WAV,
\]

assume the **full** modal pairing D is invertible and set

\[
R=VD^{-1}W,\qquad P_X=RA,\qquad P_Y=AR.
\]

Lean proves RAR=R, idempotency of both P_X and P_Y, their exact ranges
range V and range AV, and their ranks dim K when K is finite-dimensional.
R itself is a source-to-field map, not generally a projector. A one-vector
pairing does not establish the result for a larger eigenspace.

The basic algebra supplies A as a linear map. Calling it `derivative` in the
code does not identify it with the physical L'(z*). W is not a left-kernel map
until the required annihilation premise is supplied. Zero-dimensional K is
allowed, so this algebra by itself does not even assert that a pole exists.

There is a separate conditional identification theorem. For L0=L(z*), if

\[
\operatorname{range}V=\ker L_0,\quad WL_0=0,\quad
L_0R=0,\quad L_0H+AR=I_Y,
\]

then the supplied coefficient R equals VD^{-1}W. These are leading and constant
identities from a simple inverse expansion when such an expansion exists and
A is the actual derivative. Lean proves their algebraic consequence; it does
not construct that expansion or prove the general analytic semisimplicity
criterion. Full-kernel coverage and invertibility of D remain real premises.

## 2. Retain ordered higher coefficients and response-map derivatives

The proofs use actual positively oriented, zero-centered circle integrals,
normalized by 1/(2 pi i), with positive radius. For any finite Laurent sum in
the stated complete complex normed space, the integral extracts precisely its
z^-1 coefficient. Matrix examples use the finite elementwise norm.

For the supplied finite double principal part

\[
Q(z)=C_{-1}/z+C_{-2}/z^2,
\]

Lean proves the actual circle-integral identity

\[
\frac1{2\pi i}\oint Q(z)(A_0+zA_1)\,dz
=C_{-1}A_0+C_{-2}A_1.
\]

The matrix order matters. If the affine operand represents a pencil derivative,
its coefficients must be the actual derivative coefficients; that identification
is supplied by the application. The finite formula does not silently remove an
unproved holomorphic remainder from an actual inverse.

For rectangular observation O(z)=O0+zO1 and forcing B(z)=B0+zB1, the actual
finite response integral is

\[
\frac1{2\pi i}\oint O(z)Q(z)B(z)\,dz
=O_0C_{-1}B_0+O_0C_{-2}B_1+O_1C_{-2}B_0.
\]

Thus freezing the forcing or observation at a double pole can lose part of the
residue. To apply these identities to a general meromorphic physical inverse,
establish its expansion, regularity of the maps and the necessary remainder
integral identities. Those general analytic steps are not proved here.

## 3. Actual contour examples prevent the misleading shortcuts

| Pencil or response | Lean-checked conclusion | Consequence |
|---|---|---|
| Scalar L=z² | Inverse residue 0; logarithmic integral 2, whose square is 4; weighted inverse moment 1 | Zero residue does not remove the double pole, and the logarithmic count need not be a projector. |
| Scalar L=z²−1 | Zeros exactly ±1; the radius-2 logarithmic integral is 2 and is not idempotent | Enclosing multiple simple roots does not justify a global nonlinear projector claim. |
| L=zI−N, N=[[0,1],[0,0]] | Actual two-sided inverse I/z+N/z², derivative I, determinant z², one-axis kernel and explicit chain; inverse/logarithmic integral I, trace 2; nonzero higher coefficient N | A defective affine example still has a genuine state projection; its higher coefficient remains essential. |
| Faithful scalar transfer e1ᵀ(zI−N)⁻¹e2 | Equals z^-2; residue 0 but transfer value 1 at z=1 and a nonzero higher moment | Compression can hide the residue without eliminating the singular response. |

These are explicit pencils, inverse identities and actual integrals, not merely
algebraic values labeled as contours. The examples discriminate these cases;
they do not classify every nonlinear pencil or prove a general Riesz calculus.
Inverse identities have the required nonzero-frequency premises, and the
integration paths avoid the poles.

For a concrete response-map control, use the scalar transfer z^-2 with
O=2+5z and B=1+3z. Then

\[
O(z)z^{-2}B(z)=2z^{-2}+11z^{-1}+15.
\]

The actual residue is **11**. Frozen maps give residue **0**; omitting the
observation derivative gives **6**, and omitting the forcing derivative gives
**5**. All three wrong values are rejected by the paired mathematical controls.

## 4. What remains necessary for the actual S11c application

Before promoting the finite identities to a claim about the physical pencil:

1. Identify fixed field/source spaces, the operator domain, topology and sheet.
   Establish the analytic inverse framework and an admissible contour. The
   finite-matrix regularity or operator Fredholm hypotheses in
   `POLE_ASSESSMENT.md` are application/literature inputs, not NP Lean theorems.
2. At a proposed semisimple pole, use full right/left nullspaces and the actual
   full derivative pairing. Supply the inverse-expansion/coefficient hypotheses
   before calling VD^-1W the physical residue. The algebra alone is not an
   existence, isolation, normalizability or radiation-condition test.
3. At higher-order poles, retain the required principal-part coefficients and
   frequency dependence of forcing/observation. Check the actual response
   O L^-1 B; a vanishing residue alone cannot discard that channel.
4. Label each contour object by the theorem supporting it. Idempotency of a
   state projection, an algebraic multiplicity count and a physical response
   residue are different claims. General multiplicity/trace assertions require
   their own applicable analytic theorem; a general Riesz projection requires
   an identified fixed-operator realization and admissible resolvent contour.

These are the premises to supply when using the result. This handoff does not
assert that the current calculation code omitted them, request a production
rerun, or authorize changes to the pinned v10 inputs. Homotopy/displacement
theory, a general chain classification, physical pole searches and scattering
remain outside this bounded increment.

## Evidence and source entry points

The final recorded run passes seven canonical builds including the audit root,
55 standard-axiom audits, 17 mathematical rejections and 21 positive controls.
Run3 reused recorded run1 objects under full source/pin/command/input/output
guards and ran every control fresh. Each false statement produced exactly one
unsolved False diagnostic; no instrument failure counted as a rejection.

The compact native check passes 19 checks and four translation controls. It
parses selected original source assignments and reads existing synthetic
operator/map records; the older 89-check correction evidence remains historical.
This is a tested compact identification, not a kernel-certified CAS bridge or
a verification of the physical closed operator. All seven read-only native
inputs and forty VC plus five T1 compiled objects match their preservation record.

Files under `research/pde_ledger_v3/lean/s11/`:

- `S11NonlinearPole/Modal.lean`: typed pairing, projections, conditional residue.
- `S11NonlinearPole/Moments.lean`: actual finite Laurent circle integrals.
- `S11NonlinearPole/Laurent.lean`: ordered products and rectangular response maps.
- `S11NonlinearPole/Scalar.lean`: scalar contour examples and residue 11.
- `S11NonlinearPole/Jordan.lean`: actual affine inverse, chain and faithful transfer.
- `S11NonlinearPole/Controls.lean`: admissible pairing and order/frozen-map witnesses.
- `POLE_COVERAGE.md`, `POLE_FIDELITY.md`, `POLE_VERIFICATION.txt`: scope and evidence.
- `POLE_ASSESSMENT.md`: analytic application hypotheses and primary references.

The approved 31-file review packet has aggregate SHA256
`9c7e8e41dd1820647c3dc787ad171043944038c70d37b873699d414f49396e57`.
This requested handoff is outside that frozen packet. Both independent reviews
are CLEAR; their source-reading limits and optional findings are recorded in
`POLE_FIDELITY_REVIEW.md`. The zero-map control tests noninjectivity rather than
failure to construct `PairingData`, and the nonzero point-value control is
separate from the proved higher-coefficient identities. Those distinctions do
not change the general theorems. No proof or instrument changed after review.

The user subsequently authorized a `lean:` checkpoint commit after both reviews
clear and required findings are resolved; those conditions are now met. That authorization is recorded in
`research/pde_ledger_v3/_measurements/S11_lean_pole_commit_authorization.json`.
It supersedes the earlier review-launch handoff's statement that no NP commit
was requested; it does not authorize further mathematical scope.
