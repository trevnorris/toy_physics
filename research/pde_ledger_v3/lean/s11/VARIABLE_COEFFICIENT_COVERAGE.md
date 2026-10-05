# Variable-coefficient bulk identities and flat-interface traces, VC1–VC4

Authorized by the user's instruction to commit T1–T4 after successful review
and proceed to question 2. T1–T4 is complete at `ad365b5f`. Governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). This is a new bounded
increment. The recorded local proof, control and source checks pass, and both
independent non-author fidelity reviews are CLEAR. The bounded contract is
complete; see [review and closure](VARIABLE_COEFFICIENT_FIDELITY_REVIEW.md). See
[VARIABLE_COEFFICIENT_FIDELITY.md](VARIABLE_COEFFICIENT_FIDELITY.md) and
[VARIABLE_COEFFICIENT_VERIFICATION.txt](VARIABLE_COEFFICIENT_VERIFICATION.txt).

| Item | Claim and finite deliverable |
|---|---|
| VC1 | For the already classified D3 density L=-[a(div u)^2+b tr(G^2)+c norm(G)^2]/2, replace a,b,c by arbitrary smooth prescribed real profiles. Reuse the actual jet derivative for momenta, and prove the exact local EL: (a+b)grad div u+c Delta u+(grad a)div u+sum_i (partial_i b)partial_j u_i+sum_i (partial_i c)partial_i u_j. Recover constant coefficients and the residual for a=-b,c=0. |
| VC2 | Prove the generic weighted-divergence product rule on the existing spacetime, then instantiate the reviewed D3 null current and D4 odd current. For D4 L=-beta P/2 prove the actual momentum-divergence EL is (1/2)sum_i(partial_i beta)M_ij; a constant beta recovers zero. The current and momenta can remain nonzero. |
| VC3 | Prove an actual normal-slice integration-by-parts identity on two adjacent finite intervals, using separately regular fluxes and a common test trace. Retain the interface term (p_minus-p_plus)h at the interface and outer endpoint terms (or explicit zero endpoint hypotheses). Define the normal momentum/traction from the same D3/D4 momenta and prove that its pairing with every finite-dimensional test value vanishes iff the normal momenta match. |
| VC4 | Compact source/parameter identification, sign/index/factor and omitted-gradient/jump controls with admissible nonzero and constant-profile cases, sequential strict builds and standard-axiom audit; two independent non-author fidelity reviews and resolution of findings. |

Coordinate convention: G_ij=partial_i u_j, spatial derivative rows; the existing
Point D has coordinates (t,x1,...,xD). Coefficients are prescribed functions,
not varied fields. Smooth spacetime profiles are allowed in the local identities;
spatial profiles are included. These densities have zero time momentum. No
coefficient positivity, nonzero condition, transverse ansatz or stationarity
of the prescribed background is assumed. Smoothness is a sufficient stated
premise, not a minimal-regularity theorem.

VC1 and VC2 concern the exact already-reviewed local quadratic families, not
the full closed S11c operator, its nonlocal terms or its three field sectors.
The new compact native check may differentiate these same supplied families
with variable coefficients and compare the resulting identity; it must not
import a production driver or change any native source/export. Historical
D3/D4 basis/action identifications remain pinned evidence. A generic identity
is useful for S11c but does not identify all of its operator rows.

For the D3 null family define J_i=sum_j[u_i partial_j u_j-u_j partial_j u_i].
Then div J=(div u)^2-tr(G^2), and
 a div J = div(a J)-grad a dot J.
For the D4 reviewed current K, div K=P and
 beta P = div(beta K)-grad beta dot K.
Thus equalities modulo a divergence require retaining the coefficient-gradient
term. Constant-coefficient bulk equivalence alone does not grant equivalence
of traction or of an interface problem.

VC3 uses oriented finite intervals and separate flux representatives defined
on the real line, with actual derivatives on each closed interval and integrable
derivative data. These are sufficient extension hypotheses, not a weak
one-sided-trace theorem. The flux representatives need not agree at
the interface. The test has a common trace. The normal points from minus to
plus; momentum flux is n_i dL/dG_ij, while stiffness traction has the opposite
sign for the stated negative L. Zero interface contribution for every test
trace is flux continuity only when no independent surface action/source is
present. The theorem does not impose continuity of the field or assert that
all arbitrary trace values arise from solutions.

This is a normal-slice/trace identity with a finite-dimensional traction
criterion, not a general multidimensional transmission theorem. Curved
interfaces, Sobolev trace theorems, tangential Fubini lifting, discontinuous
field products, thin-interface limits, general null Lagrangian classification,
full integrated stationarity/solution existence, scattering and pole work are
excluded. Any multidimensional application must supply its trace/regularity
and integration hypotheses. A surface term cannot be erased by bulk equality.

Coverage is universal over the stated smooth profiles/fields and all finite
jet/normal/test values, with constant/zero profiles, equal and unequal traces
explicitly covered. Controls must expose omitted gradient terms, lost current
factors, incorrect derivative indices and reversed or omitted jump signs.
Instrument failures/timeouts are not mathematical rejections. Stop after
VC1–VC4, reviewed fidelity, and this explicit application boundary.

One Lean worker, -j1 -M4096, 600 seconds per process, sequential jobs with
process-group timeout cleanup and silent completion/error hooks. Preserve all
historical reviewed proof/evidence and S11c files/exports. Rebuilt unchanged
imports may have new object hashes; retain old records as historical. A new
fixed review packet was explicitly approved and reviewed. The user subsequently
authorized committing after both clearances and proceeding to question 3.
That separate authorization does not extend this contract or approve transfer
of a future packet.
