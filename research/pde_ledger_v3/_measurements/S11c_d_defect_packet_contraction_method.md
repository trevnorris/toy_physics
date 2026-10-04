# Exact finite-window contraction of the two packet pressure actions

2026-10-04. New numerical-method proposal for independent assessment. No new
algebra certificate, integral, evaluator or READY gate has been executed or
created. The formulas below are author derivations, not restored results.

The complete pressure-summand unit transport is accepted at 7d5b0d7e. It joins
544 addresses and 20 complete face/slot templates. That closes a dimensional
prerequisite; it does not make the full nested quadrature practical or ready.
The saved geometry counts 8,304,551,424 unshared A24+A48 addressed product
occurrences across its eight plans. Those are occurrences, not unique requests,
an elapsed-time forecast, a storage requirement or proof of impracticality.

The proposed change removes the repeated construction of Jwhole and Dwhole at
every outer (k,l) point by contracting packet factors first. It preserves the
continuous finite-window integral exactly, if the following source, algebra,
domain and absolute-integrability joins hold. It does not assert equality of
two different finite quadrature sums. No interpolation, approximate low-rank
compression, small-term deletion or changed physical operator is proposed.

## 1. Unchanged target and authority boundary

Retain the two complex-bilinear e_W-input/THETA-test Gaussian packets, carriers
kappa and zero, centers -5/2 and +5/2, width 8, held real omega=3, rest bulk,
LAB_HELD/RHO4_CONSTANT, saved edge momenta 1/5 and 1/10, and effective
cs=sqrt(6)/2. Here kappa=sqrt(595)/10 must still join its saved physical origin.
W=1 and L=10 retain their physical dimensions. Keep all 544 selected addresses,
including zeros, their 16 ordered source/response/consumer grade triples, both
faces and actual pressure/normal slots. The 16 completed local cells and their
two packet values are inherited without quadrature. No source, grade, profile,
response, geometry, quadrature-rule or completed numerical bank is replayed.

The original finite domains are k,l in [-K,K], t in [-T,T], with (K,T)=(27,122)
and the predeclared enlargement (29,124). The separate height-PV domain keeps
k in [-K,K], Q in [0,75], including both k+Q and k-Q branches; it is not cut
down to the k,l square. Retain all original all-real tail allocations. Reordering
an identical finite domain introduces no new tail or permission to prune it.

This proposal concerns J and the three ADDED direct terms only. Flat, first
height/contact/paired-PV, first slope and the H part of native mixed iteration
remain separate obligations under the existing method. In particular, do not
silently replace finite-T H by its analytic all-real expression. H's surviving
contact and independent physical-product check remain as previously assessed.
Native mixed iteration and the complete direct addend each enter once.

The next permitted scientific instrument would derive the missing exact
contraction/source/domain certificates and concrete dependency/readiness plan,
after its own build assessment. It would evaluate no packet integral. A full
numerical evaluator still needs a reviewed, concrete discretization, error
allocation, storage plan and responsive numerical controls. This method review
must not be treated as approval of an unspecified faster quadrature.

## 2. Actual source-bound scalar factors

For EACH original address use its own X(k)=hat[b D_j u](k) and
Y(l)=2*pi*hat[c v](-l), with derivative-before-coefficient ordering, the saved
epsilon coefficient, Fourier normalization, independent grades and physical
unit provenance. No extra complex conjugation is taken. The original jet is
(-3i)^nt (i/5)^n2 (i/10)^n3 partial_x^n1. A field name alone is not a join.

Write q(p)=sqrt(kappa^2-p^2) on the common positive-real/positive-imaginary
outgoing branch and E(p)=q(p)+beta. Keep mu=3/10,
a=1/(1-3i/10), beta=a*mu in their already established physical units; a is not
dimensionless merely because its magnitude has that displayed number.

Let N_f(l)=1 for a pressure slot, +i*q(l) for the upper normal slot and
-i*q(l) for the lower normal slot. These are the actual 20 saved template
assignments, not signs inferred from a face label. Define

    x(k) = X(k)/E(k),       y(l) = N_f(l)*Y(l)/E(l),
    A(z) = L*z/[4*sinh(pi*L*z/2)],       A(0)=1/(2*pi).

These abbreviations factor the resolvents already INSIDE the saved complete
J/D densities. They do not multiply a closed Dwhole by two more resolvents.
Their units must be obtained by transporting the accepted original operand
units through this equality. Normal depth stays at l inside y(l).

For reference, the ACTUAL saved generic middle formulas, before source/test
multiplication, are

    J(t;l,k) = Cj A(t) A(l-k-t)
      * (k+t)(2k+t) q(k)
      / [q(k+t) E(l) E(k+t) E(k) (q(k+t)+q(k))],

    D(t;l,k) = Cd A(t) A(l-k-t) / [E(k) E(l)]
      * { k(2l-t)/(q(l-t)+q(l))
          + k(t+2k)q(k)/[q(k+t)(q(k+t)+q(k))]
          + q(k)^2/q(k+t) },

    Cj = a*mu^2*W*L/4,       Cd = (-i*mu)*W*L/(4i).

The three terms in braces are reflected, height-route and quadratic, in that
order. Never introduce a product of the two root routes. Cd may be simplified
only after its original two phase factors have been joined. The future exact
certificate must extract the original `kernel_components` assignments, join
the saved runtime parameters and all 20 complete templates, and persist both
old and new full expressions before checking their residuals. No old integral
or native producer function is called for this purpose.

## 3. Finite-domain identities to be derived

Put M=K+T and I_T(m)=[max(-K,m-T), min(K,m+T)], interpreted as empty if its
lower endpoint exceeds its upper endpoint. The outer range is [-M,M]. All
integrals below are continuous integrals; their numerical evaluation is future
work. Define the following address-specific contractions:

    Y0(m) = integral_[-K,K] A(l-m) y(l) dl,
    X1(m) = integral_[-K,K] A(m-k) k x(k) dk,

    XrT(m) = integral_[I_T(m)] A(m-k) k^r x(k) dk, r=1,2,
    X02T(m) = integral_[I_T(m)] A(m-k) q(k)^2 x(k) dk,
    CrT(m) = integral_[I_T(m)]
        A(m-k) k^r q(k) x(k)/(q(m)+q(k)) dk, r=0,1,2,

    YrCT(m) = integral_[I_T(m)]
        A(l-m) l^r y(l)/(q(m)+q(l)) dl, r=0,1.

Then the proposed exact finite-window actions are

    BJ = Cj integral_[-M,M]
        Y0(m) * m/[q(m) E(m)] * [C1T(m)+m*C0T(m)] dm,

    BD_height = Cd integral_[-M,M]
        Y0(m)/q(m) * [C2T(m)+m*C1T(m)] dm,

    BD_quadratic = Cd integral_[-M,M]
        Y0(m)/q(m) * X02T(m) dm,

    BD_reflected = Cd integral_[-M,M]
        X1(m) * [Y1CT(m)+m*Y0CT(m)] dm.

For BJ, BD_height and BD_quadratic the substitution is m=k+t. Thus dt=dm,
A(t)A(l-k-t)=A(m-k)A(l-m), and |t|<=T restricts K's INPUT integration to
I_T(m). The output l integration remains the entire [-K,K]. For BD_reflected
the substitution is m=l-t. Its negative derivative is cancelled by reversing
the mapped limits, giving a positive measure; |t|<=T instead restricts the
OUTPUT l integration to I_T(m), and k remains the entire [-K,K]. These two
windows may not be identified just because both have the same endpoints as
functions of m. The variable to which each window applies is part of the proof.

In particular, m is not a new external momentum limited to [-K,K]. Retain the
whole [-K-T,K+T] interval and both clipped wings. No slab-dependent window may
be replaced by [-K,K] or by a constant inner interval. The equalities must be
certified separately for each of the four primitives, before adding them.

These identities use the original independent source and consumer functions;
no equality between two different field transforms is assumed. Identical
contractions can only be shared after all field/jet/carrier/face/slot/unit,
coordinate, window, rule, precision and route settings match exactly. Preserve
every original address and grade even when the value is shared. A grouped sum
against an addressed sum is assembly consistency, not an independent numerical
reference or evidence that a component is nonzero.

## 4. Absolute integrability, grazing and the domain map

The finite-domain Fubini argument is for ordinary J/D terms. It is not an
unsubtracted-PV Fubini argument and grants no new operation on the height/H
distributions. The fields times Gaussians have smooth bounded Fourier transforms
on each compact interval, by the accepted polynomial-field/Schwartz evidence.
q stays in the closed first quadrant and beta has strictly positive real and
imaginary parts at the saved binding. Away from simultaneous zero depths,

    |q1+q2| >= (|q1|+|q2|)/sqrt(2),
    |q1/(q1+q2)| <= sqrt(2),
    |q1+beta| >= |beta|/sqrt(2).

Consequently the J and height-route density are bounded, up to compact
polynomial factors and nonzero beta denominators, by C/|q(m)|; the quadratic
term has the same integrable majorant. The reflected term can also use
sqrt(2)/|q(m)| as an upper bound for 1/|q(m)+q(l)| when q(m) is nonzero.
The simple zeros at m=+-kappa give locally integrable inverse-square-root
majorants. Constants depend on the declared finite K,T and address envelopes;
this is not a uniform-in-cutoff error estimate. The existing all-real tail
proofs still carry that separate obligation.

The measure-zero simultaneous-depth-zero sets do not require an arbitrary
pointwise value. The certificate must keep them excluded from pointwise
equalities; a numerical worker uses open nodes and one-sided integrable bounds.
It must not assign zero to 0/0. New compact bounds need actual source/physical
joins and independent assessment; a printed inequality is not an executed
certificate or a machine proof of Fubini.

The linear substitutions have absolute Jacobian one, map the finite domains
bijectively apart from boundaries, and preserve the original dt dl dk measure.
For m=k+t the original q(k+t) becomes q(m), whereas q(l-t) would become
q(k+l-m), if that route were present. For m=l-t the original q(l-t) becomes
q(m), whereas q(k+t) would become q(k+l-m). Do not map both roots to q(m).
The splitting into three ADDED direct terms is what permits the different
substitutions. This is not a change of variables for an invented root product.

The future domain certificate must persist forward/inverse maps, inequalities,
orientation, all original labels and clipped windows. Its proposed cuts include
m=+-kappa, m=+-K+-T and the outer endpoints, input/output branch points at
their own +-kappa, and profile/packet resolution points in the actual variables.
The original collision surfaces must be transported per primitive and shown
either retained as a singular boundary or absent from that primitive by its
actual dependency. A whole-union geometric label cannot simply be forgotten.
The exact discretization and all line/window ordering intersections still need
a concrete build. This method does not approve dropping narrow intervals.

## 5. Native controls and complete numerical-route requirements

New exact certificate controls must affect the actual finite-domain/algebra
maps: apply the reflected window to k instead of l; omit a clipped wing; use
the wrong reflected substitution sign in its numerator; and omit a term from
the J polynomial m(k+m). Persist original/mutant operands and full residual or
domain-witness returns. These are algebra/domain sensitivity controls, not
packet-action values. Choose their symbolic/domain witnesses before results.

The original full numerical controls remain required. The reflected-root mutant
is NOT obtained by reusing the correct reflected contraction formula. With
q(l-t) wrongly replaced by q(k+t), use m=k+t in the mutant original term.
Writing YrC(m) as the UNCLIPPED [-K,K] version of YrCT(m), its proposed action is

    BD_reflected_wrong = Cd integral_[-M,M]
       {2*X1T(m)*Y1C(m) + [X2T(m)-m*X1T(m)]*Y0C(m)} dm.

This formula needs its own original-mutant integrand/domain join; it is not an
accepted computed control. Wrong normal q(k) versus q(l) must move the actual
normal multiplier from the Y side to the X side on an applicable saved normal
address. Flat support would make that mutation silent and cannot stand in for
the required applicable off-diagonal test. H-contact omission and Leibniz
corruption retain their original actual-address definitions. Compare complete
mutant-minus-baseline weak contributions for each carrier, with the original
10-times-summed-empirical-envelope rule. Domain failure is not numerical
responsiveness, and silence remains unestablished coverage.

A future numerical build must fix TWO complete independent implementations of
the contracted finite-domain integrals. They may share immutable sources,
accepted rule objects and analytic identities; they must not share evaluated
Fourier samples, contracted arrays, adaptive choices or action values. Keep
30-digit open squared GL24/48 and an independent 50-digit physical-variable
adaptive route as the proposed numerical families, but a concrete scheme for
the new nesting and finite error allocation is still required before launch.
An inner integral's empirical error cannot be ignored or reused as an outer
error bound without propagating its actual multiplier and outer weights.

Retain the original per-request Fourier comparisons and analytic tails, unless
a separate substantive method revision is assessed. No interpolant, transform
table interpolation or unchecked analytic H substitution is authorized here.
The old 190 Fourier requests and 38 inner points are reused only for exactly
matching settings; they are not surrogates for new arguments or uniform error
bounds. The new contraction order does not authorize rerunning those banks.

Require component- and grade-resolved A24/A48/B comparisons and the declared
K,T enlargement, with tau=1e-9+1e-7*abs(refined value); preserve positive tail
budgets and do not hide unresolved pieces by cancellation. The future review
must assess the new dependencies and costs, including unavoidable new Fourier
requests and adaptive work. This proposal makes no runtime or memory estimate.
Use the exact-storage preparation only after concrete integration assessment;
do not delete old evidence or launch an unbounded disk representation on hope.

## 6. Next bounded instrument and exclusions

If the method is substantively supported, prepare a worker for the exact
source/algebra/domain/unit joins and dependency/readiness certificate only.
Restore actual kernel ASTs, complete adapter maps, native source/test and unit
records, bounds, physical profile and finite-window operands. Derive only the
missing substitution/factorization/integrability/domain/control statements.
Do not invoke old producers, rule constructors, completed functions or numerical
integrals. No new scientific calculation has run during this preparation.

Scientific imports/restoration and the new certificate remain inside the
unchanged no-deadline shared pooled guard around the existing normalization
supervisor: 4GiB native/cgroup, 16GiB aggregate pool, zero swap, one CPU/thread,
32 tasks, 4GiB host reserve and desktop priority. Pin a fresh gate and arm the
local completion hook before launching. Preserve every failure and partial
return; no automatic retry, fallback or scheduler change.

This continues the same two bounded weak packet tests. There is no physical
current, leakage factor, scattering inverse, new mode census, finite matrix,
speed/defect sweep, production overwrite, drain or calibration claim. Preserve
all original method/review/result/failure history, Lean/S11_lean, the shared
guard and protected builder. Standing authority covers this necessary
established-Claude/Grok assessment, not author clearance or execution readiness.
