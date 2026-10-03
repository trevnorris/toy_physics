# Bounded Fourier request bank: concrete build

This implements the Fourier prerequisite of the existing wavelength-resolved
packet-action method, section 4. It evaluates no pressure response integral and
no complete packet action. The accepted preflight is restored without its
functions: 544 selected addresses, 34 complete coefficient-vector/proof sets,
actual native wave jets, fixed physical context, and K=27/U=75/T=122. Its 222
original files remain in the inputs. Additional complete published preflight
operands/returns bring this worker's saved-input count to 826 files.

The purpose is to establish a reusable numerical transform implementation on
actual operands before nesting it inside the much larger collision-aware action.
A pass at these 190 requests is not a uniform transform-error theorem and does
not waive comparisons at every future actual argument. Future requests can
reuse a complete result only with exact full argument and implementation joins.

## Actual requests and inherited source joins

`family_plan` reads the complete selected original THETA row, then visits every
live address. Exact coefficient ID and native derivative orders identify 13
source families. Three consumer fields give three test and three test-derivative
families. The bank therefore has 19 families, two packet carriers (sqrt(595)/10
and zero), and five physical transform arguments (-27, -sqrt(595)/10, 0,
+sqrt(595)/10, +27), totaling 190 requests. These are fixed before any numerical
comparison. For Y the transform argument is -l; its signs are consequently
opposite the corresponding physical output momentum labels.

Every selected address joins its complete original row entry and accepted
preflight adapter input. Every coefficient vector joins its original field,
certificate, saved numerator/denominator object, and published reconstruction
operands. Each new rational-pair transport is an exact identity; the old proof
is inherited rather than recalculated. The source derivative orders join the
actual published wave-multiplier input and literal-zero return. Fields are
finite polynomials in tanh(x/10), with coefficients in Q+iQ. Their actual degrees
in this bank range from zero through five. The all-zero families are not
numerically integrated; all 544 original addresses are retained as provenance.

The fixed physical input remains real omega=3, strict rest bulk,
LAB_HELD/RHO4_CONSTANT, tangents 1/5 and 1/10, one effective cs=sqrt(6)/2. No
source is rebound to that bulk speed. The saved packet centers -5/2,+5/2 and
width s=8 join before computation. This bank needs neither a depth evaluation
nor H/J/D, and establishes no physical profile-product, paired-height or full
summand-unit identity. Those remain obligations before the complete action.

## Fourier signs, products and derivatives

With u=exp(-(x+5/2)^2/128) exp(i p0 x) and
v=exp(-(x-5/2)^2/128) exp(-i p0 x), define
X(k)=hat[b D_j u](k), Y(l)=2*pi*hat[c v](-l), with
hat[f](k)=integral exp(-ikx) f(x) dx/(2*pi). The pairing is complex bilinear,
without conjugation. Native D_j=(-3i)^nt (i/5)^n2 (i/10)^n3 partial_x^n1.
The spatial derivative acts on u before multiplying b.

In the evaluator, the product carrier is p0 for X and -p0 for Y. The generic
transform argument is k for X and -l for Y, so its signed offset is k-p0 or
p0-l. Normalization is 1/(2*pi) for X and one for Y. The requested Y derivative
multiplies the integrand by +i*x, including the physical center. No derivative
of the source coefficient is silently added.

For w=z-x0, Q0=1 and Q(n+1)=Qn' + (i*carrier-w/64)Qn. The exact finite recurrence
is implemented on coefficient lists. At z=x0+y-i*c*sign(nu), the polynomial
Qn and the requested coordinate factors are expanded in y. The full integrand
uses exp(-w^2/128-i*nu*z), retaining the carrier/contour factors together. It
uses the original coefficient polynomial at tanh(z/10); the contour is not a
shift of the physical coefficient profile.

For constant coefficient fields, a separate reference uses
s sqrt(2*pi) exp(-s^2*(argument-carrier)^2/2-i*(argument-carrier)*x0),
multiplied by the coefficient, native constants, normalization and (i*argument)^n.
Requested argument derivatives are obtained by a separate polynomial recurrence
in the transform argument, with an additional minus sign per Y derivative.
Actual constant families include second spatial derivatives and Y-prime.

## Independent quadrature rules and absolute checks

Route A uses its own mpmath context at 30 digits and open Gauss-Legendre24/48.
Nodes and weights are constructed separately for the two orders and saved with
all monomial moment residuals through degree 2n-1. Positivity, interior nodes,
and maximum residual below 1e-25 are required before transforms.

Route B uses an independent context at 50 digits. Its embedded G7/K15 rule is
constructed without borrowing any Route A rule, sample or array. Let
P7=(429*x^7-693*x^5+315*x^3-35*x)/16. The even monic polynomial E8 is determined
by the exact rational equations integral P7 E8 x^j=0 for j=1,3,5,7 on [-1,1].
The moment formula is the exact integral of a monomial, not a response integral.
The four-by-four rational system and its solution are persisted. E8's eight
roots are bisected in the eight intervals separated by G7 nodes and +-1, to
width <=1e-48. Failed interlacing or arithmetic stagnation is fatal. Fifteen
weights are determined from moments through degree14. Every moment through23
is then checked below 1e-42, with positive weights and open distinct nodes.
These are numerical rule checks; the assessed Gauss/Kronrod construction is not
a machine proof of quadrature error.

Each route selects its own radius before quadrature. A uses c=5, B c=0. The
finite polynomial majorant is formed from absolute real/imaginary coefficient
sums of the actual shifted Q and coordinate factors. The coefficient norm is
the sum of absolute real and imaginary parts of its rational coefficients.
On the allowed strip |tanh(z/10)|<=1. The bound is
2*C0*exp(c^2/128-c*abs(nu))*sum a_n M_n(R), where M_n is the one-sided Gaussian
tail moment. M0=8 sqrt(pi/2) erfc(R/(8 sqrt(2))), M1=64 exp(-R^2/128), and
M_n=64 R^(n-1) exp(-R^2/128)+(n-1)64 M_(n-2). The smallest positive integer
m with R=8m and tail <=1e-14 is chosen. All trial radii, moments and factors
are persisted. A finite radius capacity of m=4096 is a spatial capacity, never
an elapsed-time cutoff; a miss stops before quadrature. Moment evaluation is
high-precision numerical evaluation of assessed analytic bounds, not directed
interval certification.

Every route partitions at centered y=0 and the physical profile center y=-x0.
A panel widths are <=min(10/8,8/8,pi/(4*(1+abs(nu)))); B's independently built
initial widths are <=min(10/5,8/5,pi/(3*(1+abs(nu)))). B compares K15 and G7
on each panel, splitting its own panels and dividing the allocated absolute
target between children until the sum of accepted embedded estimates is <=1e-13.
This is the prescribed adaptive algorithm, not a retry after a failed transform
comparison. The embedded estimate remains empirical. Precision stagnation stops.
No fallback algorithm or additional order campaign exists.

At every request save all three complex values, tail bounds, adaptive estimates,
and absolute pairwise differences; constant cases additionally save all three
analytic-reference differences. Every difference must be strictly below1e-12.
A failure persists its evidence and stops. Exponential damping does not establish
relative accuracy, and this finite bank provides no complete-action error bar.

## Persistence, containment and limits

The worker copies all original inputs byte-identically, verifies gate/source/
manifest/review/authority identities before importing science, and uses the
unchanged shared guard around the existing supervisor. Resource settings remain
4GiB native/cgroup in the16GiB pool, zero swap, one CPU/thread,32tasks,4GiBhost
reserve, desktop priority, RuntimeMaxUSec=infinity and Restart=no. The completion
hook targets the existing session before launch. No scientific job or READY
gate exists at build-review preparation.

`fourier-evidence.sqlite` uses FULL synchronous transactions, immutable record
keys and SHA256 over complete JSON bytes. Rule operands/returns, each transform
input, tail selection, full points/values/weights, panel sums, adaptive rejected
parents and accepted children, and comparisons are durable. Shared rule weights
are fully saved once and referenced exactly. Each completed request is also
receipted in the ordinary JSON journal. A failure closes and hashes the SQLite
file before rethrowing; all prior data remain. No payload is restored outside
containment. Later reuse must join exact values, signs, normalization, precision,
rule version and original source arguments, not only a hash or printed summary.

The numerical library's exact installed mpmath1.3.0 source files are pinned and
snapshotted. The actual loaded module path/version must join. Stdlib storage,
source and metadata tests exercise no scientific routine. No completed source,
response, endpoint, inventory or preflight calculation is replayed. This is only
a Fourier prerequisite. Complete action collision cells, local action tails,
physical h/j and paired-height joins, per-summand units, H/J/D route checks and
actual numerical controls remain future work. No current, inverse or loss claim.
