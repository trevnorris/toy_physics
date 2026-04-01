# Moving-Throat PDE — Stage 39: Transport Origin of the Lowest-Lane Source Asymmetry

## Purpose

Stage 38 left one sharp operator-level question open:

> what is the **physical** origin of the source-shape asymmetry parameter that had previously been written as the abstract exponential-bias variable `alpha`?

This stage answers that with the simplest conservative axial source-transport law on the finite throat.

The main result is that the Stage-36 exponential family is not ad hoc. It is the exact stationary zero-flux branch of a drift-diffusion transport operator,

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

on the finite throat interval `s in [0,L]`.

On the stationary recirculating branch `J=0`, the normalized source density is

`sigma_Pe(s)`
`= Pe exp(Pe s/L) / (exp(Pe)-1),`

with the single physical asymmetry parameter

`Pe := v_sigma L / D_sigma.`

So the old abstract bias parameter is exactly the axial **Peclet number** of the lowest support-source transport problem.

This converts the overlap side of Stages 36–38 from a free deformation into a concrete operator output.

---

## 1. Explicit axial source-transport operator

Let `sigma(s,t)` be the coherent source density feeding the lowest support lane along the finite throat interval `s in [0,L]`.

Take the minimal conservative axial transport law

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

with

- `D_sigma > 0` the axial spreading coefficient,
- `v_sigma` the directed recirculation / drift speed along the throat.

This is the simplest operator that can bias the support source toward one end of the throat without introducing explicit nonconservative loss.

---

## 2. Stationary recirculating branch and exact exponential profile

On a stationary branch,

`partial_t sigma = 0,`

so `J` is constant in `s`.

For the lowest closed recirculating support branch, the natural condition is zero net axial source flux,

`J = 0.`

Then the transport equation collapses exactly to

`D_sigma sigma' = v_sigma sigma,`

hence

`sigma(s) = C exp(v_sigma s / D_sigma).`

Normalizing to the same total source strength used in Stage 36,

`int_0^L ds sigma(s) = L,`

gives

`sigma_Pe(s)`
`= Pe exp(Pe s/L) / (exp(Pe)-1),`

with

`Pe := v_sigma L / D_sigma.`

So the Stage-36 constructive family is exactly the stationary zero-flux branch of the minimal transport operator.

### Physical interpretation of the sign

- `Pe = 0`: symmetric uniform source, i.e. the twin baseline;
- `Pe > 0`: constructive branch, with source weight shifted toward the D/N bottom antinode;
- `Pe < 0`: destructive branch, with source weight shifted toward the mouth where the D/N mode is smallest.

So the same operator explains both the helpful and harmful non-twin deformations.

---

## 3. Exact overlap boost on the transport branch

Keep the D/N lowest mode

`chi_0(s) = sqrt(2/L) sin(pi s/(2L)).`

Its mixed-lane uniform overlap is

`I_W = int_0^L ds chi_0(s) = 2 sqrt(2L)/pi.`

For the transport branch,

`I_Pe = int_0^L ds sigma_Pe(s) chi_0(s)`

evaluates exactly to

`I_Pe = 2 sqrt(2L) Pe (2 Pe exp(Pe) + pi)`
`       / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ].`

Therefore the physical overlap boost is

`Omega_Pe := I_Pe / I_W`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ].`

So the abstract Stage-36 factor `Omega_exp(alpha)` is now identified exactly as `Omega_Pe` on a concrete operator branch.

---

## 4. Exact monotonicity identity on the constructive branch

Introduce the probability density

`p_Pe(s) := sigma_Pe(s)/L.`

Then

`I_Pe = L E_Pe[chi_0],`

and differentiation gives the exact covariance identity

`dI_Pe/dPe = Cov_Pe(chi_0,s),`

hence

`dOmega_Pe/dPe = Cov_Pe(chi_0,s) / I_W.`

Because `chi_0(s)` is strictly increasing on `[0,L]`, the covariance is positive on the constructive branch.

Therefore

`dOmega_Pe/dPe > 0`

for `Pe >= 0`.

So the physical transport branch is not merely continuous. It is an exact **monotone** route from the twin point toward maximal overlap.

This is important later because it means the required transport bias is unique once the target overlap is known.

---

## 5. Endpoint and asymptotic structure

The exact endpoint values are

`Omega_Pe(0) = 1,`

`lim_(Pe -> +infinity) Omega_Pe = pi/2.`

So the constructive transport branch reproduces the full Stage-36 overlap window

`1 <= Omega_Pe <= pi/2.`

Useful asymptotics are:

### Small transport bias

`Omega_Pe`
`= 1 + ((4-pi)/(2 pi)) Pe + O(Pe^2).`

So the constructive branch leaves the symmetric point with a strictly positive linear slope.

### Strong transport bias

`Omega_Pe`
`= pi/2 - pi^3/(8 Pe^2) + O(Pe^-3).`

So the approach to the sharp finite-throat ceiling is algebraic, not exponential.

---

## 6. Best current theorem statement after Stage 39

### What is exact now

- the Stage-36 exponential source family is the exact stationary zero-flux solution of

  `partial_t sigma + partial_s(-D_sigma sigma' + v_sigma sigma) = 0`,

- the asymmetry parameter is the physical Peclet number

  `Pe = v_sigma L / D_sigma`,

- the overlap boost on that branch is

  `Omega_Pe`
  `= pi Pe (2 Pe exp(Pe) + pi) / [ (4 Pe^2 + pi^2)(exp(Pe)-1) ]`,

- and on the constructive branch `Pe >= 0` it is strictly increasing from

  `Omega_Pe(0)=1`

  to

  `lim_(Pe->+infinity) Omega_Pe = pi/2`.

### What this means physically

The source-shape asymmetry is no longer an abstract deformation parameter.
It is the axial transport Peclet number of the lowest support-source channel.

So one of the two Stage-38 “unknown physical inputs” has now been converted into a concrete moving-throat operator datum.

The remaining operator-level task is to combine that transport bias with the physical support-compliance ratios, so that the whole lowest-lane reachability problem is written directly in terms of real throat operator parameters rather than the old abstract pair `(alpha,x)`.
