G0_DIMCHECK_PASS

# pathA_25 G0 Dimensional Homogeneity Report

Scope: G0 freeze arithmetic only. No smectic ground state, spectrum, light-mode count, boundedness, leak, charge, cone-lock, or throat gate was solved.

## Verdict

`G0_DIMCHECK_PASS`

Engine agreement: PASS. SymPy and Mathematica both:

- recomputed the T0 freeze hash `8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064`;
- recomputed the G0 combined freeze hash `f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290`;
- verified the exact T0 freeze-action bytes occur unchanged inside the G0 combined block;
- reduced every new frozen bulk-density term to `M L^-2 T^-2`;
- reduced every layer term to `M L^-1 T^-2` and its `delta_Sigma` bulk representation to `M L^-2 T^-2`;
- verified `k_Lstar` and `k_Rstar` carry `L^-1`;
- replaced the old Family-C monomial `k->0` assertion with a derived density/P response: Cdiv shifts the directional low-`k` EOS stiffness after Goldstone integrate-out, while Cpin remains an `O(k^2)` anisotropic density Hessian and preserves `A_rho(0)`.

## Artifacts

- SymPy checker: `software/stage1_solver/tools/pathA_25_G0_dimcheck_sympy.py`
- Mathematica checker: `software/stage1_solver/tools/pathA_25_G0_dimcheck.wl`
- SymPy JSON: `software/stage1_solver/_scratch/pathA_25_G0_dimcheck_sympy.json`
- Mathematica JSON: `software/stage1_solver/_scratch/pathA_25_G0_dimcheck_mathematica.json`

Run commands from `/var/projects/toy_physics`:

```text
timeout 600 python3 software/stage1_solver/tools/pathA_25_G0_dimcheck_sympy.py
timeout 600 math -script software/stage1_solver/tools/pathA_25_G0_dimcheck.wl
```

Both commands exited `0`.

## Unit Conventions

Base dimensions are exponent triples `(M,L,T)`.

```text
[rho] = L^-4
[m] = M
[K] = M L^18 T^-2
[a] = L
[u] = L
[P^i] = 1
[delta_Sigma] = L^-1
[c_s^2(rho)=5K rho^4/m] = L^2 T^-2
```

Bulk action density target:

```text
M L^-2 T^-2
```

Layer density target before the `delta_Sigma` representation:

```text
M L^-1 T^-2
```

## New-Term Dimensional Results

Family L:

```text
c_L1 (partial_i rho)^2                  -> M L^-2 T^-2
c_L2 (partial_i partial_i rho)^2        -> M L^-2 T^-2
k_Lstar=sqrt(c_L1/(2 c_L2))             -> L^-1
```

Family R:

```text
V_R(|R|) real-space kernel              -> M L^2 T^-2
A_R Fourier-kernel amplitude            -> M L^6 T^-2
int d^4R delta_rho V_R delta_rho        -> M L^-2 T^-2
k_R, k_Rstar                            -> L^-1
```

Family C:

```text
lambda_Cdiv delta_rho partial_i P^i     -> M L^-2 T^-2
chi_Cpin (P^i partial_i rho)^2          -> M L^-2 T^-2
```

Light-sector layer terms:

```text
varrho_br (D_t u)^2                     -> M L^-1 T^-2
mu_br (curl u)^2                        -> M L^-1 T^-2
I_PSigma (D_t P_parallel)^2             -> M L^-1 T^-2
K_PSigma (partial_a P_parallel)^2       -> M L^-1 T^-2
G_PSigma (delta P_parallel)^2           -> M L^-1 T^-2
J_Pu (D_t delta P - D_t Omega_u)^2      -> M L^-1 T^-2
kappa_Pu (delta P - Omega_u)^2          -> M L^-1 T^-2
```

Multiplying each layer term by `delta_Sigma` gives the required bulk action density `M L^-2 T^-2`.

Light-sector constants/functions:

```text
varrho_br=int_layer dn m rho                         -> M L^-3
mu_br                                                -> M L^-1 T^-2
I_PSigma=int_layer dn m rho a^2                      -> M L^-1
K_PSigma=int_layer dn m rho c_s^2 a^2                -> M L T^-2
G_PSigma=int_layer dn m rho c_s^2                    -> M L^-1 T^-2
J_Pu                                                 -> M L^-1
kappa_Pu                                             -> M L^-1 T^-2
```

## k->0 Limit

The frozen baseline density driver has quadratic symbol

```text
Delta E_L(k) = c_L2 k^4 - c_L1 k^2.
```

Both engines verify:

```text
lim_{k->0} Delta E_L(k) = 0
k_Lstar = sqrt(c_L1/(2 c_L2))
```

The Family R sensitivity shape is

```text
f_R(x) = (x^4 - 2 x^2) exp(-x^2)
tilde V_R(k) = A_R f_R(k/k_R)
```

Both engines verify:

```text
tilde V_R(0) = 0
tilde V_R(k) = -2 A_R k^2/k_R^2 + O(k^4)
k_Rstar = k_R sqrt(2 - sqrt(2))
```

Family C is no longer checked by the superseded tautological monomials `lambda_Cdiv*k` and `chi_Cpin*k^2`. Both engines now derive the quadratic density response about uniform `rho0` and unit `P0`, retain the Family-L kernel, and integrate out the P fluctuations before taking the directional low-`k` limit.

```text
S_Goldstone = 5 K a^2 rho0^5
G_magnitude = 10 K rho0^5

Cdiv mixed block:
  d^2E/(d eta d pi_T) = -k lambda_Cdiv sin(theta)
  d^2E/(d eta d sigma) = -k lambda_Cdiv cos(theta)

Delta A_Cdiv =
  -lambda_Cdiv^2 sin(theta)^2/(5 K a^2 rho0^5)
  -k^2 lambda_Cdiv^2 cos(theta)^2/(5 K rho0^5 (a^2 k^2 + 2))

lim_{k->0} [A_Cdiv(k,theta)-5 K rho0^3] =
  -lambda_Cdiv^2 sin(theta)^2/(5 K a^2 rho0^5)

liminf_direction shift =
  -lambda_Cdiv^2/(5 K a^2 rho0^5)
```

Thus Cdiv is a directional low-wavelength admission failure for any nonzero `lambda_Cdiv` and becomes phase-separating if the shift overwhelms `5 K rho0^3`.

Cpin enters through a direct density Hessian, with no rho-P quadratic mixed block:

```text
d^2E_Cpin/d eta^2 = chi_Cpin k^2 cos(theta)^2
lim_{k->0} A_Cpin(k,theta) = 5 K rho0^3
```

Therefore Cpin preserves the strict EOS stiffness while adding an anisotropic `O(k^2)` shell-selection term. Family L and R conclusions above are unchanged.
