# Path-A Derivation 01: Return Sources and Static Balance

This note is target-blind. It derives only the structural return-source and
static-balance statements needed for the promoted Path-A program. No numerical
branch solve or target comparison is used.

## Grounding

- The current parent action contains GNLS matter plus localized Maxwell, with
  geometry entering through the confinement argument
  `Sigma = r - R(Omega,w,t)`; an autonomous throat PDE requires adding
  `S_Sigma[R]` to the action
  (`notes/moving_throat_pde_program_compact.md:196-218`,
  `notes/moving_throat_pde_program_compact.md:408-431`,
  `notes/moving_throat_pde_program_compact.md:480-505`).
- The promoted compact wall action is
  `S_Sigma = int L_Sigma`, with `mu_Sigma`, `T_w,Sigma`,
  `T_Omega,Sigma`, and `U_Sigma`; its quadratic stiffness is the boxed
  `K_eta` expression in
  `notes/moving_throat_pde_program_compact.md:510-528`.
- The only explicit matter confinement interaction is
  `-V_conf(X;Sigma) rho` in `L_psi`
  (`notes/moving_throat_pde_program_compact.md:556-565`).
- The localized Maxwell sector uses `Z(w)` and optional gauge-fixing weight
  `H(w)`, not a declared `Z(R,w)`
  (`notes/moving_throat_pde_program_compact.md:590-620`,
  `notes/moving_throat_pde_program_compact.md:678-685`).
- The reduced wall equation is the displayed source equation with
  `S_eta^(psi) + S_eta^(A) + f_ext` on the right
  (`notes/moving_throat_pde_program_compact.md:1259-1305`,
  `notes/moving_throat_pde_program_compact.md:1383-1451`).
- The forward confinement perturbation is fixed by
  `delta Sigma = -eta` and
  `delta V_conf = -(V_wall'(Sigma0/ell_c)/ell_c) eta`
  (`notes/moving_throat_pde_program_compact.md:1059-1089`,
  `notes/moving_throat_pde_program_compact.md:1424-1429`).
- The prior ledger gives the leading matter-return form only schematically
  (`research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md:344-356`).
- The reduced mixed-sector algebra defines
  `Delta`, `Q`, `P`, `D0`, and `N0` in
  `notes/moving_throat_pde_program_compact.md:3946-4000`, and the reduced
  audit carries the Schur-boundary stability gate in
  `research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:230-263`
  and
  `research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:497-507`.

## Conventions

Let

```text
k1 := d V_conf/d Sigma |0 = V_wall'(Sigma0/ell_c)/ell_c,
k2 := d^2 V_conf/d Sigma^2 |0 = V_wall''(Sigma0/ell_c)/ell_c^2,
rho = rho0 + delta rho,
R = R0 + eta,
Sigma = Sigma0 - eta.
```

The displayed wall-PDE source convention follows the confinement-potential
kernel used in the matter equation. Because the Lagrangian contains
`-V_conf rho`, while the matter equation contains `+V_conf`, the Lagrangian
cross-Hessian and the displayed confinement/source kernel differ by the usual
`L = kinetic - energy` sign. The symbolic verifier checks both ledgers:
the Taylor series of the declared Lagrangian, and the symmetric
Hamiltonian/equation-kernel bilinear that enters the forward and return
equations.

## D1: Matter Return Source

Taylor-expand the declared interaction:

```text
L_int = -V_conf(Sigma0 - eta) (rho0 + delta rho)
      = L0 + L1
        + k1 eta delta rho
        - (1/2) k2 rho0 eta^2
        + O(cubic).
```

The quadratic part fixed by the declared action is therefore

```text
L_int^(2) = k1 eta delta rho - (1/2) k2 rho0 eta^2.
```

Equivalently, in the confinement-energy/equation-kernel ledger,

```text
H_int^(2) = -k1 eta delta rho + (1/2) k2 rho0 eta^2.
```

The forward wall-to-matter term is

```text
delta V_conf = -k1 eta.
```

The raw linear matter-to-wall source in the displayed wall-PDE convention is

```text
S_eta^(psi),raw = -k1 delta rho + k2 rho0 eta.
```

The first term is the genuine matter-return cross term:

```text
S_eta^(psi),cross = -k1 delta rho
                  = -(V_wall'(Sigma0/ell_c)/ell_c) delta rho.
```

The second term is not a new matter-return magnitude. It is the local
self-stiffness contribution corresponding to the quadratic Lagrangian term
`-(1/2) k2 rho0 eta^2`. If it is absorbed into the left-hand wall operator,
the same displayed-source convention shifts the local coefficient by the
opposite of the right-hand self term:

```text
K_eta,eff = K_eta - k2 rho0.
```

Thus the completed quadratic-order statement is

```text
S_eta^(psi) = -k1 delta rho,
```

with the `k2 rho0 eta` piece tracked as a wall-stiffness renormalization, not
as an independent return channel. The quadratic action fixes no `delta rho^2`
term and no higher cross term. Terms such as `eta delta rho` in the source
itself would require cubic order in the interaction expansion and are not
claimed here.

Reciprocity is the single load-bearing point. The symmetric bilinear

```text
H_cross = (-k1) eta delta rho
```

has

```text
d H_cross/d(delta rho) = -k1 eta,
d H_cross/d eta       = -k1 delta rho,
d^2 H_cross/(d eta d(delta rho)) = -k1.
```

The forward and return kernels are therefore the same symmetric cross-Hessian.
Closing the return direction introduces no new coupling magnitude beyond the
already-declared confinement kernel.

## D2: Gauge Return Source

At fixed `psi` and `A`, the declared Maxwell action has no direct `eta`
variation: its localization and gauge-fixing weights are functions of `w`
only (`notes/moving_throat_pde_program_compact.md:590-620`,
`notes/moving_throat_pde_program_compact.md:678-685`). The exact matter current
is generated by varying the covariant matter action and must not be
double-counted as an independent explicit Maxwell source
(`notes/moving_throat_pde_program_compact.md:622-629`,
`notes/moving_throat_pde_program_compact.md:651-685`).

Therefore the declared action fixes no direct `eta -> Z` gauge-wall source.
The gauge return is matter-mediated:

```text
eta -> delta V_conf -> delta psi, delta rho, delta J_psi -> delta A
```

in the forward direction, and the return is the adjoint variation of that same
on-shell matter/gauge functional with respect to `eta`. A compact way to name
the fixed part is

```text
S_eta^(A) := delta Gamma_gauge,med[R] / delta eta,
```

where `Gamma_gauge,med` is the gauge-field part of the on-shell effective
action induced by the covariant matter current. This names the source without
inventing a closed local formula that the declared action has not supplied.

The admissible reduced gauge data are the covariant current and the exact mixed
gauge invariants

```text
E_w = -partial_t A_w - partial_w A_0,
C_a =  partial_a A_w - partial_w A_a,
```

which are gauge invariant under the transformations in
`notes/moving_throat_pde_program_compact.md:769-786` and
`research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:19-47`.
The verifier checks these identities symbolically.

One boundary caveat remains. The compact program records an open finite-radius
conduit with endpoint work/leakage terms
(`notes/moving_throat_pde_program_compact.md:310-325`), but the declared
Maxwell action is written on the fixed parent integration domain with `Z(w)`.
A future shape-calculus promotion of moving electromagnetic boundaries could
create additional boundary terms. Such terms are not present in the declared
compact action and are left open here.

## D3: Static Self-Consistent Throat Balance

The promoted wall action adds autonomous throat inertia, axial stiffness,
angular stiffness, and a restoring potential through
`mu_Sigma R_t^2`, `T_w,Sigma R_w^2`, `T_Omega,Sigma |nabla_Omega R|^2`, and
`U_Sigma` (`notes/moving_throat_pde_program_compact.md:510-528`). For a
stationary isotropic background `R0(w)`, the time and angular terms vanish in
the static equation, and the Euler derivative of the wall energy gives

```text
-partial_w(T_w,Sigma(R0,w) partial_w R0)
+ (1/2) T_w,Sigma,R(R0,w) (partial_w R0)^2
+ U_Sigma,R(R0,w)
= S_eta^(psi,A)|static[rho0(R0), A0(R0)].
```

The right-hand side contains the static confinement/matter force plus the
matter-mediated gauge contribution, evaluated on the self-consistent static
matter/gauge solution for the same `R0`. This replaces the previously frozen
smooth throat profile with a derived equilibrium; it is not a free post-hoc
choice.

Quadratic reduction around that same `R0` gives

```text
mu_eta    = mu_Sigma(R0,w),
T_w       = T_w,Sigma(R0,w),
T_Omega   = T_Omega,Sigma(R0,w),
K_eta     = U_Sigma,RR(R0,w)
            - partial_w(T_w,Sigma,R(R0,w) R0')
            + (1/2) T_w,Sigma,RR(R0,w) (R0')^2,
```

matching the compact stiffness formula
(`notes/moving_throat_pde_program_compact.md:523-528`) and the effective
quadratic wall action
(`notes/moving_throat_pde_program_compact.md:1269-1305`). The promotion adds
parent throat dynamics; it does not derive the constitutive functions
`mu_Sigma`, `T_w,Sigma`, `T_Omega,Sigma`, or `U_Sigma` from a deeper material
law.

Once solved, the derived `R0(w)` feeds the same downstream coefficient chain:
BdG support coefficients, wall `K,M`, gauge/mixed coefficients, and the
Schur denominator

```text
D0 = K_* - Q/Delta.
```

No new algebraic transfer channel is added by closing the return loop. This is
a consequence of the reciprocity result in D1: the return kernels are the
adjoints of the already-declared forward kernels, so closing the loop adds no
independent free magnitude.

## D4: No New Numerator Magnitude

The reduced mixed algebra defines

```text
Delta = Omega_U^2 Omega_W^2 - R_mix^2,
P     = Omega_U^2 G_W + R_mix G_U,
N0    = (P/Delta)^2,
D0    = K_* - Q/Delta.
```

The return sources derived above do not add a term to `P`. This
no-new-numerator statement is not established by differentiating newly named
return symbols out of the displayed `P`; that would only restate how `P` was
written. It follows from the genuine reciprocity check in D1. The forward
confinement/gauge kernels and their return adjoints share the same
cross-Hessian, so the return direction contributes no additive independent
coupling to the forward numerator built from `G_U`, `G_W`, `Omega_U^2`, and
`R_mix`.

The return sources can change:

1. the self-consistent background `(rho0, A0, R0)`, on which the already
   declared coefficients are evaluated;
2. the conservative wall/gauge Schur denominator through the realized
   background and wall stiffness.

They do not introduce a new free numerator magnitude. The verifier therefore
reports `no_numerator_knob` as a construction-restatement/corollary, not as an
independent physics gate. A genuinely independent numerator audit would require
the forward background-to-port-coefficient map from the mixed-Maxwell
reduction, which is outside this return-source derivation. The only structural
softening lever in this closed promoted system is the denominator/background
side, with the branch required to remain on the positive side of the Schur
boundary (`research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:230-263`,
`research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:497-507`).

## Posited or Open

- The nonlinear constitutive wall functions in `S_Sigma` remain posited. This
  derivation only shows that their quadratic reduction reproduces the effective
  wall coefficients.
- The matter return is complete only to quadratic order in the declared
  confinement interaction. Cubic and higher source terms are not fixed here.
- The gauge return is named as the matter-mediated on-shell variation fixed by
  the covariant matter/current action. A closed local formula requires the full
  closed matter/gauge solve and is not supplied by the compact action alone.
- Direct electromagnetic moving-boundary terms are not in the declared fixed
  `Z(w)` action. If a future shape-calculus boundary action is declared, those
  terms must be derived separately.
- Boundary terms from integrating the wall action by parts depend on the chosen
  endpoint class and impedance/work bookkeeping.
- This note does not solve the static balance, prove existence or uniqueness,
  or determine whether the self-consistent branch lies near a Schur boundary.
