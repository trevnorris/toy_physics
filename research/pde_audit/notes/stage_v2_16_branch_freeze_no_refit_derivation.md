# Stage V2-16 — Branch-freeze / no-refit protocol

## Purpose

This stage turns the anti-overfitting rule into an executable theorem gate.

The problem is not that the grouped `P2` algebra lacks enough variables. It has many. The problem is that, unless the branch is frozen before comparison, the same algebra contains enough freedom to make target matching look stronger than it is.

The rule for Volume 2 is therefore:

> Define the parent action, gauge convention, wall/interface action, open-exit boundary protocol, projection/source map, support-profile family, number of modes/ports, stability gates, and extraction formulas **before** evaluating any target residual.

Then, and only then, evaluate the residual packet.

---

## 1. Frozen protocol DAG

The stage encodes the protocol as a directed acyclic graph:

```text
parent_action
  -> gauge_convention
  -> wall_action_and_geometry
  -> open_boundary_protocol
  -> projection_and_source_map
  -> support_profile_family
  -> branch_solve
  -> coefficient_extraction
  -> target_residual_evaluation
  -> target_decision
```

The audit verifies three graph facts:

1. all allowed arrows point forward in the frozen order;
2. there is no path from target residuals or target decisions back into branch definitions;
3. a deliberately bad edge,
   ```text
   target_residual_evaluation -> support_profile_family
   ```
   is detected as a protocol violation.

This is the formal no-refit rule.

---

## 2. Algebraic residual packet

On the isotropic grouped `P2` branch define

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The normalized conservative response moments are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The outgoing prefactor moments are

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The residual packet used after freeze is:

### 2.1 One-pole conservative residual

\[
R_{\rm pole}
=
D_0(B_4+Z_4)
-
3(M+B_2+Z_2)^2.
\]

The script verifies

\[
(u_4-4u_2^2)D_0^2
=
R_{\rm pole}.
\]

So the one-pole condition is

\[
R_{\rm pole}=0.
\]

### 2.2 Constant-prefactor residuals

The constant-prefactor branch is

\[
P_2=0,
\qquad
P_4=0.
\]

The script verifies that these are equivalent to

\[
N_2=\frac{2D_2N_0}{D_0},
\]

and

\[
N_4=
\frac{
2D_0(D_2N_2+D_4N_0)-3D_2^2N_0
}{D_0^2}.
\]

### 2.3 Universal quadrupole normalization residual

With the port/source convention carried explicitly,

\[
P_{\rm eff}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}.
\]

The target residual is

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5}.
\]

The script verifies that

\[
R_{\rm norm}=0
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5},
\]

because

\[
\gamma_{\rm quad}^{\rm eff}
=
P_{\rm eff}\frac{a^5}{27c_s^5}.
\]

### 2.4 4PN tail-transport residual

If the tail bridge is written as

\[
C_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}
\gamma_{\rm quad}^{\rm eff},
\]

then after the 2.5PN quadrupole normalization is met, the remaining tail residual is

\[
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

The script verifies that

\[
R_{\rm tail}=0
\]

makes the toy tail coefficient match

\[
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR}.
\]

---

## 3. Incidence matrix

The residual rows are

\[
(R_{\rm pole},R_{P2},R_{P4},R_{\rm norm},R_{\rm tail}).
\]

The branch/evaluation columns are

\[
(K,M,B_0,B_2,B_4,Z_0,Z_2,Z_4,N_0,N_2,N_4,\widehat m_0,\mathcal S_{\rm port},a,c_s,\Theta_{\rm tail}).
\]

The incidence matrix printed by the script is a dependency ledger, not a fit instruction. It shows which frozen branch data are used by which residual.

The key interpretation is:

- dependencies from branch data to residuals are allowed;
- dependencies from residuals or target decisions back to branch data are forbidden.

---

## 4. Why the no-refit rule is mandatory

The script computes the Jacobian of the five residuals with respect to the post-hoc knob set

\[
(K,N_2,N_4,N_0,\Theta_{\rm tail}).
\]

It finds

\[
\det
\frac{\partial(R_{\rm pole},R_{P2},R_{P4},R_{\rm norm},R_{\rm tail})}
{\partial(K,N_2,N_4,N_0,\Theta_{\rm tail})}
=
\frac{
\mathcal S_{\rm port}\widehat m_0^{\,2}(B_4+Z_4)
}{
D_0^3
}
\left(\frac{c}{c_s}\right)^3.
\]

Equivalently, as printed in the script’s native sign convention,

\[
-\frac{
\mathcal S_{\rm port}c^3\widehat m_0^{\,2}(B_4+Z_4)
}{
c_s^3(B_0-K+Z_0)^3
}.
\]

Under generic nonzero/stable conditions,

\[
D_0\neq0,
\qquad
B_4+Z_4\neq0,
\qquad
\widehat m_0\mathcal S_{\rm port}\neq0,
\qquad
c/c_s\neq0,
\]

this determinant is nonzero.

Therefore, if post-hoc fitting were allowed, five algebraic knobs could generically tune five residuals. That is the point of this stage: the algebra is powerful enough that target comparison is not trustworthy without a freeze certificate.

---

## 5. Freeze packet

The script creates a deterministic freeze packet and SHA256 certificate.

The packet requires freezing before target evaluation:

```text
parent action and current bookkeeping
gauge-fixing convention
wall/throat action or effective interface action
open-exit impedance boundary protocol
projection/source-map convention
support profile family and number of modes/ports
stability acceptance gates
coefficient extraction formulas
```

Only after these are frozen may the branch evaluate:

```text
one-pole residual
constant-prefactor residuals
universal quadrupole normalization
tail-transport scalar
weak-axisymmetric prefactor slope Xi_1
```

Forbidden after target evaluation:

```text
changing support-cardinality
changing boundary condition class
changing gauge convention
changing port/source normalization convention
dropping dark or unstable branches only after target miss unless this was predeclared
adding compensating modes because a target residual is nonzero
```

---

## 6. Stage result

The audit passes all symbolic and protocol checks.

The important warning is:

\[
\boxed{
\text{The target surface is algebraically tuneable unless the branch is frozen first.}
}
\]

So the correct Volume 2 wording is:

> A branch may be accepted only if the frozen pre-target packet produces stable coefficients satisfying the target residuals. A branch may not be rescued by changing mode support, boundary class, gauge convention, source normalization, or port structure after the residuals are known.

---

## 7. Carry-forward to the next stage

V2-17 can now use this freeze protocol for the weak-axisymmetric branch. The next audit should derive the axisymmetric splitting line

\[
(20,21,22)\sim(1,\tfrac12,-1),
\qquad
b=3a,
\]

and evaluate the transported prefactor slope

\[
\Xi_1=\frac{P_1}{P_0}
\]

without allowing post-target changes to the branch packet.
