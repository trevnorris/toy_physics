# Stage V2-03 — Dimension Ledger and Normalization Audit

## Theorem target

This stage audits the unit consistency of the normalization packets used by the moving-throat grouped `P2` and outgoing-quadrupole bridge. The main targets are

\[
P_0^{\rm target}=\frac{54G c_s^5}{5a^5c^5},
\]

\[
\Lambda_0=\frac{27\pi^2G c_s^5}{20a^5c^5},
\]

and the invariant outgoing-quadrupole condition

\[
\widehat m_0^2 P_0\frac{a^5}{27c_s^5}=\frac{2G}{5c^5}.
\]

The audit also checks the grouped-response conversion formulas, constant-prefactor branch equations, EM zero-mode coupling dimensions, support-wave mass theorem, and the raw Maxwell/mixed one-lane transfer factor.

## Base dimensions

The script uses the exact base vector

\[
(M,L,T,Q,O),
\]

where `O` is an abstract reduced wall/operator unit. The `O` axis is kept separate because the absolute dimension of the reduced wall operator depends on the chosen normalization of the wall amplitude. Physical observables should not retain an uncancelled `O` unless they are still operator-valued.

The primitive assignments are

\[
[G]=M^{-1}L^3T^{-2},
\qquad
[c]=[c_s]=LT^{-1},
\qquad
[a]=L,
\qquad
[\omega]=T^{-1},
\]

\[
[\hbar]=ML^2T^{-1},
\qquad
[E]=ML^2T^{-2}.
\]

For the EM zero-mode ledger the script assumes a dimensionless localization profile `Z(w)`, so

\[
[Z_{\rm int}]=L.
\]

Then

\[
\mu_0^{\rm eff}=\frac{\mu_0^{(5)}}{Z_{\rm int}}
\]

is consistent if

\[
[\mu_0^{(5)}]=[\mu_0^{\rm eff}]L.
\]

Similarly,

\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}}
\]

is consistent if

\[
[q_\star]=[q_{\rm eff}]L^{1/2}.
\]

## Main target dimensions

The target prefactor has dimension

\[
\left[\frac{G c_s^5}{a^5c^5}\right]
= M^{-1}L^{-2}T^{-2}.
\]

Therefore

\[
[P_0^{\rm target}]=[\Lambda_0]=M^{-1}L^{-2}T^{-2}.
\]

The outgoing GR quadrupole coefficient has dimension

\[
\left[\frac{G}{c^5}\right]=M^{-1}L^{-2}T^3.
\]

The compact outgoing `l=2` factor has

\[
\left[\frac{a^5}{c_s^5}\right]=T^5.
\]

Therefore

\[
\left[P_0^{\rm target}\frac{a^5}{c_s^5}\right]
= M^{-1}L^{-2}T^3
=\left[\frac{G}{c^5}\right].
\]

So the central 2.5PN/4PN normalization product is dimensionally consistent.

## Outgoing fingerprint

The normalized outgoing branch is

\[
\widehat Y_2^{\rm out}(\omega)
=1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}+O(\omega^6).
\]

The script verifies

\[
\left[\frac{a^2\omega^2}{c_s^2}\right]
=
\left[\frac{a^4\omega^4}{c_s^4}\right]
=
\left[\frac{a^5\omega^5}{c_s^5}\right]
=1.
\]

So the outgoing fingerprint is dimensionless, as required.

## Operator and response moments

Let

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6).
\]

The script assigns

\[
[D_0]=O,
\qquad
[D_2]=OT^2,
\qquad
[D_4]=OT^4,
\]

so that every term in `D(omega)` has unit `O`.

The normalized response

\[
Y(\omega)=\frac{D_0}{D(\omega)}
=1+u_2\omega^2+u_4\omega^4+O(\omega^6)
\]

has

\[
[u_2]=T^2,
\qquad
[u_4]=T^4.
\]

SymPy verifies the exact coefficients

\[
u_2=-\frac{D_2}{D_0},
\]

\[
 u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The one-pole condition

\[
 u_4=4u_2^2
\]

is also dimensionally consistent because both sides have unit `T^4`.

## Outgoing prefactor moments

The Stage-5 prefactor is

\[
\mathrm{Pref}(\omega)
=\frac{D_0N(\omega)}{D(\omega)^2},
\qquad
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

The gravitationally normalized assignment is

\[
[N_0]=[P_0]O,
\qquad
[N_2]=[P_0]OT^2,
\qquad
[N_4]=[P_0]OT^4.
\]

Then

\[
[P_0]=M^{-1}L^{-2}T^{-2},
\qquad
[P_2]=M^{-1}L^{-2},
\qquad
[P_4]=M^{-1}L^{-2}T^2.
\]

SymPy verifies

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

The constant-prefactor branch equations are dimensionally consistent:

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]

## Isotropic one-pole surface

The isotropic full-bundle one-pole surface is

\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
\]

With

\[
[B_4+Z_4]=OT^4,
\qquad
[M+B_2+Z_2]=OT^2,
\]

both sides have unit

\[
O^2T^4.
\]

So the one-pole target surface is dimensionally consistent.

## Raw Maxwell/mixed transfer warning

The only nontrivial warning appears when the Stage-4 one-lane Maxwell/mixed transfer factor is interpreted as a raw canonical mechanical transfer.

For canonical internal coordinates `U,W`, the reduced kinetic terms imply

\[
[U]=[W]=E^{1/2}T.
\]

If the wall displacement coordinate has unit `L`, then the wall operator has

\[
[D_{\rm mech}]=E/L^2=M.
\]

The raw one-lane transfer factor

\[
N_0^{\rm raw}
=\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2}
\]

has dimension

\[
[N_0^{\rm raw}]=[D_{\rm mech}]T^2.
\]

Therefore

\[
P_0^{\rm raw}=\frac{N_0^{\rm raw}}{D_{\rm mech}}
\]

has unit

\[
[P_0^{\rm raw}]=T^2,
\]

not

\[
[P_0^{\rm target}]=M^{-1}L^{-2}T^{-2}.
\]

This is not a contradiction in the reduced algebra. It means the raw canonical transfer factor is not yet the gravitationally normalized `P0` used in the 2.5PN/4PN target. A port/source normalization scale is required:

\[
P_0^{\rm grav}=\mathcal S_{\rm port}P_0^{\rm raw},
\]

with

\[
[\mathcal S_{\rm port}]
=M^{-1}L^{-2}T^{-4}.
\]

Equivalently,

\[
N_0^{\rm grav}=\mathcal S_{\rm port}N_0^{\rm raw}.
\]

Volume 2 should therefore avoid identifying raw canonical `N0/D0` with the gravitational outgoing prefactor unless this port/source normalization is explicitly included.

## SymPy audit summary

The script produced:

- 29 exact dimension checks,
- 29 passes,
- 8 SymPy series-identity checks,
- 8 passes.

Key output dimensions:

```text
P0_target_dim = M^-1 L^-2 T^-2
Gamma_GR_dim = M^-1 L^-2 T^3
D0_dim = O
N0_normalized_dim = M^-1 L^-2 T^-2 O
raw_mixed_P0_dim = T^2
required_raw_to_grav_bridge_scale_dim = M^-1 L^-2 T^-4
```

## Verdict

```text
PASS for published target identities and response-moment algebra.
WARN that the raw canonical Maxwell/mixed transfer factor is a mechanical T^2 object and must be multiplied by an explicit port/source normalization scale to become the gravitational P0 used in the 2.5PN/4PN target.
```

## Carry-forward patch

Introduce an explicit normalization bridge in the Volume 2 notation:

\[
\boxed{
P_0^{\rm grav}
=\mathcal S_{\rm port}\,P_0^{\rm raw},
\qquad
[\mathcal S_{\rm port}]=M^{-1}L^{-2}T^{-4}.
}
\]

Then the actual theorem target becomes

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}P_0^{\rm raw}
=\frac{54Gc_s^5}{5a^5c^5}.
}
\]

If later branch work defines `N0` directly in the gravitationally normalized port convention, then `S_port` is implicitly already absorbed. The derivation should say which convention is being used.
