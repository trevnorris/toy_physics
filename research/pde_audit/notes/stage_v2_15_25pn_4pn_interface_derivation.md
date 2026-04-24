# Stage V2-15 — 2.5PN / 4PN Quadrupole-Normalization Interface Audit

## 0. Purpose

This stage checks whether the conservative 4PN hereditary/tail sector introduces a new quadrupole normalization problem, or whether it is controlled by the same canonically normalized outgoing STF quadrupole branch already isolated at 2.5PN.

The result is algebraic:

\[
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR},
\qquad
\gamma_{\rm GR}=\frac{2G}{5c^5},
\]

so

\[
C_{\rm tail}^{\rm GR}
=
\frac{G^2M}{5c^8}.
\]

Thus, once the model’s effective quadrupole coefficient equals the 2.5PN Burke--Thorne coefficient, the 4PN hereditary coefficient is fixed.

The only extra scalar that can remain in a toy-model tail bridge is a **tail-transport** factor, not a new grouped-\(P_2\) normalization datum.

---

## 1. Inputs from earlier Volume 2 stages

From V2-13 and V2-14, the canonically normalized moving-throat outgoing STF quadrupole coefficient is

\[
\gamma_{\rm quad}^{\rm eff}
=
P_{\rm eff}\frac{a^5}{27c_s^5},
\]

where

\[
P_{\rm eff}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}.
\]

The universal 2.5PN target is

\[
\gamma_{\rm quad}^{\rm eff}
=
\gamma_{\rm GR}
=
\frac{2G}{5c^5}.
\]

Equivalently,

\[
P_{\rm eff}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

This is the same target previously written as

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

---

## 2. GR tail coefficient bridge

The conservative 4PN tail coefficient has the GR reference value

\[
C_{\rm tail}^{\rm GR}
=
\frac{G^2M}{5c^8}.
\]

Using the 2.5PN coefficient

\[
\gamma_{\rm GR}
=
\frac{2G}{5c^5},
\]

one obtains

\[
\frac{GM}{2c^3}\gamma_{\rm GR}
=
\frac{GM}{2c^3}\frac{2G}{5c^5}
=
\frac{G^2M}{5c^8}
=
C_{\rm tail}^{\rm GR}.
\]

So the exact interface identity is

\[
\boxed{
C_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\gamma_{\rm GR}.
}
\]

---

## 3. Moving-throat / toy-model tail bridge

A conservative symbolic bridge for the model is

\[
C_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}\gamma_{\rm quad}^{\rm eff}.
\]

Here:

- \(\gamma_{\rm quad}^{\rm eff}\) is the same STF quadrupole coefficient used at 2.5PN;
- \(c_s\) is the propagation/sound speed used in the outgoing port normalization;
- \(\Theta_{\rm tail}\) is a scalar monopole-scattering/tail-transport factor.

The ratio to the GR tail coefficient is

\[
\frac{C_{\rm tail}^{\rm toy}}{C_{\rm tail}^{\rm GR}}
=
\Theta_{\rm tail}
\left(\frac{c}{c_s}\right)^3
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}.
\]

Therefore, if the 2.5PN target is already satisfied,

\[
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}=1,
\]

the remaining tail condition is

\[
\boxed{
\Theta_{\rm tail}
\left(\frac{c}{c_s}\right)^3
=
1.
}
\]

On the \(c_s=c\) branch, this reduces to

\[
\boxed{
\Theta_{\rm tail}=1.
}
\]

So 4PN does not add another grouped-quadrupole normalization. It adds, at most, a scalar tail-transport gate.

---

## 4. Constant-prefactor check

The outgoing \(l=2\) fingerprint is

\[
\widehat Y_2^{\rm out}
=
1
+
\frac{a^2\omega^2}{9c_s^2}
+
\frac{4a^4\omega^4}{81c_s^4}
+
i\frac{a^5\omega^5}{27c_s^5}
+
O(\omega^6).
\]

Let the internal prefactor be

\[
P(\omega)=P_0+P_2\omega^2+P_4\omega^4+\cdots.
\]

Then

\[
P(\omega)\widehat Y_2^{\rm out}(\omega)
\]

has leading imaginary odd coefficient

\[
i\omega^5
\left[
P_0\frac{a^5}{27c_s^5}
\right].
\]

The \(P_2\) and \(P_4\) prefactor moments first affect higher odd orders:

\[
P_2\omega^2\cdot i\omega^5=O(i\omega^7),
\qquad
P_4\omega^4\cdot i\omega^5=O(i\omega^9).
\]

So the leading 2.5PN / 4PN interface uses only the static prefactor \(P_0=N_0/D_0\), not the higher even prefactor moments.

---

## 5. Residual bookkeeping

Define fractional deviations

\[
\frac{\gamma_{\rm quad}^{\rm eff}}{\gamma_{\rm GR}}=1+\delta_Q,
\]

and

\[
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1+\delta_{\rm tail}.
\]

Then

\[
\frac{C_{\rm tail}^{\rm toy}}{C_{\rm tail}^{\rm GR}}-1
=
(1+\delta_Q)(1+\delta_{\rm tail})-1
=
\delta_Q+\delta_{\rm tail}+\delta_Q\delta_{\rm tail}.
\]

At linear order,

\[
\frac{\Delta C_{\rm tail}}{C_{\rm tail}^{\rm GR}}
\simeq
\delta_Q+\delta_{\rm tail}.
\]

This is a useful diagnostic split:

- \(\delta_Q\) is the quadrupole-normalization miss already seen at 2.5PN;
- \(\delta_{\rm tail}\) is the tail-transport/scattering miss.

They should not be conflated.

---

## 6. Dimension ledger

Use dimensions in \((M,L,T)\).

\[
[G]=M^{-1}L^3T^{-2},
\qquad
[c]=[c_s]=LT^{-1}.
\]

Then

\[
[\gamma_{\rm GR}]
=
\left[\frac{G}{c^5}\right]
=
M^{-1}L^{-2}T^3.
\]

The static outgoing prefactor target has

\[
[P_{\rm eff}]
=
\left[
\frac{Gc_s^5}{a^5c^5}
\right]
=
M^{-1}L^{-2}T^{-2}.
\]

Multiplying by

\[
\left[\frac{a^5}{c_s^5}\right]=T^5
\]

gives

\[
[\gamma_{\rm quad}^{\rm eff}]
=
M^{-1}L^{-2}T^3,
\]

matching \(\gamma_{\rm GR}\).

For the tail coefficient,

\[
[C_{\rm tail}^{\rm GR}]
=
\left[
\frac{G^2M}{c^8}
\right]
=
M^{-1}L^{-2}T^4.
\]

And

\[
\left[\frac{GM}{c^3}\right]=T,
\]

so

\[
\left[
\frac{GM}{2c^3}\gamma_{\rm GR}
\right]
=
M^{-1}L^{-2}T^4.
\]

The bridge is dimensionally consistent.

---

## 7. SymPy audit result

The audit script verifies:

1. \(P_{\rm eff}=54Gc_s^5/(5a^5c^5)\) converts \(\gamma_{\rm quad}^{\rm eff}\) into \(\gamma_{\rm GR}\).
2. \(C_{\rm tail}^{\rm GR}=(GM/2c^3)\gamma_{\rm GR}\).
3. The toy/GR tail ratio factorizes as
   \[
   \Theta_{\rm tail}(c/c_s)^3(\gamma_{\rm quad}^{\rm eff}/\gamma_{\rm GR}).
   \]
4. On \(c_s=c\), \(\Theta_{\rm tail}=1\), and the 2.5PN target, the toy tail equals the GR tail.
5. The leading \(i\omega^5\) coefficient ignores \(P_2\) and \(P_4\).
6. The tail interface ignores \(P_2\) and \(P_4\).
7. The fractional residual factorizes as
   \[
   \delta_Q+\delta_{\rm tail}+\delta_Q\delta_{\rm tail}.
   \]
8. All dimension checks pass.

The script reports:

```text
symbolic_checks_total: 9
symbolic_checks_passed: 9
dimension_checks_total: 3
dimension_checks_passed: 3
FINAL_STATUS: PASS with one explicit tail-transport gate Theta_tail*(c/c_s)^3 = 1
```

---

## 8. Carry-forward theorem statement

The clean V2-15 theorem statement is:

\[
\boxed{
C_{\rm tail}
=
\frac{GM}{2c^3}\gamma_{\rm quad}^{\rm eff}
}
\]

on the canonical \(c_s=c\), \(\Theta_{\rm tail}=1\) branch.

Therefore,

\[
\boxed{
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}
\quad\Longleftrightarrow\quad
C_{\rm tail}
=
\frac{G^2M}{5c^8}.
}
\]

So the conservative 4PN tail does not create a new grouped-\(P_2\) normalization target. It reuses the same 2.5PN STF quadrupole normalization.

The remaining PDE work is now split into two scalar gates:

\[
\boxed{
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5}
}
\]

and

\[
\boxed{
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1.
}
\]

On the \(c_s=c\) branch, the second gate is simply \(\Theta_{\rm tail}=1\).
