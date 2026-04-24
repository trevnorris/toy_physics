# Stage V2-04 — Open-Junction / Organ-Pipe Impedance Audit

## Executive verdict

This audit **accepts the architectural patch** from a capped throat to an open finite-radius conduit:

\[
R(L)>0,
\]

with the bottom of the throat opening into unconfined four-spatial-dimensional bulk. A hard cap,

\[
R(L)=0,
\]

is incompatible with a nonzero baseline/DC superfluid throughput, because the exit area vanishes and a fixed flux would require a divergent flux density.

The audit gives a split result for the D/N support ladder:

\[
\textbf{D/N survives only if the support variable is the free-end flow/displacement variable.}
\]

A sudden expansion into a large bulk reservoir is a **low-impedance open end**. For a pressure-like or potential-like scalar, that open end reflects with phase \(\pi\),

\[
R_p\to -1,
\]

which is Dirichlet-like. For the dual flow/displacement variable, the reflected amplitude has the opposite sign,

\[
R_q=-R_p\to +1,
\]

which is Neumann-like:

\[
\frac{q_w(L)}{q(L)}\to0.
\]

So the organ-pipe patch is viable, but the Volume 2 wording must be precise:

> The open throat dynamically produces the D/N ladder for the **flow/displacement support coordinate**, not for a generic pressure-like scalar field.

---

## 1. Junction model

Use the variable-area scalar wave action

\[
S_\psi^{(2)}
=
\frac12\int dt\,dx\,
A(x)\left[
\frac{1}{c_s^2}\psi_t^2-\psi_x^2
\right].
\]

The Euler-Lagrange equation is

\[
\frac{A}{c_s^2}\psi_{tt}
-
\partial_x(A\psi_x)
=0.
\]

In divided form,

\[
\frac{1}{c_s^2}\psi_{tt}
-
\frac1A\partial_x(A\psi_x)
=0.
\]

For the 1D throat tube,

\[
A(x)=A_t=\text{constant},
\]

so

\[
\psi_{tt}=c_s^2\psi_{xx}.
\]

For the open 4D spatial bulk, the radial area scales as

\[
A_{\rm bulk}(r)=S_3r^3,
\qquad
S_3=2\pi^2,
\]

so the radial wave equation is

\[
\frac{1}{c_s^2}\Psi_{tt}
-
\frac{1}{S_3r^3}\partial_r(S_3r^3\Psi_r)
=0,
\]

or

\[
\frac{1}{c_s^2}\Psi_{tt}
-
\Psi_{rr}
-
\frac3r\Psi_r
=0.
\]

The conservative pressure/potential-like scalar matching conditions are continuity of the field and continuity of the weighted normal derivative:

\[
\psi_{\rm tube}(L)=\Psi_{\rm bulk}(r_e),
\]

\[
A_t\partial_w\psi_{\rm tube}(L)
=
A_{\rm bulk}(r_e)\partial_r\Psi_{\rm bulk}(r_e).
\]

For acoustic variables it is often clearer to use effort/flow impedance:

\[
p=Z_LQ
\]

at the load.

---

## 2. Reflection and transmission coefficients

Let

\[
Z_t
\]

be the tube characteristic impedance and

\[
Z_L
\]

the effective load impedance of the open bulk. Define

\[
\epsilon=\frac{Z_L}{Z_t}.
\]

For a pressure-like amplitude,

\[
R_p
=
\frac{Z_L-Z_t}{Z_L+Z_t}
=
\frac{\epsilon-1}{\epsilon+1},
\]

\[
T_p=1+R_p
=
\frac{2\epsilon}{1+\epsilon}.
\]

For the dual flow/displacement amplitude,

\[
R_q=-R_p
=
\frac{1-\epsilon}{1+\epsilon},
\]

\[
T_q=1+R_q
=
\frac{2}{1+\epsilon}.
\]

The energy reflection and transmission coefficients are

\[
\mathcal R_E
=
\left(\frac{\epsilon-1}{\epsilon+1}\right)^2,
\]

\[
\mathcal T_E
=
\frac{4\epsilon}{(1+\epsilon)^2}.
\]

A sudden expansion has

\[
\epsilon\ll1.
\]

Then

\[
R_p=-1+2\epsilon-2\epsilon^2+O(\epsilon^3),
\]

\[
R_q=1-2\epsilon+2\epsilon^2+O(\epsilon^3),
\]

\[
\mathcal T_E=4\epsilon+O(\epsilon^2).
\]

So the AC wave is strongly reflected, with only \(O(\epsilon)\) energy leakage per encounter.

If the effective open area is written as

\[
\chi=\frac{A_{\rm eff}}{A_t}\gg1,
\]

then

\[
\epsilon\simeq\chi^{-1},
\]

and

\[
R_p=\frac{1-\chi}{1+\chi}\to -1,
\]

\[
R_q=\frac{\chi-1}{\chi+1}\to +1.
\]

---

## 3. D/N validation

For a tube mode near the exit,

\[
f(w)
=
e^{ik(w-L)}
+
R e^{-ik(w-L)}.
\]

At \(w=L\),

\[
\frac{f_w}{f}
=
ik\frac{1-R}{1+R}.
\]

### Pressure/potential-like scalar

Using

\[
R_p=\frac{\epsilon-1}{\epsilon+1},
\]

the derivative ratio becomes

\[
\frac{p_w}{p}
=
\frac{ik}{\epsilon}.
\]

As \(\epsilon\to0\),

\[
p(L)=1+R_p\to0.
\]

So the open end is Dirichlet-like for \(p\):

\[
p(L)\simeq0.
\]

If the mouth is also Dirichlet, this gives a D/D ladder,

\[
k_j=\frac{(j+1)\pi}{L}.
\]

This **does not** reproduce the desired half-shifted D/N branch.

### Flow/displacement-like support variable

For the dual support variable,

\[
R_q=\frac{1-\epsilon}{1+\epsilon},
\]

and

\[
\frac{q_w}{q}
=
ik\epsilon.
\]

As \(\epsilon\to0\),

\[
\frac{q_w(L)}{q(L)}\to0,
\]

so

\[
q_w(L)\simeq0.
\]

With the mouth condition

\[
q(0)=0,
\]

the support ladder is

\[
q_j(w)=A\sin(k_j w),
\]

\[
q_j'(L)=0,
\]

which gives

\[
\cos(k_jL)=0,
\]

and therefore

\[
\boxed{
k_j=\frac{\pi}{L}\left(j+\frac12\right).
}
\]

So the D/N support ladder is validated by the open-junction mechanism only after the support coordinate is identified as the free-end flow/displacement variable, or an equivalent conjugate variable.

---

## 4. DC leakage and the zero-mode

For a steady throughput \(\Phi\), continuity in the 4D bulk gives

\[
A_{\rm bulk}(r)J_r(r)=\Phi.
\]

Since

\[
A_{\rm bulk}(r)=S_3r^3=2\pi^2r^3,
\]

the steady flux density is

\[
J_r(r)
=
\frac{\Phi}{2\pi^2r^3}.
\]

The divergence check is exact:

\[
\frac1{A_{\rm bulk}}\partial_r(A_{\rm bulk}J_r)=0.
\]

At a finite tube exit with effective 3D-ball cross-sectional area

\[
A_{\rm exit}=\frac43\pi a_{\rm exit}^3,
\]

the tube flux density is finite:

\[
J_{\rm tube}
=
\frac{\Phi}{A_{\rm exit}}
=
\frac{3\Phi}{4\pi a_{\rm exit}^3}.
\]

But for a hard cap,

\[
A_{\rm exit}\to0,
\]

and

\[
\frac{\Phi}{A_{\rm exit}}\to\infty.
\]

Therefore a capped geometry is incompatible with a finite nonzero DC superfluid throughput.

This is the main reason the old cap language must be retired.

---

## 5. Pass/fail table

| Gate | Status | Reason |
|---|---:|---|
| Hard cap with nonzero DC flux | **Fail** | \(\Phi/A_{\rm exit}\to\infty\) as \(A_{\rm exit}\to0\). |
| Open finite-radius exit for DC flow | **Pass** | \(J_{\rm tube}\) finite and \(J_r=\Phi/(2\pi^2r^3)\) satisfies radial continuity. |
| Strong AC reflection from sudden expansion | **Pass** | \(\mathcal T_E=4\epsilon/(1+\epsilon)^2=O(\epsilon)\) for \(\epsilon\ll1\). |
| Neumann reflection for generic scalar \(\psi\) | **Fail** | Pressure/potential-like scalar has \(R_p\to-1\), Dirichlet-like. |
| Neumann reflection for flow/displacement support variable | **Pass** | Dual variable has \(R_q\to+1\), \(q_w/q=ik\epsilon\to0\). |
| D/N ladder | **Conditional pass** | Valid if the finite-throat support field is the free-end flow/displacement coordinate or an impedance-transformed equivalent. |

---

## 6. Required wording patch

Replace the old cap language with:

> The finite throat is an open conduit. At \(w=L\), the exit radius is finite, \(R(L)>0\), and the throat opens into unconfined 4D spatial bulk. The half-shifted D/N support ladder is not produced by a hard geometric cap. It is produced by the low-impedance open-junction reflection of the free-end support coordinate. The conjugate pressure/potential variable sees the usual open-end Dirichlet-like phase, while the support displacement/flow coordinate sees a Neumann-like condition. The baseline/DC superfluid current exits through the same open junction and is tracked by the existing leakage/current bookkeeping.

This preserves both requirements:

\[
\text{AC support trapping}
\]

and

\[
\text{DC mass/current leakage}.
\]

---

## 7. Required geometrical/tapering condition if using a pressure-like scalar

If the program insists that the D/N field \(\psi\) is pressure-like or potential-like, then a sudden expansion does **not** generate the required Neumann condition.

To recover Neumann for that variable, one needs an effective high-impedance input:

\[
\frac{Z_L}{Z_t}\gg1.
\]

That requires one of the following:

1. a narrow high-impedance neck,
2. a quarter-wave impedance transformer,
3. a side-branch/stub resonator whose input impedance is high at the support frequencies,
4. or a variable redefinition in which the D/N support coordinate is the dual flow/displacement variable.

The simplest patch is option 4.

---

## 8. Volume 2 consequence

The rewritten V2-04 is not a cap-regularity problem anymore. It is an open-junction impedance theorem.

The next downstream stages should use:

\[
R(L)>0,
\]

not

\[
R(L)=0,
\]

and should distinguish:

- DC throughput and leakage, handled by continuity and \(S_{\rm leak}\),
- AC support modes, handled by impedance reflection,
- pressure/potential variables, which see Dirichlet-like open-end reflection,
- flow/displacement support variables, which see Neumann-like open-end reflection.

This is the clean way to preserve the D/N ladder without violating mass continuity.
