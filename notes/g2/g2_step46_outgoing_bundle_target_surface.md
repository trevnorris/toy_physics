# Step 46 — The adiabatic constant-prefactor outgoing bundle lies on a one-parameter algebraic target surface

## Goal

Step 45 showed two things at once:

1. on the moving-throat side, the natural isotropic outgoing branch is especially simple when the prefactor is constant,
   \[
   P_2=P_4=0,
   \]
   so that
   \[
   K_0=P_0,\qquad
   K_2=\frac{P_0 a^2}{9c_s^2},\qquad
   K_4=\frac{4P_0 a^4}{81c_s^4},\qquad
   \Gamma_5=\frac{P_0 a^5}{27c_s^5},
   \]
   with the universal normalization condition
   \[
   P_0=\frac{54Gc_s^5}{5a^5c^5}
   \]
   on the natural source-map branch; and

2. on the adiabatic anomaly track, that same normalization law survives as a
   retuned target curve rather than being replaced by a new rule.

The next clean move is therefore to eliminate the intermediate variables
`(P_0,a,c_s)` and ask:

> what **purely algebraic surface** in the observable outgoing bundle
> \((K_0,K_2,K_4,\Gamma_5)\) is selected by the adiabatic constant-prefactor branch?

This step shows that the answer is extremely rigid.

---

## Step 46A — Exact bundle values on the universal target curve

Start from the constant-prefactor branch

\[
K_0=P_0,
\qquad
K_2=\frac{P_0 a^2}{9c_s^2},
\qquad
K_4=\frac{4P_0 a^4}{81c_s^4},
\qquad
\Gamma_5=\frac{P_0 a^5}{27c_s^5},
\]
together with the universal normalization target
\[
P_0=\frac{54Gc_s^5}{5a^5c^5}.
\]

Substituting immediately gives

\[
\boxed{
K_0=\frac{54Gc_s^5}{5a^5c^5},
}
\]

\[
\boxed{
K_2=\frac{6Gc_s^3}{5a^3c^5},
}
\]

\[
\boxed{
K_4=\frac{8Gc_s}{15ac^5},
}
\]

\[
\boxed{
\Gamma_5=\frac{2G}{5c^5}.
}
\]

So the whole outgoing bundle collapses to a single scale ratio
\[
s:=\frac{c_s}{a}.
\]

In terms of \(s\),

\[
\boxed{
K_0=\frac{54G}{5c^5}s^5,\qquad
K_2=\frac{6G}{5c^5}s^3,\qquad
K_4=\frac{8G}{15c^5}s,\qquad
\Gamma_5=\frac{2G}{5c^5}.
}
\]

That is already a strong result: the adiabatic constant-prefactor bundle is not a
generic 4-parameter family. It is a **one-parameter leaf** with a universal odd slot.

---

## Step 46B — Exact algebraic surface in observable bundle space

Because the bundle depends only on the single ratio \(s=c_s/a\), the observable
coefficients satisfy exact algebraic relations.

The simplest is

\[
\boxed{
K_2^2=\frac14\,K_0K_4.
}
\]

Equivalently,

\[
\boxed{
K_0=\frac{4K_2^2}{K_4}.
}
\]

Using the universal odd coefficient \(\Gamma_5=2G/(5c^5)\), the same one-parameter
leaf can also be written as

\[
\boxed{
K_2=\frac{81}{64}\frac{K_4^3}{\Gamma_5^2},
}
\]

\[
\boxed{
K_0=\frac{6561}{1024}\frac{K_4^5}{\Gamma_5^4}.
}
\]

So the adiabatic constant-prefactor branch is a **codimension-2 algebraic target
surface** inside the 4-dimensional coefficient space \((K_0,K_2,K_4,\Gamma_5)\):

- one condition fixes the odd slot universally,
- one condition ties the three even slots together.

That is exactly the kind of theorem gate a future PDE output can be tested against
without ever referring back to \((P_0,a,c_s)\) explicitly.

---

## Step 46C — Physical reading

This outgoing bundle has a very clean interpretation.

### The odd coefficient is universal

\[
\boxed{\Gamma_5=\frac{2G}{5c^5}}
\]

is independent of the adiabatic retuning parameter.
So the universal Burke–Thorne quadrupole coefficient is **not** moved by the
adiabatic anomaly once the branch is expressed on the full target curve.

### The even coefficients form a rigid hierarchy

The three even slots carry the single ratio \(s=c_s/a\):

- \(K_4\) is linear in \(s\),
- \(K_2\) is cubic in \(s\),
- \(K_0\) is quintic in \(s\).

So the anomaly can move the bundle only along this one-parameter leaf; it cannot
deform the even hierarchy arbitrarily if the constant-prefactor branch is correct.

---

## Step 46D — What is established, and what is still conditional

### Established here

Once we take both of the already isolated ingredients seriously,

1. the moving-throat constant-prefactor branch
   \[
   P_2=P_4=0,
   \]
2. the universal normalization target
   \[
   P_0=\frac{54Gc_s^5}{5a^5c^5},
   \]

the outgoing bundle is forced onto the exact surface

\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4.
\]

### Still conditional

What is **not** yet proven by the current file stack is that the true moving-throat
PDE actually selects the constant-prefactor branch rather than a more general
\((P_0,P_2,P_4)\) branch. So the exact surface derived here should be read as:

> the sharp algebraic target selected by the minimal isotropic outgoing branch.

That is still the right next target, because it is now simple enough that a future
PDE computation can fail it cleanly.

---

## Main result of the step

The adiabatic constant-prefactor outgoing bundle is a one-parameter algebraic leaf:

\[
\boxed{
K_0=\frac{54G}{5c^5}\left(\frac{c_s}{a}\right)^5,\qquad
K_2=\frac{6G}{5c^5}\left(\frac{c_s}{a}\right)^3,\qquad
K_4=\frac{8G}{15c^5}\left(\frac{c_s}{a}\right),\qquad
\Gamma_5=\frac{2G}{5c^5}.
}
\]

Equivalently, in purely observable form,

\[
\boxed{
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4,
}
\]

or

\[
\boxed{
K_0=\frac{4K_2^2}{K_4},
\qquad
K_2=\frac{81}{64}\frac{K_4^3}{\Gamma_5^2}.
}
\]

So the next PDE-facing theorem gate is no longer diffuse at all:

> compute the outgoing bundle and check whether it lands on this algebraic target surface.
