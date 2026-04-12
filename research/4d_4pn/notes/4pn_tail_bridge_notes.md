# 4PN tail / hereditary bridge audit — result ledger

## What this step accomplished

This step freezes the **nonlocal conservative 4PN tail bridge** after the local 4PN theorem.

The main result is that the tail side is much narrower than the local side was:

1. the exact GR conservative 4PN tail is one nonlocal STF-quadrupole functional,
2. its local logarithmic shadow occupies only the **degree-2 `G^4/r^4` block**,
3. the source representation is the same canonical STF quadrupole already isolated by the 2.5PN program,
4. and the remaining toy-model gap can be parameterized by **one scalar transport coefficient**.

So after the local 4PN theorem, the hereditary problem is **not** another large generic-frame lift. It is a sharply focused quadrupole-transport theorem.

---

## Exact GR tail functional now frozen

From the standard 4PN conservative literature, the nonlocal time-symmetric tail Hamiltonian has coefficient

\[
\alpha_{\rm tail}^{\rm GR}=\frac{G^2 M}{5c^8}.
\]

Equivalently,

\[
H_{4\rm PN}^{\rm tail,sym}(t)
=
-\frac{1}{5}\frac{G^2M}{c^8}
I_{ij}^{(3)}(t)
\,\mathrm{Pf}_{2s/c}\!\int_{-\infty}^{+\infty}\frac{dv}{|v|}
I_{ij}^{(3)}(t+v).
\]

The local logarithmic shadow is controlled by

\[
F_{\rm tail}(t)=\frac{2}{5}\frac{G^2M}{c^8}\bigl(I_{ij}^{(3)}(t)\bigr)^2.
\]

---

## Universal Newtonian quadrupole shadow

Using the Newtonian order-reduced STF quadrupole

\[
I_{ij}=\mu\,\mathrm{STF}(x_i x_j),
\]

the audit verified the exact order-reduced third derivative

\[
I_{ij}^{(3)}
=
-\frac{2GM\mu}{r^3}
\left(4x_{\langle i}v_{j\rangle}-3\frac{\dot r}{r}x_{\langle i}x_{j\rangle}\right),
\]

and therefore

\[
\bigl(I_{ij}^{(3)}\bigr)^2
=
\frac{8}{3}\frac{G^2M^2\mu^2}{r^4}
\left(12v^2-11\dot r^2\right).
\]

This gives the exact local logarithmic tail shadow

\[
\frac{F_{\rm tail}}{\mu}
=
\frac{16}{15}\,\nu\,\frac{U^4}{c^8}
\left(12v^2-11\dot r^2\right),
\qquad
U\equiv\frac{GM}{r}.
\]

So the hereditary shadow sits **only** in the local 4PN `U`-block (`G^4/r^4` with degree-2 kinematics).

That is the exact reason the tail problem is structurally much smaller than the local lift problem.

---

## Exact bridge to the 2.5PN quadrupole coefficient

The 2.5PN audit already froze the canonical Burke–Thorne coefficient target as

\[
\gamma_{\rm GR}=\frac{2G}{5c^5}.
\]

The present tail audit observes the exact arithmetic relation

\[
\alpha_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\,\gamma_{\rm GR}.
\]

So in GR the 4PN conservative tail coefficient is exactly **half a monopole-scattering factor** `(GM/c^3)` times the 2.5PN quadrupole coefficient.

This is the cleanest bridge between the already narrowed 2.5PN quadrupole route and the 4PN hereditary sector.

---

## Minimal toy-model tail bridge ansatz

Let

\[
\gamma_{\rm toy}
\equiv
\hat m_0^2\Gamma_5
\]

be the canonically normalized STF quadrupole coefficient already isolated by the 2.5PN notes.

Then the remaining hereditary gap can be parameterized by one scalar transport factor `Theta_tail`:

\[
\boxed{
\alpha_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}
\gamma_{\rm toy}.
}
\]

On the `c_s=c` branch,

- if `gamma_toy = gamma_GR`, and
- if `Theta_tail = 1`,

then the exact GR 4PN conservative tail coefficient is recovered automatically.

So the tail bridge has become maximally narrow:

- the **representation/source problem is already solved** by the 2.5PN quadrupole audit,
- the **local 4PN basis problem is already solved** by the local referee chain,
- and the remaining issue is only whether the moving-throat model produces the correct **monopole-scattering transport coefficient**.

---

## What remains open after this step

This step does **not** derive the full hereditary kernel from the moving-throat PDE.

What remains open is the single scalar theorem:

\[
\Theta_{\rm tail}\stackrel{?}{=}1
\]

on the same canonical STF quadrupole branch already isolated by the 2.5PN program.

Equivalently:

- no new tensor channel is needed,
- no new generic-frame existence solve is needed,
- the remaining 4PN hereditary problem is a **quadrupole tail-transport normalization theorem**.

That is the clean next target before any claim of a full conservative 4PN theorem.
