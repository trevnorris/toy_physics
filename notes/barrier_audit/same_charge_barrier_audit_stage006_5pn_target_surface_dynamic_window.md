# Same-Charge Barrier Audit — Stage 006: 5PN Isotropic Target Surface, Primitive-Branch Compatibility, and the Dynamic Survival Window

## 0. Purpose

Stage 005 showed that the linear same-charge corridor does **not** immediately die on the explicit primitive finite-throat branch. But that result was still too loose, because the primitive slice had not yet been forced onto the exact isotropic target surface already isolated by the 5PN / 2.5PN / 4PN moving-throat endgame.

So the next honest question is:

> if the same primitive one-port branch is required to satisfy the exact isotropic one-pole and outgoing-normalization conditions, does the dynamic same-charge corridor survive or die?

That is what this stage answers.

The main outputs are:

1. the exact symbolic compatibility equation between the isotropic one-pole condition and the isotropic outgoing-normalization condition,
2. its specialization to the explicit primitive finite-throat one-port family,
3. one concrete compatibility branch on the Stage-005 sample slice,
4. the resulting compatibility-branch pole census,
5. and a finite **dynamic survival window** in the branch-compatible target parameter.

So after Stage 006, the problem is no longer merely “find a good pole.” It is:

> can the same branch support a good pole **while also lying on the 5PN-compatible isotropic target surface**?

---

## 1. Frozen input carried forward

### 1.1 Primitive finite-throat one-port branch from Stage 005

Keep the same explicit finite-throat branch:

- lowest N/N zero mode for the wall and brane-like internal coordinate,
- lowest D/N half-wave for the trapped support and mixed coordinate,
- overlap constant
  \[
  \kappa = \frac{2\sqrt2}{\pi}.
  \]

With reduced couplings
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R,
\]
we still have
\[
\Delta = \Omega_U^2\Omega_W^2-R^2,
\qquad
Q = G_U^2\Omega_W^2 + 2G_UG_WR + G_W^2\Omega_U^2,
\qquad
P = \Omega_U^2G_W + RG_U.
\]

The primitive bundle moments are
\[
B_0=\frac{C^2}{\varpi^2},
\qquad
B_2=\frac{C^2}{\varpi^4},
\qquad
B_4=\frac{C^2}{\varpi^6},
\]
\[
Z_0=\frac{Q}{\Delta},
\qquad
Z_2=\frac{QS_2-H\Delta}{\Delta^2},
\qquad
Z_4=\frac{Q(S_2^2-\Delta)-S_2H\Delta}{\Delta^3},
\]
where
\[
S_2=\Omega_U^2+\Omega_W^2,
\qquad
H=G_U^2+G_W^2,
\]
and
\[
N_0=\frac{P^2}{\Delta^2}.
\]

### 1.2 Exact isotropic 5PN / 2.5PN target surface

On the isotropic one-port bundle,
\[
D_0 = K-B_0-Z_0,
\qquad
D_2 = -(M+B_2+Z_2),
\qquad
D_4 = -(B_4+Z_4),
\]
with normalized conservative response
\[
u_2 = -\frac{D_2}{D_0},
\qquad
u_4 = \frac{D_2^2-D_0D_4}{D_0^2},
\]
and outgoing prefactor
\[
P_0 = \frac{N_0}{D_0}.
\]

The exact isotropic one-pole condition is
\[
\nu_4 = 4\nu_2^2
\iff
D_0(B_4+Z_4) = 3(M+B_2+Z_2)^2.
\]
The isotropic outgoing-normalization condition is
\[
P_0 = P_{0,\mathrm{target}},
\]
where for the fully calibrated moving-throat branch
\[
P_{0,\mathrm{target}} = \frac{54Gc_s^5}{5a^5c^5\,\hat m_0^{\,2}}.
\]

Stage 006 keeps this target symbolic as \(P_{0,\mathrm{target}}\), because on the primitive reduced branch the important question is compatibility first.

---

## 2. Exact compatibility equation

The one-pole condition solves for the wall stiffness as
\[
K_{\mathrm{pole}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4} + B_0 + Z_0.
\]
The outgoing-normalization condition solves for the same wall stiffness as
\[
K_{\mathrm{norm}}
=
\frac{N_0}{P_{0,\mathrm{target}}} + B_0 + Z_0.
\]
So simultaneous isotropic one-pole success and isotropic normalization success require
\[
K_{\mathrm{pole}} = K_{\mathrm{norm}},
\]
which is equivalent to the exact compatibility equation
\[
\boxed{
\frac{N_0}{P_{0,\mathrm{target}}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]
Equivalently, the primitive branch itself induces the unique branch-compatible target
\[
\boxed{
P_{0,\mathrm{target,compat}} = \frac{N_0(B_4+Z_4)}{3(M+B_2+Z_2)^2}.
}
\]
This is the first exact place where the same-charge branch is forced to talk directly to the 5PN isotropic surface.

Two points are worth emphasizing.

First, this is **not** a generic fit condition. It is an exact consequence of trying to satisfy the same isotropic target surface from two sides.

Second, the compatibility equation does **not** determine every coupling. It tells us whether the primitive branch wants a normalization target that is even compatible with its own conservative one-pole structure.

---

## 3. Primitive specialization of the compatibility equation

Substituting the explicit primitive one-port data gives
\[
P_{0,\mathrm{target,compat}}
=
\frac{\dfrac{P^2}{\Delta^2}\left(\dfrac{C^2}{\varpi^6}+Z_4\right)}{3\left(M+\dfrac{C^2}{\varpi^4}+Z_2\right)^2},
\]
or equivalently
\[
\boxed{
\frac{P^2/\Delta^2}{P_{0,\mathrm{target}}}
=
\frac{3\left(M+\dfrac{C^2}{\varpi^4}+Z_2\right)^2}{\dfrac{C^2}{\varpi^6}+Z_4}.
}
\]
So on the primitive family the isotropic 5PN-compatible surface is already a single explicit algebraic relation in the radial/axial couplings and frequencies.

That is a substantial tightening compared with Stage 005.

---

## 4. Concrete sample compatibility branch

Now specialize to the same Stage-005 sample values
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
(0.5,0.3,0.4,0.25,1.0,1.4,2.0,1.0),
\]
with \(a=c_s=1\).

The overlap-renormalized primitive data are
\[
C \approx 0.450158158078553,
\qquad
G_U = 0.3,
\qquad
G_W \approx 0.360126526462843,
\qquad
R \approx 0.225079079039277.
\]
The static bundle quantities are
\[
\Delta \approx 1.90933940817883,
\qquad
Q \approx 0.354725283210515,
\qquad
P \approx 0.427650250174625,
\]
\[
B_0 \approx 0.0506605918211689,
\quad
B_2 \approx 0.0126651479552922,
\quad
B_4 \approx 0.00316628698882306,
\]
\[
Z_0 \approx 0.185784298847558,
\quad
Z_2 \approx 0.172955320626603,
\quad
Z_4 \approx 0.170825285860668,
\]
\[
N_0 \approx 0.0501661980249591.
\]

The exact compatibility target on this primitive slice is
\[
\boxed{
P_{0,\mathrm{target,compat}}
\approx 0.00206979231806289.
}
\]
The corresponding compatibility wall stiffness is
\[
\boxed{
K_{\mathrm{compat}}
\approx 24.4737548792910.
}
\]
So the compatibility branch is much stiffer than the loose Stage-005 sample branch. Its compatible static denominator is
\[
D_{0,\mathrm{compat}} = K_{\mathrm{compat}}-B_0-Z_0
\approx 24.2373099886222.
\]

This is already informative.

The Stage-005 sample branch had
\[
P_0 \approx 0.0181527764203329,
\]
while the same primitive family, when forced onto the exact isotropic one-pole/normalization compatibility surface, wants
\[
P_{0,\mathrm{target,compat}} \approx 0.00207.
\]
So the 5PN-compatible branch lives at a much lower static prefactor and a much higher wall stiffness on this particular primitive slice.

---

## 5. Pole census on the compatibility branch

Using \(K=K_{\mathrm{compat}}\), the conservative pole census is
\[
\omega_* \approx 0.971575315129468 \quad (\text{internal-like}),
\]
\[
\omega_* \approx 1.41651290122561 \quad (\text{internal-like}),
\]
\[
\omega_* \approx 1.99753567893361 \quad (\text{wall-like}),
\]
\[
\omega_* \approx 4.94905432364313 \quad (\text{wall-like}).
\]

The exact pure-quadrupolar residue/linewidth figure remains
\[
\mathcal R_{Q,*} = \frac{27c_s^5}{a^5\omega_*^5 N(\omega_*)}.
\]
On the compatibility branch the four values are
\[
\mathcal R_{Q,*} \approx 0.159888393135835 \quad (\text{internal-like}),
\]
\[
\mathcal R_{Q,*} \approx 0.000806281535937178 \quad (\text{internal-like}),
\]
\[
\mathcal R_{Q,*} \approx 30.1999075602499 \quad (\text{wall-like}),
\]
\[
\mathcal R_{Q,*} \approx 36.1711864832695 \quad (\text{wall-like}).
\]
So on this concrete 5PN-compatible branch the dynamic same-charge corridor is **not** carried by the internal poles at all. It is carried entirely by the wall-like poles.

This is already a nontrivial structural simplification.

---

## 6. Dynamic survival window on the compatibility surface

Carry forward the same illustrative local barrier benchmark from Stage 005 at \(x=1\):
\[
V_{\mathrm{known}}(1) \approx 1.181909222592,
\qquad
\epsilon = 0.1,
\qquad
\Delta V_{\mathrm{req}}(1) \approx 1.081909222592.
\]
Then the required residue/linewidth thresholds are
\[
\mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.1)
\approx 21.8545662963584,
\]
\[
\mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.3)
\approx 7.86187368416853.
\]

Now scan \(\lambda_W\) **along the exact compatibility surface**, i.e. always resetting \(K\) to \(K_{\mathrm{compat}}(\lambda_W)\). The resulting branch-compatible target and wall-like residue/linewidth figures are:

| \(\lambda_W\) | \(P_{0,\mathrm{target,compat}}\) | \(K_{\mathrm{compat}}\) | lower wall \(\mathcal R_Q\) | upper wall \(\mathcal R_Q\) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.000576970879843 | 29.3158464872314 | 138.814136942081 | 137.502546600713 |
| 0.4 | 0.002069792318063 | 24.4737548792910 | 30.1999075602499 | 36.1711864832695 |
| 0.6 | 0.004865681200486 | 21.1544287401845 | 12.8348600273988 | 16.7575510327116 |
| 0.8 | 0.009169913681573 | 19.0298300900561 | 7.06074242207991 | 9.69035785242054 |
| 1.0 | 0.014981190324091 | 17.7824591822917 | 4.45922850098679 | 6.30111094469551 |

This scan shows a clean monotonic tradeoff on the explicit compatibility family:

- increasing \(\lambda_W\) raises the branch-compatible static target \(P_{0,\mathrm{target,compat}}\),
- the same move lowers the required wall stiffness \(K_{\mathrm{compat}}\),
- and both wall-like dynamic figures \(\mathcal R_Q\) fall monotonically.

So the static/dynamic tension from Stage 005 survives **even after** the branch is forced onto the exact 5PN isotropic target surface.

### 6.1 Finite survival windows

At the stricter \(10\%\)-loss benchmark,

- the **lower wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00283133168555932,
  \]
- the **upper wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00359651058968466.
  \]

At the looser \(30\%\)-loss benchmark,

- the **lower wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00817339430971383,
  \]
- the **upper wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.0116633929790174.
  \]

So the dynamic corridor is not generically open on the explicit primitive family. It survives only inside a finite interval of the same branch-compatible target that the isotropic 5PN surface itself wants.

That is the sharpest same-charge compatibility statement reached so far.

---

## 7. What Stage 006 changes

Stage 006 does **not** prove the same-charge idea works.
But it does materially tighten the status.

Before this stage, the dynamic corridor lived on a primitive branch that had not yet been forced to satisfy the exact isotropic 5PN surface.

After this stage, we know four much stronger things.

### 7.1 The primitive branch can be put on the exact isotropic 5PN surface

The isotropic one-pole and normalization conditions do not conflict abstractly. They reduce to one exact compatibility equation.

### 7.2 The dynamic corridor is not killed automatically by that calibration

On the concrete sample slice, once the branch is moved to the compatibility wall stiffness, both wall-like poles still clear the stricter \(10\%\)-loss benchmark.

### 7.3 The same branch develops a genuine target window

The dynamic corridor survives only for a finite range of the branch-compatible static target. So the eventual PDE-selected normalization cannot be arbitrarily large on this primitive family if the dynamic same-charge corridor is to remain alive.

### 7.4 The wall-like poles are now the only relevant survivors

The internal poles are dynamically irrelevant on the compatibility branch. So the surviving same-charge corridor is a wall-like corridor, not a generic mixed-pole corridor.

---

## 8. Best current verdict after Stage 006

The idea is still alive.
But it is now alive in a much narrower form:

> an explicit primitive finite-throat branch can satisfy the exact isotropic 5PN compatibility surface **and** still retain a wall-like dynamic same-charge corridor, but only inside a finite branch-compatible normalization window.

So the next honest move is no longer another generic resonance scan.
It is:

1. extract the actual branch-compatible normalization target from the moving-throat PDE,
2. compare it against the finite survival window above,
3. and see whether the real branch lands inside or outside that window.

That is the next clean kill test.
