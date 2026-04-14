# Same-Charge Barrier Audit — Stage 005: Concrete Finite-Throat Primitive Branch, Pole Census, and the Residue/Linewidth Survival Test

## 0. Purpose

Stage 004 reduced the linear dynamic corridor to one exact question:

> can an admissible pole of the mixed one-port bundle produce a large enough residue-to-linewidth ratio to matter in the tunneling region before the same pole is dominated by absorptive loading?

The next honest move is therefore not another generic resonance slogan.
It is to pick one **explicit primitive branch family**, compute the first actual poles, and write the survival test in direct branch data.

That is what this stage does.

The main outputs are:

1. one concrete finite-throat primitive branch built from the lowest N/N wall/U profile and the lowest D/N half-wave support/W profile,
2. the exact quartic pole polynomial of the conservative one-port bundle,
3. the exact simple-pole cancellation
   \[
   \mathcal R_{Q,*}:=\frac{|A_{Q,*}|}{\gamma_*}
   =
   \frac{27c_s^5}{a^5\omega_*^5 N_*(\omega_*)},
   \]
   for the pure quadrupolar same-charge family,
4. the exact low-loss survival inequality against a required local barrier reduction,
5. one explicit numerical pole census,
6. and the first clean evidence of a **static/dynamic tension**: strengthening the static outgoing prefactor makes the linear dynamic residue/linewidth corridor worse on the same branch family.

So after Stage 005, the idea is still alive, but it has become much more specific.

---

## 1. Frozen input carried forward

The parent 4D/plasma ontology keeps the mixed channels
\[
A_w,\qquad F_{\mu w},\qquad J^w
\]
alive outside the strict far-field brane reduction, so the linear dynamic same-charge corridor still has to live in the mixed sector rather than in pure brane Maxwell shaping. The 5PN / moving-throat side independently says that the explicit isotropic finite-throat overlap model is already the correct place to package the grouped conservative/prefactor data, and that the surviving weak-axisymmetric grouped scalar is the outgoing-prefactor slope
\[
\Xi_1=\frac{P_1}{P_0}.
\]
So the correct next move is a concrete pole test on that same mixed/outgoing bundle. fileciteturn21file2turn24file4turn23file18

---

## 2. Concrete primitive finite-throat branch

Fix the finite throat interval
\[
s\in[0,L].
\]
Use the lowest N/N zero mode for the wall and the brane-like internal coordinate,
\[
u_0(s)=\frac1{\sqrt L},
\]
and the lowest D/N half-wave for the trapped support and mixed coordinate,
\[
f_0(s)=\sqrt{\frac2L}\sin\frac{\pi s}{2L}.
\]
Then the exact overlap constant is
\[
\kappa:=\int_0^L u_0(s)f_0(s)\,ds=\frac{2\sqrt2}{\pi}.
\]
On this branch the overlap-renormalized one-port couplings are
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R.
\]

The dynamic one-port wall/BdG/Maxwell/mixed bundle is therefore
\[
K_B(\omega)=K-M\omega^2-\frac{C^2}{\varpi^2-\omega^2},
\]
\[
A(\omega)=\Omega_U^2-\omega^2,
\qquad
W(\omega)=\Omega_W^2-\omega^2,
\]
\[
\Delta(\omega)=A(\omega)W(\omega)-R^2,
\]
\[
Q(\omega)=G_U^2W(\omega)+2G_UG_WR+G_W^2A(\omega),
\]
\[
D(\omega)=K_B(\omega)-\frac{Q(\omega)}{\Delta(\omega)}.
\]

The static admissibility conditions remain the same as earlier:
\[
\Delta_0>0,
\qquad
D_0>0,
\]
with
\[
\Delta_0=\Omega_U^2\Omega_W^2-R^2,
\qquad
D_0=K-\frac{C^2}{\varpi^2}-\frac{Q_0}{\Delta_0}.
\]

---

## 3. Exact quartic pole polynomial

Writing
\[
y=\omega^2,
\]
the conservative pole condition is exactly the quartic equation
\[
F(y)=0,
\]
with
\[
F(y)=
\Bigl((K-My)(\varpi^2-y)-C^2\Bigr)
\Bigl((\Omega_U^2-y)(\Omega_W^2-y)-R^2\Bigr)
-
(\varpi^2-y)
\Bigl(G_U^2(\Omega_W^2-y)+2G_UG_WR+G_W^2(\Omega_U^2-y)\Bigr).
\]

So the primitive finite-throat one-port branch has exactly four conservative poles in the generic admissible case.

This is already a useful narrowing:
we are no longer talking about a vague resonance landscape, but about one explicit quartic pole census.

---

## 4. Exact residue/linewidth cancellation

For a simple conservative pole \(\omega_*\) with \(\Delta(\omega_*)\neq 0\), attach the passive outgoing quadrupole port
\[
\Pi_{\rm out}(\omega)=i\Gamma_5\omega^5,
\qquad
\Gamma_5=\frac{a^5}{27c_s^5}.
\]

The exact transfer factor is
\[
N(\omega)=\frac{P(\omega)^2}{\Delta(\omega)^2},
\qquad
P(\omega)=A(\omega)G_W+RG_U.
\]

For the **wall susceptibility**
\[
\chi_{qq}(\omega)=\frac1{D(\omega)},
\]
the simple-pole residue and linewidth satisfy
\[
|A_{qq,*}|=\frac1{|D'(\omega_*)|},
\qquad
\gamma_*=rac{\Gamma_5\omega_*^5N(\omega_*)}{|D'(\omega_*)|},
\]
so the derivative cancels out exactly:
\[
\boxed{
\mathcal R_{qq,*}:=\frac{|A_{qq,*}|}{\gamma_*}
=
\frac{1}{\Gamma_5\omega_*^5N(\omega_*)}
=
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}.
}
\]

For the **pure quadrupolar same-charge family**
\[
S_Q(x)=\frac1{x^3},
\]
Stage 002/003 give
\[
\mathfrak V_Q(x,\omega)=-\frac12\chi_{qq}(\omega)\frac1{x^6},
\]
so the same exact ratio controls the useful conservative coefficient:
\[
\boxed{
\mathcal R_{Q,*}=\mathcal R_{qq,*}
=
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}.
}
\]

This is the first truly branch-level simplification of the dynamic same-charge problem.
The residue-to-linewidth figure is controlled only by

1. the pole frequency \(\omega_*\), and
2. the exact outgoing transfer factor \(N(\omega_*)\).

Everything else cancels out.

---

## 5. Exact low-loss survival inequality

Let \(\Delta V_{\rm req}(x)\) be the local barrier reduction required at radius \(x\) and impose the same low-loss condition as in Stage 004,
\[
|\Im\chi|\le \eta |\Re\chi|,
\qquad 0<\eta\le 1.
\]
Then the maximum conservative line shape allowed in that window is
\[
|\Re\chi|_{\max,\eta}
=
\frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}.
\]
For the pure quadrupolar family this gives the exact local survival criterion
\[
\frac12\,\mathcal R_{Q,*}\,\frac{\eta}{1+\eta^2}\,\frac1{x^6}
\ge
\Delta V_{\rm req}(x).
\]
Equivalently,
\[
\boxed{
\mathcal R_{Q,*}
\ge
2\,\Delta V_{\rm req}(x)
\frac{1+\eta^2}{\eta}
\,x^6.
}
\]
Substituting the exact pole ratio gives the branch-level theorem gate
\[
\boxed{
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}
\ge
2\,\Delta V_{\rm req}(x)
\frac{1+\eta^2}{\eta}
\,x^6.
}
\]
So the linear dynamic corridor lives or dies on one explicit inequality in the primitive branch data.

---

## 6. Explicit numerical primitive slice

Take the admissible sample branch
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,K,M)
=(0.5,\,0.3,\,0.4,\,0.25,\,1.0,\,1.4,\,2.0,\,3.0,\,1.0)
\]
with \(a=c_s=1\).

The overlap-renormalized couplings are
\[
C\approx 0.450158158078553,
\quad
G_U=0.3,
\quad
G_W\approx 0.360126526462843,
\quad
R\approx 0.225079079039277.
\]
The static branch data are
\[
\Delta_0\approx 1.90933940817883,
\qquad
D_0\approx 2.76355510933127,
\qquad
N_0\approx 0.0501661980249591,
\qquad
P_0\approx 0.0181527764203329.
\]
So the sample slice is statically admissible.

### 6.1 Uncoupled roots

The uncoupled wall/BdG roots are
\[
\omega_{\rm wall/BdG}^{(0)}\approx 1.68143182591478,
\qquad
2.04274007519334,
\]
and the uncoupled internal U/W roots are
\[
\omega_{UW}^{(0)}\approx 0.974601723746314,
\qquad
1.41779810977117.
\]

### 6.2 Full conservative pole census

The full quartic pole census is
\[
\omega_*\approx 0.938272741746753 \quad (\text{internal-like}),
\]
\[
\omega_*\approx 1.39141087653805 \quad (\text{internal-like}),
\]
\[
\omega_*\approx 1.72045371048003 \quad (\text{wall-like}),
\]
\[
\omega_*\approx 2.04539948783659 \quad (\text{wall-like}).
\]

The corresponding pure-Q residue/linewidth figures are
\[
\mathcal R_{Q,*}\approx 18.7069287828307,
\quad
0.380740659074003,
\quad
16.0250330226177,
\quad
32.0025481088465.
\]

So on this concrete slice the strongest candidate is the **upper wall-like pole**.

---

## 7. Comparison to the Stage-001 barrier benchmark

Use the Stage-001 illustrative reduced barrier benchmark at \(x=1\):
\[
V_{\rm known}(1)\approx 1.181909222592,
\qquad
\epsilon=0.1,
\qquad
\Delta V_{\rm req}(1)=V_{\rm known}(1)-\epsilon\approx 1.081909222592.
\]
Since \(S_Q(1)^2=1\), the exact low-loss thresholds are
\[
\mathcal R_{Q,*}^{\rm req}(\eta=0.1)
\approx 21.8545662963584,
\]
\[
\mathcal R_{Q,*}^{\rm req}(\eta=0.3)
\approx 7.86187368416853.
\]

### 7.1 Sample-slice verdict

For the lower wall-like pole,
\[
\mathcal R_{Q,*}\approx 16.0250,
\]
so it

- **fails** the \(10\%\)-loss benchmark,
- but **passes** the \(30\%\)-loss benchmark.

For the upper wall-like pole,
\[
\mathcal R_{Q,*}\approx 32.0025,
\]
so it passes both.

So the concrete primitive slice does **not** kill the dynamic corridor.
But it also does **not** make it generic.
Only some poles on the same branch are actually strong enough.

---

## 8. First static/dynamic tension

Now scan the outgoing-leg coupling \(\lambda_W\) while keeping all other sample parameters fixed.
The resulting static prefactor \(P_0\) and upper-wall residue/linewidth figure \(\mathcal R_{Q,*}\) are:

| \(\lambda_W\) | \(P_0\) | \(D_0\) | upper wall \(\omega_*\) | upper wall \(\mathcal R_{Q,*}\) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.005947405318 | 2.827234421584 | 2.044022723028 | 145.483858657863 |
| 0.4 | 0.018152776420 | 2.763555109331 | 2.045399487837 | 32.002548108846 |
| 0.6 | 0.038000163140 | 2.665913497210 | 2.047932775068 | 13.688535635681 |
| 0.8 | 0.067170782681 | 2.534309585220 | 2.051906688892 | 7.580971267466 |
| 1.0 | 0.108473308110 | 2.368743373361 | 2.057783390355 | 4.827389255647 |

So on this explicit family:

- stronger outgoing coupling raises the static prefactor \(P_0\),
- but the same move lowers the dynamic residue/linewidth figure \(\mathcal R_{Q,*}\).

That is the first clean evidence that the **static outgoing-normalization corridor** and the **linear dynamic low-loss corridor** are in real tension on the same branch family.

This is exactly the kind of tradeoff Stage 004 suggested but could not yet show concretely.

---

## 9. What Stage 005 changes

Stage 005 does not prove that the same-charge idea works.
But it does change the status materially.

Before this stage, the dynamic corridor lived or died on a generic slogan:
“find a pole with a large enough residue-to-linewidth ratio.”

After this stage, we have a much sharper statement:

1. the primitive finite-throat branch gives an exact quartic pole census,
2. the pure-Q same-charge family has the exact branch-level figure
   \[
   \mathcal R_{Q,*}=\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)},
   \]
3. the survival test is one explicit inequality in \(\omega_*\) and \(N(\omega_*)\),
4. and on a concrete admissible slice the corridor is **non-empty**, but only for part of the pole spectrum and only with a real static/dynamic tradeoff.

So the idea survives this stage, but narrowly.

---

## 10. Immediate next step

The next clean move is now even sharper than before:

1. keep the exact primitive finite-throat branch,
2. replace the sample numerical couplings by a PDE-grounded overlap extractor,
3. evaluate \(N(\omega_*)\) and the pole census on the actual moving-throat branch,
4. and compare the resulting \(\mathcal R_{Q,*}\) against the same low-loss survival threshold.

That is the first point where the surviving same-charge dynamic corridor stops being a generic mixed-sector hope and becomes a concrete branch test.
