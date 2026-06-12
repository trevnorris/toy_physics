# Moving-Throat PDE — Stage 222: Concrete Finite-Throat Primitive Branch, Pole Census, and the Residue/Linewidth Survival Test

## Status

**Exact within the carried isotropic finite-throat primitive one-port wall/BdG/Maxwell/mixed closure** once the lowest N/N wall/U profile and the lowest D/N half-wave support/W profile are fixed.

The quartic pole polynomial, the wall-like residue/linewidth cancellation, and the low-loss survival gate are exact inside that reduced closure. The explicit pole census and the outgoing-leg scan below are illustrative evaluations on one admissible sample slice of the same exact formulas.

This stage is the first post-Stage-221 insertion of an actual primitive finite-throat branch into the local mixed-port resonance ledger. It does **not** solve the full driven two-throat PDE. It computes the first concrete pole census and turns the linear same-charge dynamic question into one explicit branch inequality.

---

## Purpose

Stage 221 reduced the linear dynamic corridor to one exact question:

> can an admissible pole of the mixed one-port bundle produce a large enough residue-to-linewidth ratio to matter in the tunneling region before the same pole is dominated by absorptive loading?

The next honest move is therefore not another generic resonance slogan. It is to choose one **explicit primitive branch family**, compute the first actual poles, and write the survival test in direct branch data.

The main outputs are:

1. one concrete finite-throat primitive branch built from the lowest N/N wall/U profile and the lowest D/N half-wave support/W profile,
2. the exact quartic pole polynomial of the conservative one-port bundle,
3. the exact simple-pole cancellation
   \[
   \mathcal R_{Q,*}:=\frac{|A_{Q,*}|}{\gamma_*}
   =
   \frac{27c_s^5}{a^5\omega_*^5N(\omega_*)},
   \]
   for the pure quadrupolar same-charge family,
4. the exact low-loss survival inequality against a required local barrier reduction,
5. one explicit numerical pole census,
6. and the first clean evidence of a **static/dynamic tension**: strengthening the static outgoing prefactor makes the linear dynamic residue/linewidth corridor worse on the same branch family.

So Stage 222 keeps the dynamic same-charge idea alive, but only in a much narrower and more branch-specific form.

---

## 1. Frozen input carried forward

The parent `4+1` / plasma ontology keeps the mixed channels
\[
A_w,\qquad F_{\mu w},\qquad J^w
\]
alive beyond the strict far-field brane reduction, so the linear same-charge dynamic corridor still has to live in the mixed sector rather than in pure brane Maxwell shaping.

Stage 221 already showed that the only linear corridor left is a **resonant dispersive** one. Near a wall-like pole of the dynamic one-port bundle, the local susceptibility takes the normal form
\[
\chi_{qq}(\omega)
\approx
\frac{A_{qq,*}}{\delta-i\gamma_*},
\qquad
\delta=\omega-\omega_*,
\]
with the low-loss bound
\[
\sup_{|\Im\chi|\le \eta|\Re\chi|}|\Re\chi|
=
\frac{|A_{qq,*}|}{\gamma_*}\frac{\eta}{1+\eta^2},
\qquad 0<\eta\le 1.
\]
So the next question is no longer whether a new spatial kernel class exists. It is whether an explicit admissible pole can make the already-known quadrupolar family large enough **before** passive/outgoing widening destroys the gain.

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
f_0(s)=\sqrt{\frac2L}\,\sin\frac{\pi s}{2L}.
\]
Then the exact overlap constant is
\[
\kappa:=\int_0^L u_0(s)f_0(s)\,ds=\frac{2\sqrt2}{\pi}.
\]

On this primitive branch the overlap-renormalized one-port couplings are
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R.
\]

The conservative dynamic one-port wall/BdG/Maxwell/mixed bundle is therefore
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

The static admissibility conditions are the same as before:
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

Write
\[
y=\omega^2.
\]
Then the conservative pole condition is exactly the quartic equation
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

Equivalently,
\[
D(\omega)=\frac{F(\omega^2)}{(\varpi^2-\omega^2)\,\Delta(\omega)}.
\]
So the primitive finite-throat one-port branch has exactly four conservative poles in the generic admissible case.

This is already a major narrowing. We are no longer talking about a vague resonance landscape, but about one explicit quartic pole census.

---

## 4. Exact residue/linewidth cancellation

For a simple conservative wall-like pole \(\omega_*\) with
\[
D_0(\omega_*)=0,
\qquad
D_0'(\omega_*)\neq 0,
\qquad
\Delta(\omega_*)\neq 0,
\]
attach the passive outgoing quadrupole port
\[
\Pi_{\rm out}(\omega)=i\Gamma_5\omega^5,
\qquad
\Gamma_5=\frac{a^5}{27c_s^5}.
\]

The exact outgoing transfer factor on the primitive branch is
\[
N(\omega)=\frac{P(\omega)^2}{\Delta(\omega)^2},
\qquad
P(\omega)=A(\omega)G_W+RG_U.
\]

For the wall susceptibility
\[
\chi_{qq}(\omega)=\frac1{D_\Pi(\omega)},
\]
Stage 221 gives the local normal form
\[
\chi_{qq}(\omega)
\approx
\frac{1/D_0'(\omega_*)}{\delta-i\gamma_*},
\qquad
\gamma_*=
\frac{\Gamma_5\omega_*^5N(\omega_*)}{|D_0'(\omega_*)|}.
\]
Therefore the simple-pole residue and linewidth satisfy
\[
|A_{qq,*}|=\frac1{|D_0'(\omega_*)|},
\qquad
\gamma_*=
\frac{\Gamma_5\omega_*^5N(\omega_*)}{|D_0'(\omega_*)|},
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
the reduced barrier contribution is
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

This is the first truly branch-level simplification of the dynamic same-charge problem. The residue-to-linewidth figure is controlled only by

1. the pole frequency \(\omega_*\), and
2. the exact outgoing transfer factor \(N(\omega_*)\).

Everything else cancels out.

---

## 5. Exact low-loss survival inequality

Let \(\Delta V_{\rm req}(x)\) be the local barrier reduction required at radius \(x\), and impose the same low-loss condition as in Stage 221,
\[
|\Im\chi|\le \eta|\Re\chi|,
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
=
(0.5,\,0.3,\,0.4,\,0.25,\,1.0,\,1.4,\,2.0,\,3.0,\,1.0)
\]
with \(a=c_s=1\).

The overlap-renormalized couplings are
\[
C\approx 0.450158158078553,
\qquad
G_U=0.3,
\qquad
G_W\approx 0.360126526462842,
\qquad
R\approx 0.225079079039277.
\]
The static branch data are
\[
\Delta_0\approx 1.90933940817883,
\qquad
D_0\approx 2.76355510933127,
\]
\[
N_0\approx 0.0501661980249591,
\qquad
P_0\approx 0.0181527764203328.
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

| pole \(\omega_*\) | type | \(\mathcal R_{Q,*}\) |
|---:|:---|---:|
| 0.938272741746754 | internal-like | 18.7069287828307 |
| 1.39141087653804 | internal-like | 0.380740659074003 |
| 1.72045371048003 | wall-like | 16.0250330226177 |
| 2.04539948783659 | wall-like | 32.0025481088465 |

So on this concrete slice the strongest candidate is the **upper wall-like pole**.

---

## 7. Comparison to the illustrative barrier benchmark

Carry forward the same illustrative reduced barrier benchmark at \(x=1\):
\[
V_{\rm known}(1)\approx 1.181909222592,
\qquad
\epsilon=0.1,
\qquad
\Delta V_{\rm req}(1)=V_{\rm known}(1)-\epsilon\approx 1.081909222592.
\]

Provenance note. The \(V_{\rm known}(1)\) float is an illustrative sample-slice barrier benchmark used to compute
\(\Delta V_{\rm req}\); its genealogy is script/session-side, not an in-note derivation or external published target.

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

So the concrete primitive slice does **not** kill the dynamic corridor. But it also does **not** make it generic. Only part of the pole spectrum on the same branch is actually strong enough.

---

## 8. First static/dynamic tension

Now scan the outgoing-leg coupling \(\lambda_W\) while keeping all other sample parameters fixed. The resulting static prefactor \(P_0\) and upper-wall residue/linewidth figure \(\mathcal R_{Q,*}\) are:

| \(\lambda_W\) | \(P_0\) | \(D_0\) | upper wall \(\omega_*\) | upper wall \(\mathcal R_{Q,*}\) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.00594740531769 | 2.82723442158450 | 2.04402272302752 | 145.483858657863 |
| 0.4 | 0.01815277642033 | 2.76355510933127 | 2.04539948783659 | 32.0025481088465 |
| 0.6 | 0.03800016314041 | 2.66591349720966 | 2.04793277506821 | 13.6885356356808 |
| 0.8 | 0.06717078268091 | 2.53430958521967 | 2.05190668889211 | 7.58097126746582 |
| 1.0 | 0.10847330811048 | 2.36874337336129 | 2.05778339035510 | 4.82738925564702 |

So on this explicit family:

- stronger outgoing coupling raises the static prefactor \(P_0\),
- but the same move lowers the dynamic residue-to-linewidth figure \(\mathcal R_{Q,*}\).

That is the first clean evidence that the **static outgoing-normalization corridor** and the **linear dynamic low-loss corridor** are in real tension on the same branch family.

---

## 9. What Stage 222 achieves physically

This stage does **not** prove that the same-charge dynamic idea works. But it changes the status materially.

Before this stage, the dynamic corridor lived or died on a generic slogan: “find a pole with a large enough residue-to-linewidth ratio.”

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

## 10. What is still missing

Several things remain outside the scope of this primitive-branch theorem.

1. **Actual overlap extraction.** The sample couplings are illustrative. The real moving-throat PDE must determine the branch-level overlap data \(\lambda_B,\lambda_U,\lambda_W,\lambda_R\) on the actual branch.
2. **Compatibility with the exact isotropic target surface.** This stage does not yet force the primitive branch onto the 5PN / 2.5PN / 4PN isotropic one-pole and outgoing-normalization surface.
3. **Actual barrier data.** The benchmark \(\Delta V_{\rm req}(1)\) used here is illustrative. The true local barrier requirement must be pulled back from the actual same-charge branch.

So the right continuation is not another generic resonance story. It is to ask whether the same primitive branch can remain useful **after** it is forced onto the exact isotropic target surface.

---

## 11. SymPy-backed status

The accompanying SymPy audit verifies all of the concrete algebra used here:

- the exact primitive overlap constant
  \(\kappa=2\sqrt2/\pi\),
- the exact quartic pole polynomial of the conservative primitive bundle,
- the wall-like residue/linewidth cancellation
  \(\mathcal R_{Q,*}=27c_s^5/(a^5\omega_*^5N(\omega_*))\),
- the exact low-loss survival threshold,
- the uncoupled and coupled pole census on the sample slice,
- the residue/linewidth figures quoted in the table,
- and the monotone static/dynamic tension under the \(\lambda_W\) scan.

Supporting file:

- `moving_throat_pde_stage222_concrete_finite_throat_primitive_branch_pole_census_and_residue_linewidth_survival_test_sympy_audit.py`

---

## 12. Immediate next step

The next clean move is now even sharper than before.

1. Keep the exact primitive finite-throat branch.
2. Force it onto the exact isotropic one-pole and outgoing-normalization surface.
3. Recompute the pole census and the survival inequality on that compatible branch.
4. Then ask whether the resulting compatible branch still has a non-empty dynamic survival window.

That is exactly the Stage-223 theorem gate.
