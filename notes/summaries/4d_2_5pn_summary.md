4d_2_5pn.tex — Comprehensive summary (2.5PN theorem-audit / quadrupole-normalization package)

## 0) What this package is doing

This document is the **carry-forward 2.5PN theorem summary** for the 4D toy-model program.

It is built from the session notes in `25pn_notes.md` plus the two SymPy referee/audit scripts:

- `2_5pn_master_sympy_audit.py`
- `2_5pn_master_session_sympy_audit.py`

The package is not a finished moving-throat PDE paper. Its role is narrower and more useful:

1. freeze what the program looked like before the 2.5PN push,
2. record the decisive benchmark showing that the frozen conservative hierarchy cannot by itself produce a local 2.5PN theorem,
3. identify the **minimal retarded structure** needed for a GR-like 2.5PN sector,
4. audit the scalar, dipole, and quadrupole channels one by one,
5. reduce the surviving universal point-particle route to the **orbital/worldtube STF quadrupole**,
6. and isolate the one remaining serious theorem gap as the **final passive/outgoing quadrupole normalization** on the actual moving-throat branch.

The strongest honest verdict carried forward here is:

\[
\boxed{
\text{A full GR-like point-particle 2.5PN theorem is conditionally reachable from the present stack.}
}
\]

The remaining gap is no longer diffuse. It has collapsed to:

\[
\boxed{
\text{Does the completed moving-throat PDE realize the natural passive/outgoing quadrupole branch with the correct far-zone normalization?}
}
\]

So this summary should be read as a **reduced theorem ledger** sitting on top of the already frozen 4D / 1PN / 2PN backbone.

---

## 1) Claim taxonomy and how to read this package

### 1.1 Exact

These are exact symbolic or algebraic statements once the declared lower-order ledger, basis choice, and source/port conventions are fixed.

Examples:

- the Burke–Thorne prototype landing on the Iyer–Will family,
- the odd-channel hierarchy \(n=1,3,5\),
- the dissipation/Schott identities,
- the real-STF identification of the solved \(P_2\) ports,
- the scalarity of the low-frequency \(P_2\to\mathrm{STF}\) matching matrix under isotropy,
- the outgoing \(\ell=2\) low-frequency fingerprint,
- and the grouped-\(P_2\) inverse formulas.

### 1.2 Controlled reduction

These statements require the same sort of reduced description already used below 2.5PN:

- near-zone / local radiation-reaction reduction,
- compact outgoing partial-wave completion,
- worldtube / small-body reduction,
- projection-language interpretation of leakage,
- and grouped-channel packaging of the throat response.

This package is not a theorem of the unreduced moving-throat PDE.

### 1.3 Protocol closure

These are statements that depend on the same declared closure hierarchy already used by the conservative 1PN and 2PN papers.

Examples:

- carrying forward the 1PN and 2PN support/response data,
- using the compact passive/outgoing small-body branch as the natural retarded completion,
- and using the grouped real \(P_2\) conservative module as the higher-order front end.

### 1.4 Conditional theorem

The final 2.5PN statement is **conditional rather than unconditional**.

The package shows that once the completed PDE satisfies the natural compact/passive/outgoing, isotropic, small-body conditions isolated by the audit, the reduced dynamics has the correct 2.5PN architecture and all currently derived lower odd channels are pushed above universal point-particle 2.5PN.

What is still missing is the final PDE-side normalization of the surviving quadrupole branch.

---

## 2) What this package claims — and what it does not claim

### 2.1 What it claims

This package claims all of the following.

1. The frozen conservative 1PN/2PN hierarchy by itself cannot generate a nonzero local 2.5PN reaction theorem.
2. A minimal retarded quadrupole odd term \(i\lambda_5\omega^5\) is sufficient to reproduce a standard valid 2.5PN force.
3. The natural low-frequency odd channels are classified by
   \[
   n=1\ \text{(scalar)},\qquad n=3\ \text{(dipole)},\qquad n=5\ \text{(quadrupole)}.
   \]
4. On the natural compact passive/outgoing small-body branch, the scalar and dipole dangers are demoted above universal point-particle 2.5PN.
5. The solved 2PN \(P_2\) bundle is exactly the real STF \(\ell=2\) content needed for the universal quadrupole route.
6. The surviving universal point-particle 2.5PN branch is the **orbital/worldtube STF quadrupole**, not the compact internal support quadrupole.
7. Once the grouped low-frequency conservative quadrupole coefficients through \(O(\omega^4)\) are fixed on the natural isotropic branch, the odd Burke–Thorne coefficient follows automatically.

### 2.2 What it does **not** claim

This package does **not** claim:

- a fully solved moving-throat PDE theorem;
- an assumption-free derivation of every retarded normalization datum;
- a proof that every possible anisotropic, non-passive, or noncompact branch is excluded by first principles;
- exact projection-locking of the scalar sector from the current file stack alone;
- exact matter/geometry scalar cancellation from the current file stack alone;
- or an unconditional finished 2.5PN theorem.

So the correct reading is:

> the 2.5PN program is structurally narrowed to a single surviving universal branch and a single narrow normalization gap, but the final PDE-side normalization is still open.

---

## 3) Lower-order anchors, notation firewall, and carry-forward objects

### 3.1 Carried lower-order backbone

Frozen from the earlier program:

\[
\kappa_\rho=1,
\qquad
n=5,
\qquad
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3,
\qquad
L/a\approx 1.85.
\]

The conservative backbone already contains:

- Newtonian / 0PN inverse-square closure,
- exact conservative 1PN EIH equality within the declared hierarchy,
- exact conservative 2PN ADM equality within the same hierarchy.

### 3.2 What 2PN already forced structurally

The 2PN paper already sharpened the defect into a **small throat-response theory** with:

- a carried odd dipole wake sector,
- a genuine even \(P_0\oplus P_2\) mouth/support layer,
- and a separate geometry-closure channel.

That is why the 2.5PN push is not “just another PN coefficient solve.” It is probing the first dissipative opening of an already nontrivial response system.

### 3.3 Frozen 2PN support data used here

The 2.5PN package repeatedly uses the 2PN static even support data

\[
R_0=R_{20}=R_{21}=R_{22}=1,
\qquad
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
\qquad
\Delta_{\rm geom}=\frac{281}{80},
\]

while treating the dynamic pole scales

\[
\Omega_{20},\qquad \Omega_{21},\qquad \Omega_{22}
\]

as genuinely open observables of the unresolved throat / DtN completion.

### 3.4 Notation firewall

This package is about the dissipative gravitational sector, not electric charge bookkeeping. Still, the corrected lower-order notation firewall remains in force:

- electric charge sign is not circulation,
- historical gravity-side bare `q=1` is really
  \[
  \kappa_\rho=1,
  \]
- and the grouped \(P_2\) indices \((20,21,22)\) always label real quadrupole channels, not spacetime tensor indices.

### 3.5 Canonical orbital quadrupole and grouped \(P_2\) basis map

The canonical source is the orbital/worldtube STF quadrupole

\[
I_{ij}=\mu\,x_{\langle i}x_{j\rangle}.
\]

In the real STF basis used by the audit, the solved 2PN \(P_2\) ports map to canonical STF components by

\[
q_{20}=\sqrt{\frac23}\,\Pi_{20},
\qquad
q_{21c}=\Pi_{21c},
\qquad
q_{21s}=\Pi_{21s},
\qquad
q_{22c}=2\Pi_{22c},
\qquad
q_{22s}=2\Pi_{22s},
\]

with invertible basis-change determinant

\[
\det B=\frac{4\sqrt6}{3}\neq 0.
\]

So the solved 2PN support bundle already contains the full real STF \(\ell=2\) representation content.

---

## 4) Headline outputs / carry-forward constants

This is the fastest memory ledger for future work.

### 4.1 Decisive benchmark

The constructive Burke–Thorne prototype gives the local 2.5PN force in standard form and lands exactly on the Burke–Thorne / Iyer–Will family member

\[
\boxed{\alpha=4,\qquad \beta=5.}
\]

So a minimal retarded quadrupole odd term is sufficient to reproduce a standard valid 2.5PN force.

### 4.2 Low-frequency odd-channel hierarchy

The audit scripts make Part III explicit:

\[
K_0^{\rm odd}(\omega)\sim i\gamma_1\omega,
\qquad
K_1^{\rm odd}(\omega)\sim i\gamma_3\omega^3,
\qquad
K_2^{\rm odd}(\omega)\sim i\gamma_5\omega^5.
\]

The time-domain signs are

\[
i\omega\ \mapsto\ -\frac{d}{dt},
\qquad
i\omega^3\ \mapsto\ +\frac{d^3}{dt^3},
\qquad
i\omega^5\ \mapsto\ -\frac{d^5}{dt^5}.
\]

And the standard Schott/dissipation identities are:

\[
\dot q(-\dot q)+\dot q^2=0,
\]
\[
\dot q\,q^{(3)} = \frac{d}{dt}(\dot q\,\ddot q)-\ddot q^2,
\]
\[
\dot q(-q^{(5)})
=
-\frac{d}{dt}\bigl(\dot q\,q^{(4)}-\ddot q\,q^{(3)}\bigr)-(q^{(3)})^2.
\]

The model-specific 2PN scalar/quadrupole combinations are

\[
\gamma_1^{\rm eff}=\frac{16}{5}\delta_{01}-\frac{281}{80}\delta_{g1},
\qquad
\gamma_5^{\rm eff}=\frac{25}{16}\delta_{205}.
\]

### 4.3 Scalar final ledger

The stable scalar conclusions are:

1. a direct scalar monopole overlap generically gives
   \[
   G_R(\omega,r)=\frac{e^{i\omega r/c_s}}{4\pi r}
   =
   \frac{1}{4\pi r}+\frac{i\omega}{4\pi c_s}-\frac{\omega^2 r}{8\pi c_s^2}-\frac{i\omega^3 r^2}{24\pi c_s^3}+\cdots,
   \]
   so a linear \(i\omega\) danger is real at the naive level;
2. the continuity-derived breathing outlet is **derivative-coupled**, not direct `q`-type:
   \[
   K_{qq}^R(\omega)=\omega^2\Pi_q^R(\omega),
   \qquad
   \Im K_{qq}^R(\omega)\propto \omega^3;
   \]
3. the mouth reciprocity / radiative-order theorem kills the dangerous Ohmic scalar bleed coefficient:
   \[
   \boxed{z_2=0;}
   \]
4. compact reciprocal mouth radiation is super-Ohmic / higher-order on the natural branch;
5. a genuinely Ohmic scalar bleed would require new structure outside the current compact passive/radiative ontology.

The one exact rescue condition that would fully kill the leading scalar outlet is projection-locking:

\[
B_A(q_0)v^A=0
\quad \text{for every dynamically allowed scalar tangent } v^A,
\]

equivalently

\[
\int d^3x\,dw\;W(w)\,\Psi_{(n)}(\mathbf x,w)=0
\quad \text{for every soft scalar mode }\Psi_{(n)}.
\]

That condition is **not** proven by the current file stack.

### 4.4 Dipole/vector final ledger

The stable dipole conclusions are:

1. the carried odd wake is real and survives into the relative sector;
2. naive COM cancellation and naive vector Ward identities do **not** kill it;
3. a clean outgoing 3D \(\ell=1\) completion forces the absorptive part to start at cubic order,
   \[
   \Im Y_1^R(\omega)\propto \omega^3,
   \]
   not linearly;
4. a same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN Iyer–Will family;
5. in a strict small-body branch it is naturally demoted above universal point-particle 2.5PN.

The scaling statement used repeatedly is

\[
\epsilon\,\bigl(\sqrt\epsilon\,a/r\bigr)^3
=
\epsilon^{5/2}\left(\frac{a}{r}\right)^3,
\qquad
\frac{a}{r}\sim\epsilon^\delta
\ \Longrightarrow\ 
\epsilon^{5/2+3\delta}.
\]

So for any \(\delta>0\), the outgoing dipole wake lies above universal point-particle 2.5PN.

### 4.5 Quadrupole structural ledger

The stable quadrupole conclusions are:

1. the solved 2PN \(P_2\) bundle is exactly the real STF \(\ell=2\) content;
2. a clean outgoing compact \(\ell=2\) completion gives the correct positive parity/sign,
   \[
   \Im Y_2^R(\omega)\propto +\omega^5;
   \]
3. the compact internal \(P_2\) mouth/support route is real but finite-size suppressed,
   \[
   \epsilon^2\bigl(\sqrt\epsilon\,a/r\bigr)^5
   =
   \epsilon^{9/2}\left(\frac{a}{r}\right)^5
   =
   \epsilon^{9/2+5\delta};
   \]
4. the universal point-particle route is the **orbital/worldtube STF quadrupole**;
5. rotational invariance forces the low-frequency \(P_2\to\mathrm{STF}\) matching matrix to be scalar;
6. the feared static-overlap failure does not occur naturally:
   \[
   m_0\neq 0,
   \qquad
   \hat m_0 = 1 + O\!\left(\frac{a^2}{r^2}\right)
   \]
   on the natural source-map branch.

### 4.6 Quadrupole normalization stack

The key normalization relations are

\[
\gamma_{\rm GR}=\frac{2G}{5c^5},
\qquad
\gamma_{\rm toy}=\hat m_0^2\Gamma_5,
\qquad
\boxed{\hat m_0^2\Gamma_5=\frac{2G}{5c^5}.}
\]

In the canonical invariant packaging used by the second audit script, the target quadrupole normalization stack is

\[
N_Q^{\rm target}=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_2^{\rm target}=\frac{6Gc_s^3}{5a^3c^5},
\qquad
\overline K_4^{\rm target}=\frac{8Gc_s}{15ac^5}.
\]

On the minimal isotropic quadrupole module,

\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\overline K_0=\frac{64G\Omega_Q^5}{45c^5}.
\]

In normalized-support language,

\[
\overline K_0^{\rm target} = \frac{2G}{45c^5u_2^{5/2}}.
\]

### 4.7 Best current theorem statement

The strongest current conditional theorem is:

> If the completed moving-throat PDE introduces no new direct scalar source, keeps the scalar and dipole outlets compact/passive/outgoing, preserves rotational invariance of the isolated low-frequency defect/worldtube sector, realizes the passive/outgoing quadrupole branch, and stays on a controlled small-body expansion with \(\hat m_0=1+O(a^2/r^2)\), then the reduced nonspinning two-body dynamics has the structure
> \[
> \mathbf a
> =
> \mathbf a_{\rm N}
> +\mathbf a_{1\rm PN}
> +\mathbf a_{2\rm PN}
> +\mathbf a_{2.5\rm PN}^{\rm quad}
> +\Delta\mathbf a_{\rm fs}^{>2.5\rm PN},
> \]
> where \(\mathbf a_{2.5\rm PN}^{\rm quad}\) is the standard local quadrupolar Burke–Thorne / Iyer–Will family member and all other currently-derived odd channels are pushed above universal point-particle 2.5PN.

---

## 5) Minimal dictionary and low-frequency bookkeeping

### 5.1 Small-body PN bookkeeping

The audit repeatedly uses

\[
\epsilon\equiv\left(\frac{v}{c}\right)^2,
\qquad
\frac{a}{r}\sim\epsilon^\delta,
\qquad
\delta>0.
\]

Universal point-particle 2.5PN means the unique \(c^{-5}\) local quadrupole slot without extra positive powers of \(a/r\).

### 5.2 Canonical quadrupole action form

In canonical STF coordinates \(q_a\), the unique local odd quadratic action is

\[
L_{\rm rr}^{(2)}=-\frac{\gamma}{2}\sum_a q_{a,-}\,q_{a,+}^{(5)},
\]

which produces

\[
a_i=-\gamma\,x^j I_{ij}^{(5)}.
\]

Matching to Burke–Thorne fixes \(\gamma_{\rm GR}=2G/(5c^5)\).

### 5.3 Scalar and quadrupole effective coefficients inherited from 2PN packaging

The scalar-like 2PN support combination is

\[
S_{\rm eff}(\omega)=J_0^2Y_0(\omega)+J_{20}^2Y_{20}(\omega)-\Delta_{\rm geom}Y_g(\omega),
\]

with

\[
J_0^2=\frac{16}{5},
\qquad
J_{20}^2=\frac{25}{16},
\qquad
\Delta_{\rm geom}=\frac{281}{80}.
\]

The quadrupole side is carried by the real STF \(P_2\) bundle.

### 5.4 Grouped real \(P_2\) moments

The conservative grouped-channel target uses the low-frequency coefficients

\[
u_2^{(20)},\quad u_2^{(21)},\quad u_2^{(22)},
\qquad
u_4^{(20)},\quad u_4^{(21)},\quad u_4^{(22)}.
\]

The axisymmetric trace/anisotropy variables are

\[
\bar u_2=\frac{u_2^{(20)}+2u_2^{(21)}+2u_2^{(22)}}{5},
\qquad
a_2=\frac{2u_2^{(20)}-u_2^{(21)}-u_2^{(22)}}{10},
\qquad
b_2=\frac{u_2^{(21)}-u_2^{(22)}}{2},
\]

and similarly for \((\bar u_4,a_4,b_4)\).

The isotropy gate is

\[
a_2=b_2=a_4=b_4=0.
\]

---

## 6) Main proof chain

### 6.1 Part I–II: why 2.5PN is the decisive gate

The pre-2.5PN state was:

- conservative Newtonian + 1PN + 2PN already strong,
- particle ontology already more than a scalar monopole,
- fully solved moving-throat PDE still open,
- and dynamic retarded/open-system response clearly the next missing ingredient.

The first decisive benchmark then showed:

1. a purely real-even low-frequency kernel
   \[
   K(\omega)=K_0+K_2\omega^2+K_4\omega^4+\cdots
   \]
   cannot generate a genuine local 2.5PN reaction term;
2. the minimal admissible escape hatch is a retarded odd quadrupole term
   \[
   K^{\rm ret}(\omega)=K^{\rm even}(\omega)+i\lambda_5\omega^5+O(\omega^7);
   \]
3. a Burke–Thorne prototype built from the canonical STF quadrupole lands exactly on the standard local 2.5PN family with \((\alpha,\beta)=(4,5)\).

So after Part II the problem was no longer vague: derive the retarded quadrupole response law whose first physical odd term is \(i\omega^5\).

### 6.2 Part III: the general low-frequency selection-rule framework

The audit scripts make the missing Part III explicit.

The key framework ingredients are:

- wrong odd channels classified by \(n=1,3,5\);
- minimal retarded kernel prototypes
  \[
  K_1\sim \frac{g^2}{\Omega^2-\omega^2-i\sigma\omega},
  \qquad
  K_3\sim \frac{g^2}{\Omega^2-\omega^2-i\sigma\omega^3},
  \qquad
  K_5\sim \frac{g^2}{\Omega^2-\omega^2-i\sigma\omega^5};
  \]
- and the dissipation/Schott identities showing how odd kernels split into true dissipation plus total derivatives.

This is the formal reason scalar, dipole, and quadrupole audits can be treated as separate theorem gates.

### 6.3 Part IV: scalar sector

The scalar story proceeds in four steps.

First, the danger is real: a generic outgoing scalar channel gives a linear \(i\omega\) odd term, and the 2PN scalar packaging already organizes the dangerous combination as

\[
\gamma_1^{\rm eff}=\frac{16}{5}\delta_{01}-\frac{281}{80}\delta_{g1}.
\]

Second, naive automatic rescues fail. Exact projected continuity and the no-static-flux theorem do **not** by themselves force \(\gamma_1^{\rm eff}=0\). The Gaussian overlap counterexample makes that explicit.

Third, the breathing channel is rescued structurally: continuity forces the breathing-to-outlet vertex to be `\dot q`-like rather than direct `q`-like, so the effective breathing kernel starts at cubic odd order,

\[
\Im K_{qq}^R(\omega)\propto \omega^3,
\]

not linearly.

Fourth, the compact mouth route is rescued by reciprocity and radiative-order counting, which kill the dangerous Ohmic mouth term \(z_2\). What remains is super-Ohmic / finite-size.

The result is not an unconditional scalar theorem, but the scalar sector is no longer the default universal kill switch.

### 6.4 Part V: dipole/vector sector

The carried odd wake survives into the relative sector, so naive COM arguments and naive vector Ward identities do **not** save the model.

The crucial rescue is spectral and geometric rather than algebraic:

- a clean outgoing 3D \(\ell=1\) completion starts at \(+i k^3\), not \(ik\);
- same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN family;
- therefore strict small-body scaling is the only viable rescue.

On that branch, the outgoing dipole correction scales as

\[
\epsilon^{5/2+3\delta},
\qquad \delta>0,
\]

so it is demoted above universal point-particle 2.5PN.

### 6.5 Part VI: quadrupole sector

The quadrupole side is the surviving universal branch.

The proof chain is:

1. distinguish the compact internal \(P_2\) mouth/support sector from the universal orbital/worldtube STF quadrupole;
2. show exactly that the solved 2PN \(P_2\) ports are the full real STF \(\ell=2\) basis;
3. show that a clean outgoing compact \(\ell=2\) completion automatically produces the right parity and sign,
   \[
   +i\omega^5,
   \]
   with positive damping;
4. show that the literal compact internal \(P_2\) route is finite-size suppressed and therefore not the universal 2.5PN source;
5. construct the orbital/worldtube source map
   \[
   I_{ij}^{\rm orb}=\mu x_{\langle i}x_{j\rangle};
   \]
6. prove that rotational invariance forces the low-frequency \(P_2\to\mathrm{STF}\) matching matrix to be scalar;
7. show that the feared zero-overlap branch does not occur naturally:
   \[
   m_0\neq 0.
   \]

At that point the whole problem narrows to normalization.

### 6.6 Part VII: quadrupole normalization, extraction, and grouped \(P_2\) module

The second master audit extends the package by making the normalization stack explicit.

The load-bearing steps are:

1. canonical normalization map:
   \[
   \gamma_{\rm toy}=\hat m_0^2\Gamma_5,
   \qquad
   \hat m_0^2\Gamma_5=\frac{2G}{5c^5};
   \]
2. outgoing \(\ell=2\) PDE fingerprint:
   \[
   \Lambda_2(k)
   =
   -\frac{3}{a}
   +\frac{a}{3}k^2
   +\frac{a^3}{9}k^4
   +i\frac{a^4}{9}k^5
   -\frac{2a^5}{27}k^6+\cdots;
   \]
3. normalized outgoing branch:
   \[
   \widehat Y_2(k)=1+\frac{a^2k^2}{9}+\frac{4a^4k^4}{81}+i\frac{a^5k^5}{27}+\cdots;
   \]
4. invariant moment identities:
   \[
   m_4=\frac{4m_2^2}{m_0},
   \qquad
   \Gamma_5^{\rm PDE}=9\frac{m_2^{5/2}}{m_0^{3/2}};
   \]
5. source/port cleanup establishing \(\hat m_0=1\) on the point-particle orbital branch;
6. one-pole insufficiency and the minimal positive conservative precursor;
7. exact moment-matched minimal isotropic quadrupole module;
8. grouped \(P_2\) extraction formulas and the isotropy/branch identities.

This is the stage where the remaining theorem gap becomes purely the final passive/outgoing quadrupole normalization on the actual moving-throat branch.

---

## 7) What is effectively fixed vs what remains open

### 7.1 Effectively fixed by the current package

Treat the following as stable current results.

#### Conservative decisive benchmark

\[
(\alpha,\beta)=(4,5)
\]
for the Burke–Thorne prototype.

#### Low-frequency odd-channel ordering

\[
\text{scalar }\sim i\omega,
\qquad
\text{dipole }\sim i\omega^3,
\qquad
\text{quadrupole }\sim i\omega^5.
\]

#### Scalar rescues on the natural branch

- breathing outlet derivative-coupled,
- compact reciprocal mouth radiation super-Ohmic,
- no derived independent linear scalar source beyond continuity + reciprocal boundary motion.

#### Dipole rescue on the natural branch

- outgoing \(\ell=1\) starts at \(+ik^3\),
- same-order survival would be fatal,
- strict small-body scaling demotes it above universal point-particle 2.5PN.

#### Quadrupole structural route

- solved 2PN \(P_2\) bundle equals real STF \(\ell=2\),
- matching matrix scalar under isotropy,
- positive outgoing \(+i\omega^5\) sign,
- nonzero static overlap,
- surviving universal source is orbital/worldtube STF quadrupole.

#### Normalization/extraction module

- canonical target relation \(\hat m_0^2\Gamma_5=2G/(5c^5)\),
- outgoing \(\ell=2\) moment identities,
- grouped \(P_2\) isotropy and minimal-branch formulas.

### 7.2 Fixed only within the declared hierarchy

These are stable only within the same reduced hierarchy already used by the conservative 1PN/2PN papers.

- compact passive/outgoing outlet choice;
- strict small-body expansion;
- isotropic low-frequency quadrupole branch;
- grouped real \(P_2\) conservative module as the correct front end;
- no new direct scalar source beyond the continuity + reciprocal-boundary ontology.

### 7.3 Still genuinely open

These are the real unresolved questions.

1. the final passive/outgoing quadrupole normalization on the actual moving-throat branch;
2. whether the completed PDE stays on the natural isotropic passive/outgoing branch;
3. exclusion of every possible anisotropic, singular, or non-passive low-frequency completion;
4. exact projection-locking or exact scalar cancellation from the file stack alone;
5. fully solving the moving-throat PDE itself.

The live unresolved quantity can be stated as

\[
\Gamma_5\stackrel{?}{=}\Gamma_{5,\rm GR}
\qquad \text{or equivalently} \qquad
\hat m_0^2\Gamma_5\stackrel{?}{=}\frac{2G}{5c^5}.
\]

---

## 8) Appendix digest — quadrupole normalization and grouped \(P_2\) formulas

### 8.1 Canonical quadrupole normalization target

The canonical local odd action is

\[
L_{\rm rr}^{(2)}=-\frac{\gamma}{2}\sum_a q_{a,-}q_{a,+}^{(5)},
\]

so Burke–Thorne matching gives

\[
\gamma_{\rm GR}=\frac{2G}{5c^5}.
\]

In toy-model variables,

\[
\gamma_{\rm toy}=\hat m_0^2\Gamma_5,
\qquad
\boxed{\hat m_0^2\Gamma_5=\frac{2G}{5c^5}.}
\]

### 8.2 Outgoing \(\ell=2\) PDE fingerprint

For the compact outgoing \(\ell=2\) branch,

\[
\Lambda_2(k)
=
-\frac{3}{a}
+\frac{a}{3}k^2
+\frac{a^3}{9}k^4
+i\frac{a^4}{9}k^5
-\frac{2a^5}{27}k^6
+O(k^7),
\]

and the normalized admittance is

\[
\widehat Y_2(k)
=
1+\frac{a^2k^2}{9}+\frac{4a^4k^4}{81}+i\frac{a^5k^5}{27}+O(k^6).
\]

So

\[
m_2=\frac{a^2m_0}{9c_s^2},
\qquad
m_4=\frac{4a^4m_0}{81c_s^4},
\qquad
\Gamma_5^{\rm PDE}=\frac{a^5m_0}{27c_s^5},
\]

and therefore

\[
\boxed{m_4=\frac{4m_2^2}{m_0},}
\qquad
\boxed{\Gamma_5^{\rm PDE}=9\frac{m_2^{5/2}}{m_0^{3/2}}.}
\]

### 8.3 Source/port cleanup and invariant products

The mixed-basis static overlap is a convention artifact; after same-basis cleanup and orbital normalization,

\[
\hat m_0=1
\]

on the pure point-particle orbital branch, with worldtube corrections

\[
\hat m_0=1+O(a^2/r^2).
\]

The basis-invariant low-frequency products are

\[
\mathfrak m_0^2K_0,
\qquad
\mathfrak m_0^2\Gamma_5^{\rm port},
\]

and the target normalization is

\[
N_Q^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
\]

### 8.4 One-pole insufficiency and the minimal isotropic quadrupole module

A naive one-pole completion is not enough at \(O(\omega^4)\). The minimal positive conservative precursor is

\[
Y_{Q,\min}(\omega)
=
\frac{3(a^2\omega^2-3c_s^2)}{4a^2\omega^2-9c_s^2}.
\]

It reproduces the desired low-frequency moments and gives the exact moment-matched isotropic solution

\[
\boxed{\Omega_Q=\frac{3c_s}{2a}.}
\]

On this branch,

\[
\overline K_0 = \frac{64G\Omega_Q^5}{45c^5}.
\]

### 8.5 Grouped real \(P_2\) extraction formulas

For grouped coefficients \((x_{20},x_{21},x_{22})\), the exact inverse map is

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]

The axisymmetric grouped expansion is

\[
Y_{20}=1+(\bar u_2+4a_2)\omega^2+(\bar u_4+4a_4)\omega^4+\cdots,
\]
\[
Y_{21}=1+(\bar u_2-a_2+b_2)\omega^2+(\bar u_4-a_4+b_4)\omega^4+\cdots,
\]
\[
Y_{22}=1+(\bar u_2-a_2-b_2)\omega^2+(\bar u_4-a_4-b_4)\omega^4+\cdots.
\]

The anisotropy norms are

\[
\mathcal A_2^2=4a_2^2+\frac45 b_2^2,
\qquad
\mathcal A_4^2=4a_4^2+\frac45 b_4^2.
\]

The isotropy gate is

\[
a_2=b_2=a_4=b_4=0.
\]

If it passes, the minimal-branch identity is

\[
\boxed{\bar u_4=4\bar u_2^2,}
\]

and therefore

\[
\Omega_Q^2=\frac{1}{4\bar u_2},
\qquad
\Gamma_{5,\rm norm}=9\bar u_2^{5/2},
\qquad
\boxed{\overline K_0^{\rm target}=\frac{2G}{45c^5\bar u_2^{5/2}}.}
\]

This is the exact grouped-channel front end that future higher-order conservative work must hit.

---

## 9) Verification / SymPy / reproducibility ledger

The 2.5PN audit archive is split into two executable SymPy referee files.

### 9.1 `2_5pn_master_sympy_audit.py`

This is the first master audit. It verifies:

- the decisive Burke–Thorne benchmark,
- the low-frequency selection-rule / influence-functional framework,
- the scalar-sector no-gos and rescues,
- the dipole/vector-sector no-gos and rescues,
- the first quadrupole representation / source-map / static-overlap results,
- and the final **conditional theorem ledger**.

Its closing theorem ledger is:

- Burke–Thorne prototype lands on \((\alpha,\beta)=(4,5)\),
- scalar direct overlap generically gives \(i\omega\),
- derivative breathing outlet starts at \(i\omega^3\),
- compact reciprocal mouth radiation is super-Ohmic,
- outgoing \(\ell=1\) starts at \(+ik^3\),
- compact outgoing \(\ell=2\) starts at \(+ik^5\),
- and the remaining narrow gap lives in the final passive/outgoing quadrupole normalization.

### 9.2 `2_5pn_master_session_sympy_audit.py`

This second audit extends the first by verifying:

- the canonical quadrupole normalization map,
- the outgoing \(\ell=2\) PDE fingerprint,
- source/port normalization cleanup,
- convention-invariant normalization products,
- the canonical invariant low-frequency pair \((\overline K_0,\overline K_2)\),
- one-pole insufficiency,
- the minimal isotropic quadrupole response module,
- grouped / normalized-support \(P_2\) extraction formulas,
- the axisymmetric grouped-\(P_2\) parameterization,
- and the final **session theorem ledger**.

Its closing theorem ledger states that the full GR-like point-particle 2.5PN closure has been reduced to the **isotropic real \(P_2\) quadrupole module**, and that once the grouped low-frequency conservative coefficients through \(O(\omega^4)\) are fixed on that branch, the odd Burke–Thorne coefficient follows automatically.

### 9.3 What the archive verifies — and what it does not

The archive checks the full symbolic derivation chain **within the declared closure hierarchy**. It verifies the load-bearing identities, basis maps, sign statements, scaling laws, and normalization formulas used in the summary.

It does **not** replace the missing moving-throat PDE derivation, and it does not upgrade the package into an unconditional theorem of the full unresolved bulk dynamics.

---

## 10) Short list of the most important formulas to remember

If you only need the highest-value equations for future work, keep these.

### 10.1 Lower-order carry-forward data

\[
\kappa_\rho=1,
\qquad
n=5,
\qquad
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3,
\qquad
L/a\approx 1.85.
\]

### 10.2 Decisive benchmark

\[
\alpha=4,
\qquad
\beta=5.
\]

### 10.3 Odd-channel hierarchy

\[
\text{scalar }\sim i\omega,
\qquad
\text{dipole }\sim i\omega^3,
\qquad
\text{quadrupole }\sim i\omega^5.
\]

### 10.4 Scalar effective coefficients

\[
\gamma_1^{\rm eff}=\frac{16}{5}\delta_{01}-\frac{281}{80}\delta_{g1},
\qquad
K_{qq}^R(\omega)=\omega^2\Pi_q^R(\omega),
\qquad
\Im K_{qq}^R(\omega)\propto \omega^3,
\qquad
z_2=0.
\]

### 10.5 Small-body scaling of the demoted channels

\[
\text{scalar breathing} \sim \epsilon^{5/2+3\delta},
\qquad
\text{outgoing dipole} \sim \epsilon^{5/2+3\delta},
\qquad
\text{compact internal }P_2 \sim \epsilon^{9/2+5\delta}.
\]

### 10.6 Universal surviving branch

\[
\boxed{I_{ij}^{\rm orb}=\mu x_{\langle i}x_{j\rangle}.}
\]

### 10.7 Quadrupole normalization target

\[
\boxed{\hat m_0^2\Gamma_5=\frac{2G}{5c^5}.}
\]

### 10.8 Outgoing \(\ell=2\) moment identities

\[
\widehat Y_2(k)=1+\frac{a^2k^2}{9}+\frac{4a^4k^4}{81}+i\frac{a^5k^5}{27}+\cdots,
\]
\[
m_4=\frac{4m_2^2}{m_0},
\qquad
\Gamma_5^{\rm PDE}=9\frac{m_2^{5/2}}{m_0^{3/2}}.
\]

### 10.9 Minimal isotropic quadrupole module

\[
Y_{Q,\min}(\omega)=\frac{3(a^2\omega^2-3c_s^2)}{4a^2\omega^2-9c_s^2},
\qquad
\Omega_Q=\frac{3c_s}{2a},
\qquad
\overline K_0=\frac{64G\Omega_Q^5}{45c^5}.
\]

### 10.10 Grouped \(P_2\) minimal-branch formulas

\[
a_2=b_2=a_4=b_4=0,
\qquad
\bar u_4=4\bar u_2^2,
\qquad
\Omega_Q^2=\frac{1}{4\bar u_2},
\qquad
\overline K_0^{\rm target}=\frac{2G}{45c^5\bar u_2^{5/2}}.
\]

### 10.11 Best current theorem statement

\[
\boxed{
\text{A full GR-like point-particle 2.5PN theorem is conditionally reachable from the present stack.}
}
\]

### 10.12 Remaining narrow gap

\[
\boxed{
\text{The remaining serious gap is the final passive/outgoing quadrupole normalization on the actual moving-throat branch.}
}
\]

---

## 11) Bottom-line interpretation

This package does four things at once.

1. It proves that the frozen conservative hierarchy alone cannot produce a local 2.5PN theorem.
2. It proves that the minimal retarded quadrupole route is structurally sufficient by the Burke–Thorne benchmark.
3. It shows that the lower odd scalar and dipole channels are no longer generic universal spoilers on the natural compact passive/outgoing small-body branch.
4. It reduces the whole remaining theorem problem to the normalization of one surviving universal branch: the orbital/worldtube STF quadrupole.

So the 2.5PN program is no longer blocked by diffuse lower-order leakage scenarios.

What remains open is **not** whether the model has a viable 2.5PN architecture in principle. What remains open is whether the completed moving-throat PDE realizes the already identified passive/outgoing quadrupole branch with the **correct far-zone normalization**.

That is why this package should be treated as a major narrowing milestone:

- **not yet an unconditional theorem,**
- **no longer a broad structural crisis,**
- **and close enough that the remaining work is a sharply focused PDE/normalization problem rather than a full-model rescue.**
