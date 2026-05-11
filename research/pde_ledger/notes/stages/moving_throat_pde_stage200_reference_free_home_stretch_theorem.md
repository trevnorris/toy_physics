# Moving-Throat PDE — Stage 200: Exact Reference-Free Full Home-Stretch Theorem, Orbit-Representative Independence, and the Four-Scalar Final Verdict Packet

## Status

**Exact within the carried Stage 248 Packet-A finish-line theorem, the Stage 250 exact pairwise orbit-transport / orbit-lock theorem, and the Stage 242 two-packet home-stretch hierarchy.**

This stage does **not** introduce a new constitutive law.
It is the exact reference-free completion of the reduced home stretch.

---

## Purpose

Stage 242 had already reduced the endgame to two finite packets:
\[
\Delta_{\rm branch},
\qquad
\Delta_{\rm orbit}.
\]
Stage 248 then collapsed the full Packet-A branch-side verdict to the single scalar finish line
\[
\Delta_{\rm branch}=0
\iff
\chi_Q=1.
\]
Stage 250, in turn, removed the last orbit-side base-point privilege by proving the exact two-point theorem
\[
\mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)}
\iff
q_{\rm tr}^{(2\leftarrow1)}=q_{\rm nt}^{(2\leftarrow1)}=q_\eta^{(2\leftarrow1)}=0.
\]

So after Stages 248 and 250, the natural next question is already forced:

> can the entire reduced home stretch be written as one exact **reference-free** verdict packet with no privileged orbit representative left anywhere?

This stage answers that completely.

The main outputs are:

1. the exact proof that the Packet-B orbit packet can be attached to an **orbit** rather than to a chosen orbit representative,
2. the exact four-scalar **full reduced verdict packet**
   \[
   \Delta_{\rm full},
   \]
3. the equivalent multiplicative and mismatch versions of the same packet,
4. the exact theorem
   \[
   \boxed{\Delta_{\rm full}=0
   \iff
   \Delta_{\rm branch}=0\ \text{and}\ \mathbf x\in\mathcal O_*},
   \]
   where `\(\mathcal O_*\)` is the target similarity orbit,
5. and the sharp statement that the reduced endgame now consists of **one Packet-A scalar plus one Packet-B quotient triple and nothing else.**

So Stage 251 is the exact reference-free upgrade of the Stage 242 home-stretch theorem.

Script-backed status:
- `scripts/moving_throat_pde_stage251_reference_free_home_stretch_theorem_sympy_audit.py`
  checks the exact orbit-relative packet identities, chart conversions, and the
  final linearized four-scalar compiler.
- `mathematica/moving_throat_pde_stage251_reference_free_home_stretch_theorem_mathematica_audit.wl`
  mirrors the same checkpoint algebra in the second CAS without introducing
  extra numeric literals.

---

## 1. Carry-forward Packet-A and Packet-B sharpeners

### 1.1 Packet A after Stage 248

Within the carried isotropic grouped-real `P2` one-pole front end and the natural point-particle source-map branch, Stage 248 proved
\[
\boxed{\Delta_{\rm branch}=0 \iff \chi_Q=1.}
\]
Equivalently,
\[
\boxed{\Delta_Q:=\chi_Q-1=0,}
qquad
\boxed{N_Q=1,}
\]
with
\[
\boxed{N_Q=\chi_Q^{-1}}
\]
on the same natural source-map branch.

So Packet A no longer carries eight independent residual slots. It carries one scalar only.

### 1.2 Packet B after Stage 250

Stage 250 proved the exact two-point orbit-lock equivalences
\[
\boxed{
\mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)}
\iff
m_T^{(2\leftarrow1)}=m_K^{(2\leftarrow1)}=m_\mu^{(2\leftarrow1)}=1,
}
\]
\[
\boxed{
\mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)}
\iff
\mathfrak R_{\rm tr}^{(2\leftarrow1)}=\mathfrak R_{\rm nt}^{(2\leftarrow1)}=\mathfrak R_\eta^{(2\leftarrow1)}=1,
}
\]
\[
\boxed{
\mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)}
\iff
q_{\rm tr}^{(2\leftarrow1)}=q_{\rm nt}^{(2\leftarrow1)}=q_\eta^{(2\leftarrow1)}=0.
}
\]

So Packet B is already an exact three-scalar orbit-lock packet.

---

## 2. Exact orbit-representative independence

Let the target similarity orbit be
\[
\boxed{\mathcal O_*:=\mathcal G_*\!\cdot\mathbf x_*\subset\mathcal M_+.}
\]
Take any two representatives
\[
\mathbf x_*,\ \widetilde{\mathbf x}_*\in\mathcal O_*.
\]
Since they lie on the same exact similarity orbit, Stage 250 gives
\[
\boxed{
q_{\rm tr}^{(\widetilde *\leftarrow *)}
=
q_{\rm nt}^{(\widetilde *\leftarrow *)}
=
q_\eta^{(\widetilde *\leftarrow *)}=0.
}
\]
Using the exact Packet-B cocycle law,
\[
\mathbf q^{(x\leftarrow \widetilde *)}
=
\mathbf q^{(x\leftarrow *)}+
\mathbf q^{(*\leftarrow \widetilde *)},
\]
we immediately get
\[
\boxed{
\mathbf q^{(x\leftarrow \widetilde *)}=
\mathbf q^{(x\leftarrow *)}.
}
\]

The same argument works for the mismatch and invariant-ratio packets because Stage 250 proved the exact multiplicative composition laws
\[
\boxed{
m_\bullet^{(3\leftarrow1)}=m_\bullet^{(3\leftarrow2)}m_\bullet^{(2\leftarrow1)},
\qquad
\mathfrak R_\bullet^{(3\leftarrow1)}=
\mathfrak R_\bullet^{(3\leftarrow2)}\mathfrak R_\bullet^{(2\leftarrow1)}.
}
\]
Since
\[
m_\bullet^{(\widetilde *\leftarrow *)}=1,
\qquad
\mathfrak R_\bullet^{(\widetilde *\leftarrow *)}=1,
\]
we have
\[
\boxed{
m_\bullet^{(x\leftarrow \widetilde *)}=m_\bullet^{(x\leftarrow *)},
\qquad
\mathfrak R_\bullet^{(x\leftarrow \widetilde *)}=\mathfrak R_\bullet^{(x\leftarrow *)}.
}
\]

So the Packet-B orbit packet is now genuinely attached to the **orbit** rather than to an orbit representative.

### 2.1 Orbit-relative notation

Therefore the following are well defined for any candidate microscopic state `\(\mathbf x\)`:
\[
\boxed{
\mathbf q^{(x\leftarrow \mathcal O_*)},
\qquad
\mathbf m^{(x\leftarrow \mathcal O_*)},
\qquad
\boldsymbol{\mathfrak R}^{(x\leftarrow \mathcal O_*)}.
}
\]
By definition, one may compute them against **any** orbit representative `\(\mathbf x_*\in\mathcal O_*\)`.

This is the exact sense in which the Packet-B side is now reference free.

---

## 3. Exact four-scalar full verdict packet

### 3.1 Additive chart

Define the exact four-scalar additive full verdict packet by
\[
\boxed{
\Delta_{\rm full}^{(x\mid \mathcal O_*)}
:=
\Bigl(
\Delta_Q(x),
q_{\rm tr}^{(x\leftarrow \mathcal O_*)},
q_{\rm nt}^{(x\leftarrow \mathcal O_*)},
q_\eta^{(x\leftarrow \mathcal O_*)}
\Bigr),
}
\]
where
\[
\boxed{\Delta_Q(x):=\chi_Q(x)-1.}
\]

So the full reduced verdict is now four scalars only.

### 3.2 Multiplicative chart

Define the exact multiplicative chart of the same packet by
\[
\boxed{
\mathcal V_{\rm full}^{(x\mid \mathcal O_*)}
:=
\Bigl(
\chi_Q(x),
\mathfrak R_{\rm tr}^{(x\leftarrow \mathcal O_*)},
\mathfrak R_{\rm nt}^{(x\leftarrow \mathcal O_*)},
\mathfrak R_\eta^{(x\leftarrow \mathcal O_*)}
\Bigr).
}
\]
On the natural point-particle source-map branch one may equally use `\(N_Q\)` in the first slot because
\[
\boxed{N_Q=\chi_Q^{-1},}
qquad
\boxed{\chi_Q=1 \iff N_Q=1.}
\]

### 3.3 Mismatch chart

Define the exact mismatch version by
\[
\boxed{
\mathcal M_{\rm full}^{(x\mid \mathcal O_*)}
:=
\Bigl(
\chi_Q(x),
 m_T^{(x\leftarrow \mathcal O_*)},
 m_K^{(x\leftarrow \mathcal O_*)},
 m_\mu^{(x\leftarrow \mathcal O_*)}
\Bigr).
}
\]

### 3.4 Exact chart conversion laws

The Packet-B conversion laws from Stage 249/250 become, in orbit-relative form,
\[
\boxed{
\mathfrak R_{\rm tr}=e^{q_{\rm tr}}=(m_T)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak R_\eta=e^{q_\eta}=\frac{1}{m_K},
}
\]
\[
\boxed{
\mathfrak R_{\rm nt}=e^{q_{\rm nt}}=\frac{m_\mu}{m_K m_T^{F_*}}.
}
\]
Equivalently,
\[
\boxed{m_T=\exp\!\left(\frac{q_{\rm tr}}{1+\chi_{0,*}}\right),}
\qquad
\boxed{m_K=e^{-q_\eta},}
\]
\[
\boxed{m_\mu=\exp\!\left(q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}\right).}
\]
So the additive, multiplicative, and mismatch forms of the full four-scalar packet are exact reparameterizations of the same final verdict.

---

## 4. Exact reference-free full home-stretch theorem

\[
\boxed{\textbf{Theorem (Stage 251 exact reference-free full home-stretch theorem).}}
\]

**Within the carried Stage 248 Packet-A finish-line hierarchy and the Stage 250 exact pairwise orbit-lock hierarchy, the following are equivalent for any positive microscopic state `\(\mathbf x\)` and any target similarity orbit `\(\mathcal O_*\)`:**

1. **full reduced closure relative to the target orbit:**
   \[
   \boxed{\text{`\(\mathbf x\)` is fully reduced-admissible relative to `\(\mathcal O_*\)`};}
   \]
2. **the original Stage 242 two-packet criterion:**
   \[
   \boxed{\Delta_{\rm branch}(x)=0\quad\text{and}\quad \mathbf x\in\mathcal O_*;}
   \]
3. **the sharpened Packet-A + Packet-B scalar/quoitent form:**
   \[
   \boxed{
   \chi_Q(x)=1,
   \qquad
   q_{\rm tr}^{(x\leftarrow \mathcal O_*)}=q_{\rm nt}^{(x\leftarrow \mathcal O_*)}=q_\eta^{(x\leftarrow \mathcal O_*)}=0;
   }
   \]
4. **the multiplicative four-scalar verdict packet equals its canonical value:**
   \[
   \boxed{\mathcal V_{\rm full}^{(x\mid \mathcal O_*)}=(1,1,1,1);}
   \]
5. **the mismatch four-scalar packet equals its canonical value:**
   \[
   \boxed{\mathcal M_{\rm full}^{(x\mid \mathcal O_*)}=(1,1,1,1);}
   \]
6. **the additive four-scalar verdict packet vanishes:**
   \[
   \boxed{\Delta_{\rm full}^{(x\mid \mathcal O_*)}=0.}
   \]

### 4.1 Proof in one line

Stage 242 gave
\[
\text{full reduced closure}
\iff
\Delta_{\rm branch}=0\ \text{and}\ \Delta_{\rm orbit}=0.
\]
Stage 248 sharpened the first factor to
\[
\Delta_{\rm branch}=0\iff\chi_Q=1,
\]
and Stage 250 sharpened the second to the exact orbit-relative condition
\[
\Delta_{\rm orbit}=0
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]
with the packet already reference-independent at the orbit level. Combining the two gives the theorem.

---

## 5. Exact one-sided pairwise corollary

Suppose now that `\(\mathbf x^{(1)}\)` is already known to be fully reduced closed relative to the target orbit `\(\mathcal O_*\)`. Then
\[
\chi_Q\bigl(\mathbf x^{(1)}\bigr)=1,
\qquad
\mathbf x^{(1)}\in\mathcal O_*.
\]
By orbit-representative independence,
\[
\mathbf q^{(2\leftarrow1)}
=
\mathbf q^{(\mathbf x^{(2)}\leftarrow \mathcal O_*)}.
\]
So for any second candidate state `\(\mathbf x^{(2)}\)` the exact one-sided full verdict packet is simply
\[
\boxed{
\Delta_{\rm full}^{(2\leftarrow1)}
:=
\Bigl(
\chi_Q\bigl(\mathbf x^{(2)}\bigr)-1,
q_{\rm tr}^{(2\leftarrow1)},
q_{\rm nt}^{(2\leftarrow1)},
q_\eta^{(2\leftarrow1)}
\Bigr).
}
\]
Then
\[
\boxed{
\Delta_{\rm full}^{(2\leftarrow1)}=0
\iff
\mathbf x^{(2)}\ \text{is fully reduced closed relative to the same target orbit}.
}
\]

So once one fully closed state on the target orbit is known, every further candidate can be tested by one Packet-A scalar and one direct pairwise Packet-B quotient triple.

---

## 6. Exact two-state closed-branch corollary

For any two positive microscopic states `\(\mathbf x^{(1)},\mathbf x^{(2)}\)`,
\[
\boxed{
\chi_Q\bigl(\mathbf x^{(1)}\bigr)=\chi_Q\bigl(\mathbf x^{(2)}\bigr)=1,
\qquad
q_{\rm tr}^{(2\leftarrow1)}=q_{\rm nt}^{(2\leftarrow1)}=q_\eta^{(2\leftarrow1)}=0
}
\]
is equivalent to the statement that both states are fully reduced closed and lie on the same exact target similarity orbit.

So the completed PDE may be benchmarked either:

- against a chosen target orbit, or
- directly between any two candidate fully closed realizations.

---

## 7. What Stage 251 changes in the theorem problem

Stage 242 reduced the endgame to two packets, but not yet to its sharpest form.
Stages 248 and 250 removed the remaining redundancies separately.
Stage 251 is the exact combined endpoint.

### 7.1 The full reduced endgame is now four scalars only

The completed moving-throat PDE no longer has to be judged against an eight-slot Packet-A vector plus a separate orbit packet. After the carried sharpeners, the whole reduced verdict is
\[
\boxed{(\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta).}
\]

### 7.2 No privileged orbit representative remains

The Packet-B orbit packet is now attached to the target similarity orbit itself. Any convenient representative may be used in practice, but the verdict is mathematically independent of that choice.

### 7.3 Everything still open is a realization problem, not an algebra problem

The reduced compiler algebra is finished.
What remains is only:

1. realize the Packet-A scalar finish line `\(\chi_Q=1\)`, and
2. realize the Packet-B orbit-relative quotient triple `\(\mathbf q=0\)`.

No other reduced endgame coordinate survives.

---

## 8. Immediate next derivation step

The natural continuation after Stage 251 is now procedural rather than algebraic:

1. feed the actual completed moving-throat PDE output into the Packet-A compiler to extract `\(\chi_Q\)`,
2. feed the same output into the Packet-B orbit-relative compiler to extract
   \[
   (q_{\rm tr},q_{\rm nt},q_\eta),
   \]
3. and evaluate the exact four-scalar final verdict packet
   \[
   \Delta_{\rm full}^{(x\mid \mathcal O_*)}.
   \]

That is the sharpest possible reduced theorem gate after Stage 251.
