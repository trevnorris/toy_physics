# Moving-Throat PDE — Stage 238: Physical-Branch Transfer-Shape Compiler, Packet Factorization, and the Post-Static Dressing-Invariance Theorem

## Status

**Exact within the carried Stage-237 actual-branch dressing compiler and the later coherent-branch observable compiler already frozen in the moving-throat notes.**

This stage does **not** solve the full moving-throat PDE branch.
It does the next sharp reduction after Stage 237:

> it rewrites the actual same-charge packet directly in the physical coherent-branch variables
> \((R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta)\),
> shows that the static same-charge ceiling is exactly transfer-shape rigidity,
> and isolates the post-static barrier as pure dressing invariance.

So after this stage the same-charge chain factors into three clean gates:

1. tracking rigidity,
2. transfer-shape rigidity,
3. dressing rigidity.

And once the first two have been cleared, the barrier reduces to one scalar test only:
\[
\delta\ln\epsilon_\eta = 0.
\]

---

## Purpose

Stage 237 proved that on the rigid-mouth slice the surviving post-static obstruction is
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\mathrm{ref}}}\right),
\]
and that this scalar is exactly support-blind.

But Stage 237 still left one useful bridge implicit:

> what does the static same-charge gate `q_nt = 0` look like in the actual coherent branch variables before the dressing scalar is read off?

The later coherent-branch observable compiler already supplies the missing object. On that branch the target observable is not primitive; it is tied to one wall-normalized transfer shape
\[
\mathcal T^2.
\]
Once that identity is inserted, the same-charge packet becomes triangular in the physical variables themselves.

---

## 0. Why this stage is needed

Stage 237 showed that on the rigid-mouth slice the surviving actual-branch packet is
\[
q_{\rm nt}
=
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,\rm ref}}\right),
\qquad
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]
So after the first static gate is cleared, the same-charge obstruction is indeed the dressing coordinate `q_eta`.

But the remaining static gate `q_nt = 0` is still written there in the mixed observable pair
\((R_{\rm target},\epsilon_\eta)\).
The next clean move is therefore to use the later coherent-branch observable compiler to rewrite that gate directly in terms of the physical transfer-shape observable `\mathcal T^2`.

That step is useful because it exposes the packet factorization already present in the actual branch and turns the same-charge chain into a clean three-gate structure.

---

## 1. Exact coherent-branch observables and the transfer-shape identity

On the actual coherent branch the direct observables are
\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)},
\]
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\]
\[
R_{\rm target}
=
\Lambda_0\,\frac{\Omega_W^2(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]
with the exact selected-branch identity
\[
\boxed{
R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta).
}
\]
So the actual coherent branch itself supplies the finite physical packet
\[
(R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta),
\]
from which the quotient coordinates can be reconstructed exactly.

---

## 2. Exact finite same-charge packet in physical branch variables

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},\,\mathcal T_{\rm ref}^2,\,\epsilon_{\eta,\rm ref}),
\]
the carried finite packet is
\[
q_{\rm tr}=-C_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]
\[
q_{\rm nt}
=
B_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]
Using
\[
R_{\rm target}\,\mathcal T^2=\Lambda_0(1-\epsilon_\eta),
\qquad
R_{\rm target,ref}\,\mathcal T_{\rm ref}^2=\Lambda_0(1-\epsilon_{\eta,\rm ref}),
\]
one gets the exact finite factorization
\[
\boxed{
q_{\rm nt} + \frac{B_*}{C_*}q_{\rm tr}
=
\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right).
}
\]
So the finite nontracking packet is already the transfer-shape ratio up to the universal tracking feed-through.

On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\qquad
R_{\rm tr}=R_{\rm tr,ref},
\]
so the finite surviving packet becomes
\[
\boxed{
q_{\rm nt}=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
q_\eta=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]
Therefore the first static same-charge ceiling is exactly
\[
\boxed{
q_{\rm nt}=0
\iff
\mathcal T^2=\mathcal T_{\rm ref}^2.
}
\]
So the finite static-blind set is transfer-shape rigidity, not an additional hidden quotient condition.

---

## 3. Exact first-order physical drift compiler

Linearizing the coherent observables gives
\[
\delta\ln R_{\rm tr}
=
-
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\bigl[(1+\delta_U)\,d\ln\chi_0 + (1+\chi_0)\,d\ln\delta_U\bigr],
\]
\[
\delta\ln \mathcal T^2
=
 d\ln Z_W - d\ln \Omega_W^2
 + \frac{2\chi_0}{1+\chi_0}d\ln\chi_0
 + \frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]
\[
\delta\ln R_{\rm target}
=
 d\ln\Omega_W^2 - d\ln Z_W
 - \frac{2\chi_0}{1+\chi_0}d\ln\chi_0
 - \frac{2\epsilon}{1-\epsilon}d\ln\epsilon
 - \frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta,
\]
\[
\delta\ln\epsilon_\eta=d\ln\epsilon_\eta.
\]
So the first-order same-charge packet is
\[
\boxed{
q_{\rm tr}=-C_*\,\delta\ln R_{\rm tr},
\qquad
q_{\rm nt}+\frac{B_*}{C_*}q_{\rm tr}=\delta\ln\mathcal T^2,
\qquad
q_\eta=d\ln\epsilon_\eta.
}
\]
On the rigid-mouth branch,
\[
q_{\rm tr}=0,
\]
so the surviving first-order packet is simply
\[
\boxed{
q_{\rm nt}=\delta\ln\mathcal T^2,
\qquad
q_\eta=d\ln\epsilon_\eta.
}
\]
Thus the first-order static same-charge ceiling is
\[
\boxed{
q_{\rm nt}=0
\iff
\delta\ln\mathcal T^2=0.
}
\]
And once that ceiling has been cleared, the remaining post-static obstruction is exactly
\[
\boxed{
q_\eta=d\ln\epsilon_\eta.
}
\]

---

## 4. Exact support-blindness factorization of the physical packet

The coherent support-enhancement sector enters only through the baseline factor
\[
M_{\rm tr}=M_{\rm mix}
\left[1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}\right].
\]
But the direct same-charge observables satisfy
\[
\partial_\zeta R_{\rm tr}=0,
\qquad
\partial_\zeta \mathcal T^2=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\]
\[
\partial_{M_{\rm mix}} R_{\rm tr}=0,
\qquad
\partial_{M_{\rm mix}} \mathcal T^2=0,
\qquad
\partial_{M_{\rm mix}} \epsilon_\eta=0.
\]
Therefore the full finite and first-order packet is support-blind:
\[
\boxed{
\partial_\zeta q_{\rm tr}=\partial_\zeta q_{\rm nt}=\partial_\zeta q_\eta=0,
\qquad
\partial_{M_{\rm mix}} q_{\rm tr}=\partial_{M_{\rm mix}} q_{\rm nt}=\partial_{M_{\rm mix}} q_\eta=0.
}
\]
So support enhancement may rescue the steady normalization side of the branch, but it cannot change the actual same-charge packet at first weak-axisymmetric order.

---

## 5. Post-static dressing-invariance theorem on the actual branch

The same-charge chain on the actual coherent rigid-mouth branch now factors exactly into three gates.

### 5.1 Tracking gate
\[
q_{\rm tr}=0
\iff
R_{\rm tr}=R_{\rm tr,ref}
\iff
(1+\delta_U)d\ln\chi_0 + (1+\chi_0)d\ln\delta_U = 0
\quad\text{(at first order).}
\]

### 5.2 Static-blind transfer-shape gate
\[
q_{\rm nt}=0
\iff
\mathcal T^2 = \mathcal T_{\rm ref}^2
\iff
\delta\ln\mathcal T^2 = 0
\quad\text{(at first order).}
\]

### 5.3 Post-static dressing gate
After the first two have been cleared, the remaining obstruction is exactly
\[
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\qquad
q_\eta^{(1)}=d\ln\epsilon_\eta.
\]
So the final post-static criterion is
\[
\boxed{
q_\eta=0
\iff
\epsilon_\eta=\epsilon_{\eta,\rm ref}
\iff
 d\ln\epsilon_\eta=0
\quad\text{(at first order).}
}
\]
This is the exact post-static dressing-invariance theorem.

---

## 6. Best current summary after Stage 238

The continuation from Stage 237 is now complete.

1. The actual same-charge packet is charted directly by
   \[
   (R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta).
   \]
2. The finite static same-charge ceiling is exactly
   \[
   q_{\rm nt}=0 \iff \mathcal T^2 = \mathcal T_{\rm ref}^2.
   \]
3. The first-order static ceiling is exactly
   \[
   q_{\rm nt}=\delta\ln\mathcal T^2.
   \]
4. The whole direct packet is support-blind.
5. Therefore, once tracking rigidity and transfer-shape rigidity have both been imposed, the same-charge barrier reduces to one exact test only:

> compute `\epsilon_eta` on the actual rigid-mouth branch and check whether it is invariant.

That is the sharpest post-static physical-variable criterion reached so far.

---

## 7. SymPy-backed status

The accompanying SymPy audit verifies:

- the exact coherent-branch identity
  \[
  R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta),
  \]
- the finite packet factorization
  \[
  q_{\rm nt}+\frac{B_*}{C_*}q_{\rm tr}=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
  \]
- the rigid-mouth reduction
  \[
  q_{\rm nt}=\ln(\mathcal T^2/\mathcal T_{\rm ref}^2),
  \qquad
  q_\eta=\ln(\epsilon_\eta/\epsilon_{\eta,\rm ref}),
  \]
- the exact first-order physical drift compiler for `\delta\ln R_tr`, `\delta\ln \mathcal T^2`, and `\delta\ln R_target`,
- the first-order packet relation
  \[
  q_{\rm nt}+\frac{B_*}{C_*}q_{\rm tr}=\delta\ln\mathcal T^2,
  \qquad
  q_\eta=d\ln\epsilon_\eta,
  \]
- the support-blindness identities with respect to `\zeta` and `M_mix`,
- and the three-gate post-static dressing-invariance theorem on the actual rigid-mouth branch.

Supporting file:
- `moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
