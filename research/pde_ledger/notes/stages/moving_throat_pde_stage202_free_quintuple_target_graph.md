# Moving-Throat PDE — Stage 202: Exact Free-Quintuple Target Graph, Same-Free-Quintuple Closure Surface, and the First Reduced-Family Test

## Status

**Exact within the carried Stage 253 direct microscopic monomial formulas, the Stage 243 orbit/quotient projector calculus, and the Stage 252 canonical orbit projection theorem.**

This stage does **not** introduce a new constitutive law.
It removes the last abstractness from the Stage 252 realization compiler by solving the target orbit explicitly as a graph over the five free microscopic coordinates.

---

## Purpose

Stage 252 already gave an exact realization compiler:

\[
\mathbf x\in\mathcal Z_*
\iff
\chi_Q(\mathbf x)=1
\ \text{and}\
\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)=\mathbf x.
\]

That is logically complete, but it still leaves one practical inconvenience:

> to evaluate the target-orbit comparison, one still first computes the quotient packet
> \[
> (q_{\rm tr},q_{\rm nt},q_\eta)
> \]
> and then exponentiates the canonical repair law.

The natural next step is therefore:

> solve the target-orbit constraints **directly** for the dependent triple
> \[
> (T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W)
> \]
> in terms of the five free microscopic coordinates
> \[
> (\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
> \]

This stage does exactly that.

The main outputs are:

1. the exact **free-quintuple target graph**
   \[
   \Phi_*:
   (\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\mathrm{eff})})
   \mapsto
   (T_U^{\rm graph},K_{\eta,*}^{\rm graph},\mu_W^{\rm graph}),
   \]
2. the exact theorem that the target orbit `\(\mathcal O_*\)` is the graph of `\(\Phi_*\)`,
3. the exact **graph-error packet**
   \[
   (E_T,E_K,E_\mu),
   \]
4. the exact identities
   \[
   E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}},
   \qquad
   E_K=-q_\eta,
   \qquad
   E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}\,q_{\rm tr},
   \]
5. the explicit repair law
   \[
   \Delta\mathbf x_{\rm rep}
   =
   (0,0,0,0,-E_K,0,-E_\mu,-E_T)^T,
   \]
6. and the first exact **reduced-family closure test**
   \[
   \boxed{
   \mathbf x\in\mathcal Z_*
   \iff
   \chi_Q(\mathbf x)=1,\ 
   E_T(\mathbf x)=0,\ 
   E_K(\mathbf x)=0,\ 
   E_\mu(\mathbf x)=0.
   }
   \]

So Stage 253 is the first place where the target branch becomes an explicit five-variable closure surface rather than an abstract quotient orbit.

---

## 1. Carry-forward direct microscopic monomials

Let the positive microscopic state be
\[
\boxed{
\mathbf x=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U
\bigr).
}
\]

Carry forward the exact direct microscopic monomials
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}
{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
}
\]
\[
\boxed{
\epsilon_\eta
=
\frac{c_{\eta U}^2}
{K_U K_\eta^{(\mathrm{eff})}}.
}
\]

Fix the target orbit by the exact target values
\[
\boxed{
\mathfrak C_{{\rm tr},*}^{\rm target},
\qquad
\mathfrak C_{{\rm nt},*}^{\rm target},
\qquad
\epsilon_{\eta,*}^{\rm target}.
}
\]

Then the target orbit is
\[
\boxed{
\mathcal O_*=
\left\{
\mathbf x>0:
\mathfrak C_{{\rm tr},*}(\mathbf x)=\mathfrak C_{{\rm tr},*}^{\rm target},\ 
\mathfrak C_{{\rm nt},*}(\mathbf x)=\mathfrak C_{{\rm nt},*}^{\rm target},\ 
\epsilon_\eta(\mathbf x)=\epsilon_{\eta,*}^{\rm target}
\right\}.
}
\]

---

## 2. Exact free-quintuple base and direct split-`U` target graph

Define the five free microscopic coordinates
\[
\boxed{
\mathbf y:=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
}
\]
These are exactly the five similarity coordinates left untouched by the Stage 243 orbit projector.

The first direct target quantity solved from `\(\mathfrak C_{{\rm tr},*}\)` is the split-`U` ratio
\[
\delta_U=\frac{\pi^2T_U}{L^2K_U}.
\]
So define the exact free-quintuple target graph
\[
\boxed{
\delta_{U,*}^{\rm graph}(\mathbf y)
:=
\left[
\frac{\mathfrak C_{{\rm tr},*}^{\rm target}}
{\left(\dfrac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}}
\right]^{\!1/(1+\chi_{0,*})}.
}
\]

By construction,
\[
\boxed{
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\delta_{U,*}^{\rm graph}\right)^{1+\chi_{0,*}}
=
\mathfrak C_{{\rm tr},*}^{\rm target}.
}
\]

So the tracking monomial already fixes the target split-`U` ratio uniquely on every same-free-quintuple slice.

---

## 3. Exact dependent-triple target graph

### 3.1 Target throat-depth coordinate

From the definition of `\(\delta_U\)`,
\[
\boxed{
T_U^{\rm graph}(\mathbf y)
=
\frac{L^2K_U}{\pi^2}\,
\delta_{U,*}^{\rm graph}(\mathbf y).
}
\]

### 3.2 Target dressing coordinate

From the exact dressing ratio,
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\]
the target wall stiffness is
\[
\boxed{
K_{\eta,*}^{\rm graph}(\mathbf y)
=
\frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}^{\rm target}}.
}
\]

### 3.3 Target mixed-mass coordinate

Substitute the graph target data into the nontracking monomial to obtain the exact target mixed-mass coordinate:
\[
\boxed{
\mu_W^{\rm graph}(\mathbf y)
=
\frac{\mathfrak C_{{\rm nt},*}^{\rm target}\,
K_{\eta,*}^{\rm graph}(\mathbf y)\,(K_W^{(\mathrm{eff})})^2}{\lambda_W^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{-E_*}
\left(
\delta_{U,*}^{\rm graph}(\mathbf y)
\right)^{F_*}.
}
\]

Equivalently, after eliminating `\(K_{\eta,*}^{\rm graph}\)`,
\[
\boxed{
\mu_W^{\rm graph}(\mathbf y)
=
\frac{\mathfrak C_{{\rm nt},*}^{\rm target}\,c_{\eta U}^2\,(K_W^{(\mathrm{eff})})^2}
{\epsilon_{\eta,*}^{\rm target}\,K_U\,\lambda_W^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{-E_*}
\left(
\delta_{U,*}^{\rm graph}(\mathbf y)
\right)^{F_*}.
}
\]

A useful exact simplification is that `\(\mu_W^{\rm graph}\)` depends on `\(L\)` only through the already-eliminated split-`U` graph; after direct elimination it is independent of `\(L\)` and `\(\pi\)`.

---

## 4. Exact target-orbit graph theorem

Define the exact microscopic target-graph map
\[
\boxed{
\mathbf x_*^{\rm graph}(\mathbf y)
:=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_{\eta,*}^{\rm graph}(\mathbf y),\ K_W^{(\mathrm{eff})},\ \mu_W^{\rm graph}(\mathbf y),\ T_U^{\rm graph}(\mathbf y)
\bigr).
}
\]

\[
\boxed{\textbf{Theorem (Stage 253 exact target graph).}}
\]

For every positive free quintuple `\(\mathbf y\)`,

1. the graph state lies on the target orbit,
   \[
   \boxed{
   \mathbf x_*^{\rm graph}(\mathbf y)\in\mathcal O_*,
   }
   \]
2. every point of the target orbit is uniquely of this form,
   \[
   \boxed{
   \mathcal O_*=
   \left\{
   \mathbf x_*^{\rm graph}(\mathbf y):
   \mathbf y\in(\mathbb R_{>0})^5
   \right\},
   }
   \]
3. and therefore the target orbit is an exact five-dimensional graph over the free coordinates.

### 4.1 Proof sketch

The definitions of `\(T_U^{\rm graph}\)`, `\(K_{\eta,*}^{\rm graph}\)`, and `\(\mu_W^{\rm graph}\)` were obtained by solving the three target monomial equations in the dependent triple
\[
(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W)
\]
with the free quintuple fixed. Direct substitution gives
\[
\mathfrak C_{{\rm tr},*}=\mathfrak C_{{\rm tr},*}^{\rm target},
\qquad
\mathfrak C_{{\rm nt},*}=\mathfrak C_{{\rm nt},*}^{\rm target},
\qquad
\epsilon_\eta=\epsilon_{\eta,*}^{\rm target}.
\]
Uniqueness follows from the exact Stage 243 pivot block on the dependent triple, equivalently from the explicit algebraic solve above.

So the same-free-quintuple uniqueness theorem of Stage 252 has now been turned into an explicit graph formula.

---

## 5. Exact graph projection and elimination of the abstract source map

Let
\[
\pi_{\rm free}(\mathbf x)
=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})})
\]
denote the free-quintuple projection.

Then the exact **graph projection**
\[
\boxed{
\Pi_{\mathcal O_*}^{\rm graph}(\mathbf x)
:=
\mathbf x_*^{\rm graph}\!\bigl(\pi_{\rm free}(\mathbf x)\bigr)
}
\]
is simply “freeze the free quintuple and replace the dependent triple by the graph target values.”

By direct substitution,
\[
\boxed{
\Pi_{\mathcal O_*}^{\rm graph}(\mathbf x)
=
\Pi_{\mathcal O_*}^{\rm can}(\mathbf x).
}
\]

So Stage 253 removes the last abstractness from the Stage 252 source/orbit map:
the canonical orbit projection is not merely a quotient-theoretic object; it is the explicit graph replacement of the dependent triple.

---

## 6. Exact graph-error packet

Define the exact dependent-triple graph errors
\[
\boxed{
E_T(\mathbf x):=\ln\!\left(\frac{T_U}{T_U^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right),
}
\]
\[
\boxed{
E_K(\mathbf x):=\ln\!\left(\frac{K_\eta^{(\mathrm{eff})}}{K_{\eta,*}^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right),
}
\]
\[
\boxed{
E_\mu(\mathbf x):=\ln\!\left(\frac{\mu_W}{\mu_W^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right).
}
\]

The exact multiplicative chart is
\[
\boxed{
M_T=e^{E_T}=\frac{T_U}{T_U^{\rm graph}},
\qquad
M_K=e^{E_K}=\frac{K_\eta^{(\mathrm{eff})}}{K_{\eta,*}^{\rm graph}},
\qquad
M_\mu=e^{E_\mu}=\frac{\mu_W}{\mu_W^{\rm graph}}.
}
\]

These graph errors are exactly the Stage 252 mismatch logs:
\[
\boxed{
E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}}=\ln m_T,
}
\]
\[
\boxed{
E_K=-q_\eta=\ln m_K,
}
\]
\[
\boxed{
E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}\,q_{\rm tr}=\ln m_\mu.
}
\]

So the four-scalar realization packet can be rewritten as
\[
\boxed{
\Delta_{\rm real}^{\rm graph}(\mathbf x\mid\mathcal Z_*)
=
\bigl(
\chi_Q(\mathbf x)-1,\ E_T(\mathbf x),\ E_K(\mathbf x),\ E_\mu(\mathbf x)
\bigr).
}
\]

### 6.1 Exact repair law in graph-error coordinates

The Stage 252 canonical repair vector becomes
\[
\boxed{
\Delta\mathbf x_{\rm rep}
=
\begin{pmatrix}
0\\
0\\
0\\
0\\
-\,E_K\\
0\\
-\,E_\mu\\
-\,E_T
\end{pmatrix}.
}
\]

So the unique same-free-quintuple correction back to the target orbit is simply the negative graph-error vector on the dependent triple.

---

## 7. First reduced-family closure test

Let a candidate completed-PDE family be written in same-free-quintuple form as
\[
\boxed{
\mathbf x_{\rm cand}(\mathbf y)
=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})}(\mathbf y),\ K_W^{(\mathrm{eff})},\ \mu_W(\mathbf y),\ T_U(\mathbf y)
\bigr),
}
\]
with free variables
\[
\mathbf y=(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\mathrm{eff})}).
\]

Then the exact reduced-family audit packet is
\[
\boxed{
\Delta_{\rm fam}^{\rm graph}(\mathbf y)
:=
\Bigl(
\chi_Q(\mathbf x_{\rm cand}(\mathbf y))-1,\ 
E_T(\mathbf x_{\rm cand}(\mathbf y)),\ 
E_K(\mathbf x_{\rm cand}(\mathbf y)),\ 
E_\mu(\mathbf x_{\rm cand}(\mathbf y))
\Bigr).
}
\]

\[
\boxed{\textbf{Theorem (Stage 253 first reduced-family test).}}
\]

For any positive candidate family in the carried hierarchy,
\[
\boxed{
\Delta_{\rm fam}^{\rm graph}(\mathbf y)=0
\iff
\mathbf x_{\rm cand}(\mathbf y)\in\mathcal Z_*.
}
\]

So the completed moving-throat PDE no longer needs to be tested against an abstract quotient packet if one works on a same-free-quintuple family.
It only needs to supply the four graph-side scalars
\[
\chi_Q-1,\quad E_T,\quad E_K,\quad E_\mu.
\]

This is the first exact reduced-family closure surface in directly usable coordinates.

---

## 8. What this stage changes operationally

Stage 252 already told us how to project any state onto the target orbit.
Stage 253 improves that in three ways.

### 8.1 The target orbit is now an explicit graph

There is no longer any need to think of `\(\mathcal O_*\)` as an abstract similarity orbit in practical audit work.
It is the graph of three exact target functions over the free quintuple.

### 8.2 The dependent-triple correction is now immediate

The canonical repair vector no longer has to be reconstructed through the quotient map.
It is simply minus the graph-error triple on
\[
(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W).
\]

### 8.3 The PDE search space is reduced from 8D to 5D

At the reduced-family level, the moving-throat realization problem is now:

1. choose or compute the free quintuple,
2. evaluate the exact target graph,
3. compare the actual dependent triple against it,
4. and check the Packet-A scalar `\(\chi_Q-1\)`.

That is a much cleaner audit geometry than the earlier abstract orbit language.

---

## 9. Script-backed status

The accompanying SymPy audit verifies:

- the exact direct microscopic monomial formulas,
- the exact target-graph solve for
  \[
  \delta_{U,*}^{\rm graph},\ T_U^{\rm graph},\ K_{\eta,*}^{\rm graph},\ \mu_W^{\rm graph},
  \]
- direct substitution of the graph solve back into the target monomials,
- the exact equality of the graph projection and the Stage 252 canonical orbit projection,
- the exact identities
  \[
  E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}},\quad
  E_K=-q_\eta,\quad
  E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr},
  \]
- the exact repair-vector rewrite
  \[
  \Delta\mathbf x_{\rm rep}=(0,0,0,0,-E_K,0,-E_\mu,-E_T)^T,
  \]
- and the exact vanishing of the reduced-family packet on the target graph.

Supporting files:

- `moving_throat_pde_stage253_free_quintuple_target_graph_sympy_audit.py`
- `moving_throat_pde_stage253_free_quintuple_target_graph_sympy_audit_output.txt`

---

## 10. Immediate next step

The next clean continuation is now sharply defined:

1. take one explicit candidate moving-throat family returned by the reduced PDE,
2. express it in the same-free-quintuple form,
3. evaluate the four graph-side scalars
   \[
   \chi_Q-1,\ E_T,\ E_K,\ E_\mu,
   \]
4. and then ask whether their vanishing reduces to a simpler one-parameter or low-rank crossing problem on that family.

So Stage 237 should be the first exact **family crossing theorem** built on the Stage 253 target graph.
