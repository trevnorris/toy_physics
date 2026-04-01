# 4PN Hamiltonian-chart generic-frame lift notes

This note freezes the next exact local 4PN step after the Hamiltonian-to-ordinary translation audit.

## 1. What was actually tested

The earlier local-scaffold audit showed that the raw constant-coefficient exchange-symmetric generic-frame local 4PN interaction scaffold is already **slot-complete** after COM projection:

- \(G/r\): 5 COM interaction slots,
- \(G^2/r^2\): 4 COM interaction slots,
- \(G^3/r^3\): 3 COM interaction slots,
- \(G^4/r^4\): 2 COM interaction slots,
- \(G^5/r^5\): 1 COM interaction slot.

But the later Hamiltonian-to-ordinary audit showed that the **ordinary** target develops extra \(\nu^4\) tails in the \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks after quartic translation, whereas the fixed local **Hamiltonian** target matches the scaffold degree ceilings exactly.

So the next exact question was:

> Does the constant-coefficient generic-frame scaffold span the full fixed-chart **Hamiltonian** local 4PN comparable-mass residual, or is there still a seed/obstruction problem hiding there?

## 2. Natural Hamiltonian seed

Using the exact one-body 4PN Hamiltonian gate
\[
\frac{7}{256}p^{10},
\qquad
\frac{45}{128}\frac{p^8}{r},
\qquad
\frac{13}{8}\frac{p^6}{r^2},
\qquad
\frac{105}{32}\frac{p^4}{r^3},
\qquad
\frac{105}{32}\frac{p^2}{r^4},
\qquad
-\frac{1}{16}\frac{1}{r^5},
\]
the natural reduced COM Hamiltonian self/static seed is

\[
K^{\rm seed}=\frac{7}{256}(X_A^9+X_B^9),
\]
\[
Q_1^{\rm seed}=\frac{45}{128}(X_A^8+X_B^8),
\qquad
T_1^{\rm seed}=\frac{13}{8}(X_A^7+X_B^7),
\qquad
S_1^{\rm seed}=\frac{105}{32}(X_A^6+X_B^6),
\]
\[
U_1^{\rm seed}=\frac{105}{32}(X_A^5+X_B^5),
\qquad
W_1^{\rm seed}=-\frac{1}{16}(X_A^4+X_B^4),
\]
with \(X_A=(1+\Delta)/2\), \(X_B=(1-\Delta)/2\), \(X_A X_B=\nu\), and all subleading slots initially zero.

A first exact outcome is that this seed already captures the **entire free 4PN COM Hamiltonian slot**:
\[
\Delta K = 0.
\]

So the generic-frame local lift only has to solve the **interaction** residual.

## 3. Polynomial coefficient-space formulation

To test the lift exactly, each COM block was decomposed into explicit \(\nu\)-coefficient space:

- \(G/r\): \(5\) slots, each with \(\nu,\nu^2,\nu^3,\nu^4\) coefficients \(\Rightarrow 20\) coefficients,
- \(G^2/r^2\): \(4\) slots, each with \(\nu,\nu^2,\nu^3\) coefficients \(\Rightarrow 12\),
- \(G^3/r^3\): \(3\) slots, each with \(\nu,\nu^2,\nu^3\) coefficients \(\Rightarrow 9\),
- \(G^4/r^4\): \(2\) slots, each with \(\nu,\nu^2\) coefficients \(\Rightarrow 4\),
- \(G^5/r^5\): \(1\) slot, with \(\nu,\nu^2\) coefficients \(\Rightarrow 2\).

This is stricter than just counting COM interaction slots. It asks whether the scaffold spans the full **polynomial data** of the imported local Hamiltonian target.

## 4. Exact blockwise ranks

The constant-coefficient generic-frame Hamiltonian scaffold achieves full rank in every block:

\[
G/r:\quad 52 \to 20 \quad\text{with nullity }32,
\]
\[
G^2/r^2:\quad 46 \to 12 \quad\text{with nullity }34,
\]
\[
G^3/r^3:\quad 29 \to 9 \quad\text{with nullity }20,
\]
\[
G^4/r^4:\quad 10 \to 4 \quad\text{with nullity }6,
\]
\[
G^5/r^5:\quad 2 \to 2 \quad\text{with nullity }0.
\]

So the total interaction coefficient-space dimensions are

- rows: \(20+12+9+4+2=47\),
- columns: \(52+46+29+10+2=139\),
- total rank: \(47\),
- total nullity: \(92\).

The important point is not the large nullity — that is expected in a generic-frame scaffold — but the fact that the **rank is maximal in every block**.

## 5. Main theorem outcome

There is **no Hamiltonian-chart obstruction** analogous to the ordinary-chart \(\nu^4\)-tail problem.

More concretely:

1. the free 4PN COM Hamiltonian slot is already absorbed by the natural Hamiltonian seed;
2. every local interaction residual slot of the imported fixed-chart COM Hamiltonian target lies in the image of the constant-coefficient exchange-symmetric generic-frame scaffold;
3. no extra seed refinement is needed before the Hamiltonian-chart generic-frame local lift;
4. one exact generic-frame representative exists blockwise and reduces back to the imported COM Hamiltonian residual exactly.

So the correct local-first program is now sharply confirmed:

\[
\boxed{
\text{fixed local 4PN Hamiltonian target}
\;\longrightarrow\;
\text{generic-frame Hamiltonian lift}
\;\longrightarrow\;
\text{ordinary translation only afterward}.
}
\]

## 6. Contrast with the ordinary-chart translation audit

This is the clean structural comparison:

- In the **Hamiltonian** chart, the imported local target has blockwise degree ceilings
  \[
  G/r:4,\qquad G^2/r^2:3,\qquad G^3/r^3:3,\qquad G^4/r^4:2,\qquad G^5/r^5:2,
  \]
  exactly matching the constant-coefficient scaffold.

- In the translated **ordinary** chart, the target becomes
  \[
  G/r:4,\qquad G^2/r^2:4,\qquad G^3/r^3:4,\qquad G^4/r^4:4,\qquad G^5/r^5:2,
  \]
  which is why the ordinary-chart lift looked obstructed.

So the previous “Hamiltonian first” recommendation is no longer just a heuristic. It is now backed by an explicit exact rank test.

## 7. What remains open after this step

The local Hamiltonian-chart lift is now **structurally open but not obstructed**. The remaining tasks are:

1. fix a canonical generic-frame representative inside the \(92\)-direction null family,
2. import or reconstruct the full fixed-chart **generic-frame** local 4PN Hamiltonian target if we want uniqueness in a chosen chart rather than mere existence,
3. then translate that representative back through the quartic compiler to the ordinary chart,
4. while keeping the tail/hereditary quadrupole bridge separate.

So the local 4PN program is now in a better position than the previous ordinary-chart result suggested: the real local obstruction is no longer span, but only chart-fixing / canonical-choice bookkeeping.
