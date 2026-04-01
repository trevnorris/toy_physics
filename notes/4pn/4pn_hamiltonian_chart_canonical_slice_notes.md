
# 4PN Hamiltonian-chart canonical-slice notes

This note freezes the next exact local 4PN step after the Hamiltonian-chart generic-frame lift audit.

## 1. What was actually fixed

The previous lift audit proved **existence**: every local interaction block of the fixed-chart COM 4PN Hamiltonian residual lies in the image of the constant-coefficient exchange-symmetric generic-frame scaffold, with total nullity \(92\).

The unresolved local question was therefore no longer span. It was:

> What is a clean canonical generic-frame Hamiltonian representative inside that 92-direction COM-blind family?

The present audit answers that.

## 2. The COM-null ideals at 4PN

The generic-frame invariants satisfy the same COM identities as at 3PN:

\[
\mathcal C_1 = p a + q c,\qquad
\mathcal C_2 = p c + q b,\qquad
\mathcal C_3 = p d + q e,
\]
\[
\mathcal C_4 = ab-c^2,\qquad
\mathcal C_5 = ae-cd,\qquad
\mathcal C_6 = bd-ce.
\]

The 4PN local Hamiltonian-chart nullities are

\[
Q:32,\qquad T:34,\qquad S:20,\qquad U:6,\qquad W:0.
\]

The exact ideal-membership result is:

- all \(32\) \(Q\)-block null vectors lie in the **full** COM ideal
  \[
  \langle \mathcal C_1,\mathcal C_2,\mathcal C_3,\mathcal C_4,\mathcal C_5,\mathcal C_6\rangle;
  \]
- all \(34\) \(T\)-block null vectors lie already in the **linear-momentum** ideal
  \[
  \langle \mathcal C_1,\mathcal C_2,\mathcal C_3\rangle;
  \]
- all \(20\) \(S\)-block null vectors lie already in the same linear-momentum ideal;
- all \(6\) \(U\)-block null vectors also lie in that linear-momentum ideal;
- the \(W\)-block has **zero nullity**.

So the whole 92-direction ambiguity is purely COM-blind algebraic freedom, and the top static \(G^5/r^5\) block is already fully fixed once the COM target is fixed.

## 3. Canonical quotient slice

The canonical slice used in the audit is the exact Gauss–Jordan lift with **all null coordinates set to zero**.

This produces a sparse explicit generic-frame Hamiltonian representative with

\[
Q:11,\qquad T:12,\qquad S:9,\qquad U:4,\qquad W:2
\]
nonzero scaffold directions.

## 4. Explicit canonical representative

Write the local comparable-mass Hamiltonian residual as

\[
\Delta H_{4,\mathrm{loc}}^{\mathrm{can}}
=
\frac{Gpq}{r}Q_{\mathrm{can}}
+
\frac{G^2pq}{r^2}T_{\mathrm{can}}
+
\frac{G^3pq}{r^3}S_{\mathrm{can}}
+
\frac{G^4}{r^4}U_{\mathrm{can}}
+
\frac{G^5}{r^5}W_{\mathrm{can}}.
\]

Then one exact canonical representative is:

### \(G/r\) block
\[
\begin{aligned}
Q_{\mathrm{can}}={}&
-\frac{11}{64}a^2b^2
+\frac{5}{256}(a^2bc+ab^2c)
-\frac{3}{32}(a^2bd^2+ab^2e^2)
+\frac{1}{64}(a^2bde+ab^2de)
\\
&-\frac{27}{64}(a^2c^2+b^2c^2)
-\frac{9}{64}(a^2d^2e^2+b^2d^2e^2)
+\frac{3}{128}(a^2de^3+b^2d^3e)
+\frac{3}{64}(a^2e^4+b^2d^4)
\\
&-\frac{5}{32}(ad^2e^4+bd^4e^2)
+\frac{5}{64}(ad^3e^3+bd^3e^3)
-\frac{35}{256}(d^5e^3+d^3e^5).
\end{aligned}
\]

### \(G^2/r^2\) block
\[
\begin{aligned}
T_{\mathrm{can}}={}&
-\frac{87}{128}(a^2bp+ab^2q)
-\frac{2227}{256}(a^2bq+ab^2p)
+\frac{63}{64}(a^2cq+b^2cp)
\\
&+\frac{49}{16}(a^2d^2p+b^2e^2q)
-\frac{435}{64}(a^2dep+b^2deq)
+\frac{2435}{256}(a^2e^2p+b^2d^2q)
\\
&-\frac{183}{16}(ad^2e^2p+bd^2e^2q)
-\frac{8305}{768}(ad^2e^2q+bd^2e^2p)
+\frac{889}{192}(ad^3eq+bde^3p)
\\
&-\frac{5343}{1280}(d^3e^3p+d^3e^3q)
+\frac{325}{128}(d^4e^2q+d^2e^4p)
-\frac{369}{160}(d^5eq+de^5p).
\end{aligned}
\]

### \(G^3/r^3\) block
\[
\begin{aligned}
S_{\mathrm{can}}={}&
\left(-\frac{3307423}{28800}+\frac{40483\pi^2}{16384}\right)(a^2p^2+b^2q^2)
+\left(-\frac{211189}{19200}+\frac{2749\pi^2}{8192}\right)(a^2pq+b^2pq)
\\
&+\left(-\frac{184897}{900}+\frac{34985\pi^2}{8192}\right)abpq
+\left(\frac{139241}{1200}-\frac{12507\pi^2}{2048}\right)(ad^2p^2+be^2q^2)
\\
&+\left(\frac{63347}{1600}-\frac{1059\pi^2}{1024}\right)(ad^2pq+be^2pq)
+\left(-\frac{716971}{9600}+\frac{10389\pi^2}{2048}\right)(adep^2+bdeq^2)
\\
&+\left(-\frac{16727}{384}-\frac{35655\pi^2}{16384}\right)(d^2e^2p^2+d^2e^2q^2)
+\left(-\frac{51193}{960}-\frac{36405\pi^2}{8192}\right)d^2e^2pq
\\
&+\left(\frac{23533}{1280}-\frac{375\pi^2}{8192}\right)(d^3eq^2+de^3p^2).
\end{aligned}
\]

### \(G^4/r^4\) block
\[
\begin{aligned}
U_{\mathrm{can}}={}&
\left(\frac{930047}{9600}-\frac{551243\pi^2}{49152}\right)(ap^2q+bpq^2)
+\left(\frac{500761}{19200}-\frac{21837\pi^2}{8192}\right)(apq^2+bp^2q)
\\
&+\left(\frac{274387}{1600}-\frac{62047\pi^2}{49152}\right)(d^2p^2q+e^2pq^2)
+\left(\frac{3401779}{57600}-\frac{28691\pi^2}{24576}\right)(d^2pq^2+e^2p^2q).
\end{aligned}
\]

### \(G^5/r^5\) block
\[
W_{\mathrm{can}}=
\left(-\frac{169799}{2400}+\frac{6237\pi^2}{1024}\right)(p^3q+pq^3)
+\left(-\frac{609427}{3600}+\frac{44825\pi^2}{3072}\right)p^2q^2.
\]

## 5. Exact verification

The audit verifies directly that this canonical representative reduces back to the imported fixed-chart local COM Hamiltonian residual in every slot:

- \(Q\): all 5 slots match exactly,
- \(T\): all 4 slots match exactly,
- \(S\): all 3 slots match exactly,
- \(U\): both slots match exactly,
- \(W\): the static slot matches exactly.

So this is not just one heuristic sparse guess. It is an exact canonical quotient-slice representative.

## 6. What this does and does not settle

This step **does** settle the local algebraic existence problem at the generic-frame Hamiltonian level:

- the 92-direction ambiguity is purely COM-blind;
- the relevant null ideals are now identified;
- and one exact sparse canonical representative is frozen.

This step does **not** yet prove fixed-chart uniqueness against the true full generic-frame ADM Hamiltonian target. For that we would still need either:

1. the fixed-chart generic-frame local 4PN Hamiltonian target itself, or
2. an exact generic-frame Hamiltonian compiler uniqueness theorem strong enough to eliminate every remaining chart ambiguity directly.

So the local-first program is now:

\[
\boxed{
\text{canonical Hamiltonian generic-frame representative}
\;\longrightarrow\;
\text{generic-frame quartic translation to the ordinary chart}
\;\longrightarrow\;
\text{tail bridge kept separate}.
}
\]

That is the next clean move.
