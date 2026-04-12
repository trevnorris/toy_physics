# 5PN continuation notes — stages 242 through 246

These stages pick up exactly from the Stage-241 equilibrium-alignment result and carry the support/source theorem through the first explicit parent confinement branch.

The structural change is that the support theorem is no longer phrased in terms of the abstract gain
`G_eq = g_phi^2 I_1 / K_X`
with free parent data `(g_phi, I_1, K_X)`.
It is now organized in five increasingly concrete layers:

1. the first explicit thin-wall confinement family,
2. the one-number wall figure of merit,
3. the explicit GNLS wall-shell reduction of `T_X`, `K_X`, and `kappa`,
4. the canonical `tanh` wall / local mouth closure branch,
5. the final explicit branch-placement problem in the three parent dimensionless controls
   `(chi_s, Lambda_ell, Upsilon_w)`.

So after Stage 246 the next theorem gate is no longer “compute `g_phi` somehow.”
It is:

> compute the actual moving-throat branch placement in `(chi_s, Lambda_ell, Upsilon_w)` and compare it directly to the exact support-threshold surfaces.

## Stage 242 — explicit thin-wall confinement branch

Take the first concrete parent wall family
\[
V_{\rm conf}(r;a) = V_0\, f\!\left(\frac{r-a}{\ell}\right),
\qquad
\xi := \frac{r-a}{\ell}.
\]

A support displacement `a -> a + phi(s)` gives
\[
\delta V_{\rm conf}
=
-\partial_a V_{\rm conf}\,\phi(s)
=
+\frac{V_0}{\ell} f'(\xi)\,\phi(s),
\]
so the parent support-loading amplitude is fixed exactly:
\[
g_\phi = \frac{V_0}{\ell}.
\]

Using the shell measure
\[
d^3y = 4\pi r^2 dr = 4\pi \ell (a+\ell \xi)^2 d\xi,
\]
the exact shell integral entering the equilibrium gain becomes
\[
I_1
=
4\pi \ell
\left[
a^2 J_1 + 2 a \ell J_2 + \ell^2 J_3
\right],
\]
where
\[
J_1 = \int d\xi \frac{f'(\xi)^2}{H(\xi)},
\qquad
J_2 = \int d\xi \frac{\xi f'(\xi)^2}{H(\xi)},
\qquad
J_3 = \int d\xi \frac{\xi^2 f'(\xi)^2}{H(\xi)}.
\]

For a centered symmetric wall layer, `J_2 = 0`, so
\[
I_1 = 4\pi \ell\left[a^2 J_1 + \ell^2 J_3\right].
\]

Therefore the exact equilibrium-aligned gain is
\[
G_{\rm eq}
=
\frac{g_\phi^2 I_1}{K_X}
=
\frac{4\pi V_0^2}{K_X}
\left[
\frac{a^2 J_1}{\ell}
+ 2 a J_2
+ \ell J_3
\right].
\]

On the centered branch, the thin-wall leading term is
\[
G_{\rm eq}^{(\rm tw)}
=
\frac{4\pi a^2 V_0^2 J_1}{K_X \ell}.
\]

So the parent support/source theorem is now a direct wall-amplitude theorem on the first explicit wall family.

## Stage 243 — wall figure of merit and direct fail/succeed thresholds

Introduce the support geometry parameter
\[
\kappa = \frac{K_X L^2}{T_X}.
\]
Then the natural dimensionless wall control variable is
\[
W_{\rm wall}
:=
\kappa G_{\rm eq}^{(\rm tw)}
=
\frac{4\pi a^2 L^2 J_1 V_0^2}{T_X \ell}.
\]

Using the operator-selected threshold pair
\[
G_{\rm fail} = \frac{Pe_{\rm req}}{\kappa \Delta_\infty(\kappa,\eta)},
\qquad
G_{\rm suff} = \frac{Pe_{\rm req}}{\kappa \Delta_0(\kappa,\eta)},
\]
the whole branch collapses to the exact support theorem
\[
W_{\rm wall} \le \frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)}
\quad\Rightarrow\quad
\text{fail},
\]
\[
W_{\rm wall} \ge \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)}
\quad\Rightarrow\quad
\text{succeed}.
\]

The physical wall-amplitude thresholds are therefore
\[
V_{0,\rm fail}^2
=
\frac{T_X \ell\, Pe_{\rm req}}{4\pi a^2 L^2 J_1 \Delta_\infty},
\qquad
V_{0,\rm suff}^2
=
\frac{T_X \ell\, Pe_{\rm req}}{4\pi a^2 L^2 J_1 \Delta_0}.
\]

So the support/source gate is no longer an abstract gain comparison. It is one dimensionless wall figure-of-merit comparison.

## Stage 244 — explicit GNLS wall-shell reduction and exact `W_wall = Xi`

Now derive the support coefficients directly from the parent GNLS energy on a thin active shell around the throat wall.

With the shell support mode
\[
\delta\rho(s,y)=q(s)\chi_\phi(y),
\]
the quadratic reduced support energy is
\[
E^{(2)}[q]
=
\frac12 \int_0^L ds \left[T_X q'(s)^2 + K_X q(s)^2\right],
\]
where
\[
T_X = \frac{\hbar^2}{4m\rho_w} N_{\phi\phi},
\qquad
K_X = H_w N_{\phi\phi} + \frac{\hbar^2}{4m\rho_w} G_{\phi\phi},
\]
and
\[
H_w = \frac{m c_{s,w}^2}{\rho_w}.
\]

For the explicit thin shell with `chi_phi = f'(\xi)`,
\[
N_{\phi\phi} = 4\pi a^2 \ell I_f,
\qquad
G_{\phi\phi} = \frac{4\pi a^2 I_g}{\ell},
\]
so
\[
T_X = \frac{\pi a^2 \ell I_f \hbar^2}{m\rho_w},
\]
\[
K_X = 4\pi a^2 \ell I_f H_w + \frac{\pi a^2 I_g \hbar^2}{m\rho_w \ell}.
\]

Therefore
\[
\kappa = \frac{K_X L^2}{T_X}
=
4\left(\frac{m c_{s,w} L}{\hbar}\right)^2
+
\frac{I_g}{I_f}\left(\frac{L}{\ell}\right)^2.
\]

The thin-wall matched-layer wall figure collapses exactly to
\[
W_{\rm wall}
=
\frac{4\rho_w^2 V_0^2 L^2}{\hbar^2 c_{s,w}^2 \ell^2}.
\]

And the Stage-41/42 fixed-point coupling is
\[
\Xi
=
\frac{g_\phi^2 I_1 L^2}{T_X}
=
W_{\rm wall}
\]
exactly on this explicit matched thin-wall branch.

So the wall figure of merit is not merely analogous to the support/source fixed-point coupling. It is the same object.

## Stage 245 — canonical `tanh` wall branch and natural local mouth closure

Choose the canonical smooth wall
\[
f(\xi) = \frac{1+\tanh \xi}{2},
\qquad
f'(\xi)=\frac12 \operatorname{sech}^2\xi,
\qquad
f''(\xi)=-\operatorname{sech}^2\xi \tanh\xi.
\]

Its exact shell moments are
\[
I_f = \int_{-\infty}^{+\infty} d\xi\, f'(\xi)^2 = \frac13,
\qquad
I_g = \int_{-\infty}^{+\infty} d\xi\, f''(\xi)^2 = \frac{4}{15},
\]
so
\[
\frac{I_g}{I_f} = \frac45.
\]

Therefore the explicit branch coefficients become
\[
T_X = \frac{\pi a^2 \ell \hbar^2}{3m\rho_w},
\]
\[
K_X = \frac{4\pi a^2}{15\ell m\rho_w}\left(5m^2 c_{s,w}^2 \ell^2 + \hbar^2\right),
\]
\[
J_1 = \frac{1}{3H_w},
\]
and
\[
\kappa
=
4\left(\frac{m c_{s,w} L}{\hbar}\right)^2
+
\frac45\left(\frac{L}{\ell}\right)^2.
\]

Using the natural local mouth closure
\[
K_m = \frac{T_X}{\ell},
\]
the Robin ratio is
\[
\eta = \frac{K_m L}{T_X} = \frac{L}{\ell}.
\]

So the first canonical explicit parent branch is controlled by the three dimensionless variables
\[
\chi_s := \frac{m c_{s,w} L}{\hbar},
\qquad
\Lambda_\ell := \frac{L}{\ell},
\qquad
\Upsilon_w := \frac{4\rho_w^2 V_0^2}{\hbar^2 c_{s,w}^2},
\]
with exact branch map
\[
\kappa = 4\chi_s^2 + \frac45 \Lambda_\ell^2,
\qquad
\eta = \Lambda_\ell,
\qquad
W_{\rm wall} = \Upsilon_w \Lambda_\ell^2.
\]

So the first explicit parent branch is no longer described by a long symbolic list. It is a three-parameter branch-placement problem.

## Stage 246 — explicit branch threshold surfaces

Insert the canonical branch formulas into the universal matched-branch theorem:
\[
W_{\rm wall} \le \frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)}
\quad\Rightarrow\quad
\text{fail},
\qquad
W_{\rm wall} \ge \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)}
\quad\Rightarrow\quad
\text{succeed}.
\]

With
\[
\kappa = 4\chi_s^2 + \frac45\Lambda_\ell^2,
\qquad
\eta = \Lambda_\ell,
\qquad
W_{\rm wall} = \Upsilon_w \Lambda_\ell^2,
\]
this becomes the exact explicit-branch surfaces
\[
\Upsilon_w \le \Upsilon_{\rm fail}(\chi_s,\Lambda_\ell)
\quad\Rightarrow\quad
\text{fail},
\]
\[
\Upsilon_w \ge \Upsilon_{\rm suff}(\chi_s,\Lambda_\ell)
\quad\Rightarrow\quad
\text{succeed},
\]
where
\[
\Upsilon_{\rm fail}
=
\frac{Pe_{\rm req}}{\Lambda_\ell^2\,
\Delta_\infty\!\left(4\chi_s^2+\frac45\Lambda_\ell^2,\Lambda_\ell\right)},
\]
\[
\Upsilon_{\rm suff}
=
\frac{Pe_{\rm req}}{\Lambda_\ell^2\,
\Delta_0\!\left(4\chi_s^2+\frac45\Lambda_\ell^2,\Lambda_\ell\right)}.
\]

Two asymptotic regimes are already exact:

### Shell-gradient dominated
If
\[
\frac45 \Lambda_\ell^2 \gg 4\chi_s^2,
\]
then
\[
\Upsilon_{\rm fail} \sim \frac{2 Pe_{\rm req}}{\sqrt5\,\Lambda_\ell},
\qquad
\Upsilon_{\rm suff} \to \frac45\left(1+\frac{2}{\sqrt5}\right) Pe_{\rm req}.
\]

### Compression dominated
If
\[
4\chi_s^2 \gg \frac45 \Lambda_\ell^2,
\]
then
\[
\Upsilon_{\rm fail} \sim \frac{2 Pe_{\rm req}\chi_s}{\Lambda_\ell^2},
\]
\[
\Upsilon_{\rm suff} \sim \frac{4 Pe_{\rm req}\chi_s^2(\Lambda_\ell+2\chi_s)}{\Lambda_\ell^3},
\]
and for `chi_s >> Lambda_ell`,
\[
\Upsilon_{\rm suff} \sim \frac{8 Pe_{\rm req}\chi_s^3}{\Lambda_\ell^3}.
\]

So compression-dominated branches are much harder to drive across the universal success threshold than shell-gradient dominated branches.

## Net result after Stage 246

The support/source theorem gap is now much smaller than it was after Stage 241.

1. The abstract loading amplitude `g_phi` has been replaced by the explicit thin-wall relation `g_phi = V0/ell`.
2. The abstract gain `G_eq` has collapsed to the wall figure of merit `W_wall`.
3. The matched thin-wall GNLS shell reduction fixes `T_X`, `K_X`, `kappa`, and shows `W_wall = Xi` exactly.
4. On the canonical `tanh` wall with the natural local mouth closure, the whole branch is controlled by the three parent dimensionless variables
   `(chi_s, Lambda_ell, Upsilon_w)`.
5. The remaining PDE-side theorem gate is now

> compute the actual moving-throat branch placement in `(chi_s, Lambda_ell, Upsilon_w)` and compare it directly to the explicit threshold surfaces `Upsilon_fail`, `Upsilon_suff`.

That is the cleanest continuation point after Stage 241.


## Stage 247 — Family-1 reference-branch geometry map

The next reduction is to stop treating the wall-geometry ratio `Lambda_ell = L/ell` as free and evaluate it on the first carried constructive throat branch.

Use the balanced thin-layer-consistent Family-1 wall width
\[
\epsilon_r = \frac{\ell}{a} = \frac1{20},
\]
together with the carried preferred throat aspect ratio
\[
\frac{L}{a} = \frac{37}{20}.
\]

Then the explicit branch fixes
\[
\Lambda_\ell = \frac{L}{\ell} = \frac{(L/a)}{(\ell/a)} = 37.
\]

On the natural local mouth closure from Stage 245,
\[
\eta = \frac{K_m L}{T_X} = \frac{L}{\ell} = \Lambda_\ell,
\]
so the same reference branch fixes
\[
\eta = 37.
\]

So the first explicit moving-throat support branch is no longer free in the Robin geometry variable. It already sits at one concrete large-`eta` point.

## Stage 248 — healing-length lock and exact support-scale reduction

Now identify the active shell width with the local GNLS healing width,
\[
\ell = \ell_h = \frac{\hbar}{2 m c_{s,w}}.
\]

Then the support-scale variable from Stage 245 becomes
\[
\chi_s = \frac{m c_{s,w} L}{\hbar} = \frac{L}{2\ell} = \frac{\Lambda_\ell}{2}.
\]

Since Stage 247 fixed `Lambda_ell = 37`, the same reference branch fixes
\[
\chi_s = \frac{37}{2} = 18.5.
\]

Therefore
\[
\kappa
=
4\chi_s^2 + \frac45 \Lambda_\ell^2
=
\frac95 \Lambda_\ell^2
=
\frac{12321}{5}
=
2464.2,
\]
and
\[
\alpha := \sqrt{\kappa}
=
\frac{111}{\sqrt5}
\approx 49.6407091.
\]

So after the healing lock the explicit branch has fixed
\[
\Lambda_\ell = 37,\qquad
\eta = 37,\qquad
\chi_s = \frac{37}{2},\qquad
\kappa = \frac{12321}{5}.
\]

The only remaining Stage-245 control is now the wall-loading amplitude `Upsilon_w`.

## Stage 249 — explicit Family-1 threshold window and reduction to one wall-depth datum

With
\[
\Lambda_\ell = 37,\qquad
\eta = 37,\qquad
\kappa = \frac{12321}{5},
\]
the exact operator threshold functions evaluate numerically to
\[
\Delta_0\!\left(\frac{12321}{5},37\right)
\approx
1.73302079021525\times 10^{-4},
\]
\[
\Delta_\infty\!\left(\frac{12321}{5},37\right)
\approx
2.01447565540522\times 10^{-2}.
\]

Therefore the explicit branch thresholds become
\[
\Upsilon_{\rm fail}
\approx
0.0362605617972939\,Pe_{\rm req},
\]
\[
\Upsilon_{\rm suff}
\approx
4.21495341569977\,Pe_{\rm req}.
\]

Equivalently, the fixed-point coupling window is
\[
\Xi_{\rm fail}
=
\frac{Pe_{\rm req}}{\Delta_\infty}
\approx
49.6407091004953\,Pe_{\rm req},
\]
\[
\Xi_{\rm suff}
=
\frac{Pe_{\rm req}}{\Delta_0}
\approx
5770.27122609299\,Pe_{\rm req}.
\]

Now write the Family-1 wall depth as
\[
V_0 = \alpha_r \mu_*,
\qquad
\alpha_r = 10,
\]
so
\[
\Upsilon_w = \alpha_r^2 \Theta_w = 100 \Theta_w,
\]
with the one remaining microscopic datum
\[
\Theta_w := \frac{4\rho_w^2 \mu_*^2}{\hbar^2 c_{s,w}^2}.
\]

The explicit threshold window is then
\[
\Theta_{\rm fail}
\approx
3.62605617972939\times 10^{-4}\,Pe_{\rm req},
\]
\[
\Theta_{\rm suff}
\approx
4.21495341569977\times 10^{-2}\,Pe_{\rm req}.
\]

So on the explicit Family-1 / healing-locked reference branch the reduced support/source placement problem is now almost finished:

- the geometry ratio is fixed,
- the support/healing ratio is fixed,
- the operator scales are fixed,
- and the only remaining explicit branch datum is the wall-depth amplitude `Theta_w`.

## Net result after Stage 249

The continuation past Stage 241 is now substantially sharper.

1. The abstract parent loading amplitude `g_phi` is replaced by the explicit thin-wall law `g_phi = V0/ell`.
2. The support theorem collapses to the dimensionless wall figure of merit `W_wall`.
3. The matched thin-wall GNLS shell reduction fixes `T_X`, `K_X`, `kappa`, and proves `W_wall = Xi`.
4. On the canonical `tanh` wall with the natural local mouth closure, the whole branch reduces to the three control variables
   \[
   (\chi_s,\Lambda_\ell,\Upsilon_w).
   \]
5. On the explicit Family-1 / healing-locked reference branch, even those collapse further:
   \[
   \Lambda_\ell = 37,\quad
   \eta = 37,\quad
   \chi_s = 37/2,\quad
   \kappa = 12321/5,
   \]
   so the only remaining branch datum is the wall-depth amplitude `Theta_w`.

That means the next honest theorem gate is no longer to derive more reduced support/source algebra.

It is:

> extract the actual Family-1 wall-depth datum `Theta_w` (or equivalently `Upsilon_w`) from the real moving-throat branch and compare it directly to the explicit threshold window above.
