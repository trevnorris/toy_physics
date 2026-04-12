# 5PN / Moving-Throat continuation — Stages 336–339

This pass takes the exact support-regime classifier from Stages 333–335 and ties it back into the physical selected-branch softening variable `xi` and the coherent-kernel placement map. That was the cleanest open step left on the support/source side.

The continuation point before this session was:

> express the actual support regime directly on the physical selected branch, instead of only in the abstract ratio `rho_alpha = Pi_tr / C_mix`.

The result is that the support/source side is now controlled by one exact scalar branch ratio
\[
\rho_\alpha^{(\mathrm{phys})}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{G_{\rm tr}}{M_{\rm mix}},
\]
and every regime boundary is an exact algebraic surface in the same selected-branch variable `xi`.

---

## Stage 336 — exact selected-branch loading ratio

On the tracking branch,
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{[9\delta+(9+2R^2)\xi]^2[9\delta+(9+2R)\xi]^2}
{81(1-\xi)\,[9\delta^2+18\delta\xi+(9+2R^2)\xi^2]^2},
\]
so the exact branch product is
\[
\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}.
\]

On the physical target branch, where
\[
R_{\rm target}=F_{\rm tr},
\qquad
C_{\rm mix}=R_{\rm target}M_{\rm mix},
\]
the loading ratio becomes exactly
\[
\boxed{
\rho_\alpha^{(\mathrm{phys})}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{G_{\rm tr}}{M_{\rm mix}}.
}
\]

So the support regime no longer needs the full product language once the branch is on-shell. It is just the ratio of the required tracking load to the mixed-only baseline.

The regime split is therefore

\[
G_{\rm tr}\le M_{\rm mix}
\quad\Longrightarrow\quad
\text{mixed-only enough},
\]

\[
M_{\rm mix}<G_{\rm tr}\le 2M_{\rm mix}
\quad\Longrightarrow\quad
\text{lowest symmetric twin enough},
\]

\[
G_{\rm tr}>2M_{\rm mix}
\quad\Longrightarrow\quad
\text{non-twin asymmetry required}.
\]

Because
\[
\frac{dG_{\rm tr}}{d\xi}>0,
\]
one also gets
\[
\frac{d\rho_\alpha^{(\mathrm{phys})}}{d\xi}>0.
\]
So the required support regime becomes strictly harder as the selected branch softens more deeply.

---

## Stage 337 — exact saturation depths `xi_(1x)` and `xi_(2x)`

The mixed-only threshold is the exact positive root of
\[
G_{\rm tr}=M_{\rm mix},
\]
namely
\[
\boxed{
\xi_{(1x)}
=
\frac{
M_{\rm mix}(9+2R^2)-9\delta
+
\sqrt{[M_{\rm mix}(9+2R^2)-9\delta]^2+324M_{\rm mix}\delta}
}{18}.
}
\]

The lowest-twin threshold is the exact positive root of
\[
G_{\rm tr}=2M_{\rm mix},
\]
namely
\[
\boxed{
\xi_{(2x)}
=
\frac{
2M_{\rm mix}(9+2R^2)-9\delta
+
\sqrt{[2M_{\rm mix}(9+2R^2)-9\delta]^2+648M_{\rm mix}\delta}
}{18}.
}
\]

Because `G_tr` is strictly increasing,
\[
\xi_{(1x)}<\xi_{(2x)}.
\]

So the support/source regime can now be read directly from the selected-branch depth:

- `xi_phys <= xi_(1x)` means mixed-only is enough,
- `xi_(1x) < xi_phys <= xi_(2x)` means the lowest symmetric twin is enough,
- `xi_phys > xi_(2x)` means non-twin asymmetry is required.

That is the cleanest direct classifier we have had so far on the support side.

---

## Stage 338 — coherent-kernel regime map

On the coherent local D/N branch,
\[
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)}.
\]

Substituting this into the selected-branch loading ratio gives
\[
\boxed{
\rho_\alpha^{(\mathrm{coh})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{8Z_W(1+\chi_0)^2}\,
G_{\rm tr}(\xi,\delta;R).
}
\]

So at fixed selected point, increasing `Z_W` lowers the required regime ratio and makes support success easier.

This gives two exact overlap thresholds:

### Mixed-only threshold
\[
\boxed{
Z_W^{(\mathrm{mix})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{8(1+\chi_0)^2}
\,G_{\rm tr}(\xi,\delta;R).
}
\]

### Lowest-twin threshold
\[
\boxed{
Z_W^{(\mathrm{twin})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{16(1+\chi_0)^2}
\,G_{\rm tr}(\xi,\delta;R).
}
\]

Equivalently, using
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\]
the exact radiative-demand thresholds are

\[
\boxed{
\Lambda_{\rm mix}
=
\frac{\pi^2}{8(1-\epsilon)}\,
\Pi_{\rm tr},
}
\qquad
\boxed{
\Lambda_{\rm twin}
=
\frac{\pi^2}{16(1-\epsilon)}\,
\Pi_{\rm tr}.
}
\]

So the support theorem is now directly tied to the same coherent-kernel variables that place the actual moving-throat branch.

---

## Stage 339 — exact non-twin asymmetry requirement

The support-ratio demand is still
\[
\zeta_{\rm req}
=
\frac{\rho_\alpha-1}{1-\epsilon(2-\rho_\alpha)}.
\]

The exact excess over the symmetric-twin value is
\[
\boxed{
\zeta_{\rm req}-1
=
\frac{(1-\epsilon)(\rho_\alpha-2)}{1+\epsilon(\rho_\alpha-2)}.
}
\]

So the three support regimes become completely explicit:

- `rho_alpha <= 1`: mixed-only enough,
- `1 < rho_alpha <= 2`: lowest symmetric twin enough, with `0 < zeta_req <= 1`,
- `rho_alpha > 2`: non-twin asymmetry required, with `zeta_req > 1`.

On the physical selected branch,
\[
\rho_\alpha^{(\mathrm{phys})}=\frac{G_{\rm tr}}{M_{\rm mix}},
\]
so the non-twin condition is exactly
\[
G_{\rm tr}(\xi,\delta;R) > 2 M_{\rm mix},
\]
equivalently
\[
\xi_{\rm phys} > \xi_{(2x)}.
\]

That means the first true asymmetry requirement is no longer vague. It is the exact amount by which the selected branch has crossed above the twin threshold.

---

## Net result after Stages 336–339

The support/source side is now effectively solved at the reduced-theorem level.

1. The exact regime classifier has been pulled back to the physical selected-branch variable `xi`.
2. The mixed-only and lowest-twin thresholds are exact algebraic depths `xi_(1x)` and `xi_(2x)`.
3. The coherent-kernel placement map rewrites the classifier directly in terms of the actual branch variables `(chi_0, eps_eta, eps, Z_W)`.
4. The non-twin asymmetry requirement is now an exact excess formula rather than a qualitative statement.

So the remaining reduced question is no longer “what support regime might the branch be in?”
It is narrower:

> once the actual PDE gives the physical branch point `(xi_phys, delta, R, chi_0, eps_eta, eps, Z_W)`, the support regime and any required asymmetry are read off immediately from these exact formulas.  
> The only thing still open is where the completed moving-throat branch actually lands.
