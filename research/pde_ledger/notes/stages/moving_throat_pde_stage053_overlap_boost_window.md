# Moving-Throat PDE — Stage 053: Exact Overlap-Boost Window for the Lowest Support Lane

## Purpose

Stage 35 reduced the non-twin rescue problem to two concrete lowest-lane resources:

- overlap boost,
- support softening.

The next honest step is to ask whether the overlap side can be made explicit from a finite-throat operator instead of being left as an abstract factor
`Omega_0 = I_(phi,0)/I_W`.

This stage does that.

The main results are:

1. for any nonnegative lowest-lane source density `sigma_phi(s)` on the finite throat with the same total source strength as the mixed lane,
   the overlap boost satisfies the exact bound

   `0 <= Omega_0 <= pi/2,`

   so

   `A_I := Omega_0^2 <= pi^2/4;`

2. the symmetric uniform source is the baseline point `Omega_0 = 1`;
3. an explicit exponentially bottom-biased source family continuously deforms the overlap from `1` up to the sharp upper limit `pi/2`;
4. therefore **pure overlap asymmetry alone** can beat the Stage-35 support threshold only if

   `zeta_req <= pi^2/4.`

So this stage turns the phrase “maybe the support lane overlaps better” into an exact finite-throat theorem window.

---

## 1. Exact lowest-mode overlap for a general coherent source density

Keep the same D/N lowest mode on the finite throat interval `s in [0,L]`:

`chi_0(s) = sqrt(2/L) sin(pi s / (2L)).`

The mixed lane in Stage 32 used the uniform source density

`sigma_W(s) = 1,`

so its lowest-mode source overlap is

`I_W = int_0^L ds chi_0(s) = 2 sqrt(2L) / pi.`

Now let the physical lowest support lane couple through a general nonnegative source density `sigma_phi(s)` with the same total source strength,

`sigma_phi(s) >= 0,`

`int_0^L ds sigma_phi(s) = L.`

Define the support overlap

`I_(phi,0) := int_0^L ds sigma_phi(s) chi_0(s),`

and the overlap boost

`Omega_0 := I_(phi,0) / I_W.`

Then the lowest-lane coherent asymmetry factor from Stage 35 is

`A_I := Omega_0^2.`

---

## 2. Exact finite-throat overlap bound

Because `chi_0(s)` is nonnegative on `[0,L]` and satisfies

`0 <= chi_0(s) <= max chi_0 = sqrt(2/L),`

one has the exact bound

`I_(phi,0) <= (max chi_0) int_0^L ds sigma_phi(s) = L sqrt(2/L) = sqrt(2L).`

Dividing by the uniform mixed overlap gives

`Omega_0 <= sqrt(2L) / (2 sqrt(2L)/pi) = pi/2.`

So the overlap boost window is exactly

`0 <= Omega_0 <= pi/2,`

and therefore

`0 <= A_I <= pi^2/4.`

The upper bound is sharp: it is approached by a source density that concentrates at the antinode of `chi_0`, i.e. near the D/N bottom end of the finite throat.

---

## 3. Uniform branch and explicit constructive family

The symmetric mixed branch uses the uniform source density `sigma_W=1`, so

`Omega_0 = 1.`

To exhibit a constructive non-twin family, choose the bottom-biased exponential density

`sigma_alpha(s) = alpha exp(alpha s/L) / (exp(alpha) - 1),`

with `alpha >= 0`.

It has the same total source strength,

`int_0^L ds sigma_alpha(s) = L,`

and its exact overlap is

`I_alpha = int_0^L ds sigma_alpha(s) chi_0(s)`
`       = 2 sqrt(2L) alpha (2 alpha exp(alpha) + pi)`
`         / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ].`

Therefore the exact overlap boost is

`Omega_exp(alpha)`
`= I_alpha / I_W`
`= pi alpha (2 alpha exp(alpha) + pi)`
`  / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ].`

Its exact endpoint values are

`Omega_exp(0) = 1,`

`lim_(alpha -> +infinity) Omega_exp(alpha) = pi/2.`

So this explicit family interpolates continuously from the symmetric twin value to the sharp finite-throat upper bound.

For small asymmetry,

`Omega_exp(alpha)`
`= 1 + (2/pi - 1/2) alpha + O(alpha^2),`

and since `2/pi - 1/2 = (4-pi)/(2pi) > 0`, the constructive branch immediately moves upward from the symmetric point.

---

## 4. Exact pure-overlap rescue criterion

Stage 35 showed that, at equal stiffness,

`A_I = Omega_0^2 >= zeta_req`

is the exact rescue condition.

Stage 36 now bounds the left-hand side by

`A_I <= pi^2/4.`

Therefore **pure overlap asymmetry alone** is possible only if

`zeta_req <= pi^2/4.`

Equivalently:

- if `zeta_req <= pi^2/4`, a non-twin lowest lane can in principle beat the threshold using overlap enhancement alone;
- if `zeta_req > pi^2/4`, then overlap enhancement alone is impossible, and support softening must contribute.

This is the first exact no-go/sufficiency split for the overlap resource by itself.

---

## 5. Best current theorem statement after Stage 36

### What is exact now

- the finite-throat lowest-mode overlap boost satisfies

  `0 <= Omega_0 <= pi/2`,

- the corresponding asymmetry factor satisfies

  `0 <= A_I <= pi^2/4`,

- the symmetric uniform source gives `Omega_0 = 1`,
- the explicit exponential bottom-bias family gives

  `Omega_exp(alpha)`
  `= pi alpha (2 alpha exp(alpha) + pi) / [ (4 alpha^2 + pi^2)(exp(alpha)-1) ]`,

  with endpoints `1` and `pi/2`,
- and pure overlap rescue is possible only if

  `zeta_req <= pi^2/4.`

### What this means physically

The overlap side of the Stage-35 non-twin budget is no longer vague.

The moving-throat operator may indeed produce `Omega_0 > 1`, but even the most favorable finite-throat concentration can supply only a finite factor `pi^2/4` in `A_I`.

So the remaining question is already sharper:

> if the physical branch demands `zeta_req > pi^2/4`, then the lowest support lane must become softer as well; overlap enhancement alone cannot rescue it.
