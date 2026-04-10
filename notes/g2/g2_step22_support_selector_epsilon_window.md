# Step 22 — Support-demand selector and the exact `\epsilon_*` window

## Goal

Step 21 left one real microscopic ambiguity in the quartic layer:

```math
q_W \text{ versus } q_\Lambda.
```

The threshold is exact:

```math
q_W = q_\Lambda
\iff
\epsilon_* = \frac{1}{2+\beta^2},
\qquad
\beta = \frac{2}{11+9\delta_{U,*}}.
```

So the obvious next question is no longer “what are the primitive weights?”
It is:

> **what does the moving-throat support-selection side allow the coherent branch
> to do in `\epsilon_*` at all?**

That question is now answerable, because the later moving-throat support notes
already sharpen the coherent branch to an exact support-demand law:

- the coherent support enhancement is
  ```math
  S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
  ```
- the mixed-only product scale is
  ```math
  C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2},
  ```
- and the branch demand is carried by the selected-branch product `\Pi_{\rm tr}`. 

So Step 22 turns the support theorem into an exact selector on `\epsilon_*`, and
then intersects that selector with the Step-21 `q_W/q_\Lambda` crossover.

The main result is:

```math
\boxed{
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
}
```

collapses the whole support-selection problem to one scalar demand parameter, and

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho.
\end{aligned}
}
```

So the support-selection side already constrains the last live quartic dominance
ambiguity much more strongly than Step 21 by itself.

---

## Inputs from the moving-throat support theorem

From the coherent support-compensation stages,

```math
S_{\rm req} = \frac{\Pi_{\rm tr}}{C_{\rm mix}},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

The exact regime split is

```math
\Pi_{\rm tr} \le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only branch already enough,}
```

```math
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin support already enough,}
```

```math
\Pi_{\rm tr} > 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry is required.}
```

Those are the exact reduced-level support-selection facts carried forward from the
moving-throat notes.

---

## Step 22A — One scalar support-demand selector

Define

```math
\boxed{
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
}
```

Then

```math
\Pi_{\rm tr} = \frac{16\Lambda\varrho}{\pi^2},
```

and the exact required support enhancement becomes

```math
S_{\rm req}
=\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=\frac{2\varrho}{1-\epsilon_*}.
```

So the support theorem is now phrased in the simplest possible way: one scalar
`\varrho` competes against the split-blocking variable `\epsilon_*`.

---

## Step 22B — Exact `\epsilon_*` windows for the three support regimes

### Mixed-only enough

The condition

```math
\Pi_{\rm tr} \le C_{\rm mix}
```

becomes

```math
\frac{16\Lambda\varrho}{\pi^2}
\le
\frac{8\Lambda(1-\epsilon_*)}{\pi^2},
```

hence

```math
\boxed{
\epsilon_* \le 1 - 2\varrho.
}
```

So the mixed lane alone is sufficient only when the branch demand is low enough
that `\epsilon_*` lies below that ceiling.

### Symmetric lowest twin enough

The condition

```math
\Pi_{\rm tr} \le 2 C_{\rm mix}
```

becomes

```math
\boxed{
\epsilon_* \le 1 - \varrho.
}
```

Therefore the exact twin-sufficient window is

```math
\boxed{
1-2\varrho < \epsilon_* \le 1-\varrho.
}
```

### Non-twin asymmetry required

As soon as

```math
\Pi_{\rm tr} > 2 C_{\rm mix},
```

one has

```math
\boxed{
\epsilon_* > 1-\varrho,
}
```

and the symmetric lowest-twin support branch is no longer enough.

So the exact support selector is

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho.
\end{aligned}
}
```

### Existence conditions

Two immediate corollaries are exact.

- The mixed-only window is nonempty only if
  ```math
  \varrho < \frac12.
  ```
- The twin-sufficient window is nonempty only if
  ```math
  \varrho < 1.
  ```
- If
  ```math
  \varrho \ge 1,
  ```
  then every constructive branch point already lies in the non-twin-required
  regime.

So `\varrho` is already a branch classifier.

---

## Step 22C — Exact `\sigma` windows

The coherent quartic microledger does not use `\epsilon_*` directly. It uses

```math
\sigma = \frac{2\epsilon_*}{1-\epsilon_*}.
```

So the support selector is most useful after conversion to `\sigma`.

### Mixed-only ceiling

Insert `\epsilon_* = 1-2\varrho`:

```math
\sigma_{\rm mix,max}
=\frac{2(1-2\varrho)}{2\varrho}
=\boxed{\frac{1}{\varrho}-2}.
```

### Twin-sufficient ceiling

Insert `\epsilon_* = 1-\varrho`:

```math
\sigma_{\rm twin,max}
=\frac{2(1-\varrho)}{\varrho}
=\boxed{\frac{2}{\varrho}-2}.
```

So the exact `\sigma` windows are

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\sigma \le \frac{1}{\varrho}-2,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& \frac{1}{\varrho}-2 < \sigma \le \frac{2}{\varrho}-2,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \sigma > \frac{2}{\varrho}-2.
\end{aligned}
}
```

That is the cleanest Step-22 bridge back into the quartic anomaly ledger.

---

## Step 22D — Intersecting the support selector with the Step-21 `q_W/q_\Lambda` crossover

Step 21 proved

```math
q_W = q_\Lambda
\iff
\epsilon_* = \epsilon_\times := \frac{1}{2+\beta^2},
```

equivalently

```math
\sigma_\times = \frac{2}{1+\beta^2}.
```

So the support-selection question is now simple:

> does the support-compatible `\epsilon_*` interval ever rise above
> `\epsilon_\times`?

### Mixed-only branch

The mixed-only branch can reach `q_\Lambda > q_W` only if its upper ceiling lies
above the crossover:

```math
1-2\varrho > \frac{1}{2+\beta^2}.
```

Solving gives the exact mixed-only threshold

```math
\boxed{
\varrho < \varrho_{\rm mix}^\times
:= \frac{1+\beta^2}{2(2+\beta^2)}.
}
```

So if

```math
\varrho \ge \varrho_{\rm mix}^\times,
```

then every mixed-only-compatible branch point is forced into

```math
q_W > q_\Lambda.
```

### Symmetric-twin-compatible branch

The twin-sufficient branch can reach `q_\Lambda > q_W` only if its upper ceiling
lies above the same crossover:

```math
1-\varrho > \frac{1}{2+\beta^2}.
```

So the exact twin threshold is

```math
\boxed{
\varrho < \varrho_{\rm twin}^\times
:= \frac{1+\beta^2}{2+\beta^2}.
}
```

Hence if

```math
\varrho \ge \varrho_{\rm twin}^\times,
```

then even the whole mixed-only-plus-twin-compatible support sector lies below the
Step-21 crossover, and the quartic branch is forced into

```math
q_W > q_\Lambda.
```

---

## Step 22E — Exact `(\varrho,\beta)` phase diagram for the last dominance ambiguity

The support selector and the Step-21 crossover now combine into three exact
phases.

### Phase A — strong enough support demand that `q_W` always wins on every support-compatible branch

If

```math
\boxed{
\varrho \ge \varrho_{\rm twin}^\times,
}
```

then even the twin ceiling

```math
\epsilon_* = 1-\varrho
```

sits below `\epsilon_\times`, so every mixed-only or symmetric-twin-compatible
branch obeys

```math
\boxed{q_W > q_\Lambda.}
```

### Phase B — mixed-only still wall-dominated, but the twin window can cross into `q_Λ` dominance

If

```math
\boxed{
\varrho_{\rm mix}^\times \le \varrho < \varrho_{\rm twin}^\times,
}
```

then

- the mixed-only branch is forced into `q_W > q_\Lambda`, but
- the upper part of the symmetric-twin window can cross into `q_\Lambda > q_W`.

So in this band, outgoing-scale dominance is possible only with support help or
on a still-harder non-twin branch.

### Phase C — low enough support demand that even the mixed-only branch can reach `q_Λ > q_W`

If

```math
\boxed{
\varrho < \varrho_{\rm mix}^\times,
}
```

then even the mixed-only branch can enter the Step-21 strong-blocking regime
where

```math
q_\Lambda > q_W.
```

So the last live quartic dominance ambiguity is now sharply tied to one support
selector.

---

## Step 22F — Numerical size of the thresholds on the constructive coherent branch

Because

```math
0 < \beta < \frac{2}{11},
```

these support-demand thresholds live in very narrow windows:

```math
\boxed{
\frac14 < \varrho_{\rm mix}^\times < \frac{125}{492} \approx 0.25407,
}
```

```math
\boxed{
\frac12 < \varrho_{\rm twin}^\times < \frac{125}{246} \approx 0.50813.
}
```

So the physics is clean.

- Once the normalized support demand is even moderately large,
  ```math
  \varrho \gtrsim 0.5,
  ```
  all mixed-only and symmetric-twin-compatible branches are forced into the
  `q_W`-dominant phase.
- Only at relatively low demand,
  ```math
  \varrho \lesssim 0.25,
  ```
  can the mixed-only branch itself reach `q_\Lambda > q_W`.

That is much sharper than the purely parametric Step-21 phase diagram.

---

## Main result of the step

The branch-selection side now contributes a real selector, not just more symbols.

The exact support-demand variable

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
```

fixes the coherent support regimes by

```math
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho,
\end{aligned}
```

or equivalently

```math
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\sigma \le \frac{1}{\varrho}-2,\\
&\text{symmetric lowest twin enough} &&\iff&& \frac{1}{\varrho}-2 < \sigma \le \frac{2}{\varrho}-2,\\
&\text{non-twin asymmetry required} &&\iff&& \sigma > \frac{2}{\varrho}-2.
\end{aligned}
```

Intersecting that with the Step-21 crossover shows:

```math
\boxed{
\begin{aligned}
&\text{mixed-only can reach } q_\Lambda > q_W
&&\iff&&
\varrho < \frac{1+\beta^2}{2(2+\beta^2)},\\[1mm]
&\text{twin-compatible branch can reach } q_\Lambda > q_W
&&\iff&&
\varrho < \frac{1+\beta^2}{2+\beta^2}.
\end{aligned}
}
```

So the branch-selection side has now turned the last live quartic dominance
ambiguity into one exact support-demand phase diagram.

For moderate or large support demand, every support-compatible coherent branch is
forced into the

```math
q_W > q_\Lambda
```

phase. The `q_\Lambda > q_W` regime is then a low-demand corner or a genuinely
non-twin effect.

---

## Continuation point

The next clean move is now very focused:

> derive or estimate the actual selected-branch demand parameter
> ```math
> \varrho_* = \frac{\pi^2\Pi_{{\rm tr},*}}{16\Lambda_*}
> ```
> from the moving-throat normalization side.

Once that is known, the support regime and the `q_W/q_\Lambda` dominance phase
are no longer qualitative at all; they become a direct branch verdict.
