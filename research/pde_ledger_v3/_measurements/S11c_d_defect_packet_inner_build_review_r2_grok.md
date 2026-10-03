I'll review the inner-kernel build against the assessed packet method and the saved operands, starting with the named sources and keeping peer reports out of scope.The review is still open. I am re-reading the control, tail, and join sites so the verdict cites the saved source rather than the session summary.**NEEDS REVISION**

No science has run, no payload was restored, and no READY gate exists. This is a source review of the inner bank only. It is not result acceptance, and it does not authorize a run.

The integrand, the 38-point set, the collision cuts, the global-leaf adaptive rule, the J/D tail algebra, and the failure journal are in place. Two declared gates can still succeed without performing the check they name.

## Blockers

**1. The H-contact control never omits the contact.**

`numerical-library.py` adds the contact only on the momentum routes, then gates a typed formula:

```232:250:packet/numerical-library.py
contact=5*self.profile(c,Q)/4
...
a24=[contact+self.gauss(...)[0]]
a48=[contact+self.gauss(...)[0]]
...
omission=-5*self.profile(c,Q)/4
```

```280:282:packet/numerical-library.py
Q=l-k;omission=-5*self.profile(c,Q)/4
self.emit(c,'controls/H-contact',{'point':point,'baseline':H,'mutated':H+omission,'movement':omission})
require(c.isfinite(omission) and abs(omission)>c.mpf('1e-12'),'responsive H-contact control')
```

`H` here is the cached physical route-B value (`self.cache[...] = (b[0], receipt)` at line 251). That integrand is `h*j*exp(-iQx)/(2π)` and has no separate contact addend. `movement` is `-5*A(Q)/4` at 50 digits. The require reads that formula alone. Baseline and mutant are not inputs to the gate.

`build.md` lines 144–148 ask to omit the surviving H contact at `(κ/5, 2κ/5)`, persist baseline, mutation, and movement, and require finite movement above `1e-12`. Method §5 (`method.md` lines 381–387) says a formal symbolic movement is not a numerical success. This gate passes for any nonzero profile value at that `Q`, including when the momentum routes never add the contact. Route comparison can still fail a one-sided omission later; the control record itself cannot.

Required correction: on the control point only, recompute the momentum H on the same panels, rules, and `T`, once with the contact addend and once without it. Persist those two integrals and their difference. Gate on that measured difference. Keep the physical product as the independent H reference; do not subtract the typed contact from it and call that an omission.

**2. The primitive comparison accepts a short vector.**

```207:216:packet/numerical-library.py
def compare(self,key,reference,candidates):
    ...
    for label,vals in candidates.items():
        for i,(a,b) in enumerate(zip(ref,vals)):
```

`zip` stops at the shorter sequence. A candidate that drops `D_quadratic` or the sum can pass when the shared prefix matches. `build.md` lines 89–94 require every primitive and the sum to meet `1e-9 + 1e-7|A48|`. Current routes all return five entries (four kernel terms plus the sum) or one H value, so this is latent on the present call path.

Required correction: require equal component counts for every candidate, then compare every index. A length mismatch must fail closed with the same no-retry preservation as a tolerance miss.

## Fixed 38 points

Yes. They exercise the declared internal lines and the two-sided external approaches.

`point_plan` (`numerical-library.py` lines 13–28) builds 38 distinct pairs: four families `(1/3,-1/3)`, `(1/3,5/3)`, `(-1/3,-5/3)`, `(1/3,1/3)`, each with `l` shifts `0, ±1/64, ±1/128` (20); each of `k` and `l` approaching `+1` and `-1` from both sides at `±1/64` and `±1/128`, other coordinate `1/5` (16); plus `(0,0)` and the control `(1/5,2/5)`. No sample has exact `|k|=κ` or `|l|=κ`.

`exact_cuts` (lines 31–42) keeps exact `aκ+b` keys and merges only identical keys, retaining every label. On `k+l=0` the height roots and reflected roots coincide. On `k+l=±2κ` one height root coincides with one reflected root. On `l=k` the `t=0` and `t=Q` labels coincide. The `±1/128` shifts do not share those keys. Distinct resolved cuts that compare equal refuse (line 125). Endpoints stay the original objects (`partition_interval`, lines 45–49).

At a merged root, `D` stays a sum of three terms (`kernel_components`, lines 53–60). There is no `1/(qh qs)` product. `qs+qo=0` would need both outgoing roots zero, which is external grazing and is not in the set. `(0,0)` leaves the reflected denominator `qs+qo` nonzero because `q(0)=κ`.

## What can silently pass

- The H-contact gate in blocker 1.
- A short or long candidate in `compare`, blocker 2.
- A permutation of the three `D` addends. The runtime lock is the commutative sum (`worker.py` line 143). The current source order matches the three summands in `saved/preflight/new-D-numerical-adapter-input.json`: reflected `k(2l-t)/(qo+qs)`, height `k qi(2k+t)/(qh(qh+qi))`, quadratic `qi²/qh`, with the shared prefactor equal to `-15 t (-10k+10l-10t)/(32 (qi+β)(qo+β) sinh(5πt) sinh(5π(l-k-t)))` at `β=(30+9I)/109`. A later swap would still zero the sum, and both routes call the same function.
- A coordinated retype, inside `numerical-library.py` only, of `h`, `j`, `A`, the contact, the analytic product, `1/(2π)`, and `sqrt(595)/10`. Symbolic joins in `worker.py` lines 146–177 check the lemma and adapters against themselves. They do not evaluate `profile`, `fa`, or `fb`. One-sided disagreement among A24, A48, physical B, and the analytic product still fails `compare`. The present literals match the saved operands: `h=(1+tanh(x/10))/4`, `j=sech²(x/10)/4`, `A(t)=5t/(2 sinh(5π t))`, contact `5 A(Q)/4`, product `ĵ(Q)(1/4-10 i Q/8)`, Fourier factor `1/(2π)`.
- The H tail coefficient. `worker.py` line 186 types `Rational(55,3)*2**(-T)`. The require at line 188 checks that each typed tail is positive and `<1e-11`. With `T=122`, a much larger coefficient still passes. `2^{-T}` is a valid overestimate of `e^{-T}`.

These do fail closed: a changed saved-input hash, an inexact MP tuple, a non-open or non-positive restored weight, a Gauss subset that is not seven Kronrod nodes, `qi=qo=qh=0`, float-aliased distinct kernel cuts, `a<mid<b` failing, and an equal-length comparison miss. There is no order retry.

## What holds on the current source

Native scale is a runtime join. `native_profile_scale` (`worker.py` lines 83–89) requires the saved rule `L_W**number_of_native_spatial_indices`, the saved source fragment, `functionExecuted` false, and length 10 on the physical, declared, and saved operands. One spatial index supplies one factor of `L`. At `σ=1`, the lower outward slope joins lemma `j`. At `η=1`, the plus and minus lab heights join `+h` and `-h` (lines 154–167). The display grades `η=1/100` and `σ=1/1000` stay outside this kernel. Live frequency 3, `cs=sqrt(6)/2`, and `κ=sqrt(595)/10` are checked from the physical plan (line 103).

Both-face assembly matches the eight saved mixed and direct templates. Pressure is `-3 I H k qo / (10 (qi+β)(qo+β)) + J` or `D`. Plus normal multiplies by `+I qo`; minus normal by `-I qo` (`worker.py` lines 191–197 and `numerical-library.py` lines 267–269). The mapped templates already use one `packet_H`, one `packet_J`, and one `packet_D`; the originals’ `Jwhole_plus` / `Jwhole_minus` and `Dwhole_plus` / `Dwhole_minus` names were identified with those symbols by the inherited zero returns. `q` is one outgoing sheet: positive real inside the cut, positive imaginary outside (lines 95–98). Reflected `q(l-t)` stays distinct from `q(k+t)` except in the reflected-node mutant.

The paired subtracted density joins the half-line integrand. `positive(t)+positive(-t)` from `new-H-subtracted-adapter-input.json` equals `5 A(t) (A(Q-t)-A(Q+t)) / (2 I t)`; the `A(Q)` pieces cancel (`worker.py` lines 171–172). The contact join equals `5 A(l-k)/4` (line 170). Lemma `hTransformContact` is the constant `1/4`, a different object from this product contact.

Route A halves each panel and uses `x=a+(c-a)z²` and `x=b-(b-c)z²`, Jacobian `2*length*z`, then divides by 2, which is `dx` for `z=(n+1)/2`. The `z` factor cancels a square-root endpoint, and open nodes stay off the cut (`gauss`, lines 134–149). Route B integrates in physical `t` with the restored G7/K15 rule and does not read route A’s nodes.

The physical adaptive routine does what this build now claims (`adaptive`, lines 172–205). The target is one global budget per component, `1e-11/(4*38)`. The heap key is the largest maximum component error. The parent is deleted and replaced by two physical halves that share one midpoint object. Acceptance breaks only after the error vector is recomputed from the current leaves and every component is within budget. The interval 64 recomputes the sum; it is not a refinement cap. A negative incremental sum also forces that recompute. If `a<mid<b` fails, the panel prefix is emitted and the run raises. The stored criterion says the error is empirical. An integrable root can meet a fixed budget as the leaf shrinks; if the estimate does not fall, the run stops. That behavior was not executed here. The synthetic oracle in `tooling-tests.py` is not a kernel result.

J and D tails match the kernel under the stated hypotheses `|t|>T≥K+4`, `κ<3`, depth moduli above 1, quadrant sums, `|q+β|≥b=3000/11101`, `|a|≤1`, and `A(t)A(Q-t)≤121 exp(-|t|)` (`worker.py` lines 183–188). Direct is the sum `2 D_reflected + D_quadratic`. These are absolute reference-coordinate majorants. The same lines set `evaluatedIntegral` false. Nothing in the worker return claims a packet action, a current, a loss, or a proved quadrature bound: `packetActionEvaluated` false, `currentOrLoss` null, `fullActionAccuracyClaim` false, and `scientificAcceptance` forced false in `finally`.

Journal writes are durable. `DurableStore.put` uses `PRAGMA synchronous=FULL` and one committed `INSERT` per record. A failed panel is written before the exception leaves `gauss` or `adaptive_panel`. `worker.py` `finally` closes the store, hashes the SQLite file, writes `posthashes.json` and `checks.json`, and forces exit code 1 on an integrity miss. An existing journal path is refused. There is no automatic retry.

The reflected control does recompute the density: at `t=1/10` on the control point, `middle(..., mutate=True)` replaces `qs` by `qh`, and the gate uses `wrong[1]-base[1]` (lines 274–278). The lower-normal record negates the stored factor `-I qo D` on the completed direct sum. That factor sits outside the integral, so the algebraic flip matches the narrow check in `build.md` lines 146–147. It is not a second integration.

## Optional improvements

- Zero each D addend against its saved summand, not only the sum.
- Feed `profile`, `fa`, `fb`, and `κ` from the joined lemma, contact, subtracted density, and physical plan, so a self-consistent numerical retype cannot pass route comparison.
- Give the H certificate an equality join to a saved expression. The prose line in `tail-bound-derivation.json` says `H middle tail ≤(55/3) exp(-T)`. `pressure/H-bound.json` field `tailBound` is the different expression `10*exp(-5*pi)/pi**2`. Whether `55/3` follows from those envelope fields is unresolved here; only the sentence match and the `2^{-T}` overestimate were checked.
- Gate the lower-normal movement on the stored mutant minus the stored baseline.
- On the H momentum half-line, add the profile offsets `±1/10, ±1/5, ±2/5, ±4/5` about `|Q|`, keep labels, and refuse a float alias. Line 228 currently places `0`, `0.1` through `0.8`, `|Q|`, and `T` into an unlabeled set. The `A(Q-t)` peak sits on the existing `|Q|` cut, and a gross miss still has to survive the analytic and physical comparisons.
- `settingsIdentical: True` at line 221 is a constant. Inside one run the pins do not change; the flag does not compare them.

## Coverage limitations

Method §5’s `q(k)`-instead-of-`q(l)` normal control, the Leibniz corruption, addressed ancestry, and the ten-times-envelope bar are outside this bank. `build.md` lines 149 and 168–172 defer them. The reflected mutant is one node, not the integrated primitive.

Twelve of the twenty templates (flat, height, slope) are copied, ancestry-checked, and not integrated. `reference_height_hat` and `reference_slope_hat` are not this `H`. Normal-height zeros stay in the inherited templates.

`J` and `D` are also integrated at `T=124`. `H` is not. That enlargement is not the outer `K→K+2`, `T→T+2` check. `fullSummandUnitsEstablished` is false. The chi identity at `worker.py` lines 174–175 uses free symbols `A_even`, `f_plus`, `f_minus`, `f_zero`, and `chi_even`. The even cancellation that is actually tied to a saved operand is the paired-density join above. No exact external-grazing integral is in the point set. Route agreement does not certify the shared kernel convention. Adaptive errors and tail inequalities are not a global quadrature bound.

The `(55/3)` derivation from the saved H envelopes remains unresolved. No integral in this review was executed, so no numerical value, comparison margin, or tail integer is a measured result.