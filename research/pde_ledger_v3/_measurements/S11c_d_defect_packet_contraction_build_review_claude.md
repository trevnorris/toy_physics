CLEAR FOR THIS BOUNDED FINITE-WINDOW CONTRACTION BUILD

This is a source-only assessment. I ran nothing, and no runtime evidence exists. I found no equation, source-binding, domain, unit or executable blocker. The items below are non-blocking limits and risks. Any one of them would make the worker refuse and preserve its evidence. None could produce a silent false acceptance.

## What I read
- **Packet documents:** `build.md`, `evidence-guide.md`, `worker.py` (all 371 lines), `library.py` (all), `launcher.py` (all), `execution-authority.json`, `source-contracts.json` (all fragments), `method.md` (all) and `runtime-source/containment-helper.py` (all).
- **Manifest:** `input-manifest.json` lines 1–120 and 7735–7781, plus greps. The grep counted exactly 1,534 `savedInputs` byte records, matching the build's stated 1,534 files. I did not hash or recount the files on disk.
- **Saved evidence:**
  - `accepted-units/complete-template-plus-normal-NATIVE_MIXED_ITERATION.json` (all).
  - The `minus-normal` and `minus-pressure` DIRECT templates, via grep: `-I*packet_D*packet_qo`, `normalSign` = null for pressure.
  - `accepted-units/summand-9009-return.json`, `new-a-units-return.json` and `new-beta-units-return.json`.
  - `preflight/Fourier-envelope-constants.json`: the first two entries and the proof block.
  - `accepted-units/posthashes.json`: grep only, which showed `inner_lib` and `fourier_lib` entries with `intact: true`.
  - `tooling-tests.log`: grep only, "34 tests OK".
  - `pressure-addresses.json`: one grep for the applicable-component count (64 matches).
  - Greps over all summand returns found no `"argumentDerivative": 1`.

I did not open `runtime-source/shared-guard.py`, `supervisor.py` or `completion-hook.py`. I did not open the other 543 summand records, the 17 factor operand sets, the 34 field/polynomial files, the eight geometry plans, `physical-plan.json`, `binding.json`, the `whole-origin` files or the full 2,652-entry THETA row. I did not open `tooling-tests.py`, `AGENTS.md` or `packet-action-method.md`.

I could not verify any SHA pin, and no gate or build-review record exists. Runtime-only joins are therefore not proof: the 544-address join, the full THETA selection, the 20-template equality and the unit `==` checks. They are written as fail-closed `require`s.

## Findings by requested item

**1. Kernel AST, middle call, routes and branch (`worker.py:158-168`, `source-contracts.json`).**
- The worker compares `middle` assignments by `ast.dump` against: `qi,qo=q(k),q(l)`, `qh=q(k+t)`, `qs=q(l-t)`, `used_qs=qh if mutate else qs`, `a=1/(1-3j/10)`, `mu=3/10`, `beta=a*mu`, `aa,bb=profile(t),profile(l-k-t)` and the `kernel_components(...)` call.
- `q` is checked as `595/100-p*p` with `sqrt(d) if d>=0 else j*sqrt(-d)`. That is the positive-real/positive-imaginary outgoing branch.
- Nothing is called. `Arithmetic.statements` interprets only the pure assignment AST on fresh symbols. It rejects calls, attributes and unknown syntax, and requires the formal-argument set to equal the environment (15 names).
- `qi`, `qo`, `qh` and `qs` are distinct symbols. `qother` stands for the opposite root and is checked absent from each primitive's left side (`worker.py:210`). This check is structural and weak, but it does what the method asks.

**2. Four separate identities (`worker.py:183-210`).** I re-derived each by hand from the `kernel_components` text and found no discrepancy.
- **J:** `Cj·Y0·m/(qm(qm+β))·[C1T+m·C0T]`. This uses (k+t)=m and (2k+t)=k+m. The `(qo+β)`, `(qi+β)` factors go into y and x, and the `(qh+β)` factor stays inside.
- **Dr:** with t=l−m, 2l−t=l+m and qs→qm. This gives `X1·[Y1CT+m·Y0CT]`, with A1=A(l−m) and A2=A(m−k).
- **Dh:** `Y0/qm·[C2T+m·C1T]`.
- **Dq:** `Y0/qm·X02T`.
- `Cd=(-I·mu)·W·L/(4I)` matches `pref`, and `Cj` matches the J prefactor.
- Only the existing external resolvents are removed. No whole D is multiplied by another middle integral. Each identity is checked by `exact()`, which saves operands and raw residual, then requires `cancel(together(raw)) is S.Zero`.

**3. Window orientation.**
- The input route (m=k+t) clips k. The output route (m=l−t) clips l.
- The slack identities check out. `worker.py:289-293` is redundant but correct: input is `T−(m−k)=k−(m−T)` and `T+(m−k)=(m+T)−k`; output is `T+(l−m)=l−(m−T)` and `T−(l−m)=(m+T)−l`. The forward/inverse identities are `k+(m−k)=m` and `l−(l−m)=m`.
- `library.finite_window` and `mapped_contains` give the right clipped interval, full wings and empty cases. The witnesses cover |m|≤M, M−1, M+1, ±(T−K) and 0 for both (27,122) and (29,124).
- This is a continuous identity. No mesh or quadrature is created.

**4. Wrong-root mutant (`worker.py:213-215`).** I hand-checked it. Substituting qs→qh, then t=m−k, gives numerator k(2l−m+k) over (qm+qo), with k clipped and l unclipped. That equals `2·X1T·Y1C+(X2T−m·X1T)·Y0C`, so the formula in `method.md` is correctly reproduced from its own source substitution.

**5. Controls (`worker.py:302-318`).**
- Wrong window: m=T+K−½, k=−K, l=K gives good=True, bad=False. The point lies in the rectangle (t=l−m=−121.5, inside ±T).
- Wing omission: k=K, l=0, m=T+K−½ is retained, and the omission mutation drops it.
- Algebra mutants: `k(3l−m)` (sign mutation of 2l−t with t=l−m) and `m·k` (J polynomial missing m²). Each is built by factoring the actual complete expression. At the predetermined point (k,l,m)=(1,2,3) with all other symbols 1, the movements are nonzero and rational, so there is no accidental zero denominator.
- All of these are algebra/domain sensitivity only, as stated. The wrong-window and wing controls are tests of `library` predicates on predeclared witnesses. Those predicates are short and are the same logic the certificate relies on. The checks are therefore weakly independent, but they fit the bounded purpose.

**6. Templates, 544 addresses, units.**
- The injection test (`worker.py:228-236`) matches the saved templates I read. For plus/normal/NATIVE_MIXED the template is `I*packet_qo*(… packet_H … + packet_J)`. The injected increment gives `I·packet_qo·δ`, which equals `normal·δ` with normalSign=1. For minus/normal/DIRECT it is `−I·packet_D·packet_qo`. The `H`, flat, height and slope slots are untouched.
- Selection is `selected==[v for v in row if row=='THETA_BALANCE' and channel=='e_W']` with 2,652 entries. This is a real full-row join if it runs.
- Address-level joins are strong: address == `sourceTransportInput.address`, unit return == `complete-summand-dimensions` record, adapter, normal and wave operands, factor operands, field identity, `addressJoins[index]`. Explicit-zero addresses are kept, with `epsilonCount` 0.
- New units are derived from the saved `ur['X']`, `ur['Y']`, `a`, `mu`, `W_0`, `L_W` and include each dk, dl and dm. They are compared to `ur['total']`. For 9009, the saved total of (−2,−1,1) is X(2,−1,0)+Y(1,1,0)+kernel(−3,−1,1)+measure(−2,0,0), which is consistent. A mismatch would refuse.

**7. Majorants.** Hand-checking the constants in `worker.py:275-281`:
- `|Cj|≤(130/109)(9/100)WL/4`. Since |a|=1/√1.09≈0.958, 130/109≈1.193 is an upper bound.
- `|Cd|=µWL/4`.
- `|A|≤1/(2π)<1` from |sinh y|≥|y|.
- `|q|≤|p|+3` (κ≈2.44).
- `|q+β|≥Re β=30/109` holds because Re q≥0. This is even stronger than the √2 bound the method cites.
- `|q1/(q1+q2)|≤√2≤2` and `1/|qm+qo|≤√2/|qm|` hold for first-quadrant values.
- The factors M(M+K) for J, K(M+K) for Dh and Dr, and (K+3)² for Dq all match the polynomial factors. The extra `/b` for the (qm+β) denominator appears only in J.
- The envelope X bound is a global supremum: the proof says "CX·exp(−5|k−p0|)", so the constant is valid on [−K,K]. The saved bounds are rational constants.
- The quadrant polynomial is `2(ux·vx+uy·vy)`, and the root-sum identity holds by hand expansion. The Fubini application is correctly labeled as analytic, not a machine proof.

**8. Launch and integrity.**
- `verify_gate` requires pins on worker, manifest, library, launcher, guard, supervisor and authority. It also requires a build record with both literal CLEAR verdicts, the method record, `argv`/output equality and `durationLimits is None`. No gate exists, so nothing can run today.
- `containment()` is a single inert function, exec'd alone. It checks memory.max=4 GiB, swap 0, pids 32, one CPU, one thread, no CPU rlimit, and sets RLIMIT_AS.
- The launcher hook-first handshake is correct: the coordinator waits for the watcher's "waiting" state and a "1" byte, and the sympy import happens only after the gate.
- `Journal.emit` does an exclusive `'x'` open, `fsync`, then a chained hash log. The `finally` block always writes posthashes and `checks.json`, and the traceback goes to stderr.
- The authority `scope` is identical to the manifest `scope` as I read both, with `scienceExecutionsAuthorized` 1 and `noDeadline` true.

## Non-blocking limits and risks (hypotheses, not established failures)
1. **Window and profile labels are not tied to the integrands.** `Amk`/`Alm` are free symbols. `windows`/`variables` are labels on the integrands, not enforced. The profile-argument identification and the window assignment rest on the separate affine and interval identities, plus the explicit A1/A2 mapping. The evenness of A is not needed. I judge this adequate because the identities are explicit, but the symbolic check cannot detect a label error.
2. **Label transport is shallow** (`worker.py:321-328`). It copies each original line into four per-primitive records with `numericalMeshReady: False`. It adds an "actualDepthDependencies" record and the new branch/window label lists, but does not map individual old cut labels into m-coordinates. `method.md` defers the exact discretization, so this is acceptable here. A numerical build must redo it.
3. **Units for β and `qother`.** β is declared MOMENTUM by assumption, not checked against unit(a)+unit(mu) (a is (3,1,−1)). `new-beta-units-return` is momentum and the final totals would expose an inconsistency, but only indirectly.
4. **Envelope scope.** I could not confirm that the saved `X` and `Y` constants include `nativeJetMagnitude` and the Gaussian/Fourier normalization. The proof text says they do. I also confirmed that no address has a differentiated test factor (no `argumentDerivative: 1`), so using `bound['Y']` rather than `Yprime` is sound.
5. **`cancel(together(…))` over complex rationals** might not return a literal 0, for example for `a·mu − (30+9I)/109`. If so, the worker refuses and preserves its evidence. That would be a spurious refusal, not a false acceptance.
6. **Memory and disk.** All 1,534 JSON files (~67 MB) are parsed in memory under RLIMIT_AS = 4 GiB. Per-address journal records re-emit the full address and factor operands 544 times, with no stated storage forecast. A failure here would be preserved but would waste the single authorized run.
7. **The tooling tests are stdlib/synthetic only.** No scientific expression has been restored or checked.
8. **The old sinhc profile approximation** is only a documented numerical obligation. The certificate uses the exact analytic A, as the build states.

## Inherited assumptions
- The accepted unit bridge (7d5b0d7e) and the old cancellation proofs for fields, factors and templates.
- The envelope constants and their original proof.
- The physical binding (a, mu, β=(30+9i)/109, κ²=119/20, W=1, L=10).
- The measure-theory facts (Fubini and local integrability of s^(−½)).
- The geometry plans, preserved but not re-derived.

## Deliberately deferred (not required here)
The numerical evaluator, quadrature, nested error propagation, positive tail budgets, new Fourier requests, two independent routes, storage and runtime forecasts, numerical controls, H/contact/height-PV/slope/flat terms, and any integral, action, current, leakage factor, loss or inverse.

## Permitted runtime scope
Exactly one guarded run of this worker, after the gate, the build-review record and the pins exist. It runs under the shared pooled guard with 4 GiB, a 16 GiB pool, zero swap, one CPU/thread, 32 tasks, a 4 GiB host reserve and no deadline, and the completion hook is armed first. It certifies only the J/Dr/Dh/Dq algebra, the wrong-root mutant, the domain and unit joins, the majorants and the label/dependency records. There is no automatic retry. Runtime evidence must be inspected before the certificate counts as accepted. A successful run would permit only numerical-evaluator preparation. It would accept no integral or leakage factor.