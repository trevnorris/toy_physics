I'll review this bounded contracted numerical build from the stated inputs, source, and saved operands only—no execution and no runtime inference. First I'll read the build contract, method, and evidence guide, then the worker, geometry, library, persistence, launcher, and native operands.The review is source-and-operand only. I’ll locate the contract, method, evidence guide, and the worker sources those documents name.`input/build.md` is present; the method and evidence guide are not at the work root. I’ll look through the packet layout for those files and the source they name.The packet sources are under `input/`. I’ll read the worker, preamble, geometry, numerical library, persistence, and launcher in full, then the native operands they join.The review is still open. I am checking the whole-definition joins, selector gates, and tail ancestry in the source and saved operands before giving one verdict.The preamble restores inherited zeros and writes new joins. I am reading the factor operands and the numerical routes those joins actually feed.**NEEDS REVISION**

This instrument does not yet establish that the values it would integrate are only the ordinary part of `Jwhole` and the three added parts of `Dwhole` for the twenty eligible addresses. The quadrature, clipping, and route machinery are aimed at the saved contraction skeletons. The joins that would attach those skeletons to the address wholes accept identical sides.

## Blockers

1. The eligible whole-definition joins are structural self-identities, and the numerical target never names the wholes.

`prepare.py` lines 77–82 accept every `address-full-factor-0` … `16` return solely because `cancelled` is `Integer(0)`. The refusal text says this is an inherited full-expression proof and that the sides are not structurally equal. The inputs for the two eligible templates are equal on both sides:

- Address 8346 names `address-full-factor-3`. Both sides are the mixed expression `Hwhole(...) + Jwhole_plus(...)` (`address-full-factor-3-full-mapped-residual-input.json` lines 2–8). The return is raw `0` and cancelled `0`.
- Address 8347 names `address-full-factor-4`. Both sides are `Dwhole_plus(...)` (that residual input, lines 2–8; `address-8347-input.json` lines 227–283).

The same pattern is what the preamble treats as the native injection proof. `prepare.py` lines 127–128 require only `residual == Integer(0)`. These inputs have both sides equal to the symbol `contract_injected_whole`:

- Plus mixed `packet_J`, first address 8346. The context template is the `packet_H` term plus `packet_J`; the original is `Hwhole + Jwhole_plus`.
- Minus mixed `packet_J`, first address 9672. The original is `Hwhole + Jwhole_minus`.
- Plus direct `packet_D`. The original is `Dwhole_plus`.

`numeric.py` `combine()` (lines 110–116) and the inner vector (lines 227–230) evaluate `J`, `Dr`, `Dh`, and `Dq` only. Those expressions match `new-contraction-definitions.json` outer integrands, which `prepare.py` lines 150–151 cancel against `combine()`. They contain no `Hwhole`, `Jwhole_plus`, `Jwhole_minus`, or `Dwhole_plus`. `saved/pressure/whole-tags.json` holds the saved `Jwhole` and `Dwhole` definitions and is not read by `prepare.py`.

The distinct-side records are the contraction factorizations, and they are inherited only as a stored zero:

- `new-full-factorization-J`, `Dr`, and `Dh` have different left and right kernel sides. Their `originalDensity` values match the ordinary and three-piece shapes in the `kernel_components` fragment of `numeric-source-contracts.json`. They do not name `Jwhole` or `Dwhole`.
- `new-full-factorization-Dq-input.json` lines 2–7 have the same text on both sides, so that stored zero is also a self-identity.

Address 8346’s response coefficient is still the `Hwhole` term plus `Jwhole_plus` (`address-8346-input.json` lines 84–86). The return field `contractedOrdinaryComponent: true` is present (`address-8346-return.json` line 4) and is not required by `prepare.py`. Dropping `H` by calling `combine()['J']` is therefore unjoined to these address wholes.

2. The original `argumentDerivative` selector is written as 0 and never read.

`prepare.py` line 119 stores `argumentDerivative: 0` in every constructed X/Y spec. `numeric.py` line 178 repeats that constant in the transform settings. Neither reads `matchingXFamilyInterfaces` or `matchingYFamilyInterfaces`. Method §2 requires the original selector to be 0 and requires every other selector to be rejected. The `constant_reference` fragment applies extra derivative powers when that field is nonzero. I read selector 0 on the X interface of 8346 (return lines 93–96) and 8347 (return lines 62–66). The other eighteen unit returns were not opened in this pass, so their stored selectors remain unresolved. The missing rejection stands either way.

3. The tail preamble does not perform the required ancestry joins, and its budget is one grand total.

`prepare.wt` (lines 164–166) inlines the same rational expression as `exponential_moment` and `weighted_tail` in `S11c_d_defect_packet_preflight.py` lines 130–139. Method §5 requires an AST join of `exponential_moment`, `weighted_tail`, and `contributions`. There is no `ast.dump` comparison. The `tail-contributions` fragment is stored as journal context.

`contributions()` adds, for `NATIVE_MIXED_ITERATION`, the positive H middle term `(4/b)*cx*cy*E30*F[2]**2*(55/3)*2**(-T)` (`numeric-source-contracts.json` lines 57–58). The K29 branch (`prepare.py` lines 176–178) uses the displayed `middle_J` and `middle_D` and omits that term. `includesExtraHOvercount` is set only for `K==27` and primitive `J` (line 182), so the enlarged J rows are not labeled as the larger mixed budget.

The direct-piece ancestry is the prose lemma at `prepare.py` line 161, plus a nonnegative-coefficient check on the four numerator files (lines 155–158). There is no exact residual joining `2 P^2(1+|t|)` or `18 P^2(1+|t|)` to those polynomials.

Lines 183–185 require only the sum of every address and primitive tail at one K to be `< 1/10**11`. Method §5 requires a positive per-address and per-primitive check at `epsilon = 1e-11/80`. One row can exceed that primitive cap while the grand total still passes. Duplicating the full positive D envelope on `Dr`, `Dh`, and `Dq` is present and is the declared overcount.

4. The outer B leaf choice does not rank the weighted nested guard that can stop the run.

`numeric.py` lines 387–390 stop when every addressed component satisfies `max(error/eps, nested/(eps/4)) <= 1`. The leaf score is `part[key].error/eps` only. A component whose binding quantity is `nested/(eps/4)` can lose the rank to a leaf with a larger Gauss–Kronrod gap and a smaller nested contribution. A result that returns has still met both guards. The declared global choice, the leaf with the largest contribution to the current worst normalized component, is not what this score selects, so the loop can keep splitting a non-binding leaf.

## What the source does implement

These parts match the contracted method on the records and functions read. They do not close the gaps above.

- New exact joins exist for `beta`, `Cj`, and `Cd` on the exact stand-ins; for the centered moment identity at the actual `n` in `{0,2}`; for the Y phase `(p0-l)*(-5/2) = (l-p0)*(5/2)`; for the outgoing `q` radicand AST; and for `combine()` against the five saved outer integrands, including `Dr_wrong_root`.
- Constant quotients for an eligible entry join the degree-zero polynomial, the address field, and `constantValue`. Address 8346/8347 consumer originals contain one `epsilon_shape`, and the saved consumer fields are the epsilon-free values. The preamble sets `epsilonAlreadyExtracted: True` without a residual `consumerOriginal/epsilon_shape - consumerField`.
- Eligible census gates in source are pressure, normal `Integer(1)`, grades source/consumer `00` and target/response `11`, `epsilonCount == 1`, `n in {0,2}`, 64 applicable, 20 live, 10 ordinary J, 12 `EXACT_ZERO_SOURCE_JET`, and 32 `EXACT_ZERO_CONSUMER`. Address 8348’s status is `EXACT_ZERO_SOURCE_JET` (input line 196). The other explicit zeros were not opened one by one.
- X uses `nu = k-p0`, center `-5/2`, and `p^n/(2π)`. Y uses `nu = p0-l`, center `+5/2`, and no extra `1/(2π)`. `i^n` sits in `alpha` once. Route B’s moment sum is divided by `i^n`, so `(ik)^n` is not applied twice. The derivative mutant sets route B’s coefficient list to the single value `(i*p0)^n`, and the recurrence loop is empty (`numeric.py` lines 100–106). Route A replaces `p^n` by `carrier^n`. Y is unchanged.
- Clipping matches the definitions: `J`/`Dh`/`Dq` and the wrong-root control zero the input-`k` pieces outside `I(m)`; correct `Dr` clips output `l`; `Y0C`/`Y1C` stay on `[-K,K]`. The reflected denominator in the Dr factorization is `qm+qo` with the `l-t` root. The wrong-root formula is the saved `Dr_wrong_root` integrand.
- `geometry.py` splits at `{-M,-(T-K),-kappa,0,kappa,T-K,M}`, keeps six open slabs under `T-K>kappa>0`, intersects non-parallel graphs inside each slab, coalesces exact line keys in `canonical_lines` while retaining every label, and audits `true_window` at both endpoints and the midpoint. `wing_control` forces central `[-K,K]` clips onto each outer wing and requires the same audit to raise `actual max/min clipping window`. No sample is used to drop a root or an interval. `arrangement()` is not called.
- A24 and A48 use separate 30-digit contexts, open squared maps, positive Jacobians `2*half*u` and weight factor `1/2`, and their own inner GL24/GL48 pair, with absolute panel gaps accumulated before the contracted product. B50 uses an independent 50-digit context and physical `mid+half*node`. `Estimate` multiplication is `|a|e_b+|b|e_a+e_a e_b`. Addressed scaling uses one `alpha` per address after the shared normalized integral, so the largest `|alpha|` tightens that shared budget. Comparisons cover individual address/primitive keys and the face, grade `11`, J, D, and J+D groups, for A24–A48, A48–B50, and both window enlargements, at `1e-9+1e-7*|A48|`. Controls are 8347/`Dr` and 8350/`J`, all three routes, both carriers and both windows, with movement required to exceed ten times the summed finite-window indicators, cross-route gaps, and enlargement movements. Mutant rows set `allRealMutantClaim` false and do not receive the baseline tails.
- `make_mpf` (`ctx_mp_python.py` lines 597–600) stores the supplied tuple without normalizing it. `decode` therefore keeps the encoded mantissa, exponent, and bit count. Opposite-formula arithmetic still runs at the other context’s working precision. Rule nodes are restored through that codec; endpoints of the saved open rules are required to lie in `(-1,1)`.
- Requests use full descriptors, separate route/purpose/precision namespaces, `PREPARING` then `PENDING` then `COMPLETE`, and lookup refuses anything except `COMPLETE`. Cache size is 0. Panels require 1–48 aligned columns. Disk reserve is 20 GiB plus twice the record and 131072 bytes. Records above 8 MiB refuse. The worker journals a failure traceback, writes the final chain, active-leaf count, and source/copy posthashes in `finally`, and does not retry. The launcher command is one pooled `s11c_guarded_run.py` job at `--memory-gib 4 --tasks-max 32` around the existing supervisor, with the hook armed first and `automaticRetry` false.

## Mandatory runtime checks

Source clearance is not a run, and this packet has no science or READY gate. These remain unchecked because nothing was executed:

- The live census over all 544 selected rows, all 20 eligible coefficient and unit records, and all 44 exact-zero reasons.
- The preamble’s exact SymPy residuals, the six-slab meshes, and both wing refusals on all four plans.
- Formula-check gates, per-panel inner gaps, B’s actual leaf sequence, pointwise and weighted budgets, route and window comparisons, and the two finite-window control movements.
- SQLite state counts, disk admission, active leaves, and original/copy posthashes after a real launch.

Unknown numerical values stay unknown. No completed calculation is replayed here.

## Scope limits

H, flat, height/contact/PV, and slope stay pending. Ordinary J is not full native mixed. J plus the three direct primitives is not a packet action, current, inverse, or leakage factor. Both action routes share the assessed analytic Gaussian identity; the 60-digit formula comparison is an algebra and roundoff check. The native-sinh profile branch records error 0 with `empiricalRoundoffNotCertified` (`numeric.py` lines 209–215). Those indicators are empirical. The finite-window mutants do not establish all-real mutant tails. The planned worker is the single no-deadline pooled job in `launcher.py` lines 36–41: 4 GiB native/cgroup inside the 16 GiB pool, zero swap, one CPU/thread, 32 tasks, and a 4 GiB host reserve.

## Coverage

Read in this pass: `build.md`, `method.md`, `evidence-guide.md`, `prepare.py`, `geometry.py`, `numeric.py`, `worker.py`, `launcher.py`, `request-index.py`, `evidence-store.py` (header), `input-manifest.json` through the source pins, `new-contraction-definitions.json`, `numeric-source-contracts.json`, whole-tags, factor residual inputs 0–9 and factor-3’s operands and return, factorization inputs for J, Dr, Dh, and Dq, native injections for plus `packet_J`, minus `packet_J`, and plus `packet_D`, address inputs/returns for 8346, 8347, and the status region of 8348, `geometry_lib.py` through `canonical_lines`, `ctx_mp_python.py` `make_mpf` and tuple construction, and preflight `exponential_moment` / `weighted_tail`.

Not fully read: the other eighteen unit interfaces, the other explicit-zero address bodies, factor inputs 10–16, the rest of `arrangement()`, the guard body, the helper `containment` function, A-GL48 node tuples, and the manifest’s `savedInputs` list after the first entries. No shell and no scientific execution.