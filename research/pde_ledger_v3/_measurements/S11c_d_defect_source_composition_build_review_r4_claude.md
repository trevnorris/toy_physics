**NEEDS REVISION**

I read the method, guide, `worker.py`, the launcher and helpers, the manifest, and the saved retained-response, sheet, closed-density, factor and consumer files that the checks depend on. Nothing was executed. Apart from the fixes below, the build is a faithful limited-inventory emitter. The single scientific/validation fix is small.

## Blockers: scientific / validation

**S1. The cs inventory is asserted, not computed.**
- `worker.py:273` emits the literal `'csOnlyInDepth':True`. This happens before any source join, and nothing computes it.
- Method §1 says any actual cs occurrence in the source/consumer coefficients "must be inventoried rather than assumed absent".
- The only related check is `worker.py:400`. It tests that `c_s0` is absent from the source after `bind()` has already substituted it. That check is vacuous.
- A text search of the saved consumer operands finds `c_s0` (=10) only in `binding-context.json`, and no other cs-like names. The claim is probably true, but the instrument does not show it.
- **Minimum correction:**
  - Before `bind()`, collect the free-symbol names of the raw source, chemical and consumer operands that match cs-like names (`c_s*`, `*speed*`).
  - Emit the list and require it empty, or record each hit as native-bound and not tuned.
  - Replace the literal with that computed result.

**S2. The L=10 profile scaling is not joined to the native convention (same patch as S1).**
- `profile()` at `worker.py:363` hard-codes `10**n`, and line 356 hard-codes `tanh(x/10)`.
- Native `Inputs.at_source` multiplies the all-1 jets by `values['L_W']**n` (`native-c2.py:265`). Method §5 asks for this to be checked. The worker only emits the native text.
- `physical-input.json` has `L_W=10`, so the value is consistent today.
- **Minimum correction:** add `require(context['numeric']['L_W']==10)` and derive the scale from it.

## Tooling only

**T1. The worker does not join its own argv to the gate.**
- `verify_gate` (`worker.py:63-90`) and `main` (`worker.py:987-991`) never compare `--out`, `--inputs` or `--gate` with `gate['command']` or `gate['outputDirectory']`. Only the launcher does (`launcher.py:31,36-41`).
- A direct invocation under the guard with a fresh `--out` would bypass the hook-first and one-run constraints.
- **Minimum correction:** `require(args.out.resolve()==Path(gate['outputDirectory']))`, and require the argv tail to equal the gate command tail.
- Also require `gate['buildReviewRecord']==manifest['reviewRecordWillBe']`.

**T2. The launcher's hook-arming window is 5 seconds.** Raising `RuntimeError` after 5 seconds leaves the watcher running. This is optional.

## Verified as sound

- **Typed factors and double resolvent:**
  - The three typed objects are kept separate, and `Rprod = qi·qo·E` is a pure algebraic check (`worker.py:442-446`).
  - Isolated-factor, closed-reference and normal-jet pairs are inherited zero returns joined to current operands, never compared as "density − Rprod" (`worker.py:448-471`).
  - The bare linear factor is a single `R(qi)R(qo)`.
  - `Dwhole` enters once, unmultiplied, with signature `(l,k,3,cs,1/5,1/10,1,10)`.
- **Independent-D adapter:**
  - It is declared and runs only the five pinned `build_face` assignments.
  - Its argument order and refusals match the native `kernel_apply` call at `native-c2.py:586`.
  - The normal factor is `i·f·qo·P` for both faces, and the excluded (2,1) height term equals `−H·i·f·qo·P`.
  - The whole trace matrix is injected and `(T·Δ)[1,2]=0`.
  - The native reference location is overwritten by the actual assignment, then compared.
- **Source census and joins:**
  - The broad `delta_p`/`d_w_` census is conservative: it cross-checks substring counts against literal AST hits and refuses unsupported constructors.
  - All four slot ablations and the affine reconstruction are present.
  - Per-face native chemical, density and velocity-normalization joins feed a simultaneous stage-2 substitution that must equal the saved combined source.
  - The success flag at line 391 follows those joins.
- **Grades and coverage:**
  - The quotient recurrence is correct for the bivariate rectangle, and the unsplit remainder and excluded pure grades are kept.
  - Exactly 16 triples are enumerated per row, face and slot, with zeros recorded (`worker.py:763`).
- **Controls:**
  - Supports are right: flat `k=l`, constant-source `k=p`, constant-consumer `r=l`.
  - The depth points `q = 2, 3/2, 25/26` at `cs²=10/7` check out and sit on the positive branch.
  - Source-10 and consumer-10 fields are certified nonconstant by an exact tanh polynomial, with no transform-value claim.
  - Unknown applicability never counts as applicable (`select_certified_candidate`).
- **Write-before-guard and strict failure:** `Journal.zero` and `exact_nonzero_number` emit before they require. `except BaseException` preserves the failure and posthashes.

## Coverage limitations

- **Actual-row substitution:** it is an exact grade-series identity with opaque response tokens, so it cannot detect a wrong response operand. The per-address full-factor residual (`worker.py:701-723`) carries that identity. Together they are adequate for an inventory.
- **c1 velocity and shifts:** these are not text-joined to the native constructor. Only the chemical and density leaves are, so velocity relies on the saved `nativeVelocity` and the inherited normalization.
- **Off-delta support:** it shows the transform is not supported at {0}, not that it is nonzero at 1/2 or 4/13. The worker correctly claims nothing more.
- **Zero-profile reductions:** they are formal tag-zero checks. The constant-height reduction is `NOT_COMPUTED`. The constant-source map at lines 811-815 is emitted but never checked. This is sufficient for the bounded inventory and says nothing about global operator applicability.
- **Control movements:** they are formal tag-coefficient movements. The direct omit/double movement is `±cp·sourceval` and excludes the tag.
- **Fourier normalization:** the 1/(2π) factor is definitional. Only the saved powers `(-3,-3,0)` and the plain edge measure are checked.
- **Scope limits:** source units are inferred only. Composition is limited to certified `k,l∈[-3,3]` and `cs∈[1,2]`. Global composition, test space and the grazing limit stay UNRESOLVED.

## Runtime evidence still needed

I could not verify these without running. All of them fail closed, so they cost a run but do not corrupt results:

- `sympify` of native constructor text (`worker.py:305`)
- SymPy's derivative form of tanh for the tau substitution (`worker.py:859`)
- structural equality of the Piecewise branch conditions (`worker.py:586-589`)
- sign decisions on the √3 constants (`worker.py:347-353`)
- the saved-operand symbol-name censuses (`worker.py:596,607`)
- namespace completeness for the exec'd assignments (`worker.py:499`)
- `sp.cancel` on the 4 GiB `RLIMIT_AS` token rows (`worker.py:783`)
- constancy of the consumer00 pressure coefficient and uniqueness of the control address (`worker.py:937,916`)
- pooled-guard containment values
- the hook arming handshake