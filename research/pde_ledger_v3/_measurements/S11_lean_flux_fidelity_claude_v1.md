# Claude independent fidelity report — fixed F1–F4 packet

Packet: `5945c5b4b9d980613e7ff1c31261547a50f559799c14250173e325d0e02584ae`

Session: `618fcd0e-3ab3-4126-a752-308ca46acbe5`

Model: claude-opus-5-5

Raw-response SHA256: `91bcac4b234d5be47d5e579c996518abbf19f10684c0fd933a79eaf5c65d8adb`

The following is the verbatim reviewer result. Local execution was not independently repeated by the reviewer.

---

# Independent fidelity review: S11 scattering current/flux F1–F4

**Verdict: CLEAR for F1–F4.** I found no blocking findings. There are six optional improvements, listed at the end. This is one of the two non-author reviews the policy requires; the author should record its disposition.

**How this review was done.** I only read files. I did not run the compiler, the instruments or any Python, and I had no shell, so I could not recompute any hash. What I checked instead is that the hashes agree with each other across `MANIFEST.json`, `S11_lean_flux_contract_checks.json` and `S11_lean_flux_validation.json`. They do for all ten contract inputs, the four compiled objects, the four native sources and the evidence files. The build and control results are the author's recorded run8 evidence; I did not reproduce them.

## 1. F1: the current form and flux (`Current.lean`)

- **Definitions** (`Current.lean:13-15`):
  - `pair J x y = star x ⬝ᵥ (J *ᵥ y)`. `star` acts entry by entry, so this is exactly conj(x)ᵀ J y.
  - `flux = (pair J x x).re`.
  - `pullback J C = Cᴴ * J * C`.
  - The index types `n, m, k` can be any finite type. No positivity, diagonal, Hermitian, transverse-only or nondegenerate condition appears in any variable or type.
  - `empty_flux` (`Controls.lean:45`) covers the zero-dimensional case, and it has a positive control.
- **`pair_coordinates`** (`:17`) is Σⱼ Σᵢ conj(xᵢ) Jᵢⱼ yⱼ. It covers every matrix entry.
- **`flux_add`** (`:50`) keeps both cross contractions, Re(pair x y + pair y x), and holds for every J, including non-Hermitian ones. `flux_add_iff_cross_zero` makes additivity exactly equivalent to that real sum being zero. Nothing silently assumes the cross terms vanish.
- **Reality.** `hermitian_flux_real` (`:44`) proves the raw pairing has zero imaginary part, and only under `Jᴴ = J`. `nonhermitian_imaginary_witness` (J = i, x = 1, so Im = 1) shows the premise is needed, and `current_reality_premise_mutant` (claiming Im = 0) is rejected.
- **What `flux` does with a non-Hermitian J.** Taking the real part keeps only the Hermitian part of J. The docs say this is not a certificate for arbitrary imaginary defects (`SCATTERING_FLUX_FIDELITY.md:10-15`), which is correct as written. Optional item O1 would state it as a theorem.

The signs, conjugations and cross terms all match shared physics §3a and §3c. Section 3a says "cross terms inside a non-diagonal current block are retained"; §3c defines Jᴴ⁽¹⁾ and Jᴴ⁽²⁾ as sums of B[a,b]+B[b,a] pairs, which is exactly the structure `flux_add` keeps.

## 2. F2: pullback and covariance

- **`pair_pullback` / `flux_pullback`** (`:60`, `:64`) hold for all complex, possibly rectangular, C. That pins the definition to Cᴴ: a transpose version would not be provable for complex C.
- **`pullback_comp`** (`:69`) is correct: composing C then D gives the pullback by C·D.
- **`scattering_flux_covariant`** (`:75`):
  - The quantifiers are right: for every J, S, Cout, Dout, Cin, the hypothesis `Cout * Dout = 1` implies the result for every x.
  - `Cout` and `Dout` are typed as square n×n matrices. Over ℂ, a right inverse of a square matrix is also a left inverse, so the output side really is a complete invertible basis change.
  - `Cin` may be rectangular (m×k), which models representing only some inputs. The docs are accurate that this is "for the represented input" and not a channel census (`FIDELITY.md:24-30`).
  - The map S ↦ Dout·S·Cin is the correct transformation when old coordinates = Cout·new on the output side and old input = Cin·new on the input side.
  - Incoming-flux covariance is simply `flux_pullback` applied with Cin. There is no combined statement that the ratio is invariant (optional O2).

The maps and quantifiers are faithful.

## 3. F3: orientation, fraction and case coverage (`Balance.lean`)

**Signs.** `outward` = −1 for left and +1 for right, and `incident e j = −outward e · j`. `outgoing_eq` gives −l + r, a signed sum rather than a sum of absolute values. This matches §3a (s₋ = −1, s₊ = +1, J_in = −s·𝓙ₙ, J_out = Σₑ sₑ·𝓙ₙ,ₑ) exactly. It also matches the engine (`S11c_d_mixing_scattering_sympy_audit.py:4738`, where `orientation` is LEFT −1 / RIGHT +1 and `outward = orientation*current`) and `open_metrics` (`S11c_d_continuum_currents.py:123-124`, outgoing `s*a`, incoming `-s*a`).

**Where the sign is applied.** The native code multiplies the sign into the metric blocks and then block-diagonalises them (`diagonal_series`). Lean applies the sign to scalar fluxes afterwards. These agree by linearity, and because a block-diagonal metric has no cross-end terms. That identity is not proved in Lean and not executed natively; it rests on source inspection, as declared (optional O3).

**Fraction.**
- `fraction` returns `none` exactly when the denominator is 0 (`fraction_undefined_iff`).
- The zero-denominator mutant, `fraction 7 0 = some 0`, is precisely Lean's default x/0 = 0 behaviour, so rejecting it is meaningful.
- `null_flux_nonzero_amplitude` keeps zero-flux nonzero vectors in the domain.
- `flux_sign_coverage` is proved trichotomy. Disjointness of the three cases is asserted in the docs, not proved (trivial; optional O4).

**Bounds.**
- `fraction_nonneg` needs both 0 ≤ num and 0 < den.
- `fraction_le_one_iff` needs 0 < den.
- Both are stated on `num / den` rather than on the `Option` value. Combined with `fraction_defined`, that is adequate.
- The witnesses `negative_fraction_witness` and `fraction_above_one_witness` show that neither bound holds without its premises.

## 4. `conditional_balance` (`Balance.lean:33`)

The premises are the explicit balance equation and a nonzero incoming flux. The result keeps the defect term. `conservation_requires_zero_defect` shows that a normalised total of 1 is equivalent to defect = 0. Neither S, unitarity, a physical conservation law, omitted channels nor bound-state capture appears anywhere; the docstring and `FIDELITY.md:41-45` say so. This matches §3b and §3c: bound spectral overlap is never added into a continuum current, and no capture probability is formed.

## 5. The paired controls (`S11_lean_flux_contract_check.py`, `verify.adjudicate`)

**Adjudication is strict.** A rejection only counts if all of these hold:
- the compiler exits with status 1;
- the output contains none of the bad-instrument patterns (unknown identifier or module, heartbeats, memory, sorry, failed to synthesize);
- there is exactly one error, it is in `contract_control`, and its message contains `⊢ False` exactly once;
- there are no warnings.

A positive requires exit 0 with no `error:` or `warning:` text. The recorded outputs satisfy this: every mutant is exactly `error: unsolved goals\n⊢ False`, and every positive log is empty. I cross-checked the counts: 11 rejections, 15 positives (11 paired + 4 extra), 4 builds and 1 native check, for 31 records.

**Why a rejection is meaningful.** Each mutant is closed with `norm_num only [W]`, where W is the corresponding kernel-checked witness. Simp steps are equivalences, so reaching `False` shows the mutant statement is actually false, not merely unproved. Each wrong value encodes a specific wrong formula:

| Control | Rejected value | Wrong formula it represents |
|---|---|---|
| interference | 2 | cross terms dropped (same value as `diagonal_flux`) |
| conjugation | −1 | no complex conjugate |
| basis_metric | 1 | stale metric (basis change ignored) |
| orientation | 4 | unsigned sum |
| signed_incident | 3 | absolute value |
| zero_denominator | `some 0` | Lean's total division, 7/0 = 0 |
| normalization | `some 2` | inverted ratio |
| negative_domain | `some (1/2)` | sign of the denominator dropped |
| no_unconditional_unit_bound | `some 1` | ratio clamped to 1 |
| no_unconditional_conservation | 1 | conservation assumed |
| current_reality_premise | Im = 0 | reality assumed without Hermiticity |

The balance and zero-defect positives show that the premises can be satisfied, so the theorems are not vacuous.

**Limits of the controls:**
- They are instance-level value controls that depend on the canonical witnesses. They are not independent re-derivations, and they are not source-replacement mutations (which the docs correctly do not claim). The phrase "independent fresh true/false statements" (`SCATTERING_OBSERVABLE_COVERAGE.md:41`) overstates this (O5).
- The `basis_metric` control uses the real scalar C = 2, so in Lean it cannot tell Cᴴ from Cᵀ. That distinction is covered by the universally quantified `pair_pullback` proof and, natively, by the complex nonunitary fixture. It is not blocking, but a complex-C Lean witness would close it directly (O1).
- No control drops `hInv` from the covariance theorem or `hd` from the balance theorem (O2).

For the compact claims of this contract, the controls are meaningful.

## 6. Native identification (`S11_lean_flux_source_check.py`)

**Only the intended functions run.** The instrument extracts `multiply`/`quadratic`, `adjoint`/`gram` and the engine's `RectangularModeJets.multiply` from the AST, removes decorators, and runs them without importing their modules. The sources are exactly what the Lean definitions claim to mirror:
- `adjoint` is `v.conj().T`;
- `gram` is adjoint·metric·series;
- `quadratic` is adjoint(a)·b·a.

**The fixture values check out.** I recomputed them by hand:
- q(a, j) = 3.
- The wrong formulas all differ from it: transpose = −1, diagonal-only = 5, jᵀ = 7, stale metric = 10 vs 3.
- The complex nonunitary c reproduces the congruence cᴴ·j·c.

**The adjoint mutation is live.** Removing `.conj()` from `adjoint` is checked with `reject`. If the transformer had silently done nothing, that check would have failed, so its pass shows the mutation actually took effect. All values are Gaussian integers, so exact array equality is legitimate.

**Declared limits are honest.** Orientation is confirmed only by a string check on `open_metrics`. That check does not confirm that `ends[...]['orientation']` in the saved packet is the engine's −1/+1. Nothing in the packet claims a scientific output, numerical current reality, a channel census, a scattering solve or higher grades.

**Observation for the later bookkeeping increment (not F1–F4 evidence).** The native `quotient` (`S11c_d_continuum_currents.py:38-45`) divides by `np.diag(denominator[0])` without a zero guard, and it does not take a real part. Lean's `fraction` therefore has no checked native counterpart. The F1–F4 docs never claim one, so this does not block (O6).

## 7. Provenance, build evidence, axioms and limits

- **Axioms:** all 34 `#print axioms` lists contain only `propext`, `Classical.choice` and `Quot.sound`. The four sources contain no `axiom` or `sorry`, and builds ran with `-DwarningAsError=true` and empty logs.
- **Theorems not in the audit list:** `pair_add_matrix` and the helper lemmas were not audited directly. They are either used by audited theorems or compiled cleanly under warnings-as-errors.
- **Pins:** Lean 4.33.0, 15 clean package pins, and 7 direct Mathlib source/object pairs. The transitive Mathlib cache is a pinned baseline, not a rebuild; that is declared and acceptable, since the kernel check of these small files does not depend on rebuilding Mathlib.
- **Resource receipt:** peak about 250 MB against a 2 GiB limit, no swap, no out-of-memory events.
- **Failure handling:** failed earlier runs (runs 6–7, the missing cache) are recorded as failures, not counted as rejections.
- **Python/NumPy versions** were not captured. That weakens reproducibility of the native leg only. It does not invalidate the claim, because the fixtures are exact integer arithmetic.
- **D5B registration** is correctly described as catalog-only (`INSTALL_D5_BULK_REGISTRATION.json`), with no new proof replay claimed.

None of these limits invalidates the bounded claim, and none is silently upgraded.

## Findings

**Blocking: none.**

**Optional improvements:**
- **O1.** Add a Lean witness with a complex basis map (e.g. C = i or 1+i) to `basis_flux`, so a transpose-for-adjoint error is distinguished inside Lean. Optionally also add a lemma stating that `flux J x` equals the flux of J's Hermitian part, to make explicit what taking the real part discards.
- **O2.** Add a combined theorem that the fraction is invariant under both coordinate maps (Cout and Cin). Add premise-drop controls: a non-inverse pair (e.g. Cout = Dout = 2) and incoming = 0 in the balance theorem.
- **O3.** Add a one-line Lean lemma: the flux of the block-diagonal sign-weighted metric (−J_L ⊕ J_R) equals `outgoing (flux J_L x_L) (flux J_R x_R)`. This turns the matrix-level vs scalar-level sign link from source inspection into a proof and shows where cross-end terms vanish.
- **O4.** Prove the three sign cases are pairwise disjoint (currently asserted in `COVERAGE.md:44-45`).
- **O5.** Reword "independent fresh true/false statements" to say the controls are fresh statements decided by rewriting with the canonical witnesses.
- **O6.** In the F1–F4 fidelity record, say explicitly that `fraction` has no native counterpart checked, and carry the unguarded, complex-valued, diagonal-only `quotient` into the bookkeeping increment's contract.

**F1: CLEAR. F2: CLEAR. F3: CLEAR. F4: CLEAR.**
