I'll review only the indexed packet files, starting with the index and the listed source documents.The packet is indexed. Next I’ll read the required first files, then the remaining packet sources.The evidence file is large. I’ll keep reading it in sections while loading the implementation directive and response script.I’ll inspect the remaining evidence in smaller slices and start reading the proposed worker.The worker’s first half is loaded. I’ll continue through the transform, domain checks, and remaining packet files.The worker is fully loaded. Next I’ll read the remaining indexed reports, checkpoints, and supporting source files.The transform rules are in view. I’ll now inspect the remaining evidence sections and the accepted checkpoints.I’ll pull the rest of the evidence keys and the scripts that define term structure, poles, and joins.The local/nonlocal split still needs confirmation. I’ll read the case records and the scripts that define term values and the PV/delta action.The native cells end around here. Next I’ll read the case records and how forcing terms are assembled.The PV/delta assembly and the rest of the case records are the remaining load-bearing pieces.One nonlocal operand shows the opposite plane factor. I’ll confirm whether that is the same case and inspect the seed/lift records.**NEEDS REVISION**

The local dictionary, Fourier mass, and the two distributional identities are the right objects for these thickness sources. The worker still peels one seed-plane factor that some actual substituted nonlocal integrands do not carry. That is enough to stop the stage from claiming a source-specific outgoing action.

## Substantive findings

### 1. One seed-plane factor is assumed; some saved nonlocal terms carry the opposite factor

`construct_case` in `S11c_d_localized_thickness_response.py` takes `k0` from `blocks[index]['k']` and checks `f['seed']==exp(I*k0*z)*residue`. For block 17 that seed is uniformly

`exp(-sqrt(555)*I*s11cdNormalPosition/30)`

so `k0=-sqrt(555)/30`. Local forcing in the same case uses that same factor (`localForcingPart` in `saved-source-evidence.json`).

`transform_nonlocal` always peels

- rule 1: `I*kout*z + I*(k0-kout)*zs`
- rule 2: `I*kout*z + I*(k0-p)*zs`

That is the kernel character times **this** `k0`.

The substituted inspection operands are not uniform. Example, still in `block17-grade01`:

- address `(row 3, sourceColumn 0, frameColumn 0, term 4)`: `exp(-sqrt(555)*I*zs/30)` — matches the seed
- address `(row 3, sourceColumn 2, frameColumn 0, term 4)`: `exp(+sqrt(555)*I*zs/30)`

The second is `exp(-i k0 zs)`. The derivative prefactors in that integrand match a plane factor `exp(+i|k0|zs)` with amplitude `-1/10` (seed entry `(2,0)`), and they do not match `d/dz` of the printed seed `exp(i k0 z)`.

On those terms the stated peel leaves a leftover `exp(±2 i|k0| zs)`. `profile_transform` then hits `unrecognized localized source factor` (or rule 2 hits `plane source leaves a multiplier only`). The worker never transforms the complete source.

`source_checks` joins the seed matrix and `direct == -forcing['total']`. It never joins each `term['value']` leftover `zs`-plane factor to that seed.

**Smallest repair:** After removing only the kernel output character (`exp(i k_out(z-zs))` or `exp(i k_out z - i p zs)`), save the leftover `zs`-plane factor of every substituted term and require it to equal the saved seed factor `exp(i k0 zs)`. If a term disagrees, stop with that operand and do not emit an action candidate. Do not add a second dictionary until that join holds.

Until that join exists, the two native rules do not implement the supplied nonlocal operands.

### 2. The rest of the two rules is right for terms that do carry the seed factor

Independent identities, same Fourier convention as the implementation (`Fhat(k)=∫ e^{-ikz} F dz`, `mass=2π`):

- Local: `F(z)=P(tanh(z/L)) e^{i k0 z}` with `P(T)=(1-T^2)Q(T)` gives `L B0(L(k-k0))` times the recurrence `p_{n+1}=(n p_{n-1}-i q p_n)/(n+2)`, `p_0=1`, `B0(0)=2`. That is `basis_polynomial` / `profile_transform` / `transform_local`. Sign is `F=-forcing['local']` and `-term['value']`, matching `direct==-forcing['total']` and zero lift.
- Rule 1: `c ∫ dk_out dzs e^{i k_out(z-zs)} M(k_out) H(zs)` → `2π c M(k) Hhat(k)`. Bound-variable order, whole-line limits, and `mass*coefficient` match the saved `(k_out, zs)` carriers.
- Rule 2: `c ∫ dk_out dp dzs e^{i k_out z-i p zs} M(k_out,p) e^{i k0 zs} hhat(k_out-p)` → `(2π)^2 c M(k,k0) hhat(k-k0)`. Nested `∫ 10(1/2-tanh^2 ξ/2) e^{-i 10(k_out-p)ξ} dξ` is already dimensionless; `scale=1` with transfer `L(k-k0)` keeps the saved `L` in the exponent. Collapse is `p→k0`, `k_out→k`.

`forcing['terms']` are nonlocal only (`apply_action` in `S11c_d_two_asymptote_forcing_continue.py`). Local matrix plus those terms is the total. No double count.

### 3. Branch/denominator checks justify smooth multiplication and rapid decay, on the seed-plane terms

`amplitude_domain` is the right local test: one radicand `-a-b p^2` with `a,b>0` (`-4-100 p^2` in these cells), `sqrt=i r`, `r=sqrt(a+b p^2)≥sqrt(a)>0`, integer inverse powers, same-sign real or imaginary coefficient lists, fail-closed on mixed signs or explicit `k` denominators.

That covers:

- collapse at real `p=k0` (`native-input-multiplier-domain-*` on `M(k_out,p)` before substituting `p=k0`)
- `Fhat` at both inherited poles, including removable `B0(0)=2` when `kj=k0`
- Schwartz decay of `B0(L(k-k0))` times polynomial growth in `(k,r)`

No extra global operator theorem is required. The smoothness claim applies to the complete `Fhat` only after every nonlocal term actually enters that class (finding 1).

### 4. Smooth-source PV/delta action needs no extra diagonal kernel

The accepted kernel is a Fourier multiplier

`G(z,z') = PV ∫ e^{ik(z-z')} R(k) dk/(2π) + Σ_j e^{i kj(z-z')} D_j/(2π)`

with `Ne(z,z')`, decaying tails, and zero polynomial contact (`S11c_d_outgoing_prescription_unlimited.py`, `prescription-domain`). For a Schwartz `Fhat`,

`V(z) = PV ∫ e^{ikz} R(k) Fhat(k) dk/(2π) + Σ_j e^{i kj z} D_j Fhat(kj)/(2π)`

is the same multiplier acting on these sources. `construct_case` copies `deltaCoefficient`, `principalValueExclusionIntervals`, and `fourierMass=2π`. No pointwise `G(z,z)` is used.

That extension is justified once `Fhat` is the complete smooth rapidly decreasing transform. It is not justified while finding 1 remains.

## Controls, joins, persistence

Sufficient for this limited ingredient, after finding 1 is closed:

- source/domain/unit joins in `construct_case`
- remainder and endpoint checks in `profile_transform`
- independent `profile_probe` quadrature of used `sech^2 tanh^n` orders
- actual native-term omission (`transformed_terms[0]`), not the old density mutation
- new-source `P(R Fhat)-Fhat` at `k=0,±1/L` and pole null residuals (implementation checks, not observable accuracy)
- `inverseEntryUnits+forcingUnits=endFieldUnits` (equivalent to the Fourier unit identity because spectral measure cancels)
- persist-before-guard, failure stack, no automatic retry, ordinary 900s/840s containment

`equation_probe` tests the saved inverse against whatever `Fhat` was built. It does not certify the transform.

## Optional (not a further review cycle)

- `'forcingSmoothAndRapidlyDecreasingFromSavedCertificates': True` is a typed flag after `require`; the certificates already live in the `amplitude_domain` returns.
- One quadrature of an actual local envelope against the assembled closed form would catch a `remove_phase`/`Poly.div` error the basis probes miss.
- Unmeasured cost: 36 nonlocal terms per case, nested `expand_complex`/`Poly(EX)`, and 80-digit `fixedInverse*Fhat` at three points. Timeout fails closed; that is acceptable, not a duration exception.

## Scope left open

Full response/Green/FORM/A11/A12, density/mixed grades, physical current normalization, retarded equivalence, radiating coverage, and the historical end-lift literals remain open on the unchanged evanescent slice.