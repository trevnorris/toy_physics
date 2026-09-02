---
schema_version: 2
id: source-paper-4d-3pn
title: Paper 4D 3Pn
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/4d_3pn/paper/4d_3pn.tex
content_owner: ai_generated
last_updated: '2026-09-01'
generated_from_commit: 5d3b78e8e39c0ea9dcc6c37b8f940d36d1183c0b
source_kind: paper
source_unit:
  id: paper-4d-3pn
  shape: file
  entrypoint: research/4d_3pn/paper/4d_3pn.tex
  unit_digest_sha256: 5d62367d03c9d28dd80f5eba442a7bf65ac91cdf6db692a1b1ef970bae4a49dd
  members:
  - path: research/4d_3pn/paper/4d_3pn.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 4f67e082f79ce08c10547e9c7dba209a652b8d48
    blob_size: 320909
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-4d-3pn--conservative-3pn-scope — Conservative two-body 3PN assembly

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper assembles the nonspinning conservative two-body sector through order \(c^{-6}\) in a fixed ADM chart, using the reduction and closure hierarchy carried from the program’s lower-order work. It presents an exact algebraic assembly after that hierarchy, its targets, and its chart are fixed—not an assumption-free derivation from a solved moving-throat PDE. Spin, radiation reaction, dissipative normalization, many-body response, nonperturbative completion, and higher-PN sectors are excluded.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-nonclaims}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-taxonomy}`

## Source-unit map

- Entry point and sole member: `research/4d_3pn/paper/4d_3pn.tex`
- Role: `primary`
- Read mode: `semantic`
- Git mode: `100644`
- Object type: `blob`
- Blob ID: `4f67e082f79ce08c10547e9c7dba209a652b8d48`
- Blob size: `320909` bytes

## Key statements

### source-paper-4d-3pn--one-body-repair-ledger — Strict one-body gate

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Expanding the isotropic Schwarzschild worldline action through 3PN and matching the paper’s denominator, static, and new self-invariant slots fixes
\[
(\mu_{\rho,3},d_3,s_{24})
=
\left(\frac14,-\frac{45}{4},-\frac1{16}\right).
\]
This establishes that the 3PN one-body sector is not a one-parameter continuation of the 2PN denominator-style closure. The algebraic matching is exact within the paper’s basis, while its one-body repair packaging is classified as protocol closure rather than a moving-throat PDE derivation.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:onebody-target}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:onebody-repair}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-taxonomy}`

### source-paper-4d-3pn--cubic-legendre-compiler — Exact cubic compiler

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For \(L=L_0+\varepsilon L_1+\varepsilon^2L_2+\varepsilon^3L_3\) with quadratic \(L_0\), the paper derives
\[
H_3=-L_3(v_0)+A_0^TM^{-1}B_0
-\frac12A_0^TM^{-1}C_0M^{-1}A_0.
\]
With the lower-order ledger fixed, a new generic-frame residual obeys \(\Delta H_3=-\Delta L_3(v_0)\); in the fixed momentum basis its compiler is \(C_{\rm gen}=-I_{50}\). COM-blind directions therefore cease to be Hamiltonian-blind after the ADM chart and full Hamiltonian target are fixed.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:onebody-legendre}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:compmass-uniqueness}`

### source-paper-4d-3pn--hamiltonian-first-reduction — Hamiltonian-first COM recovery

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Direct COM substitution into the imported generic-frame ordinary Lagrangian does not reproduce the exact reduced ordinary target. Applying the cubic Legendre lift before COM reduction does reproduce the target Hamiltonian, selecting the generic-frame representative through the Hamiltonian-first route rather than naive ordinary-Lagrangian reduction.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:compmass-naive-com}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:compmass-hamiltonian-first}`

### source-paper-4d-3pn--grouped-p2-middle-block — Grouped real \(P_2\) closure

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The minimally demoted grouped-\(P_2\) front end has rank three and imposes six slot relations violated by the target. Enlarging it to
\[
(T_{20},T_{21},T_{22},S_{20},S_{21},S_{22},V_{20},V_{21},V_{22})
\]
gives a nine-dimensional compiler with \(\det M_{\rm mid}=-4/27\), and the paper supplies the unique coefficients reproducing all nine residual slots. The compiler result is exact for the stated family; choosing that richer local constitutive family belongs to the paper’s protocol-closure layer.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:grouped-minimal-failure}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:grouped-richer-family}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:grouped-gr-coeffs}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-taxonomy}`

### source-paper-4d-3pn--geometry-static-completion — Unique geometry-side static completion

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

After grouped-\(P_2\) closure, the remaining static coefficient is
\[
\Delta l_{15}^{(g)}
=
\frac{\nu(408\nu^2+1232\nu-2080+63\pi^2)}{96}.
\]
Here \(p\equiv m_A\) and \(q\equiv m_B\) are mass shorthands, not electric charges. The identity
\[
\nu(p^3+q^3)=(1-3\nu)(p^2q+pq^2)
\]
removes the apparent body-local \(P_0\)/pair split parameter and yields the unique generic-frame counterterm
\[
\Delta L_{3,\rm ct}^{\rm scalar}
=
\frac{G^4pq}{96r^4}
(408\nu^2+1232\nu-2080+63\pi^2)(p^2q+pq^2).
\]
The sigma-collapse is exact mass algebra; interpreting the surviving carrier as geometry-side is protocol closure.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:dict-firewall}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:geometry-u-block}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:geometry-residual}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:geometry-sigma-collapse}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-taxonomy}`

### source-paper-4d-3pn--pure-kinetic-collapse — Free-kinetic compiler image

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The generic-frame free 3PN Lagrangian has separate-body octic terms but no mixed comparable-mass kinetic interaction. Comparison of its naive COM seed with the exact free two-body COM Hamiltonian through the compiler gives
\[
\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16},
\qquad
\Delta l_2=\cdots=\Delta l_5=0.
\]
The kinetic remainder is therefore compiler compensation for universal free kinematics, not an additional throat-response datum.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:kinetic-generic}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:kinetic-com-target}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:kinetic-collapse}`

### source-paper-4d-3pn--three-lane-assembly — Final three-lane theorem

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Within the declared closure hierarchy and fixed ADM chart, the comparable-mass residual is
\[
\Delta\hat L_3^{\rm GR}
=
\Delta l_1v^8+L_{P_2}^{\rm mid}+\Delta l_{15}^{(g)}U^4.
\]
The three lanes are the free-COM compiler image, grouped real \(P_2\) middle-block closure, and geometry-side static completion. Once the protocol-closure ingredients are inserted, the paper’s full-assembly consequence is exact algebra and satisfies \(H_3^{\rm derived}=H_3^{\rm ADM}\); this does not promote the constitutive ingredients to PDE-derived results.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:final-com}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:final-adm}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:final-main}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-taxonomy}`

## Computed evidence represented by the source

### source-paper-4d-3pn--computed-evidence-boundary — Named referee archive absent from unit

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit

The paper names a dual-CAS referee ledger: `3pn_referee_master_sympy_audit.py`, 18 staged SymPy audits with pinned SHA256 digests, and `mathematica/wl_notes.txt` with standalone Wolfram Language audits and embedded `Output:` transcripts. None is a member of this prepared unit, so their commands, literal outputs, integrity records, and comparison chain cannot be assessed here. The paper also limits the archive to an algebraic backstop within its hierarchy, not a solved moving-throat PDE theorem.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:verification}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:verification-master}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:verification-stages}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:verification-does-not}`

## Assumptions, exclusions, and open questions

### source-paper-4d-3pn--lower-order-no-retune — Frozen lower-order viability guard

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The construction freezes
\[
\kappa_\rho=1,\quad n=5,\quad
\kappa_{\rm add}=\frac12,\quad
\kappa_{\rm PV}=\frac32,\quad
\beta_{\rm 1PN}=3
\]
as imported Newtonian–2PN viability data. A 3PN mismatch may not be repaired by retuning them; inherited parity-even wake data are likewise treated as frozen branch data rather than higher-order fit parameters.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:dict-carry}`

### source-paper-4d-3pn--moving-throat-origin-open — Three distinct moving-throat tasks remain open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper leaves three distinct tasks for the true moving-throat branch: derive the low-frequency constitutive origin of the conservative grouped-\(P_2\) data and geometry completion; independently test conservative isotropy first at \(O(\omega^2)\) through \(a_2=b_2=0\), then at \(O(\omega^4)\) through \(a_4=b_4=0\) and the stated single-pole relation; and, only after those gates, realize the passive/outgoing quadrupole normalization
\[
\hat m_0^{\,2}\Gamma_5=\frac{2G}{5c^5}.
\]
The isotropy tests are conservative branch tests, not part of the unresolved dissipative normalization itself.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:interface-isotropy}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{eq:interface-first-isotropy-gate}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{eq:interface-second-isotropy-gate}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:interface-open}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{eq:interface-final-normalization}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:fixed-open-open}`

## Revision and supersession relationships

### source-paper-4d-3pn--lower-order-carry-forward — Extension of the carried program

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper carries forward the earlier 4D reduction backbone and frozen Newtonian, 1PN, and 2PN conservative ledger. It turns the earlier 2.5PN narrowing to the grouped real quadrupole route into a conservative 3PN carrier test, while explicitly stating that its conservative theorem does not replace the 2.5PN program or close that program’s outgoing-normalization condition. It does not declare wholesale supersession of the earlier treatments.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:intro-motivation}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:interface-25pn}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:interface-open}`

## Related topics and scripts

No related topic pages or script domains were supplied by the task.