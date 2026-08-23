# Measurements behind S11b_card_repoint_directive.md (rule 2)

Every artifact claim in the directive is a reading of the committed card, step record, decision list, or
macros. Commands run from `research/pde_ledger_v3/`.

## Claim: the current card is the A/B version needing re-point (source records, transverse, verification)
```
$ sed -n 's/.*//;p' paper/steps/S11b_interface_coupling_law.tex   # (full read; key lines:)
# "The only source records for this card are \StageFile{...S11bA_interface_response.md} and
#  \StageFile{...S11bB_interface_assembly.md}."                         <- re-point target (change 1)
# "\rho_{\mathrm{br}}^0\omega^2=\mu_R k^2, \operatorname{Im}\omega=0"    <- unconditional (change 4)
# "\paragraph{Verification scope.} ... their \texttt{VERDICT: PASS} strings ..."  <- old framing (change 5)
# "yields \textbf{ten} independent quadratic invariants" + the 5 new terms K_L(∇·u)²,...  <- (change 6)
# admissibility region uses \Lambda_A^0, \Lambda_V^0, \Lambda_X^0 with Λ_A^0 never defined <- G12(c)
# "\paragraph{Departure --- the finite-memory velocity channel...}" with no Λ_p⁰=0 qualifier <- G12(c)
```
(card baseline: `paper/steps/S11b_interface_coupling_law.tex`, 12334 bytes, dated Aug 4.)

## Claim: the suppressed field is `Verification` (keep the plain-paragraph handling)
```
$ sed -n '18,24p' paper/macros.tex
\newcommand{\stagefield}[2]{%
  \ifnum\pdfstrcmp{#1}{Verification}=0\relax
    \ifshowstageverification
      \paragraph{#1.} #2
    \fi
  \else
    \paragraph{#1.} #2
```
⇒ `\stagefield{Verification}{...}` prints only when `\ifshowstageverification` is true (off by default); the
card uses a plain `\paragraph{Verification scope.}` to stay visible. Preserve that.

## Claim: G12(c)/(d) are the owed card items
```
$ sed -n '122,127p' directives/S11b_unified_decisions.md
- **(c)** The two owed card items (`Λ_A⁰` used undefined; the dropped `Λ_p⁰=0` qualifier on the
  finite-memory departure) are fixed when the card is re-pointed.
- **(d)** B's **uncarried** background-flow (convective) correction is recorded as a standing scope
  limit, ⛔ not silently dropped again.
```

## Claim: the step record's corrected physics (the card must match these)
Source of truth = committed `steps/S11b_interface_coupling_law.md` (`8ddccb74`):
- Transverse: "stable real-frequency mode (Im ω=0) only where μ⊥=μ_R+μ_S/2 ≥ 0 … μ⊥<0 gives a growing root;
  only DISSIPATION is unconditionally 0." (change 4)
- Slice map: both engines emit `Λ_p⁰ = −Λ_A⁰/ρ_m` (`ZPERM_SLICE_MAP`; step record §"cross-engine"). (change 2)
- Adjudication: "No compared object is a physics contradiction … format/coefficient-basis/naming/convention +
  coverage gaps (DEGENERATE_LOCI_SOLUTION, ONSAGER_DETERMINABLE)." (change 5)
- Background-flow: emitted relative term `−2 q v₀/ω`; failure `|q v₀/ω| ≳ 1` (kc_s0≫ω necessary not
  sufficient). (change 3)
- Energy basis: 10 under §5's total-divergence quotient; representative non-unique (PY `st²` / WL `(∇·u)²`);
  SymPy over-counted to 11 (X-1), corrected. (change 6)
```
$ grep -c 'only where the transverse\|no reciprocal traction\|Coverage gaps\|B_div + 2μ_S/3\|q v₀ / ω' steps/S11b_interface_coupling_law.md
```
(the corrected phrases are present in the committed step record — see commit `8ddccb74`.)
