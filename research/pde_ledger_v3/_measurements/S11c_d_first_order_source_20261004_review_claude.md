**CLEAR FOR THIS FIRST-ORDER TRANSVERSE SOURCE METHOD**

This clears source-instrument preparation only. It does not clear a worker, a receiving inverse, a physical plane-wave action, a regular grazing expansion, or any loss or leakage number. I found no incorrect hand identity and no missing essential term. I read the saved operands and did the algebra by hand; I ran no code.

**What I checked**

1. **Lift and chart.**
   - `incident/LEFT-invariant-lift.json` has uniform columns (−iN, 0, i/5) and (i/10, −i/5, 0). Both are exactly orthogonal to k = (1/5, 1/10, N) for a symbolic N. The θ and e_W rows are exactly zero.
   - `actual-source-chart.json` maps weak (x, e1, e2) to uniform axes (3, 1, 2).
   - Column 2 has U_x = 0 and equals i(0, h2, −h1) in the weak chart. That is the plan's test polarization, fixed by geometry.
   - Column 1 has U_x = i/5.
2. **S00 = S10 = 0.**
   - In `sources/{plus,minus}-source-jets-00.json`, the only u dependence is div u, and everything else is θ or e_W. The 10 records have the same structure, with u entering only as −w·div u.
   - With k·U = 0 and θ = e_W = 0 these vanish identically in x, not just on-shell.
   - Any linear R, even a distributional one, acting on an identically zero function gives zero. So the C10·R00·S00, C00·R10·S00 and C01·R00·S00 terms drop, and C00·R00·S01 is the whole pressure forcing.
   - The same holds for the chemical zero-grade part (div u only, `consumer/chemical-amplitude-grade-split.json`) and the velocity.
3. **S01 hand identity.**
   - I re-derived it term by term from the 01 source text, with profile jets in x only, u_j,k → i k_k U_j, and h1U2 + h2U3 = −pU_x.
   - The m_1 bracket is 400 + 2100p² + 600(h1²+h2²) − 1500p² = 400 + 600K². The w_1 bracket is 1900 + 6600p² + 2100(h1²+h2²) − 4500p² = 1900 + 2100K².
   - The m_2 coefficient is 210 − 80 = 130, and the w_2 coefficient is 660 − 230 = 430. The prefactor c matches the source prefactor.
   - Only the incident p enters; l does not.
   - The identity is correct as written, including K² = 5.95 + 0.05 = 6.
   - L powers are not an error. Profile jets have dimension 0 (`native-wave-profile-contract.json`) and the native rule multiplies them by `L_W**n` (`native/S11c_d_defect_source_composition.py:137`). So m_j and w_j are the ξ-derivatives. The instrument must not add another 1/L^j, but it must carry the factor L from dx = L dξ in any transform.
4. **Chemical source.** The 01 part of the chemical amplitude (`consumer/chemical-amplitude-grade-split.json`) is σ_W/10100 times the same U_x bracket. I checked every listed coefficient ratio. The U_x = 0 selection rule therefore applies to the chemical channel too. Keep the two channels separate, because their face maps differ.
5. **Velocity.** `consumer/plus-source-input.json` and the minus file give nativeVelocity = W0·e_W,t·ε/2 on both faces, which is zero for e_W = 0. There is no convected-sheet term. The normal jets are +iq (plus) and −iq (minus) in `consumer/reference-response-census.json`.
6. **R00.** The census flat factor is ω/(10(q + ω/(1000(−iω/1000 + 1/100)))). At ω = 3 this is 3/(10(q + β3)) with β3 = 3/(10−3i). That matches the plan and is not q/(q+β). Re β3 = 30/109 > 0, so R00 has no pole on either branch.
7. **Grazing.** q(p)² = 6 − 5.95 − 0.05 = 0, so the incident wave sits exactly at the exterior threshold. Using q(p) would be degenerate, as the plan says.
8. **Endpoints.** With w = (1+tanh)/2 and m = (1−tanh²)/3, all of w_1, w_2, m_1, m_2 decay exponentially at both ends. The S01 transform is therefore regular, and its Q = 0 moment is set only by Δw = 1 (the m_1, m_2 and w_2 integrals vanish). By hand that is c·U_x·i·14500·L, nonzero. The plan was right not to impose a zero there.
9. **k_out·U.** k_out·U = (l−p)U_x is exact. The diagnostic is fine as labeled. It vanishes for column 2 only for the part of f_U proportional to U, so it is not a statement about all of f_U.

**Non-blocking obligations for the instrument**
- Decide S01 against the full 01 operands for both faces. Keep ξ-scaled jets and track the L from dx.
- The first-order pressure and chemical drive is only part of the face data. Still required are the local-row x-derivative forces with their step terms, the matched-end work, the opposite-sign normal maps per face, and the θ/e_W rows.
- A localized S01 says nothing about the receiving side. R00 is finite, but q(l) has a square-root branch at |l| = p, and the power weight is unresolved.

**Coverage**
- Fully read: `plan.txt`, `guide.txt`, `review-prompt.md`, the lift and chart, the plus 00, 10 and 01 source texts, the chemical grade split, the velocity records, and the flat and normal-jet census.
- The minus 01 source was checked only by matching key segments (the full-line comparison failed) and the minus 00 prefix.
- Local rows: one return file and the 400-cell count only. I did not inspect the pressure splits, flat-proofs, route inspection, `governing/`, or the other `native/` files. I also did not check the RIGHT end record.