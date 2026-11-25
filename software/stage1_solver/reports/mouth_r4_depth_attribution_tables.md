# Step-108 log-depth diagnostic · Cell tables

Cells are ranked by the absolute wide-probe force discrepancy at ε = 10⁻⁸. Indices are radial, depth and polar, starting at zero. The coordinates are physical cell centers. The signed force error is the negative of the energy-derivative error. All ten cells lie on the collar side, clear of both the core–collar weld and the outer boundaries.

| Rank | Cell | Radius | Depth | θ (rad) | Vertex χ range | Wide force error, ε = 10⁻⁸ | Wide force error, ε = 10⁻¹⁰ |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | (23, 2, 5) | 3.9317491 | 0.06468283 | 0.2754652 | 0.1190941–0.2076513 | +2.087570055e-04 | +2.087570050e-04 |
| 2 | (23, 2, 1) | 3.9317909 | 0.06468283 | 0.0488450 | 0.1195383–0.2080838 | +7.039279345e-05 | +7.039279327e-05 |
| 3 | (22, 2, 1) | 3.8965477 | 0.06468283 | 0.0488450 | 0.1191378–0.2075689 | +5.690346798e-06 | +5.690346784e-06 |
| 4 | (24, 2, 5) | 3.9552439 | 0.06468283 | 0.2754652 | 0.1194098–0.2079958 | +5.293403280e-06 | +5.293403267e-06 |
| 5 | (23, 2, 9) | 3.9316192 | 0.06468283 | 0.5753610 | 0.1179182–0.2062266 | -7.431023346e-10 | -7.431023264e-10 |
| 6 | (23, 2, 13) | 3.9314198 | 0.06468283 | 0.9243165 | 0.1159892–0.2035716 | -2.656939518e-10 | -2.656939456e-10 |
| 7 | (23, 2, 8) | 3.9316603 | 0.06468283 | 0.4952220 | 0.1182956–0.2067097 | -1.954271378e-10 | -1.954271367e-10 |
| 8 | (23, 2, 14) | 3.9313689 | 0.06468283 | 1.0178838 | 0.1155173–0.2027815 | -4.291274338e-11 | -4.291274578e-11 |
| 9 | (23, 2, 4) | 3.9317669 | 0.06468283 | 0.2102967 | 0.1192546–0.2078408 | -2.668848461e-11 | -2.668848354e-11 |
| 10 | (23, 8, 15) | 3.8434154 | 0.21992161 | 1.1137403 | 0.0000000–0.0028706 | +2.441910633e-12 | +2.441910261e-12 |

## Orientation and conditioning

Entries are the minima over the native physical quadrature points evaluated on each map, including the azimuthal direction in the singular-value calculation. The two minima in a column can occur at different points. The signed endpoints use h = 2⁻¹⁸. Only ranks 1–4 change reference triangle, all on the negative arm.

| Rank | Baseline min Jinv | Minus min Jinv | Plus min Jinv | Baseline min σ | Minus min σ | Plus min σ |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.897497793 | 0.897501535 | 0.897494050 | 0.733310377 | 0.726474346 | 0.733275768 |
| 2 | 0.898405892 | 0.898409624 | 0.898402160 | 0.761316726 | 0.761318797 | 0.761314655 |
| 3 | 0.962889891 | 0.962893668 | 0.962886115 | 0.803255423 | 0.803257003 | 0.803253842 |
| 4 | 0.731856839 | 0.731859906 | 0.731853772 | 0.683289072 | 0.683290464 | 0.683287679 |
| 5 | 0.903556907 | 0.903560633 | 0.903553181 | 0.726645781 | 0.726703066 | 0.726588498 |
| 6 | 0.903161356 | 0.903165076 | 0.903157636 | 0.745688251 | 0.745689869 | 0.745686633 |
| 7 | 0.903265188 | 0.903268913 | 0.903261462 | 0.759394575 | 0.759396480 | 0.759392671 |
| 8 | 0.904823083 | 0.904826816 | 0.904819350 | 0.736674103 | 0.736675652 | 0.736672553 |
| 9 | 0.907349001 | 0.907352740 | 0.907345261 | 0.731013419 | 0.731023667 | 0.731003171 |
| 10 | 1.341502467 | 1.341509062 | 1.341495872 | 0.803444778 | 0.803448542 | 0.803441014 |

## Material-face trace

The four leading cells form two pairs sharing radial endpoints. Each crossing adds a subinterval in the neighboring triangle of the collar reference mesh. The reference depth is a material label. The physical depth of these points is 0.06468283.

| Rank | Depth/angle Gauss indices | Signed log-depth change at crossing | Added fraction of radial line at wide negative probe | Order-weighted energy density before → after, at face |
|---:|---|---:|---:|---:|
| 1 | (1, 1) | -3.5998495649e-06 | 1.509442954e-04 | 0.00688738316 → 0.0264637018 |
| 2 | (1, 2) | -2.1169319724e-06 | 1.127474077e-03 | 0.000859831593 → 0.0103099643 |
| 3 | (1, 2) | -2.1169319724e-06 | 5.530937224e-05 | 0.000908212831 → 0.0120379808 |
| 4 | (1, 1) | -3.5998495649e-06 | 8.784169056e-06 | 0.0117241709 → 0.0259400336 |

At the dominant face, the two axial components of the reference-gradient jump are −0.2475117536 and −0.02506941688. The radial cut builder inserts the crossed depth plane, then the constitutive evaluation selects the adjacent stored reference triangle. [Radial cuts](../tools/mouth_r4_eulerian_solid.cpp) · [Piecewise reference embedding](../tools/mouth_free_shape_solid.cpp) · [Per-point witnesses and exact sided values](../../../research/pde_ledger_v3/_measurements/mouth_r4_depth_cell_geometry.json)

## Signed crystal-energy derivatives

The plus column uses [E(+h) − E(0)]/h and the minus column uses [E(0) − E(−h)]/h. Each subtracts the analytic derivative on the baseline branch. The ratio divides the central error by that at half the width. A ratio whose denominator is near the numerical precision of the slope has limited diagnostic value. E″ is a second derivative, measured from the difference of the two one-sided slopes divided by h.

| Rank | h | Central error | Plus error | Minus error | Halving ratio | E″ |
|---|---:|---:|---:|---:|---:|---:|
| All cells | 2⁻¹⁸ | -2.90132315e-04 | +8.28040618e-05 | -6.63068691e-04 | -908147 | 195.526067 |
| All cells | 2⁻¹⁹ | +3.19477309e-10 | +4.14017147e-05 | -4.14010758e-05 | 4.11434 | 43.4125094 |
| All cells | 2⁻²⁰ | +7.76496763e-11 | +2.07007947e-05 | -2.07006394e-05 | — | 43.4125501 |
| 1 | 2⁻¹⁸ | -2.08757005e-04 | -2.85232503e-06 | -4.14661686e-04 | -1.50867e+07 | 107.953353 |
| 1 | 2⁻¹⁹ | +1.38371638e-11 | -1.42617630e-06 | +1.42620397e-06 | 4.0972 | -1.49546875 |
| 1 | 2⁻²⁰ | +3.37722058e-12 | -7.13091569e-07 | +7.13098323e-07 | — | -1.49546849 |
| 2 | 2⁻¹⁸ | -7.03927935e-05 | +2.51747771e-10 | -1.40785839e-04 | -3.87027e+13 | 36.9062289 |
| 2 | 2⁻¹⁹ | +1.81881052e-18 | +1.25873886e-10 | -1.25873882e-10 | -0.00506857 | 0.000131988334 |
| 2 | 2⁻²⁰ | -3.58841042e-16 | +6.29363294e-11 | -6.29370471e-11 | — | 0.000131987802 |
| 3 | 2⁻¹⁸ | -5.69034680e-06 | +3.94579818e-10 | -1.13810882e-05 | 4.57183e+10 | 2.98358742 |
| 3 | 2⁻¹⁹ | -1.24465330e-16 | +1.97290030e-10 | -1.97290279e-10 | 0.247652 | 0.000206873721 |
| 3 | 2⁻²⁰ | -5.02580838e-16 | +9.86441302e-11 | -9.86451354e-11 | — | 0.000206872789 |
| 4 | 2⁻¹⁸ | -5.29340328e-06 | +1.30408422e-09 | -1.05881106e-05 | -3.2279e+09 | 2.77595153 |
| 4 | 2⁻¹⁹ | +1.63989126e-15 | +6.52042038e-10 | -6.52038758e-10 | -0.244815 | 0.000683713912 |
| 4 | 2⁻²⁰ | -6.69849081e-15 | +3.26012288e-10 | -3.26025685e-10 | — | 0.000683711371 |
| 5 | 2⁻¹⁸ | +7.43102335e-10 | -1.89095763e-05 | +1.89110625e-05 | 4.00077 | -9.91445355 |
| 5 | 2⁻¹⁹ | +1.85739657e-10 | -9.45497352e-06 | +9.45534500e-06 | 4.04738 | -9.91445308 |
| 5 | 2⁻²⁰ | +4.58913135e-11 | -4.72753289e-06 | +4.72762468e-06 | — | -9.91445131 |
| 6 | 2⁻¹⁸ | +2.65693952e-10 | +2.38044230e-05 | -2.38038916e-05 | 4.00297 | 12.480234 |
| 6 | 2⁻¹⁹ | +6.63742210e-11 | +1.19021456e-05 | -1.19020129e-05 | 4.20028 | 12.4802347 |
| 6 | 2⁻²⁰ | +1.58023444e-11 | +5.95105662e-06 | -5.95102502e-06 | — | 12.4802372 |
| 7 | 2⁻¹⁸ | +1.95427138e-10 | +1.57704400e-05 | -1.57700491e-05 | 4.00192 | 8.26814997 |
| 7 | 2⁻¹⁹ | +4.88333083e-11 | +7.88517137e-06 | -7.88507371e-06 | 4.12119 | 8.26815025 |
| 7 | 2⁻²⁰ | +1.18493203e-11 | +3.94257367e-06 | -3.94254997e-06 | — | 8.2681514 |
| 8 | 2⁻¹⁸ | +4.29127434e-11 | +5.26221701e-06 | -5.26213119e-06 | 4.00461 | 2.75889473 |
| 8 | 2⁻¹⁹ | +1.07158489e-11 | +2.63109794e-06 | -2.63107651e-06 | 4.4371 | 2.75889492 |
| 8 | 2⁻²⁰ | +2.41505484e-12 | +1.31554639e-06 | -1.31554156e-06 | — | 2.75889568 |
| 9 | 2⁻¹⁸ | +2.66884846e-11 | -2.06365266e-06 | +2.06370603e-06 | 4.00421 | -1.08196232 |
| 9 | 2⁻¹⁹ | +6.66509975e-12 | -1.03183292e-06 | +1.03184625e-06 | 4.28761 | -1.08196223 |
| 9 | 2⁻²⁰ | +1.55450205e-12 | -5.15918069e-07 | +5.15921178e-07 | — | -1.08196187 |
| 10 | 2⁻¹⁸ | -2.44191063e-12 | +7.09908986e-09 | -7.10397368e-09 | 3.99517 | 0.00372324789 |
| 10 | 2⁻¹⁹ | -6.11215532e-13 | +3.55015748e-09 | -3.55137991e-09 | 4.04758 | 0.00372325084 |
| 10 | 2⁻²⁰ | -1.51007489e-13 | +1.77522977e-09 | -1.77553178e-09 | — | 0.00372324334 |

The machine-readable banks include unrounded values for both gauges, individual one-sided slopes, the baseline crystal energies, and each cell’s analytic work from geometry, label transport and movement of quadrature cuts. [Both gauge ladders and closure](../../../research/pde_ledger_v3/_measurements/mouth_r4_depth_cell_work.json)

## Certificate witness in cell (25, 11, 10)

The physical kernel reports the lower corner of the last uncertified subbox. The exact quadratic edge minimum is nearby. Chart singular values refer to normalized label coordinates. Physical singular values use the orthonormal metrics on both sides of the map and include the fourth, hoop direction.

| Measurement | Reported subbox corner | Minimum of chart determinant along edge |
|---|---:|---:|
| Physical (r, w, θ) | 3.85106331, 0.284604431, 0.649924286 | 3.85106107, 0.284604431, 0.649990891 |
| Chart determinant | 1.5686515812e-07 | 1.1973704821e-07 |
| Chart singular values | 1.19701586e-07, 0.922110231, 1.42116252 | 9.14448012e-08, 0.92203418, 1.4201117 |
| Physical Jinv | 2.2683711757e-07 | 1.7314834492e-07 |
| Physical singular values of F⁻¹ | 1.79551186e-07, 0.759999516, 1.03707756, 1.60288103 | 1.3711008e-07, 0.760152926, 1.03708043, 1.60190055 |
| Smallest eigenvalue of C = FᵀF | 0.389222038 | 0.389698646 |
| det(C) | 1.94344279e+13 | 3.33551817e+13 |
| Smallest squared inverse stretch | 3.22386284e-14 | 1.8799174e-14 |

At depth 8 the native test first checks all 27 point values, then the Bernstein coefficients. If the lower coefficient is nonpositive and the subdivision cap is reached, it returns code 3. A sampled nonpositive determinant would instead return code 2. The alternative chart certificate allows depth 12 and subtracts a roundoff margin of 9.2818 × 10⁻¹⁴ on this cell. At depth 9 the margin-adjusted bound agrees to within 5.3 × 10⁻¹⁹. [Native orientation predicate](../tools/mouth_r4_eulerian_solid.cpp) · [Whole-cell chart check](../tools/mouth_r4_eulerian_relax_labels.py) · [Exact-rational bounds and physical spectra](../../../research/pde_ledger_v3/_measurements/mouth_r4_depth_cert.json)
