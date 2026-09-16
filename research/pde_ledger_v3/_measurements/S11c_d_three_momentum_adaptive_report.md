# S11c-d independent adaptive outer quadrature

Production passed in 2h12m with exit zero and empty stderr in all four workers.
Adaptive GK21 outer quadrature agrees with accepted Gauss-144 raw integrals
within 1.74e-16 (field 0) and 1.33e-16 (field 1). Complete actions round equal
at working precision; every raw integral and native term difference is retained.
Summed half-interval error estimates are 7.29e-12 and 4.33e-12. Each half targets
5e-11 absolute in the declared numerical unit frame, with zero relative tolerance.

All ten native three-momentum rows and both approved fields are covered. Inner
orders remain 24/24, source/profile orders 256/256, momentum/source/profile
bounds 2/32/10 and regulator 0.2. Single/pair contributions remain unchanged.
This is independent outer quadrature on held inner/source/profile operands;
it establishes no uniform or independent-grade convergence or physical limit.

All 420 conditional points, 6,764 partials, four worker results and four records
pass source/field/profile/limit/provenance and conditional-unit joins. Finite
conditional masses agree within 5.69e-14 and every actual measure mutation
responds. All 67 current/frozen sources, exact read-only caches and pre/post
packet hashes pass. Replay covers 3,906 tags, 1,951 fresh keys and 20,678 metadata
paths. Peak worker/coordinator RSS is 134.8/194.1 MiB; workspace budgets are
separate. The native engine and accepted constructor/emitter remain unchanged.

The 1,184,850-byte transcript is prepared for DataLad/git-annex publication at
`scripts/out/S11c_d_three_momentum_adaptive.out`, SHA256
`7247f932018e700746c7ebda0dd09826c492c7255bb29735584fe7b7ed57b1dc`.
Preflight was accepted at 19cfa05c; its coarse smoke was an instrument check.
Next: source/profile domain and tail tests, followed by momentum-domain and
Abel checks as supported by the evidence. Scattering and poles remain work.
Preserve all saved operands, approved inputs and the retained solver contract.
