# Numerical interior coefficient matrices

All four independent grades 00,01,10,11 now have complete 645-by-645 interior
operators, with local and nonlocal parts retained separately. The constructor
reuses all 80 accepted momentum rows and contracts all 160 native cell terms;
no numerical integration was repeated. It completed in 18.68 seconds with peak
RSS 322628 KiB, exit zero and empty stderr.

Recombination at the approved point agrees with the accepted total/local
operator to 4.79e-16/3.58e-16 in the scaled reference-frame comparison. Original
unsplit binding at zero and a second formal point agrees to 0/3.90e-16.
Omitting the computed mixed term gives a nonzero 9.17e-7 diagnostic. These
arithmetic controls are not new physical input instances. Full residual arrays,
12 component matrices and source bindings remain in durable packets.

All 858 tags, 427 write keys and 4179 metadata paths replay. Source, original
quadrature and pre/post-emission packet hashes pass. The 399729-byte transcript
is published at 9471a10f; its actual MD5E backend key, symlink, MD5 and independent
SHA256 are verified in S11c_d_continuum_matrix_checkpoint.json. The matrix packet
SHA256 is 1db95fb176c970298b33005a94e25f8ada3ba1423fad605231be66d8a711954e.
Floating weighted projections supplement full-array hashes; they are not exact
prime-field proofs.

These are finite interior matrices before modal boundary replacement, at
interval/source64, momentum4, profile14 and regulator0.1. Boundary/channel/current
coefficients and the recursive inverse response remain next. No continuum
scattering expansion, flux coefficient or pole follows from matrix assembly
alone. The approximate modal-boundary and exploratory acceptance scope remain.
