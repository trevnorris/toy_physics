# Principal-power reader repair

The original helper SHA c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec
is immutable. The failed guard was the extra positive-real requirement in its
new numerical reader. The accepted numerical engine, SHA
5bf77729682960bb59f10c672a5fe3f67a98b0a9533da7c6189a37183596050f,
already evaluates noninteger powers on the principal complex branch and retains
any explicit Piecewise predicates. This continuation removes only that extra
restriction. Unsupported nodes still fail closed; no predicate is discarded.

The incomplete full coefficient call has no returned value or completion receipt.
All completed source coefficients/actions/comparison are restored from exact bytes.
No statement of the remaining native-integrand and quadrature body is altered.
The recovery compiles and reverses its exact AST changes before calling the suffix.
