# SparseLWE-estimator

A Sage module provides several functions that estimates bit-security of Learning With Errors (LWE) problem with "small and sparse" secrets.

## Requirement

This script is checked with SageMath version 9.0 (using Python 3.7.3).

Lower version SageMath would be problematic because they are based on Python 2.x.

## How to Run Example

For Hybrid-Dual Attack,

    > sage: from estimator_sparseLWE import *
    > sage: n = 8192; q = next_prime(2^125); alpha = 8/q; h = 64;
    > sage: hybrid_dual(n, alpha, q, secret_distribution=((-1,1),h), reduction_cost_model = BKZ.sieve)
   
For Hybrid-Primal Attack,

    > sage: from estimator_sparseLWE import *
    > sage: n = 8192; q = next_prime(2^125); alpha = 8/q; h = 64;
    > sage: hybrid_primal(n, alpha, q, secret_distribution=((-1,1),h), reduction_cost_model = BKZ.sieve)
   
## Reference

This module is based on "LWE-estimator" : Albrecht et al, https://bitbucket.org/malb/lwe-estimator/src/master/),
which provides estimation of "normal" LWE problem (non-sparse).

Hybrid-Dual Attack: https://eprint.iacr.org/2019/1114

Hybrid-Primal Attack: https://eprint.iacr.org/2019/1019
