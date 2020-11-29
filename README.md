# SparseLWE-estimator

A Sage module provides several functions that estimates bit-security of Learning With Errors (LWE) problem with "small and sparse" secrets.

This module is based on "LWE-estimator" : Albrecht et al, https://bitbucket.org/malb/lwe-estimator/src/master/),
which provides estimation of "normal" LWE problem (non-sparse).

### How to Run

    > sage: from estimator_hybrid import *
    > sage: n = 8192; q = next_prime(2^125); alpha = 8/q; h = 64;
    > sage: MITM_estimator(n, alpha, q, h, reduction_cost_model=BKZ.sieve)
    
    Chosen Parameters :
             n =  8192, log(q) = 125.0, stddev =  3.19, HW(s) =   64
     
    Start Estimation . . .

    Optimizing with beta =  240 . . .
    Optimizing with beta =  247 . . .
    Optimizing with beta =  243 . . .
    Optimizing with beta =  244 . . .

    == Bit-security : 129.9 with optimal parameters
         k =  4861, h1 =  7, h2 =  9, beta =  240, mem =  80.8
                 (For simplicity, we set k1 = k2 = k/2)

    == Details
             rop:  2^129.9
               m:   2^12.8
             red:  2^102.3
         delta_0: 1.005603
            beta:      240
          repeat:  170.419
               d:   2^12.8
               c:   23.025
            post:  2^97
        prob_inv:   2^27.6
               k:   2^12.2
              h1:        7
              h2:       9
             mem:   2^80.8
   
