from estimator_with_hyb_dual import *
from sage.functions.error import Function_erf
import logging
import sage.crypto.lwe

oo = PlusInfinity()

def my_binary_search(f, start, stop, param, predicate=lambda x, best: x<=best, *arg, **kwds):
    kwds[param] = stop
    D = {}
    D[stop] = f(*arg, **kwds)
    best = D[stop]
    b = ceil((start+stop)/2)
    direction = 0
    while True:
        if b not in D:
            kwds[param] = b
            D[b] = f(*arg, **kwds)
        if b == start:
            best = D[start]
            break
        if not predicate(D[b], best):
            if direction == 0:
                start = b
                b = ceil((stop+b)/2)
            else:
                stop = b
                b = floor((start+b)/2)
        else:
            best = D[b]
            if b-1 not in D:
                kwds[param] = b-1
                D[b-1] = f(*arg, **kwds)
            if predicate(D[b-1], best):
                stop = b
                b = floor((b+start)/2)
                direction = 0
            else:
                if b == stop:
                    break
                elif b+1 not in D:
                    kwds[param] = b+1
                    D[b+1] = f(*arg, **kwds)
                if not predicate(D[b+1], best):
                    break
                else:
                    start = b
                    b = ceil((stop+b)/2)
                    direction = 1
    return best

def prob_M(h_M, h, n, r):
    result = 0
    for i in range(h_M+1):
        result += binomial(h, 2*i) * binomial(n - h, r - 2*i)
    return result / binomial(n,r)
    
def prob_c(h_M, h, n, r):
    result = 0
    denom = 0    
    for i in range(h_M+1):
        tmp = binomial(h, 2*i) * binomial(n-h, r-2*i)
        denom += tmp
        tmp *= binomial(2*i, i) * binomial(r - 2*i, h_M - i) / 2**i
        result += tmp
    return result / (denom * binomial(r, h_M))

def prob_s(R, alpha, q, m):    
    scaling_factor = RR(sqrt(pi) / (2*alpha*q))    
    probs = [RR(2 * R[i] * scaling_factor).erf() for i in range(len(R))]
    results = [probs[i] + (exp(-(2 * R[i] *scaling_factor)**2) - 1) / (2 * sqrt(pi) * R[i] * scaling_factor) for i in range(len(R))]
    return prod(results)

def sizeW(h_M, r):
    return 2**h_M * binomial(r, h_M)

def prob_NP_and_GSA(m, n, alpha, q, r, delta_0, opt_GSA, h, h_M = False):
    dim = m + n-r + 1    
    scaling_factor = RR(sqrt(pi) / (2*alpha*q)) 
    if h_M is False:
        h_ = h * (n-r) / n + 1
        scale = RR(sqrt((n+1-r)/h_)*alpha*q/sqrt(2*pi))
    else:
        scale = RR(sqrt((n+1-r)/(h-2*h_M+1))*alpha*q/sqrt(2*pi))

    if opt_GSA is True:
        R = Mod_GSA(m, q, dim, delta_0, scale)
    else:
        R = GSA(m, q, dim, delta_0, scale)        
     
    probs = [RR(R[i] * scaling_factor).erf() for i in range(len(R))]
    return prod(probs), m, dim, R

def time_M(dim, L, flag = False):
    return dim**2 * L / 2**(1.06)

def GSA(m, q, dim, delta_0, scale):
    det = RR((q**m * scale**(dim-m))**(1/dim))
    b = [delta_0**(dim-2*i) * det for i in range(dim)]
    return b

def Mod_GSA(m, q, dim, delta_0, scale):
    
    # Modified GSA for q-ary lattices 
    #       proposed by [Wun16]

    k = min(int(sqrt((dim - m)/log(delta_0, q/scale)).n()), dim)
    b_1 = [q for i in range(dim-k)]
    b_2 = [RR(delta_0**(k - 2*(i-(dim-k)-1)) * (q/scale)**((k-dim+m)/k) * scale) for i in range(dim - k, dim)]
    return b_1 + b_2

def cost_function_delta(n, q, alpha, h, delta_P, m_max,
                        reduction_cost_model, opt_GSA, verbose, file):
    
    line = 'Current delta_P = %f \n' % delta_P
    if verbose is True:
        print(line)
    # if file is not None:
    #     file = open(str(n)+ "_" + str(int(log(q,2))) + "_" + str(h) + ".txt", "a")
    #     file.write(line)
    #     file.close()
    
    kwds = {"n":n, "q":q, "alpha":alpha, "h":h, "delta_P":delta_P, "m_max":m_max,
            "reduction_cost_model": reduction_cost_model,
            "opt_GSA":opt_GSA, "verbose":verbose}

    best_delta = my_binary_search(cost_function_r, start = 0, stop = n, param="r", 
                                predicate=lambda x, best: RR(x["rop"])<=RR(best["rop"]), **kwds)

             
    line = '    * So far Best rop : %f with (delta_0, r, mitm, post, m) = (%f, %d, %d, %d, %d)\n' % (log(best_delta["rop"], 2), best_delta["delta_0"], best_delta["r"], best_delta["mitm"], best_delta["post"], best_delta["m"]) 
    if verbose is True:   
        print(line)
    # if file is not None:
    #     file = open(str(n)+ "_" + str(int(log(q,2))) + "_" + str(h) + ".txt", "a")
    #     file.write(line)
    #     file.close()
    line = '                            log(p_M, p_NP, p_s, p_c) = (%f, %f, %f, %f)\n' % (RR(log(best_delta["p_M"],2)), RR(log(best_delta["p_NP"],2)), RR(log(best_delta["p_s"],2)), RR(log(best_delta["p_c"],2)))
    if verbose is True:   
        print(line)  
    # if file is not None:
    #     file = open(str(n)+ "_" + str(int(log(q,2))) + "_" + str(h) + ".txt", "a")
    #     file.write(line)
    #     file.close()       

    return best_delta

def cost_function_r(n, q, alpha, h, r, delta_P, m_max,
                 reduction_cost_model, opt_GSA, verbose):
    

    kwds = {"n":n, "alpha":alpha, "q":q, "r":r, "delta_0":delta_P, "opt_GSA":opt_GSA, "h":h}

    p_NP, m, dim, R = my_binary_search(prob_NP_and_GSA, start = 50, stop = 2*n, param="m", 
                                predicate=lambda x, best: RR(x[0])>=RR(best[0]), **kwds)

    current = lattice_reduction_cost(reduction_cost_model, delta_P, dim, B=log(q, 2))
    
    p_s = 1.0
    p_c = 1.0
    p_M = 1.0
    h_M = max(0, int((h-n+r+1)/2))
    mitm_flag = True
    
    if RR(p_NP).is_NaN() or n == r:
        current["rop"] = oo

    else:
        p_s = prob_s(R, alpha, q, m)
        log_p_s = RR(log(p_s, 2))
        time_lat = current["rop"]
        time_post = 0

        thres = RR(log(time_lat,2))
        Wsize = sizeW(h_M, r)
        # Approximate computation for p_c, which is always smaller than the actual comptuation
        while h_M < min(int(r/2), int(h/2)):
            log_p_c = h_M * log(RR(h_M / r), 2)
            if RR(log(Wsize, 2)) < RR((log(Wsize, 2)) - log_p_c - log_p_s) /2:
                log_time_post = RR(log(Wsize,2)) + 15.1
            else:
                log_time_post = RR((log(Wsize, 2) - log_p_c - log_p_s) / 2) + 15.1
            if log_time_post > thres:
                break
            h_M = min(h_M+1, min(int(r/2),int(h/2)))
            Wsize = sizeW(h_M, r)
    
        p_c = prob_c(h_M, h, n, r)
        scale = RR(sqrt((n+1-r)/(h-2*h_M+1))*alpha*q/sqrt(2*pi))
        if opt_GSA is True:
            R = Mod_GSA(m, q, dim, delta_P, scale)
        else:
            R = GSA(m, q, dim, delta_P, scale)   
        p_s = prob_s(R, alpha, q, m)
        if Wsize * p_c * p_s < 1:
            L = Wsize
            mitm_flag = False
        else:
            L = sqrt(Wsize / (p_c * p_s))
            mitm_flag = True
        time_post = time_M(dim, L)

        # Compensate with Actual value
        while h_M > max(0, int((h-n+r+1)/2)): 
            p_c = prob_c(h_M-1, h, n, r)
            scale = RR(sqrt((n+1-r)/(h-2*(h_M-1)+1))*alpha*q/sqrt(2*pi))
            if opt_GSA is True:
                R = Mod_GSA(m, q, dim, delta_P, scale)
            else:
                R = GSA(m, q, dim, delta_P, scale)   
            p_s = prob_s(R, alpha, q, m)
            Wsize = sizeW(h_M-1, r)
            if Wsize * p_c * p_s < 1:
                L = Wsize # Wsize is too small -> Exhaustive search
                mitm_flag = False
            else:
                L = sqrt(Wsize / (p_c * p_s))
                mitm_flag = True
            time_post = time_M(dim, L)
            if RR(log(time_post,2)) < thres:
                break
            h_M -= 1

        p_M = prob_M(h_M, h, n, r)

        #if verbose is True:
        #    print '     * Current r = %d, h_M = %d, log(|W|) = %f, m = %d' % (r, h_M, RR(log(Wsize,2)), m)
        #    print '         log(p_M, p_NP, p_s, p_c) = (%f, %f, %f, %f)' % (RR(log(p_M,2)), RR(log(p_NP, 2)), RR(log(p_s, 2)), RR(log(p_c, 2)))
        #    print '         time_lat = %f, time_post = %f' % (thres, RR(log(time_post,2)))

        probability = p_M * p_NP
        current["rop"] = time_lat + time_post
        current = current.repeat(1/probability, select={"m": False, "red": False})  

    if mitm_flag == True:
        current["post"] = 2*h_M
    else:
        current["post"] = h_M
    current["mitm"] = mitm_flag
    current["p_M"] = p_M
    current["p_NP"] = p_NP
    current["p_s"] = p_s
    current["p_c"] = p_c
    current["r"] = r
    current["m"] = m
    current = current.reorder(["rop"])
    return current

def primal_hyb(n, alpha, q, secret_distribution, 
               m = oo, h = None, success_probability=0.99, 
               reduction_cost_model=reduction_default_cost, 
               init_delta_P = None, opt_GSA = False, verbose=True):

    m_max = m
    delta_step = min(0.005, (init_delta_P - 1)/2)
    f = None

    best = cost_function_delta(n, q, alpha, h, init_delta_P, m_max,
                               reduction_cost_model=reduction_cost_model,
                               opt_GSA=opt_GSA, verbose=verbose, file = f)
        
    while delta_step > 0.00005:
        
        current_pre = None
        if best["delta_0"] != init_delta_P:
            current_pre = cost_function_delta(n, q, alpha, h, best["delta_0"] - delta_step, m_max,
                                              reduction_cost_model=reduction_cost_model, 
                                              opt_GSA=opt_GSA, verbose=verbose, file = f)
            
        current_post = cost_function_delta(n, q, alpha, h, best["delta_0"] + delta_step, m_max,
                                          reduction_cost_model=reduction_cost_model, 
                                          opt_GSA=opt_GSA, verbose=verbose, file = f)

        delta_step /= 2
        min_cost = 0
        
        cost_pre = oo
        cost_old = RR(log(best["rop"],2))
        cost_post = RR(log(current_post["rop"],2))
        
        min_cost = min(cost_old, cost_post)
        if current_pre is not None:
            cost_pre = RR(log(current_pre["rop"],2))
            min_cost = min(min_cost, cost_pre)

        if min_cost == cost_post:
            best = current_post
        elif min_cost == cost_pre:
            best = current_pre
        
        if 0 < cost_old - min_cost < 0.5:
            break

        # f = open(str(n)+ "_" + str(int(log(q,2))) + "_" + str(h) + ".txt", "a")
        line = '** So far Best rop : %f with (delta_0, r, mitm, post, m) = (%f, %d, %d, %d, %d)\n' % (float(log(best["rop"], 2)), best["delta_0"], best["r"], best["mitm"], best["post"], best["m"]) 
        if verbose is True:
            print(line)
        # f.write(line)
        line = '                        log(p_M, p_NP, p_s, p_c) = (%f, %f, %f, %f)\n' % (RR(log(best["p_M"],2)), RR(log(best["p_NP"],2)), RR(log(best["p_s"],2)), RR(log(best["p_c"],2))) 
        if verbose is True:
            print(line)
        # f.write(line)
        # f.close()

    # f = open(str(n)+ "_" + str(int(log(q,2))) + "_" + str(h) + ".txt", "a")
    line = '***** Final Best rop : %f with (delta_0, r, mitm, post, m) = (%f, %d, %d, %d, %d)\n' % (float(log(best["rop"], 2)), best["delta_0"], best["r"], best["mitm"], best["post"], best["m"]) 
    if verbose is True:
        print(line)
    # f.write(line)
    line = '                        log(p_M, p_NP, p_s, p_c) = (%f, %f, %f, %f)\n' % (RR(log(best["p_M"],2)), RR(log(best["p_NP"],2)), RR(log(best["p_s"],2)), RR(log(best["p_c"],2))) 
    if verbose is True:
        print(line)
    # f.write(line)
    # f.close()

    return best