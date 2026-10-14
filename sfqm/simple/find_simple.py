from sfqm.simple.find_simple_c import _find_simple_anisotropic
from sage.all import Integer, is_squarefree, is_even, is_odd, squarefree_part, walltime, RR, uniq
from sage.parallel.decorate import parallel
from sfqm.fqm.genus_symbol import (GenusSymbol, gamma0_N_genus_symbol,
                                   gamma1_genus_symbol, t1_genus_symbol)

@parallel(ncpus=4)
def find_simple_anisotropic_wrapper(lower, upper, max_order=None, sig=None, k=None, weights=list(list(range(2, 14))), test_dim=True, dynamic=True, only_2_n=False, reduction=False, bound=0):
    return _find_simple_anisotropic(lower, upper, max_order, sig, k, weights, test_dim, dynamic, only_2_n, reduction, bound)

def find_simple_anisotropic(num, s, weights=list(list(range(2, 14))), dynamic=True, lower=1, only_2_n=False, reduction=False, bound=0):
    r"""
      Test for anisotropic $k$-simple fqm's for $k$ in weights.

      INPUT::
        - num: upper bound for the order
        - s: number of iterations
        - lower: lower bound for the level

      TODO: don't mix level and order!
    """
    return find_simple_anisotropic_parallel(num, s, weights, dynamic, lower, only_2_n, reduction, bound)

def find_simple_anisotropic_parallel(num, s, weights=list(list(range(2, 14))), dynamic=True, lower=1, only_2_n=False, reduction=False, bound=0):
    r"""
      Test for anisotropic $k$-simple fqm's for $k$ in weights.

      INPUT::
        - num: upper bound
        - s: number of iterations
        - lower: lower bound
    """
    simple = dict()
    for kk in weights:
        simple[kk] = list()
    from .search_ranges import inclusive_chunks
    args = [(lo, hi, lower + num - 1, None, None, weights, True, dynamic,
             only_2_n, reduction, bound)
            for lo, hi in inclusive_chunks(lower, num, s)]
    tests = find_simple_anisotropic_wrapper(args)
    done = 0

    #starttime = walltime()
    for test in tests:
        done += 1
        if not isinstance(test[1], tuple) or len(test[1]) != 2:
            raise RuntimeError("Anisotropic search worker failed: {0}".format(test))
        simple_part = test[1][0]
        if not isinstance(simple_part, dict):
            raise RuntimeError("Anisotropic search worker failed: {0}".format(test))
        tt = test[1][1]
        for kk in list(simple_part.keys()):
            print("Result from process: ", simple_part[kk])
            simple[kk] = list(uniq(simple[kk] + simple_part[kk]))
        #tt = walltime(starttime)
        sl = [len(sp) for sp in list(simple.values())]
        num_simple = sum(sl)
        timeest = (s - done) * RR(tt) / RR(60)
        if timeest > 1:
            if timeest > 60:
                print(("%g%% done, ETA: %d hour" + ("s" if timeest / 60 >= 2 else "") + ", %d minutes, simple lattices: %d") % (
                    RR(done) / RR(s) * 100, timeest / 60, (timeest / 60).frac() * 60, num_simple))
            else:
                print(("%g%% done, ETA: %d minute" + ("s" if timeest >= 2 else "") +
                       ", simple lattices: %d") % (RR(done) / RR(s) * 100, timeest, num_simple))
        else:
            print("%g%% done, ETA: %d seconds, simple lattices: %d" % (RR(done) / RR(s) * 100, timeest * 60, num_simple))
    return simple

def simple_gamma0_genus_symbols(r=list(range(1,500)), precomputed=None):
    if not isinstance(precomputed, list):
        raise ValueError("precomputed must be a list of genus symbols")
    l = precomputed
    comp = False
    res = []
    for s in l:
        if Integer(4).divides(s.level()):
            if s.defines_isomorphic_module(gamma0_N_genus_symbol(s.level()/Integer(4))):
                res.append(s)
                print(s)
    
    return res, [x for x in l if x not in res]
    

def simple_gamma1_genus_symbols(r=list(range(1,500)), precomputed=None):
    if not isinstance(precomputed, list):
        raise ValueError("precomputed must be a list of genus symbols")
    l = precomputed
    comp = False
    res = []
    
    if not comp:
        for s in l:
            if Integer(4).divides(s.level()):
                if s.defines_isomorphic_module(gamma1_genus_symbol(s.level()//4)):
                    res.append(s)
                    print(s)
                if (s.level()//2) % 2 == 0:
                    if s.defines_isomorphic_module(gamma1_genus_symbol(s.level()//2)):
                        res.append(s)
                        print(s)
                    if s.defines_isomorphic_module(gamma1_genus_symbol(s.level())):
                        res.append(s)
                        print(s)
    else:
        raise NotImplementedError
    
    return res, [x for x in l if x not in res]

def simple_t1_genus_symbols(Nr=list(range(1,100)), ar=list(range(1,100)), precomputed=None):
    if not isinstance(precomputed, list):
        raise ValueError("precomputed must be a list of genus symbols")
    l = precomputed
    comp = False
    res = []
    for s in l:
        added = False
        for a in s.level().divisors():
            for N in s.level().divisors():
                if s.defines_isomorphic_module(t1_genus_symbol(a,N)):
                    print(s)
                    added = True
                    res.append(((a,N), s))
                    break
            if added: break
        if not added:
            print("{0} not of type t1".format(s))
    return res, [x for x in l if x not in [symbol for _, symbol in res]]

def contains_isomorphic_module(l, s):
    for t in l:
        if s.defines_isomorphic_module(GenusSymbol(t._symbol_dict)):
            print(t)
            return t
    return None


def theoretical_upper_bound_anisotropic(k, lower_bound = 0):
    # find a theoretical upper bound B for the order of
    # an anisotropic(!) finite quadratic module A, such that the dimension
    # of the space of cusp forms of weight k for
    # the (dual) Weil representation attached to A
    # is <= lower_bound, i.e., if |A| > B,
    # then dim S(k,A) > lower_bound.
    # We use the estimate from Corollary 4.13 in [BEF].
    from sage.all import RR, log, ceil, find_root
    k = RR(k)
    if k <= 1 or lower_bound < 0:
        raise ValueError("require k > 1 and lower_bound >= 0")
    def bf(x):
        x = RR(x)
        return (x+1)*(k-1)/24 - 3 - RR('0.86')*x**(RR(5)/8)*log(2*x) - lower_bound
    # Seek the eventual crossing, also when lower_bound is zero (the default).
    left = RR(1)
    right = RR(2)
    while bf(right) <= 0:
        left, right = right, 2*right
    if bf(left) > 0:
        return Integer(ceil(right))
    result = Integer(ceil(find_root(bf, left, right))) + 1
    while bf(result) <= 0:
        result += 1
    return result
