from simple_modules_graph import SimpleModulesGraph
from ..fqm.genus_symbol import anisotropic_symbols
from sage.all import Integer, is_squarefree, is_even, is_odd, squarefree_part, walltime, RR
from sage.parallel.decorate import *

#########################################
# Testing functions
#########################################

@parallel
def find_simple_anisotropic(lower, upper, max_order=None, sig=None, k=None, weights=range(4, 27), test_dim=True, dynamic=True, only_2_n=False, reduction=False):
    r"""
      ```find_simple_anisotropic``` computes a list of all anisotropic finite quadratic modules
      up to a given maximal order ```max_order``` that are $k$-simple for given weights.

      INPUTS:
      - lower: an integer, lower bound for the level
      - upper: an integer, upper bound for the level
      - max_order: an integer, upper bound for the order, if None: take upper
      - sig: an integer, the signature of the modules to be found (modulo 8)
      - k: a half-integer, only check for $k$-simple modules
      - weights: a list of half-integers, check for $k$-simple modules for $k$ in weights
      - test_dim: boolean,
      - dynamic: boolean,
      - only_2_n: boolean,
    """

    # type checking
    lower = Integer(lower)
    upper = Integer(upper)
    max_order = Integer(max_order) if max_order is not None else Integer(upper)
    k = Integer(2*k)/2 if k is not None else None
    if isinstance(weights, list):
        weights = map(lambda x: Integer(x), weights)
    elif k is None:
        raise ValueError("Either k or a list of half-integers for weights has to be provided as input.")
        
    t = walltime()
    if k != None:
        weights = [k]
    simple = dict()
    for kk in weights:
        simple[kk] = list()
        
    cdef double d = 0 # the last computed dimension
    numpos = 0 # the number of consecutive positive dimensions
    usedyn = True # dynamic checking
    # N is the level
    cdef long long N = 0
    for N in xrange(lower, upper + 1):
        # check that N has the correct type, otherwise continue
        if not is_squarefree(N):
            if is_even(N):
                if not Integer(N) / squarefree_part(N) in [4, 8]:
                    continue
        if d >= 0:
            # if the last computed dimension was non-negative,
            # we first ignore the quadratic form and
            # do the simplest possible check
            # that does only depend on the structure of the
            # finite abelian group
            syms = anisotropic_symbols(N, None, True)
            all_syms_computed = False
        else:
            if sig != None:
                syms = anisotropic_symbols(N, sig)
            else:
                syms = anisotropic_symbols(N)
            all_syms_computed = True
        check = False
        # first, we only check the dimension estimate
        # for anisotropic symbols
        # as soon as a non-positive estimate occurs,
        # we compute the dimensions exactly.
        # we need to do this in two steps because
        # the exact dimension depends on the quadratic form
        # and not only on the structure of the finite abelian group.
        for s in syms:
            # print numpos
            if s.order() > max_order:
                continue
            for kk in weights:
                # if we saw a few non-negative dimensions,
                # we try to apply the simpler check
                if numpos >= 10 * len(weights):
                    usedyn = False
                d = s.dimension_estimate_for_anisotropic(
                    Integer(kk), usedyn and dynamic)
                if d <= 0:
                    check = True # check exact dimension
                    numpos = 0
                    usedyn = True
                    break
                else:
                    numpos = numpos + 1
            if check:
                break
        if check and test_dim:
            if not all_syms_computed:
                # compute all the symbols if not done so before
                if not sig == None:
                    syms = anisotropic_symbols(N, sig)
                else:
                    syms = anisotropic_symbols(N)
            for ss in syms:
                for kk in weights:
                    # print ss.signature(), kk, ss
                    if only_2_n:
                        t1 = ss.signature() % 8 == (4 - 2*kk) % 8
                    else:
                        t1 = ss.signature() % 4 == (-2*kk) % 4
                    if t1 and s.dimension_estimate_for_anisotropic(kk, usedyn and dynamic) <= 0:
                        if ss.is_simple(kk, no_inv=True, aniso_formula=True, reduction=reduction):
                            simple[kk].append(ss)
        else:
            if check:
                print "check but not test_dim", N
    if test_dim:
        return simple, walltime(t)


cpdef test_simple_parallel(num, s, weights=range(4, 27), dynamic=True, lower=1, only_2_n=False, reduction=False):
    r"""
      Test for anisotripic $k$-simple fqm's for $2k$ in weights.

      INPUT::
        - num: upper bound
        - s: number of iterations
        - lower: lower bound
    """
    simple = dict()
    for kk in weights:
        simple[Integer(kk) / 2] = list()
    m = round(RR(num) / s)
    args = [(m * a + lower, min(m * (a + lower), lower + num - 1),
                          lower + num - 1, None, None, weights, True, dynamic,
                          only_2_n, reduction) for a in range(s)]
    tests = test_simple(args)
    done = 0

    #starttime = walltime()
    for test in tests:
        done += 1
        simple_part = test[1][0]
        if isinstance(simple_part, str):
            print "Result: ", test
            continue
        tt = test[1][1]
        for kk in simple_part.keys():
            print "Result from process: ", simple_part[kk]
            simple[kk] = uniq(simple[kk] + simple_part[kk])
        #tt = walltime(starttime)
        sl = [len(sp) for sp in simple.values()]
        num_simple = sum(sl)
        timeest = (s - done) * RR(tt) / RR(60)
        if timeest > 1:
            if timeest > 60:
                print ("%g%% done, ETA: %d hour" + ("s" if timeest / 60 >= 2 else "") + ", %d minutes, simple lattices: %d") % (
                    RR(done) / RR(s) * 100, timeest / 60, (timeest / 60).frac() * 60, num_simple)
            else:
                print ("%g%% done, ETA: %d minute" + ("s" if timeest >= 2 else "") +
                       ", simple lattices: %d") % (RR(done) / RR(s) * 100, timeest, num_simple)
        else:
            print "%g%% done, ETA: %d seconds, simple lattices: %d" % (RR(done) / RR(s) * 100, timeest * 60, num_simple)
    return simple
