from simple_modules_graph import SimpleModulesGraph
from ..fqm.genus_symbol import anisotropic_symbols
from sage.all import Integer, is_squarefree, is_even, is_odd, squarefree_part, walltime, RR
from sage.parallel.decorate import *

#########################################
# Testing functions
#########################################

def _find_simple_anisotropic(lower, upper, max_order=None, sig=None, k=None, weights=range(2, Integer(27)/2), test_dim=True, dynamic=True, only_2_n=False, reduction=False, bound=0):
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
      - bound: consider the module to be simple if dimension is less or equal to this bound
    """

    # type checking
    lower = Integer(lower)
    upper = Integer(upper)
    max_order = Integer(max_order) if max_order is not None else Integer(upper)
    k = Integer(2*k)/2 if k is not None else None
    if isinstance(weights, list):
        weights = map(lambda x: Integer(2*x)/2, weights)
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
        if d >= bound:
            # if the last computed dimension was at least our bound,
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
                d = s.dimension_estimate_for_anisotropic(kk, usedyn and dynamic)
                if d <= bound:
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
                    if t1 and s.dimension_estimate_for_anisotropic(kk, usedyn and dynamic) <= bound:
                        if ss.is_simple(kk, no_inv=True, aniso_formula=True, reduction=reduction, bound=bound):
                            simple[kk].append(ss)
        else:
            if check:
                print "check but not test_dim", N
    if test_dim:
        return simple, walltime(t)
