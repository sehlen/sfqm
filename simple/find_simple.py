from find_simple_c import _find_simple_anisotropic
from sage.all import Integer, is_squarefree, is_even, is_odd, squarefree_part, walltime, RR, uniq
from sage.parallel.decorate import *

@parallel
def find_simple_anisotropic_wrapper(lower, upper, max_order=None, sig=None, k=None, weights=range(2, Integer(27)/2), test_dim=True, dynamic=True, only_2_n=False, reduction=False):
    return _find_simple_anisotropic(lower, upper, max_order, sig, k, weights, test_dim, dynamic, only_2_n, reduction)

def find_simple_anisotropic(num, s, weights=range(2, Integer(27)/2), dynamic=True, lower=1, only_2_n=False, reduction=False):
    r"""
      Test for anisotropic $k$-simple fqm's for $k$ in weights.

      INPUT::
        - num: upper bound for the order
        - s: number of iterations
        - lower: lower bound for the level

      TODO: don't mix level and order!
    """
    return find_simple_anisotropic_parallel(num, s, weights, dynamic, lower, only_2_n, reduction)

def find_simple_anisotropic_parallel(num, s, weights=range(2, Integer(27)/2), dynamic=True, lower=1, only_2_n=False, reduction=False):
    r"""
      Test for anisotropic $k$-simple fqm's for $k$ in weights.

      INPUT::
        - num: upper bound
        - s: number of iterations
        - lower: lower bound
    """
    simple = dict()
    for kk in weights:
        simple[Integer(kk)] = list()
    m = round(RR(num) / s)
    args = [(m * a + lower, min(m * (a + lower), lower + num - 1),
                          lower + num - 1, None, None, weights, True, dynamic,
                          only_2_n, reduction) for a in range(s)]
    tests = find_simple_anisotropic_wrapper(args)
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
    
