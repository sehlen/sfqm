from find_simple_c import *

def find_simple_anisotropic(num, s, weights=range(4, 27), dynamic=True, lower=1, only_2_n=False, reduction=False):
    r"""
      Test for anisotropic $k$-simple fqm's for $k$ in weights.

      INPUT::
        - num: upper bound for the order
        - s: number of iterations
        - lower: lower bound for the level

      TODO: don't mix level and order!
    """
    return _find_simple_anisotropic_parallel(num, s, weights, dynamic, lower, only_2_n, reduction)
    
