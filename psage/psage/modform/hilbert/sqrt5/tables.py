import sqrt5_fast
import sqrt5
from sage.misc.all import cputime
from sage.rings.all import Integer, ZZ

F = sqrt5.F

def ideals_of_bounded_norm(B):
    return sum([v for n, v in F.ideals_of_bdd_norm(B).iteritems() if n != 1], [])
    
def ideals_of_norm(v):
    try:
        v = list(v)
    except TypeError:
        v = [Integer(v)]
    z = F.ideals_of_bdd_norm(max(v))
    return sum([z[n] for n in v if n>1],[])

def canonical_gen(I):
    """
    Return a canonical choice of generator for this ideal I (of a
    PID).

    The implementation computes the Hermite normal form (HNF) basis of
    the ideal, which is canonical, then finds a generator for the
    ideal it defines.

    EXAMPLES::

        sage: import psage.modform.hilbert.sqrt5.tables as tables
        sage: a = tables.F.gen()
        sage: z = a^30  * (-45*a+28); z
        -37284985*a - 23043388
        sage: tables.canonical_gen(tables.F.ideal(z))
        -45*a + 28
    """
    try:
        g =  I.ring().ideal(I.basis()).gens_reduced()
        assert len(g) == 1
        return g[0]
    except AttributeError:
        # sometimes I is just a usuaul integer
        return Integer(I)


def canonical_rep(z):
    """
    Return canonical generator for the ideal generated by z.  See the
    documentation for the canonical_gen function.
    """
    return canonical_gen(z.parent().ideal(z))
    
def test_canonical_gen(B=50):
    a = F.gen()
    z = -45*a + 28
    v = [canonical_rep(a**i * z)  for i in range(-B,B)]
    assert len(set(v)) == 1
        

def no_space(s):
    return str(s).replace(' ', '')

def dimensions(v, filename=None):
    """
    Compute dimensions of spaces of Hilbert modular forms for all the levels in v.
    The format is:

        Norm   dimension  generator  time
    """
    F = open(filename,'a') if filename else None
    for N in ideals_of_norm(v):
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        tm = cputime(t)
        s = '%s %s %s %s'%(N.norm(), H.cardinality(), no_space(canonical_gen(N)), tm)
        print s
        if F:
            F.write(s+'\n')
            F.flush()

def charpolys(v, B, filename=None):
    """
    Compute characteristic polynomials of T_P for primes P with norm <= B
    for spaces of Hilbert modular forms for all the levels in v.
    """
    F = open(filename,'a') if filename else None
    P = [p for p in ideals_of_bounded_norm(B) if p.is_prime()]
    for N in ideals_of_norm(v):
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        T = [(p.smallest_integer(),H.hecke_matrix(p).fcp()) for p in P if
             gcd(Integer(p.norm()), Integer(N.norm())) == 1]
        tm = cputime(t)
        s = '%s %s %s %s'%(N.norm(), no_space(canonical_gen(N)), tm, no_space(T))
        print s
        if F:
            F.write(s+'\n')
            F.flush()

def one_charpoly(v, filename=None):
    """
    Compute and factor one characteristic polynomials for all the
    levels in v.  Always compute the charpoly of T_P where P is the
    smallest prime not dividing the level.
    """
    F = open(filename,'a') if filename else None
    P = [p for p in ideals_of_bounded_norm(100) if p.is_prime()]
    for N in ideals_of_norm(v):
        NN = Integer(N.norm())
        t = cputime()
        H = sqrt5_fast.IcosiansModP1ModN(N)
        t0 = cputime(t)
        for p in P:
            if Integer(p.norm()).gcd(NN) == 1:
                break
        t = cputime()            
        T = H.hecke_matrix(p)
        t1 = cputime(t)
        t = cputime()
        f = T.fcp()
        t2 = cputime(t)
        s = '%s\t%s\t%s\t%s\t%s\t(%.1f,%.1f,%.1f)'%(N.norm(), no_space(canonical_gen(N)), 
                                                    p.smallest_integer(), no_space(canonical_gen(p)), no_space(f),
                                                    t0, t1, t2,)
        print s
        if F:
            F.write(s+'\n')
            F.flush()


def elliptic_curves(v, B=100, filename=None):
    from hmf import HilbertModularForms
    F = open(filename,'a') if filename else None
    for N in ideals_of_norm(v):
        H = HilbertModularForms(N)
        for i, E in enumerate(H.elliptic_curve_factors()):
            v = E.aplist(B)
            s = '%s\t%s\t%s\t%s'%(N.norm(), no_space(canonical_gen(N)), i, ' '.join([no_space(x) for x in v]))
            print s
            if F:
                F.write(s+'\n')
                F.flush()
                              
def elliptic_curves_parallel(v, B, dir, ncpu=16):
    from hmf import HilbertModularForms
    from sage.all import parallel, set_random_seed
    import os
    @parallel(ncpu)
    def f(N):
        set_random_seed(0)  # to replicate any linear algebra issues?
        level = no_space(canonical_gen(N)).replace('*','').replace('-','_')
        F = open(os.path.join(dir,'%s.txt'%level),'w')
        H = HilbertModularForms(N)
        level = no_space(canonical_gen(N))
        try:
            D = H.elliptic_curve_factors()
            F.write('count %s %s %s\n'%(N.norm(), level, len(D)))
            F.flush()
            for i, E in enumerate(D):
                v = E.aplist(B)
                s = '%s\t%s\t%s\t%s'%(N.norm(), level, i, ' '.join([no_space(x) for x in v]))
                print s
                F.write(s+'\n')
                F.flush()
        except Exception, msg:
            F.write('exception %s %s "%s"\n'%(N.norm(), level, msg))
        F.close()

    for X in f(ideals_of_norm(v)):
        print X
        
#################################################################

from sage.libs.all import pari
from sage.rings.all import primes

def primes_of_bounded_norm(B):
    """
    Return the prime ideals of the integers of the field Q(sqrt(5)) of
    norm at most B, ordered first by norm, then by the image of the
    golden ratio mod the prime in GF(p)={0,1,...,p-1}.

    INPUT:

        - B -- integer

    OUTPUT:

        - list of prime ideals

    EXAMPLES::

        sage: import psage
        sage: psage.modform.hilbert.sqrt5.primes_of_bounded_norm(4)
        [Fractional ideal (2)]
        sage: len(psage.modform.hilbert.sqrt5.primes_of_bounded_norm(10^4))
        1233
        sage: v = psage.modform.hilbert.sqrt5.primes_of_bounded_norm(11); v
        [Fractional ideal (2), Fractional ideal (2*a - 1), Fractional ideal (3), Fractional ideal (3*a - 1), Fractional ideal (3*a - 2)]

    Check that the sort order is as claimed::
    
        sage: P0 = v[-2]; P1 = v[-1]
        sage: K = P0.number_field(); K
        Number Field in a with defining polynomial x^2 - x - 1
        sage: P0.residue_field()(K.gen())
        4
        sage: P1.residue_field()(K.gen())   # yep, 4 < 8
        8
    """
    v = []
    g = F.gen()
    for p in primes(B+1):
        if p == 5:
            v.append((5, F.ideal(2*g-1)))
        elif p%5 in [2,3]:
            Norm = p*p
            if Norm <= B:
                v.append((Norm, F.ideal(p)))
        else:
            s = pari(5).Mod(p).sqrt()
            a = int(((1 + s)/2).lift()); b = int(((1 - s)/2).lift())
            v.append((p, a, F.ideal([p, g - a])))
            v.append((p, b, F.ideal([p, g - b])))
    v.sort()    
    return [z[-1] for z in v]