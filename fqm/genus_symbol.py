from sage.all import ZZ, Zmod, sys, parallel, is_prime, colors, cached_function, Integer, Partitions, Set, QQ, RR, is_prime_power, next_prime, prime_range, is_squarefree, uniq, MatrixSpace, kronecker, CC, exp, walltime, RealField, floor, pari, pi, ComplexField, sqrt, text, arrow, is_even, squarefree_part, polygon2d, CyclotomicField, is_odd, is_even, is_prime, CartesianProduct, prod, log, gcd, sign, valuation, binomial
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from psage.modform.weilrep_tools.dimension import VectorValuedModularForms
from sage.misc.decorators import options
from sage.misc.flatten import flatten
from sage.misc.latex import LatexExpr
from sage.modular.dims import dimension_cusp_forms
from copy import copy, deepcopy
from sage.parallel.decorate import *
from sage.misc.cachefunc import *
import itertools
from sfqm.tools import BB, h

class GenusSymbol(object):

    _K = CyclotomicField(8)
    _z = _K.gen()

    def __init__(self, s='1^+1', reduce_symbol=True):
        if isinstance(s, dict):
            self._symbol_dict = s
            # print "in init: ", s
            if not self._is_valid():
                raise ValueError
        else:
            if s == '' or s == '1' or s == '1^1' or s == '1^+1':
                self._symbol_dict = {}
            else:
                self._from_string(s)
        #self.__reset_sigma = False
        self.__sigma = dict()
        if reduce_symbol:
            self._reduce()
        self._B = dict()

    def canonical_symbol(self):
        raise NotImplementedError

    def defines_isomorphic_module(self, other):
        if isinstance(other, str):
            other = GenusSymbol(other)
        elif not isinstance(other, GenusSymbol):
            raise TypeError
        if self == other:
            return True
        if self.group_structure() != other.group_structure():
            return False
        divs = self.level().divisors()
        for t in divs:
            if self.char_invariant(t) != other.char_invariant(t):
                return False
        return True

    def finite_quadratic_module(self):
        r"""
          Returns the finite quadratic module corresponding to this genus symbol.
        """
        return FiniteQuadraticModule(str(self))

    def p_rank(self, p):
        if not is_prime(p):
            return Integer(0)
        if not self._symbol_dict.has_key(p):
            return Integer(0)
        return sum([s[1] for s in self._symbol_dict[p]])

    def order(self):
        if self._symbol_dict == {}:
            return Integer(1)
        else:
            return Integer(prod([p ** (s[0] * s[1]) for p, l in self._symbol_dict.iteritems() for s in l]))

    @cached_method
    def group_structure(self):
        l = []
        # print "type(level): ", type(self.level())
        # print self.level().prime_factors()
        for p in self.level().prime_factors():
            # print p, type(p)
            for c in self.jordan_components(p):
                # print c, type(c)
                for r in range(c.p_rank(p)):
                    # print r, c.valuation(p)
                    l.append(p ** c.valuation(p))
        return l

    def torsion(self, m):
        fac = self._torsion_factors(m)
        if fac == None:
            return 1
        else:
            return prod([p ** r for p, r in fac])

    def _torsion_factors(self, m):
        if not is_prime_power(m):
            raise NotImplementedError
        p, r = Integer(m).factor()[0]
        J = self.jordan_components(p)
        if J == [GenusSymbol('1^+1')]:
            return None
        nt = list()
        for s in J:
            n = s.valuation(p)
            nt.append([p, s.p_rank(p) * min(n, r)])
        return nt

    def sigma(self, s):
        res = Integer(1)
        for p in self.level().prime_factors():
            res = res * (Integer(1) + sum([p ** (r * s + self._torsion_factors(p ** r)[0][1] / Integer(2))
                                  for r in range(1, self.level().valuation(p) + Integer(1))]))
        return res

    def oddity(self):
        if not self._symbol_dict.has_key(2):
            return 0
        else:
            return (sum([_[4] for _ in self._symbol_dict[2]]) + 4 * len([_ for _ in self._symbol_dict[2] if is_odd(_[0]) and _[2] == -1])) % 8

    @cached_method
    def dimension_estimate_for_anisotropic(self, k, use_simple_est=False):
        N = self.level()
        A2 = self.torsion(2)
        A3 = self.torsion(3)
        A = self.order()
        d = RR(A + A2) / 2
        sqA2 = RR(sqrt(A2))
        sqA3 = RR(sqrt(A3))
        sqA = RR(sqrt(A))
        s = self.sigma(-1)
        #a5est = sqA/(2*pi)*(3/2+ln(N))*(s-sqA/N)
        # print "d = %d, A2 = %d, A3 = %d, s = %d, a5est = %d"%(d,A2,A3,s,
        # a5est)
        return d * RR(k - 1) / RR(12) - sqA2 / 4 - RR(1 + sqA3) / RR(3 * sqrt(3)) - A2 / RR(8) - RR(1) / RR(2) \
               - self.beta_est(use_simple_est) / RR(2)

    def beta(self):
        x = [BB(a) * m for a, m in self.values().iteritems()]
        return sum(x)

    def gaussum(self, s, p=None, F=CC):
        o = self.order() if p == None else p ** self.order().valuation(p)
        return F(self.char_invariant(s, p)[0] * self.char_invariant(s, p)[1] * o)

    def eps(self, p=None):
        a = self.gaussum(1, p)
        return a / abs(a)

    def chi(self, n, p=None):
        return self.gaussum(n, p) / self.gaussum(1, p)

    def beta_formula(self, debug=0):
        prec = max(log(self.order(), 2), 53)
        RF = RealField(prec)
        A = self.order()
        q = 1
        N = Integer(self.level())
        if debug > 0:
            print "N={0}".format(N)
        for p in N.prime_factors():
            c = self.jordan_component(p)._symbol_dict[p][0]
            if len(c) == 1:
                if len(c[0]) == 3:
                    if c[0][1] == 2:
                        q = q * p
                elif c[0][3] == 0 and c[0][1] == 2:
                    q = q * 2
        if debug > 2:
            print A, N, q
        r1 = lambda d: (-1) ** (len(Integer(gcd(d, q)).prime_factors()))

        R2 = self.p_rank(2)
        if debug > 0:
            print "R2 = {0}".format(R2)

        def r2(d):
            if R2 == 0:
                return 1
            elif is_odd(d):
                return 2 ** (R2 - 2)
            else:
                return 2 ** (floor(R2 / 2) - 1)

        def r3(d):
            R3 = Integer(
                len([_ for p in Integer(d / gcd(d, q)).prime_factors() if p % 4 == 3]))
            if is_odd(R3):
                return kronecker(-4, R3)
            else:
                return (-1) ** (R3 / 2)

        def eps_odd(d):
            dd = odd_part(d / gcd(d, q))
            a = kronecker(N / d, odd_part(dd))
            sig = lambda p: self.jordan_component(p)._symbol_dict[p][0][2]
            a *= prod(sig(p) * kronecker(8, p) for p in dd.prime_factors())
            return a

        delta = kronecker(-4, R2)

        def eps2(d):
            d = Integer(d)
            if N.valuation(2) == 1:
                dd = d / (gcd(d, q))
                if d % 4 == 2:
                    if dd % 4 == 3:
                        return 1
                    else:
                        return 0
                elif d % 4 == 0:
                    return 2

            if is_odd(d) or d % 4 == 2:
                return 1

            N2 = self.level() / odd_part(self.level())
            dd = d / ((gcd(N2, d)) * gcd(d, q))
            t = (self.oddity() / R2) % 8
            s = len(self.jordan_components(2))
            if s == 2:
                if debug > 0:
                    print "s = 2"
                t1 = self.jordan_component(2)._symbol_dict[2][0][4]
                t2 = self.jordan_component(4)._symbol_dict[2][0][4]
                t = t2 / t1 % 8
                if debug > 1:
                    print "t = {0}".format(t)
                if d.valuation(2) == 1 or d.valuation(2) == 2:
                    return 0
                dd = d / (8 * gcd(d, q))
                if debug > 1:
                    print "d = {0}, dd =  {1}".format(d, dd)
                if dd % 4 == 3:
                    if (t + 1) % 8 in [2, 6]:  # d.valuation(2) == 0:
                        return 0
                    a = kronecker(8, N / d) * kronecker(2, t)
                else:
                    if (t + 1) % 4 == 0 and dd.valuation(2) == 0:
                        if debug > 1:
                            print "returning 0 for t = {0} and d = {1}".format(t, d)
                        return 0
                    a = kronecker(-8, N / d) * kronecker(2, t)
                    if debug > 1:
                        print "a = {0}".format(a)
            else:
                if dd % 4 == 3:
                    a = delta * \
                        kronecker(
                            N2 / gcd(N2, N / d), (N / d) / gcd(N2, N / d) * t)
                else:
                    a = kronecker(-N2 / gcd(N2, N / d),
                                  (N / d) / gcd(N2, N / d) * t)
            if debug > 1:
                print "d = {0}, eps2(d) = {1}, delta = {2}".format(d, a, delta)
            return a

        eps = lambda d: eps2(d) * eps_odd(d) * r1(d) * r2(d) * r3(d)

        if debug > 0:
            for d in N.divisors():
                if d * gcd(d, q) % 4 in [0, 3]:
                    print "d = {0}, r1(d) = {1}, r2(d) = {2}, r3(d) = {3}, eps_odd(d) = {4}, eps2(d) = {5}, H(-d(q,d)) = {6}".format(d, r1(d), r2(d), r3(d), eps_odd(d), eps2(d), h(d, prec))
        return -sum(q / gcd(d, q) * eps(d) * h(d * gcd(d, q), prec) for d in N.divisors() if d * gcd(d, q) % 4 in [0, 3])

    def beta_est(self, simple_est=False, debug=0):
        prec = max(log(self.order(), 2), 53)
        RF = RealField(prec)
        A = self.order()
        q = Integer(1)
        N = self.level()
        if debug > 0:
            print "N={0}".format(N)
        for p in N.prime_factors():
            c = self.jordan_component(p)._symbol_dict[p]
            if len(c) == 1:
                if len(c[0]) == 3:
                    if c[0][1] == 2:
                        q = q * p
                elif c[0][3] == 0 and c[0][1] == 2:
                    q = q * 2

        R2 = self.p_rank(2)
        if debug > 0:
            print "R2 = {0}".format(R2)

        def r2(d):
            if R2 == 0:
                return 1
            elif is_odd(d):
                return 2 ** (R2 - 2)
            else:
                return 2 ** (floor(R2 / 2) - 1)
        if simple_est:
            if debug > 0:
                print "simple estimate for {}".format(self)
            hh = lambda d: log(d) * sqrt(d) / RF.pi()
        else:
            hh = h
        return sum(q / gcd(d, q) * r2(d) * hh(d * gcd(d, q)) for d in N.divisors() if d * gcd(d, q) % 4 in [0, 3])

    #@cached_method
    def sigmag(self):
        N = self.level()
        res = 1
        for p in self.level().prime_factors():
            resp = 0
            for c in self.jordan_components(p):
                for r in range(0, c.level().valuation(p) + 1):
                    g = c.char_invariant(p ** r)
                    # print g
                    resp = resp + abs(imag_part(CC(g[0])) / RR(p ** r))
            res = res * resp
        return res * RR(sqrt(self.order()))

    def _component_values(self, p):
        vals = []
        J = self.jordan_components(p)
        for c in J:
            if p == 2:
                M = c.finite_quadratic_module()
                vals.append(M.values())
                continue
            n = c.valuation(p)
            q = p ** n
            if c.jordan_component(q)._symbol_dict[p][0][2] * kronecker(2, p) ** c.p_rank(p) == -1:
                if p % 4 == 3:
                    a = -1
                else:
                    a = primitive_root(p)
            else:
                a = 1
            # print k
            #vc[0] = p**k
            # print vc
            #vc[Integer(a % p**n)/p**n] = 1
            vc = self._rank_one_values(p, n, a)
            r = c.p_rank(p)
            # print r
            if r > 1:
                # print p, n, c.jordan_component(q)[0][2]
                if p != 2 and is_even(r) and c.jordan_component(q)._symbol_dict[p][0][2] != 1:
                    vc2 = self._rank_one_values(p, n, 1)
                    r -= 1
                else:
                    vc2 = None
                values_comb = itertools.combinations_with_replacement(vc, r)
                vc = self._combine_values(vc, values_comb, r, q)
                # print vc, vc2
                if vc2 != None:
                    vcnew = {}
                    for v1 in vc:
                        for v2 in vc2:
                            # print v1,v2
                            val = Integer((v1 + v2) * q % q) / q
                            mult = vc[v1] * vc2[v2]
                            if not vcnew.has_key(val):
                                vcnew[val] = mult
                            else:
                                vcnew[val] += mult
                    vc = vcnew
            vals.append(vc)
        return vals

    def _combine_values(self, vc, values_comb, r, q):
        vcnew = {}
        for v in values_comb:
            val = Integer(sum(vv * q for vv in v) % q) / q
            mult = prod(vc[vv] for vv in v)
            R = r
            for vv in uniq(v):
                a = v.count(vv)
                mult *= binomial(R, a)
                R -= a
            # print val, mult
            if not vcnew.has_key(val):
                vcnew[val] = mult
            else:
                vcnew[val] += mult
        return vcnew

    @cached_method
    def values(self):
        ps = self.level().prime_factors()
        levels = []
        for p in ps:
            for c in self.jordan_components(p):
                levels.append(c.level())
        vc = flatten([self._component_values(p) for p in ps], max_level=1)
        # print "vc = ", vc
        only_vals = map(lambda x: x.keys(), vc)
        vals = itertools.product(*only_vals)
        # print list(vals)
        vccomb = {}
        N = self.level()
        # print vals
        for v in vals:
            val = sum(v[i] for i in range(len(levels)))
            val = (val * N % N) / N
            mult = prod(vc[i][v[i]] for i in range(len(v)))
            if not vccomb.has_key(val):
                vccomb[val] = mult
            else:
                vccomb[val] += mult
        return vccomb

    def _rank_one_values(self, p, n, a):
        # print "a=", a
        q = p ** n
        vc = {}
        if is_even(n):
            k = Integer(n) / 2
            kk = k
        else:
            k = Integer(n - 1) / 2
            kk = k + 1
        if p != 2:
            Q1 = [(x * x) % p for x in range(1, p / 2 + 1)]
            Q = [Integer(x) for x in xrange(q) if x % p in Q1]
        else:
            Q1 = [1, 7]
            Q = [Integer(x) for x in xrange(q) if x % 8 in Q1]
        Q = Q + uniq([(m * p ** (2 * r)) % q for r in range(1, kk) for m in Q])
        Q.append(0)
        # print "Q=", Q
        # print Q1, Q
        for x in Q:
            nom = Integer((x * a) % q)
            val = Integer((x * a) % q) / q
            if not p.divides(nom):
                mult = 2
            else:
                if nom == 0:
                    mult = p ** k
                else:
                    r = nom.valuation(p)
                    mult = 2 * p ** (r / 2)
            # print val
            if not vc.has_key(val):
                vc[val] = mult
            else:
                vc[val] += mult
        return vc

    def char_invariant(self, s, p=None):
        r"""
        If this quadratic module equals $A = (M,Q)$, return
        the characteristic function of $A$ (or $A(p)$ if $p$ is a prime)
        at $s$, i.e. return
        $$\chi_A (s)= |M|^{-1}\sum_{x\in M} \exp(2\pi i s Q(x))).$$

        NOTE
            We apply the formula in [Sko, Second Proof of Theorem 1.4.1].

        EXAMPLES NONE
        """
        # manual caching
        s = s % self.level()
        if self.__sigma.has_key((s, p)):
            return self.__sigma[(s, p)]

        #t = cputime()
        if s == 0:
            return 1, 1
        if not p is None and not is_prime(p):
            raise TypeError
        if p and 0 != self.order() % p:
            return 1, 1
        _p = p

        K = self._K
        z = self._z

        jd = self.jordan_decomposition(flat = True)

        ci = ci1 = 1
        for p, c in jd:
            # c: basis, ( prime p,  valuation of p-power n, dimension r,
            # determinant d over p [, oddity o])
            n = c.valuation(p)
            r = c.p_rank(p)
            d = c._eps(p**r)
            #print p,n,r,d
            if _p and p != _p:
                continue
            o = None if p != 2 else (0 if c.is_even() else c._symbol_dict[2][0][4])
            # print o
            odd = c.is_odd()
            k = valuation(s, p)
            s1 = Integer(s / p ** k)
            h = max(n - k, 0)
            # print h
            q = p ** h
            if p != 2:
                lci = z ** ((r * (1 - q)) % 8) * d ** (h % 2) if h > 0 else 1
                lci1 = q ** (-r) if h > 0 else 1
            elif k == n and odd:
                # print "t!"
                return 0, 0
            else:
                f = z ** o if odd else 1
                lci = f * d ** (h % 2) if h > 0 else 1
                lci1 = q ** (-r) if h > 0 else 1
                # print f, d, lci
            if 2 == p:
                lci = lci ** s1
            # print "lci=",lci
            ci *= lci * kronecker(s1, 1 / lci1)
            ci1 *= lci1
        v = (ci, QQ(ci1).sqrt())
        self.__sigma[(s, p)] = v
        return v

    @cached_method
    def two_torsion_values(self):
        J = self.jordan_components(2)
        vals = dict()
        for s in J:
            svals = {}
            if s.valuation(2) >= 3 or (s.is_even() and s.valuation(2) >= 2):
                svals = {0: 2 ** s.p_rank(2)}
            else:
                if s.is_even():
                    if s.valuation(2) == 1:
                        n = s.p_rank(2) / 2
                        ll = [3 ** (n - 2 * k) * binomial(n, 2 * k)
                              for k in range(0, n / 2 + 1)]
                        oh = sum(ll)
                        if s.oddity() == 4:
                            svals = {
                                0: 2 ** s.p_rank(2) - oh, QQ(1) / QQ(2): oh}
                        else:
                            svals = {
                                0: oh, QQ(1) / QQ(2): 2 ** s.p_rank(2) - oh}
                else:
                    # TODO: improve this part as well
                    M = s.finite_quadratic_module()
                    svals = M.kernel_subgroup(2).as_ambient()[0].values()

            if vals == {}:
                vals = svals
            else:
                vals_new = dict()
                for k, m in svals.iteritems():
                    for kk, mm in vals.iteritems():
                        den = lcm(QQ(k).denominator(), QQ(kk).denominator())
                        v = (((k + kk) * den) % den) / Integer(den)
                        if not vals_new.has_key(v):
                            vals_new[v] = m * mm
                        else:
                            vals_new[v] += m * mm
                vals = vals_new
        return vals

    def valuation(self, p):
        r'''Return the p-valuation of the level of self
        '''
        if not self._symbol_dict.has_key(p):
            return 0
        if not isinstance(self._symbol_dict[p], list):
            return 0
        if len(self._symbol_dict[p]) == 0:
            return 0
        return max(_[0] for _ in self._symbol_dict[p])

    @cached_method
    def level(self):
        l = Integer(prod([p ** self.valuation(p)
                          for p in self._symbol_dict.keys()]))
        # print l
        if self.p_rank(2) > 0:
            ll = [_[0] for _ in self._symbol_dict[2] if _[3] == 1]
            if ll == []:
                return l
            if max(ll) == self.valuation(2):
                l = l * 2
        return l

    @cached_method
    def signature(self):
        sig = 0
        if len(self._symbol_dict) == 0:
            return 0
        for p, l in self._symbol_dict.iteritems():
            if p == 2:
                for s in l:
                    sig = sig + s[4] + \
                        (4 if is_odd(s[0]) and s[2] == -1 else 0)
                    sig = sig % 8
            else:
                for s in l:
                    eps = s[2]
                    n = s[1]
                    r = s[0]  # symbol= q^eps*n, where q = p^r
                    # print r,n,eps
                    if r % 2 == 0:
                        continue
                    sig += n * (1 - p ** r) % 8 + \
                        (4 * r if eps == -1 else 0) % 8
        return Integer(sig % 8)

    @cached_method
    def is_simple(self, k, no_inv=False, aniso_formula=False, reduction=True):
        d = self.dimension_cusp_forms(
            k, no_inv, aniso_formula, test_positive=True, reduction=reduction)
        if d < 0:
            raise ValueError("Negative dimension for {0}".format(self))
        if d == 0:
            return True
        else:
            return False

    def is_global(self, r, s):
        r""" Checks if this symbol can be realized as
             the genus symbol of an even integral lattice
             of type (r,s).
        """
        D = self.order() * (-1) ** s
        # print D
        if (r - s) % 8 != self.signature():
            return False
        for p in self.level().prime_factors():
            if self.p_rank(p) > r + s:
                return False
            elif self.p_rank(p) == r + s:
                eps = self._eps(p)
                a = D / p ** (self.order().valuation(p))
                if p == 2:
                    A2 = self.jordan_component(2)
                    if A2.is_odd():
                        continue
                if not eps == kronecker(a, p):
                    return False
        return True

    def is_global_unique(self, r, s):
        r"""
          Returns True if there is only one lattice in the genus of a
          global realization of self. We apply Theorem 1.13.2 of [Ni].
          This is a sufficient but not a necessary condition.
          Therefore, we return None if the condition is not satisfied.
        """
        if not self.is_global(r,s):
            raise RuntimeError("The finite quadratic module is not global!")
        r = Integer(r)
        s = Integer(s)

        nonstr = "The criterion is not satisfied. We cannot decide if the finite quadratic module is unique."

        if not (r >= 1 and s >= 1 and r+s >= 3):
            print nonstr
            return None
        
        for p in self.level().prime_factors():
            satisfied = True
            if not (r+s >= 2 + self.p_rank(p)):
                satisfied = False
            if p > 2 and not satisfied:
                found = False
                for j in self.jordan_components(p):
                    if j.p_rank(p) >= 2:
                        found = True
                        break
                if not found:
                    print nonstr
                    return None
            if p == 2 and not satisfied:
                found = False
                for j in self.jordan_components(2):
                    if j.is_odd():
                        jj = self.jordan_component(2**(j.valuation(2)+1))
                        if jj.p_rank(p) > 0 and jj.is_odd():
                            found = True
                            break
                    if j.p_rank(2) >= 2:
                        try:
                            j - U(2**(j.valuation(2)))
                            found = True
                            break
                        except ValueError:
                            try:
                                j - V(2**(j.valuation(2)))
                                found = True
                                break
                            except ValueError:
                                continue
                if not found:
                    print nonstr
                    return None
        return True

    def _eps(self, p):
        if not self._symbol_dict.has_key(p):
            return 1
        l = self._symbol_dict[p]
        if len(l) == 0:
            return 1
        return Integer(prod(s[2] for s in l))

    @cached_method
    def dimension_cusp_forms(self, k, no_inv=False, aniso_formula=False, test_positive=False, reduction=False):
        r"""
          Computes the dimension of the space of cusp forms of weight $k$
          for the Weil representation associated with this finite quadratic module.

          INPUT:
          - k: a half-integer, the weight
          - no_inv: bool, if True: we do not compute invariants in the case of weight 2 (so assume they are 0)
          - aniso_formula: bool, if True use the formula involving class numbers to compute the dimension.
                           works only for anisotropic modules. The formula is given in Theorem XX of [BEF].
          - test_positive: bool, if True, only test if the dimension is positive
                           returns a positive number if this is the case,
                           which is not necessarily equal to the dimension.
          - reduction: bool, if True, we use reduction mod p to compute the invariants
                             of the Weil representation, whose dimension is needed for $k = 3/2$ and $k=2$.

          OUTPUT:
          - an Integer, the dimension of the space of cusp forms of weight $k$
          for the Weil representation associated with this finite quadratic module

          SEE ALSO:
          - For more details, see the implementation in PSAGE.
        """
        s = str(self)
        if not s == '1^+1':
            #M = FiniteQuadraticModule(s)
            V = VectorValuedModularForms(
                s, True, aniso_formula=aniso_formula, use_reduction=reduction)
            d = V.dimension_cusp_forms(
                k, no_inv=no_inv, test_positive=test_positive)
        else:
            d = dimension_cusp_forms(1, k)
        return d

    def is_even(self):
        r'''
         Returns True if there is no odd 2-adic Jordan component.
        '''
        if self.p_rank(2) == 0:
            return True
        else:
            for s in self._symbol_dict[2]:
                if s[3] == 1:
                    return False
        return True

    def is_odd(self):
        return not self.is_even()

    def jordan_components(self, p):
        r'''
        Returns all p-adic Jordan components of self as a list GenusSymbols.
        '''
        if not self._symbol_dict.has_key(p):
            return [GenusSymbol()]
        self._reduce()
        return [GenusSymbol({p: [x]}) for x in self._symbol_dict[p]]

    def jordan_decomposition(self, flat = False):
        l = {}
        for p in self.level().prime_factors():
            l[p] = [c for c in self.jordan_components(p)]
        if not flat:
            return l
        else:
            return [(p,c) for c in l[p] for p in l.keys() ]

    def jordan_component(self, q):
        if not is_prime_power(q):
            raise ValueError
        else:
            p, n = list(Integer(q).factor())[0]

        if not self._symbol_dict.has_key(p):
            return GenusSymbol('1^+1')
        else:
            r = []
            for s in self._symbol_dict[p]:
                if s[0] == n:
                    r.append(s)
            if len(r) == 0:
                return GenusSymbol('1^+1')
            else:
                return GenusSymbol({p: r})

    def rank_of_jordan_component(self, q):
        c = self.jordan_component(q)
        if c == []:
            return 0
        return sum([s[1] for s in c])

    def max_rank(self):
        r = 0
        for p in self.level().prime_factors():
            r = max(r, self.p_rank(p))
        return r

    def B(self, m):
        r'''
        Computes the B-Set of this symbol by brute force.
        '''
        return map(lambda x: GenusSymbol(x), Bbf(self.finite_quadratic_module(), m))

    #@cached_method
    def C(self, m, unique=True):
        r'''
          Computes the C-Set of this genus symbol.

          INPUTS:
          - m: an integer, currently assumed to be prime.
          - unique: boolean, return only non-isomorphic modules.

          ALGORITHM:

            See [BEF].
        '''
        CS = C(self, m)
        if unique:
            CU = []
            for s in CS:
                cont = False
                for t in CU:
                    if t.defines_isomorphic_module(s):
                        cont = True
                if not cont:
                    CU.append(s)
            CS = CU
        return CS

    def _test_C(self, m, all=False):
        r"""
          Test the implementation of the C-sets


          INPUTS:
          - m: an integer
          - all: bool,
                 - if False, we check if all modules in C(A,m)
                   are actually contained in B(A,m),
                   by checking if they are obtained as a quotient.
                 - if True, we check if $C(A,m) \subset B(A,m)$
                   using the brute-force implementation of B(A,m).
                   We also check correctness for $p=2$
                   and level $2, 4, 8$.
                   

          OUTPUT:
          
            Returns True, if the test succeeded, otherwise False.
        """
        CS = self.C(m)
        M = self.finite_quadratic_module()
        if not all:
            for s in CS:
                N = s.finite_quadratic_module()
                q = False
                for G in N.subgroups():
                    if G.is_isotropic() and G.order() == m:
                        if (N / G).is_isomorphic(M):
                            q = True
                if not q:
                    print "Ooops, this one seems to be wrong: ", s
                return False
        else:
            BS = self.B(m)
            print BS
            print len(BS), len(CS)
            if self.level() in [2, 4, 8] and self.p_rank(2) <= 3:
                if len(BS) > len(CS):
                    print 'Missing symbols!'
                    return False
            for c in CS:
                found_isomorphic = False
                for b in BS:
                    if b.defines_isomorphic_module(c):
                        found_isomorphic = True
                        break
                if not found_isomorphic:
                    print c, " is not isomorphic to any module in s.B(" + str(m) + ")!"
                    return False

        return True

    def _from_string(self, s):
        r'''
           Initializes a GenusSymbol from a string.
           Most parts are copied from finite_quadratic_module.py
           in psage. This should be made work nicely together instead.
        '''
        sl = s.split('.')
        d = dict()
        for s in sl:
            L1 = s.split('^')
            if len(L1) > 2:
                raise ValueError()
            elif len(L1) == 1:
                nn = 1
            else:
                nn = L1[1]
            n = Integer(nn)
            L1 = L1[0].split("_")
            q = Integer(L1[0])
            if len(L1) > 2:  # more than one _ in item
                raise ValueError
            elif len(L1) == 2:
                if Integer(L1[1]) in range(8):
                    t = Integer(L1[1])
                else:
                    raise ValueError, "Type given, which ist not in 0..7: %s" % (
                        L1[1])
            else:
                t = None
            if not (n != 0 and q != 1 and q.is_prime_power()
                    and (None == t or (is_even(q) and t % 2 == n % 2))
                    and (not (None == t and is_even(q)) or 0 == n % 2)
                    ):
                raise ValueError, "{0} is not a valid signature!".format(s)
            p = q.prime_factors()[0]
            r = q.factor()[0][1]
            eps = sign(n)
            n = abs(n)
            if not d.has_key(p):
                d[p] = list()

            c = None
            if p == 2:
                c = None
                if t == None:
                    d[p].append([r, n, eps, 0, 0])
                else:
                    d[p].append([r, n, eps, 1, t % 8])
            else:
                d[p].append([r, n, eps])
        self._symbol_dict = d
        if not self._is_valid():
            raise ValueError

    def _is_valid(self):
        r"""
          Determines if the _symbol_dict dictionary actually defines
          a finite quadratic module.
        """
        for p, l in self._symbol_dict.iteritems():
            if not is_prime(p):
                return False
            if p != 2:
                if not isinstance(l, list):
                    return False
                for s in l:
                    if not isinstance(s, list) or len(s) != 3:
                        return False
                    r, n, eps = s
                    if not eps in [1, -1]:
                        return False
        if not self._symbol_dict.has_key(2):
            return True
        for r, n, eps, tp, t in self._symbol_dict[2]:
            # print r,n,eps,tp,t
            if tp == None:
                tp = 0
            if tp == 0 and n % 2 == 1:
                raise ValueError("Not a valid signature!")
            c = None
            if not t in range(8):
                return False
            if t == None or t == 0:
                if not is_even(n):
                    return False
            if tp == 0 and t != 0:
                return False
            if tp != 0:
                if 1 == n:
                    if not eps * n == kronecker(t, 2):
                        return False
                if n > 1:
                    c = None
                    CP = eval(
                        "CartesianProduct(" + "[1,3,5,7]," * (n - 1) + ")")
                    # TODO: find better algorithm
                    for x in CP:
                        s = sum(x) % 8
                        if kronecker(prod(x) * (t - s), 2) == eps:
                            x.append(t - s)
                            c = x
                            break
                    # print "c=", c
                    if not c:
                        return False
        return True

    def _to_string(self):
        if self._symbol_dict == {}:
            return '1^+1'
        symstr = ''
        for p, lsym in sorted(self._symbol_dict.iteritems()):
            for s in lsym:
                if s[0] == 0:
                    continue
                if len(symstr) != 0:
                    symstr = symstr + '.'
                symstr = symstr + str(p ** s[0])
                if p == 2:
                    sgn = '+' if s[2] == 1 else '-'
                    if s[3] == 1:
                        symstr = symstr + '_' + str(s[4])
                    symstr = symstr + '^' + sgn + str(s[1])
                    # else:
                    #    symstr = symstr + '^' + str(s[1])
                else:
                    sgn = '+' if (s[2] == 1) else '-'
                    symstr = symstr + '^' + sgn + str(s[1])
        return symstr

    def _reduce(self):
        sym = self._symbol_dict
        for p, sl in sym.iteritems():
            for s in sl:
                if s[1] == 0:
                    sl.remove(s)
                if s[0] == 0:
                    sl.remove(s)
            if p != 2:
                sl.sort()
                for s in sl:
                    # print s
                    ssl = copy(sl)
                    for j in range(sl.index(s) + 1, len(ssl)):
                        t = ssl[j]
                        if s[0] == t[0]:
                            s[1] = s[1] + t[1]
                            s[2] = s[2] * t[2]
                            sl.remove(t)
            else:
                sl.sort()
                for s in sl:
                    # print s
                    ssl = copy(sl)
                    for j in range(sl.index(s) + 1, len(ssl)):
                        t = ssl[j]
                        if s[0] == t[0]:
                            # print s, t
                            if s[3] == None:
                                s[3] = 0
                            if t[3] == None:
                                t[3] = 0
                            s[1] = s[1] + t[1]
                            s[2] = s[2] * t[2]
                            if s[3] != t[3]:
                                s[3] = 1
                                s[4] = t[4] if t[3] == 1 else s[4]
                            else:
                                s[4] = (s[4] + t[4]) % 8
                            # print s
                            sl.remove(t)
                # print sym

    def __repr__(self):
        return "Genus symbol " + self._to_string()

    def __str__(self):
        return self._to_string()

    def __add__(self, other):
        d = deepcopy(self._symbol_dict)
        e = other._symbol_dict
        for p, l in e.iteritems():
            if not d.has_key(p):
                d[p] = deepcopy(e[p])
            else:
                d[p] = deepcopy(d[p] + e[p])
        # print "in add: ", d_new
        s = GenusSymbol(d)
        s._reduce()
        return s

    def __sub__(self, other):
        r"""
          Try to formally 'subtract' ```other``` (B) from ```self``` (A).
          That is, check if
          \[
            A = B \oplus C
          \]
          holds for some finite quadratic module $C$ and if this is the case,
          returns $C$.

          Otherwise, returns an error.
        """
        debug = 0
        # print "In sub"
        err = ValueError("Result does not define a genus symbol")
        # print "in sub"
        # print "self = ", self, " other = ", other
        s = GenusSymbol(deepcopy(self._symbol_dict))
        # print "s = ", s
        s._reduce()
        d = s._symbol_dict
        e = other._symbol_dict
        # print d, e
        for p, l in e.iteritems():
            if not d.has_key(p):
                raise err
            else:
                for c in e[p]:
                    # print c
                    try:
                        j = s.jordan_component(p ** c[0])._symbol_dict[p]
                        # print j
                    except:
                        raise err
                    if not len(j) == 1:
                        raise ValueError("Strange error...")
                    else:
                        j = j[0]
                    if not j[1] - c[1] >= 0:
                        raise err
                    else:
                        if j[1] == c[1] and not j[2] == c[2]:
                            # same multiplicity, different eps
                            if debug > 0:
                                print 'same multiplicity, different eps'
                            if not (p == 2 and j[0] == 1 and j[4] == (c[4] + 4) % 8 and j[3] == c[3]):
                                raise err
                        if p == 2:
                            if j[1] == c[1]:
                                # same multiplicity
                                if not j[3] == c[3]:
                                    # different types
                                    raise err
                                if not j[4] == c[4]:
                                    # different oddity
                                    if debug > 0:
                                        print 'same multiplicity, different oddity'
                                    if not (j[0] == 1 and j[4] == (c[4] + 4) % 8 and j[2] != c[2] and j[3] == c[3]):
                                        # q=2, oddity differs by 4 and sign has
                                        # changed is ok.
                                        if debug > 0:
                                            print j[0], (j[4] == (c[4] + 4) % 8), (j[2] != c[2], j[2], c[2]), (j[3] == c[3])
                                        raise err
                            # print j[4], c[4], j,c
                            j[4] = (j[4] - c[4]) % 8
                            # print j[4]
                        j[2] = j[2] * c[2]
                        j[1] = j[1] - c[1]
        # print s
        s._reduce()
        if not s._is_valid():
            raise err
        return s

    def __eq__(self, o):
        r"""
          Check for equality by using the string representation.

          NOTE:
          This does not check for isomorphy.
        """
        if str(self) == str(o):
            return True
        else:
            return False

    def __hash__(self):
        r"""
          We use the string representation for hasing.
        """
        # if not hasattr(self,'_hash') or self._hash == None:
        self._reduce
        #l = map(lambda x: (x[0],tuple(map(lambda y: tuple(y),x[1]))),list(self._symbol_dict.iteritems()))
        #self._hash = hash(tuple(l))
        # return hash(tuple(l))
        return hash(str(self))

    def latex(self):
        r"""
          Returns a latex representation of ```self```.
        """
        o = r'$'
        l = str(self).split('.')
        for s in l:
            ss = s.split('^')
            if len(ss) == 2:
                s1, s2 = ss
            else:
                s1 = ss[0]
                s2 = '+1'
            o = o + s1 + r'^{' + s2 + r'}'
        o = o + r'$'
        return LatexExpr(o)


@cached_function
def C(genus_symbol, m):
    r"""
      Return the set C(genus_symbol, m) as defined in [BEF].
    """
    m = Integer(m)
    
    Cs = list()
    if not is_prime(m):
        p = m.prime_factors([0])
        rec = True
    else:
        rec = False
        p = m

    # print "in C", genus_symbol
    Cs = Cs + trivial_rule(genus_symbol, p)
    if genus_symbol.p_rank(p) == 0:
        if not rec:
            return Cs
        else:
            sum(C(s,n) for s in Cs)

    # print Cs
    if p == 2:
        Cs = Cs + two_power_up_rules(genus_symbol)
        Cs = Cs + two_level_4_rules(genus_symbol)
    else:
        Cs = Cs + odd_power_up_rules(genus_symbol, p)
    if not rec:
        return Cs
    else:
        return sum(C(s,n) for s in Bs)


def two_level_4_rules(genus_symbol):
    if genus_symbol.p_rank(2) == 0:
        return []
    J = genus_symbol.jordan_components(2)
    Bs = []
    for s in J:
        if s.valuation(2) == 1:
            try:
                # E4
                # print "E4"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, 1, 0, 0]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, 1, 1, 0]]})
                Bs.append(sym_new)
            except ValueError:
                pass
            try:
                # E5
                # print "E5"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, -1, 0, 0]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, -1, 1, 4]]})
                Bs.append(sym_new)
            except ValueError:
                pass
            try:
                # E6.1
                # print "E6.1"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, 1, 1, 2]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, 1, 1, 2]]})
                Bs.append(sym_new)
            except ValueError:
                pass
            try:
                # E6.2
                # print "E6.2"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, 1, 1, 6]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, 1, 1, 6]]})
                Bs.append(sym_new)
            except ValueError:
                pass
            try:
                # E7.1
                # print "E7.1"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, 1, 1, 2]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, -1, 1, 2]]})
                Bs.append(sym_new)
            except ValueError:
                pass
            try:
                # E7.2
                # print "E7.2"
                sym_new = genus_symbol - GenusSymbol({2: [[1, 2, 1, 1, 6]]})
                sym_new = sym_new + GenusSymbol({2: [[2, 2, -1, 1, 6]]})
                Bs.append(sym_new)
            except ValueError:
                pass
    return Bs


def two_power_up_rules(genus_symbol):
    if genus_symbol.p_rank(2) == 0:
        return []
    J = genus_symbol.jordan_components(2)
    Bs = []
    for s in J:
        # print s
        try:
            # E1
            # print "E1"
            sym_new = genus_symbol - \
                GenusSymbol({2: [[s.valuation(2), 2, 1, 0, 0]]})
            # print 'sym_new', sym_new
            sym_new1 = sym_new + \
                GenusSymbol({2: [[s.valuation(2) + 1, 2, 1, 0, 0]]})
            Bs.append(sym_new1)
            # print sym_new, J
        except ValueError:
            pass
        try:
            # E2
            # print "E2"
            sym_new = genus_symbol - \
                GenusSymbol({2: [[s.valuation(2), 2, -1, 1, 4]]})
            # print 'sym_new', sym_new
            sym_new1 = sym_new + \
                GenusSymbol({2: [[s.valuation(2) + 1, 2, -1, 0, 0]]})
            Bs.append(sym_new1)
            # print sym_new, J
        except ValueError:
            pass
        # And now E3
        for o in [1, 3, 5, 7]:
            # print "trying o = %d for s = %s"%(o,s)
            try:
                # E3
                # print "E3:", genus_symbol
                s1 = GenusSymbol(
                    {2: [[s.valuation(2), 1, kronecker(o, 2), 1, o]]})
                # print s1
                sym_new = genus_symbol - s1
                # print genus_symbol
                # print 'sym_new', sym_new
                sym_new1 = sym_new + \
                    GenusSymbol(
                        {2: [[s.valuation(2) + 2, 1, kronecker(o, 2), 1, o]]})
                # print genus_symbol
                Bs.append(sym_new1)
                # print sym_new, J
            except ValueError:
                pass
            # print "Done, first part"
    return Bs


def odd_power_up_rules(genus_symbol, p):
    if not is_odd(p) or not is_prime(p):
        raise ValueError
    Bs = []
    if genus_symbol.p_rank(2) == 0:
        return Bs
    J = genus_symbol.jordan_components(p)
    for s in J:
        # print s
        try:
            sym_new = genus_symbol - GenusSymbol({p: [[s.valuation(p), 1, 1]]})
            sym_new1 = sym_new + GenusSymbol({p: [[s.valuation(p) + 2, 1, 1]]})
            Bs.append(sym_new1)
        except:
            pass
        try:
            sym_new = genus_symbol - \
                GenusSymbol({p: [[s.valuation(p), 1, -1]]})
            sym_new1 = sym_new + \
                GenusSymbol({p: [[s.valuation(p) + 2, 1, -1]]})
            Bs.append(sym_new1)
        except:
            pass
        try:
            sym_new = genus_symbol - GenusSymbol({p: [[s.valuation(p), 2, 1]]})
            sym_new1 = sym_new + \
                GenusSymbol({p: [[s.valuation(p) + 1, 2, ep]]})
            Bs.append(sym_new1)
        except:
            pass
        try:
            sym_new = genus_symbol - \
                GenusSymbol({p: [[s.valuation(p), 2, -1]]})
            sym_new1 = sym_new + \
                GenusSymbol({p: [[s.valuation(p) + 1, 2, -ep]]})
            Bs.append(sym_new1)
        except:
            pass
    return Bs


def trivial_rule(s, p):
    r"""
      Apply the 'trivial' rules to $s$ for the prime $p$,
      that is the rules that are of the form $1^+1 \# mapsto A$
      for a finite quadratic module $A$ of order $p^2$.
    """
    p = Integer(p)
    if not is_prime(p):
        raise ValueError("p={0} has to be a prime number.".format(p))
    if p == 2:
        sp2 = GenusSymbol("2^+2")
        so2 = GenusSymbol("2_0^+2")
        s1 = s + sp2
        s2 = s + so2
        if not s1 == s2:
            return [s1, s2]
        else:
            return [s1]
    else:
        ep = -1 if p % 4 == 3 else 1
        sp = []
        sp.append(GenusSymbol(str(p) + "^" + str(ep * 2)))
        sp.append(GenusSymbol(str(p ** 2) + "^+1"))
        sp.append(GenusSymbol(str(p ** 2) + "^-1"))
        return [s + _ for _ in sp]


@cached_function
def prime_anisotropic_symbols(p, fake=False):
    r"""
      Return a list of all anisotropic symbols of prime level $p$.

      INPUTS:
       - p: A prime number or equal to 1
       - fake: Do not return all symbols,
               but just $p^+1$ and $p^+2$ for odd $p$,
               $1^+1$ for $p = 1$
               and similar for $p = 4$ and $p = 8$.
               This is used to just reflect the possible structures
               of the corresponding finite abelian groups.
    """
    p = Integer(p)
    
    if not (p == 1 or is_prime(p) or p in [4, 8]):
        return []
    else:
        if is_odd(p) and is_prime(p):
            if fake:
                l = [str(p) + '^+1', str(p) + '^+2']
            else:
                if p % 4 == 3:
                    l = [str(p) + '^+1', str(p) + '^-1', str(p) + '^+2']
                else:
                    l = [str(p) + '^+1', str(p) + '^-1', str(p) + '^-2']
        else:
            if p == 1:
                l = ['1^+1']
            elif p == 2:
                l = ['2^-2']
            elif p == 4:
                if fake:
                    l = ['2_1^+1', '2_2^+2', '2_3^+3']
                else:
                    # TODO: are these correct?? (2_6^+2, 2_5^+3 sollte hier
                    # nicht stehen!?!?)
                    l = [
                        '2_1^+1', '2_7^+1', '2_2^+2', '2_6^+2', '2_3^+3', '2_5^+3']
            elif p == 8:
                if fake:
                    l = ['4_1^+1', '2_1^+1.4_1^+1']
                else:
                    l = ['4_1^+1', '4_3^-1', '4_5^-1', '4_7^+1', '2_1^+1.4_1^+1',
                         '2_1^+1.4_3^-1', '2_1^+1.4_5^-1', '2_1^+1.4_7^+1']
        # print l
        return map(lambda x: GenusSymbol(x), l)


@cached_function
def anisotropic_symbols(N, sig=None, fake=False):
    r"""
      Return a list of all anisotropic symbols of given level and signature.
      If ```fake``` is True, then only return all possible structures
      of finite abelian groups.
    """
    N = Integer(N)
    if sig != None and fake == True:
        raise ValueError(
            "You provided a signature and requested fake symbols. This does not make any sense.")
    syms = []
    if N == 1 and (sig == None or sig % 8 == 0):
        return prime_anisotropic_symbols(1)
    if is_prime_power(N):
        syms = prime_anisotropic_symbols(N, fake=fake)
    else:
        p = max(N.prime_factors())
        M = N / p
        q = p ** N.valuation(p)
        syms = anisotropic_symbols(M)
        syms = syms = map(lambda r: r[
                          0] + r[1], itertools.product(syms, prime_anisotropic_symbols(q, fake=fake)))
    if not sig == None:
        syms = filter(lambda s: s.signature() == sig, syms)
    have = []
    syms_new = []
    return syms


def gamma0_N_genus_symbol(N):
    s = GenusSymbol()
    # print s
    N = Integer(2 * N)
    for p in N.prime_factors():
        # print s, p
        if p == 2:
            v = N.valuation(2)
            # print v
            gs = str(2 ** (v)) + "_" + str((N / 2 ** v) % 8) + "^" + \
                ("+1" if kronecker(N / 2 ** v, 2) == 1 else "-1")
            # print gs
            g = GenusSymbol(gs)
            s = s + g
        else:
            v = N.valuation(p)
            Np = p ** v
            np = N / Np
            gs = str(Np) + "^" + ("+1" if kronecker(np, p) == 1 else "-1")
            # print gs
            s = s + GenusSymbol(gs)
    return s

def gamma1_genus_symbol(N):
    N = Integer(N)
    return GenusSymbol(FiniteQuadraticModule([2,N,N],[1/Integer(4),0,0,0,1/N,0]).jordan_decomposition().genus_symbol())

def t1_genus_symbol(a,N):
    a = Integer(a)
    N = Integer(N)
    return GenusSymbol(FiniteQuadraticModule([2*a,N,N],[1/(4*a),0,0,0,1/N,0]).jordan_decomposition().genus_symbol())
