r"""
Class representing a space of vector valued modular forms
transforming with the Weil representation of a discriminant module.

THIS IS JUST A TEST.
...

AUTHORS:

- Stephan Ehlen (2012-11-12): initial version
...

"""

#*****************************************************************************
#       Copyright (C) 2012 Stephan Ehlen <ehlen@mathematik.tu-darmstadt.de>
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
#  as published by the Free Software Foundation
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from psage.modules import *
from sage.all import SageObject, Integer, RR, is_odd, next_prime, floor, \
                     RealField, ZZ, ceil, log, ComplexField, real, sqrt, exp, round, imag
#import sys
from .weight_one_half import *

try:
    from psage.modules.finite_quadratic_module import FiniteQuadraticModule
    from psage.external.weil_invariants.weil_invariants import cython_invariants_dim
except ImportError:
    raise

from psage.modules.finite_quadratic_module import FiniteQuadraticModule

def BB(x):
    RF=RealField(100)
    x = RF(x)
    onehalf = RF(1)/2
    return x - onehalf*(floor(x)-floor(-x))

class VectorValuedModularForms(SageObject):
    r"""
    Class representing a space of vector valued modular forms
    transforming with the Weil representation of a discriminant module.


    EXAMPLES:

        As input we use any input that one can use for a finite quadratic module.

        (For instance a genus symbol.)

        ::

            sage: V=VectorValuedModularForms('2_7^+1.3^-2')
            sage: V
            Vector valued modular forms for the Weil representation corresponding to: 
            Finite quadratic module in 3 generators:
            gens: e0, e1, e2
            form: 3/4*x0^2 + 2/3*x1^2 + 1/3*x2^2

    SEE ALSO:

        :class:`psage.modules.finite_quadratic_modules.FiniteQuadraticModule`
    """

    def __init__(self, A, use_genus_symbols = False, aniso_formula = False, use_reduction = False):
        try:
            from sfqm.fqm.genus_symbol import GenusSymbol
        except ImportError, e:
            print e
            use_genus_symbols = False
            print "Not using genus symbols."
        self._use_reduction = use_reduction
        if use_genus_symbols:
            if isinstance(A, str):
                g = GenusSymbol(A)
            else:
                try:
                    g = GenusSymbol(A.jordan_decomposition().genus_symbol())
                except:
                    raise ValueError
            self._g = g
            n2 = self._n2 = g.torsion(2)
            self._v2 = g.two_torsion_values()
            self._M = None
            self._aniso_formula = aniso_formula
        else:
            self._M = FiniteQuadraticModule(A)
            self._g = None
            self._level = self._M.level()
            self._aniso_formula = False
            if is_odd(self._M.order()):
                self._n2 = n2 = 1
                self._v2 = {0: 1}
            else:
                self._M2 = M2 = self._M.kernel_subgroup(2).as_ambient()[0]
                self._n2 = n2 = self._M2.order()
                self._v2 = self._M2.values()

        if use_genus_symbols:
            self._signature = g.signature()
            m = g.order()
        else:
            self._signature = self._M.signature()
            m = self._M.order()
            
        self._m = m
        self._alpha3 = None
        self._alpha4 = None

    def __repr__(self):
        return "Vector valued modular forms for the Weil representation corresponding to: \n" + self._M.__repr__()

    def finite_quadratic_module(self):
        return self._M

    def dimension(self,k,ignore=False, debug = 0):
        if k < 2 and not ignore:
            raise NotImplementedError("k has to >= 2")
        s = self._signature
        if not (2*k in ZZ):
            raise ValueError("k has to be integral or half-integral")
        if (2*k+s)%2 != 0:
            return 0
        m = self._m
        n2 = self._n2
        if self._v2.has_key(0):
            v2 = self._v2[0]
        else:
            v2 = 1

        if self._g != None:
            if not self._aniso_formula:
                vals = self._g.values()
            #else:
                #print "using aniso_formula"
            M = self._g
        else:
            vals = self._M.values()
            M = self._M

        if (2*k+s)%4 == 0:
            d = Integer(1)/Integer(2)*(m+n2) # |dimension of the Weil representation on even functions|
            self._d = d
            if not self._aniso_formula: self._alpha4 = 1/Integer(2)*(vals[0]+v2) # the codimension of SkL in MkL
        else:
            d = Integer(1)/Integer(2)*(m-n2) # |dimension of the Weil representation on odd functions|
            self._d = d
            if not self._aniso_formula: self._alpha4 = 1/Integer(2)*(vals[0]-v2) # the codimension of SkL in MkL
            
        prec = ceil(max(log(M.order(),2),52)+1)+17
        #print prec
        RR = RealField(prec)
        CC = ComplexField(prec)
        
        if debug > 0: print "d, m = {0}, {1}".format(d,m)
        eps = exp( 2 * CC.pi() * CC(0,1) * (s + 2*k) / Integer(4) )
        eps = round(real(eps))        
        if self._alpha3 is None or self._last_eps != eps:
            self._last_eps = eps
            if self._aniso_formula:
                self._alpha4 = 1 if eps == 1 else 0
                self._alpha3 = -sum([BB(a)*mm for a,mm in self._v2.iteritems() if a != 0])
                #print self._alpha3
                self._alpha3 += Integer(d) - Integer(1) - self._g.beta_formula()
                #print self._alpha3, self._g.a5prime_formula()
                self._alpha3 = self._alpha3/RR(2)
            else:
                self._alpha3 = eps*sum([(1-a)*mm for a,mm in self._v2.iteritems() if a != 0])
                if debug>0: print "alpha3t = ", self._alpha3
                self._alpha3 += sum([(1-a)*mm for a,mm in vals.iteritems() if a != 0])
                #print self._alpha3
                self._alpha3 = self._alpha3 / Integer(2)
        alpha3 = self._alpha3
        alpha4 = self._alpha4
        if debug > 0: print alpha3, alpha4
        g1=M.char_invariant(1)
        g1=CC(g1[0]*g1[1])
        #print g1
        g2=M.char_invariant(2)
        g2=CC(g2[0]*g2[1])
        if debug > 0: print g2, g2.parent()
        g3=M.char_invariant(-3)
        g3=CC(g3[0]*g3[1])
        if debug > 0: print "eps = {0}".format(eps)
        if debug > 0: print "d/4 = {0}, m/4 = {1}, e^(2pi i (2k+s)/8) = {2}".format(RR(d) / RR(4), sqrt(RR(m)) / RR(4), CC(exp(2 * CC.pi() * CC(0,1) * (2 * k + s) / Integer(8))))
        if eps == 1:
            g2_2 = real(g2)
        else:
            g2_2 = imag(g2)*CC(0,1)
        alpha1 = RR(d) / RR(4) - sqrt(RR(m)) / RR(4)  * CC(exp(2 * CC.pi() * CC(0,1) * (2 * k + s) / Integer(8)) * g2_2)
        if debug > 0: print alpha1
        alpha2 = RR(d) / RR(3) + sqrt(RR(m)) / (3 * sqrt(RR(3))) * real(exp(CC(2 * CC.pi() * CC(0,1) * (4 * k + 3 * s - 10) / 24)) * (g1 + eps*g3))
        if debug > 0: print "alpha1 = {0}, alpha2 = {1}, alpha3 = {2}, g1 = {3}, g2 = {4}, g3 = {5}, d = {6}, k = {7}, s = {8}".format(alpha1, alpha2, alpha3, g1, g2, g3, d, k, s)
        dim = real(d + (d * k / Integer(12)) - alpha1 - alpha2 - alpha3)
        if debug > 0:
            print "dimension:", dim
        if abs(dim-round(dim)) > 1e-6:
            raise RuntimeError("Error ({0}) too large in dimension formula for {1} and k={2}".format(abs(dim-round(dim)), self._M if self._M is not None else self._g, k))
        dimr = dim
        dim = Integer(round(dim))
        if k >=2 and dim < 0:
            raise RuntimeError("Negative dimension (= {0}, {1})!".format(dim, dimr))
        return dim

    def dimension_cusp_forms(self, k, ignore=False, no_inv = False, test_positive = False, proof = False, debug=0):
        if debug>0:
            if self._g is not None:
                print "Computing dimension for {}".format(self._g)
            else:
                print "Computing dimension for {}".format(self._M)
        if (2*k+self._signature)%2 != 0:
            return 0
        if k == Integer(3)/2:
            dim = self.dimension(k, True, debug=debug) - self._alpha4
            if not test_positive or dim <= 0:
                if self._M == None:
                    self._M = self._g.finite_quadratic_module()
                corr = weight_one_half_dim(self._M, self._use_reduction, proof = proof)
                if debug > 0: print "weight one half: {0}".format(corr)
                dim += corr
        else:
            dim = self.dimension(k, ignore, debug=debug) - self._alpha4
        if k == Integer(2) and not no_inv:
            if test_positive and dim > 0:
                return dim
            if self._M == None:
                self._M = self._g.finite_quadratic_module()
            if self._M.level() == 1:
                return dim + 1
            dinv = cython_invariants_dim(self._M,self._use_reduction)
            dim = dim + dinv
        if dim < 0:
            raise RuntimeError("Negative dimension (= {0}, alpha4 = {1})!".format(dim, self._alpha4))
        return dim
        
def test_real_quadratic(minp=1,maxp=100,minwt=2,maxwt=1000):
    for p in prime_range(minp,maxp):
        if p%4==1:
            print "p = ", p
            gram=Matrix(ZZ,2,2,[2,1,1,(1-p)/2])
            M=VectorValuedModularForms(gram)
            if is_odd(minwt):
                minwt=minwt+1
            for kk in range(minwt,round(maxwt/2-minwt)):
                k = minwt+2*kk
                if M.dimension_cusp_forms(k)-dimension_cusp_forms(kronecker_character(p),k)/2 != 0:
                    print "ERROR: ", k, M.dimension_cusp_forms(k), dimension_cusp_forms(kronecker_character(p),k)/2
                    return False
    return true

#sys.path.append('/home/stroemberg/Programming/Sage/sage-add-ons3/nils')
#from jacobiforms.all import *

def test_jacobi(index=1,minwt=4,maxwt=100,eps=-1):
    m=index
    gram=Matrix(ZZ,1,1,[-eps*2*m])
    M=VectorValuedModularForms(gram)
    if is_odd(minwt):
        minwt=minwt+1
    for kk in range(0,round(maxwt/2-minwt)+2):
        k = minwt+2*kk+(1+eps)/2
        print "Testing k = ", k
        if eps==-1:
            dimJ=dimension_jac_forms(k,m,-1)
            dimV=M.dimension(k-1/2)
        else:
            dimJ=dimension_jac_cusp_forms(k,m,1)
            dimV=M.dimension_cusp_forms(k-1/2)
        if dimJ-dimV != 0:
            print "ERROR: ", k, dimJ, dimV
            return False
    return true
