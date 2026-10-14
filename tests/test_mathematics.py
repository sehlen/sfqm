"""Independent enumeration checks; run with sage -python -m unittest discover -s tests."""
from collections import Counter
from itertools import product
import unittest
from sage.all import ZZ, QQ, Integer, sqrt, divisors, prod, CyclotomicField
from sfqm.fqm.genus_symbol import GenusSymbol, anisotropic_symbols, C
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from psage.modform.weilrep_tools.dimension import VectorValuedModularForms


def direct_values(A, torsion=None):
    return dict(Counter(x.norm() for x in A if torsion is None or torsion*x == A(0)))


class MathematicalTests(unittest.TestCase):
    def test_representation_numbers(self):
        symbols = ['1^+1', '2^+2', '2^-2', '2_1^+1', '2_7^+1',
                   '2_0^+2', '2_2^+2', '2_3^+3', '4_3^-1',
                   '4_0^+2', '4^-2', '8_1^+1', '3^+1', '3^-2',
                   '9^-1', '3^+1.9^-1', '2_1^+1.4_3^-1',
                   '2^+2.3^-1', '2_7^+1.3^+1.5^-1']
        for symbol in symbols:
            with self.subTest(symbol=symbol):
                g = GenusSymbol(symbol)
                A = g.finite_quadratic_module()
                self.assertEqual(g.values(), direct_values(A))
                self.assertEqual(g.two_torsion_values(), direct_values(A, 2))
                self.assertEqual(sum(g.values().values()), g.order())

    def test_sigma_against_definition(self):
        for symbol in ['3^+1.9^-1', '2_1^+1.4_1^+1', '3^-2.5^+1']:
            g = GenusSymbol(symbol)
            A = g.finite_quadratic_module()
            for s in [-1, 0, 1]:
                expected = sum(Integer(d)**s * sqrt(sum(d*x == A(0) for x in A))
                               for d in divisors(g.level()))
                self.assertEqual(g.sigma(s), expected)

    def test_anisotropic_formula(self):
        # Independent sum of the sawtooth function over actual module elements.
        for N in range(1, 61):
            for g in anisotropic_symbols(N):
                A = g.finite_quadratic_module()
                values = direct_values(A)
                self.assertEqual(values[0], 1, str(g))
                beta = sum((a - QQ(1)/2)*n for a, n in values.items() if a)
                self.assertLess(abs(g.beta_formula()-beta), 1e-10, str(g))
                for k in [Integer(2), Integer(3), Integer(4), QQ(5)/2]:
                    if (k > 2 or N <= 12) and (2*k + g.signature()) % 4 == 0:
                        self.assertEqual(g.dimension_cusp_forms(k, aniso_formula=True),
                                         g.dimension_cusp_forms(k), (str(g), k))

    def test_scalar_dimensions(self):
        from sage.modular.dims import dimension_modular_forms, dimension_cusp_forms
        for k in [2, 4, 6, 8, 10, 12, 14, 24]:
            g = GenusSymbol()
            self.assertEqual(g.dimension_modular_forms(k), dimension_modular_forms(1,k))
            V = VectorValuedModularForms(FiniteQuadraticModule())
            self.assertEqual(V.dimension(k), dimension_modular_forms(1,k))
            self.assertEqual(V.dimension_cusp_forms(k), dimension_cusp_forms(1,k))

    def test_modular_dimension_api(self):
        self.assertEqual(GenusSymbol('3^+1').dimension_cusp_forms(3), 0)
        g = GenusSymbol('3^+1')
        self.assertEqual(g.dimension_modular_forms(3), 1)
        self.assertEqual(g.dimension_modular_forms(QQ(5)/2), 0)

    def test_constructor_does_not_alias(self):
        data = {3: [[1, 1, 1], [1, 1, -1]]}
        g = GenusSymbol(data)
        self.assertEqual(data, {3: [[1, 1, 1], [1, 1, -1]]})
        data[3][0][1] = 42
        self.assertEqual(g.order(), 9)
        g.jordan_component(3)._symbol_dict[3][0][1] = 40
        self.assertEqual(g.order(), 9)

    def test_addition_subtraction_and_order(self):
        a, b = GenusSymbol('3^+1'), GenusSymbol('9^-1')
        self.assertEqual((a+b)-b, a)
        self.assertEqual((a+b)-(a+b), GenusSymbol())
        self.assertFalse(a < a)
        self.assertTrue(a <= a)
        self.assertEqual(sorted([b,a,a]), [a,a,b])
        self.assertEqual(hash(a), hash(GenusSymbol(str(a))))

    def test_C_paper_example_5_1(self):
        A = GenusSymbol('3^+1.9^-1.7^-3')
        expected3 = ['3^-3.9^-1.7^-3', '3^+1.9^-2.7^-3', '3^+1.9^+2.7^-3',
                     '9^-1.27^+1.7^-3', '3^+1.81^-1.7^-3']
        expected7 = ['3^+1.9^-1.7^+5', '3^+1.9^-1.7^-3.49^+1',
                     '3^+1.9^-1.7^-3.49^-1', '3^+1.9^-1.7^-2.343^+1',
                     '3^+1.9^-1.7^+2.343^-1']
        self.assertEqual(set(C(A,3)), set(map(GenusSymbol,expected3)))
        self.assertEqual(set(C(A,7)), set(map(GenusSymbol,expected7)))
        self.assertEqual(C(A,1), [A])
        self.assertEqual(set(C(A,9)), {t for B in C(A,3) for t in C(B,3)})
        with self.assertRaises(ValueError): C(A,0)

    def test_low_weight_corrections(self):
        from psage.modform.weilrep_tools.weight_one_half import weight_one_half_dim
        expected = {'2_1^+1': 1, '2_7^+1': 0, '4_1^+1': 1, '4_7^+1': 0}
        for symbol, dimension in expected.items():
            A = GenusSymbol(symbol).finite_quadratic_module()
            for local in [False, True]:
                self.assertEqual(weight_one_half_dim(A, False, local=local), dimension)
        self.assertEqual(GenusSymbol('2_1^+1').dimension_cusp_forms(QQ(3)/2), 0)

    def test_large_integer_multiplicities(self):
        # Dict inputs use native Python ints; division must remain exact.
        g = GenusSymbol({2: [[1, 60, 1, 0, 0]]})
        expected = {0: (2**60 + 2**30)//2, QQ(1)/2: (2**60 - 2**30)//2}
        self.assertEqual(g.two_torsion_values(), expected)

    def test_C_rules_have_isotropic_quotients(self):
        for symbol in ['1^+1', '2^-2', '2_0^+2', '2_2^+2', '3^+1']:
            g = GenusSymbol(symbol)
            A = g.finite_quadratic_module()
            for p in [2, 3]:
                for extension in C(g, p):
                    B = extension.finite_quadratic_module()
                    self.assertEqual(B.order(), A.order()*p*p)
                    self.assertEqual(B.signature(), A.signature())
                    self.assertTrue(any(U.is_isotropic() and U.quotient().is_isomorphic(A)
                                        for U in B.subgroups(p)), str(extension))

    def test_invariants_against_exact_S_matrix(self):
        from psage.external.weil_invariants.weil_invariants import cython_invariants_dim
        from sage.all import matrix
        for symbol in ['1^+1','2^+2','2^-2','3^+2','3^-2','5^+2','5^-2']:
            A = FiniteQuadraticModule() if symbol == '1^+1' else FiniteQuadraticModule(symbol)
            elements = list(A)
            isotropic = [x for x in elements if not x.norm()]
            K = CyclotomicField(8*A.level())
            z = K.zeta(A.level())
            phase = K.zeta(8)**(-A.signature())
            order_sqrt = K(sqrt(A.order()))
            SminusI = matrix(K, [[phase/order_sqrt * z**Integer(-A.level()*x.dot(y))
                                 - (x == y) for y in isotropic] for x in elements])
            expected = SminusI.right_kernel().dimension()
            for reduction in [False, True]:
                self.assertEqual(cython_invariants_dim(A, reduction), expected,
                                 (symbol, reduction))
            self.assertEqual(cython_invariants_dim(A, True, proof=True), expected)


if __name__ == '__main__':
    unittest.main()
