"""Search regressions, including failure propagation and pruning options."""
import unittest
from unittest.mock import patch
from sage.all import Integer
from sfqm.fqm.genus_symbol import GenusSymbol, anisotropic_symbols
from sfqm.simple.find_simple_c import _find_simple_anisotropic
from sfqm.simple import find_simple as search
from sfqm.simple.simple_modules_graph import SimpleModulesGraph


class SearchTests(unittest.TestCase):
    def test_worker_matches_unpruned_enumeration(self):
        weights = [Integer(3), Integer(4)]
        expected = {k: set() for k in weights}
        for N in range(1, 25):
            for g in anisotropic_symbols(N):
                if g.order() > 24:
                    continue
                for k in weights:
                    if (2*k+g.signature()) % 4 == 0 and g.dimension_cusp_forms(k) == 0:
                        expected[k].add(g)
        for dynamic in [False, True]:
            actual, _ = _find_simple_anisotropic(1, 24, 24, weights=weights, dynamic=dynamic)
            self.assertEqual({k:set(v) for k,v in actual.items()}, expected)

    def test_parallel_matches_serial_nondefault_lower(self):
        serial, _ = _find_simple_anisotropic(7, 24, 24, weights=[Integer(4)])
        parallel = search.find_simple_anisotropic_parallel(18, 4, [Integer(4)], lower=7)
        self.assertEqual({k:set(v) for k,v in parallel.items()},
                         {k:set(v) for k,v in serial.items()})

    def test_worker_failure_is_not_a_partial_result(self):
        with patch.object(search, 'find_simple_anisotropic_wrapper',
                          return_value=[(None, 'NO DATA')]):
            with self.assertRaises(RuntimeError):
                search.find_simple_anisotropic_parallel(4, 1, [Integer(4)])

    def test_graph_expands_without_preliminary_pruning(self):
        graph = SimpleModulesGraph(signature=0, weight=10, level_limit=2, rank_limit=2)
        # Exercise actual parallel graph serialization, including Sage's uniq generator.
        graph.compute_from_startpoints([GenusSymbol()], p=[2], cut_nonsimple_aniso=False)
        self.assertGreater(graph.num_verts(), 1)
        self.assertEqual(graph._simple, [GenusSymbol()])

    def test_helpers_require_precomputed_input(self):
        for function in [search.simple_gamma0_genus_symbols,
                         search.simple_gamma1_genus_symbols, search.simple_t1_genus_symbols]:
            with self.assertRaises(ValueError):
                function()
        g = GenusSymbol('2_1^+1')
        matched, rest = search.simple_t1_genus_symbols(precomputed=[g])
        self.assertTrue(matched)
        self.assertEqual(rest, [])

    def test_default_bound_brackets_a_root(self):
        self.assertGreater(search.theoretical_upper_bound_anisotropic(2), 1)
