"""These tests intentionally need no Sage installation."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    'search_ranges', Path(__file__).parents[1] / 'sfqm/simple/search_ranges.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SearchRangesTests(unittest.TestCase):
    def test_exhaustive_coverage(self):
        for lower in (1, 2, 7, 100):
            for count in range(21):
                for parts in range(1, 25):
                    chunks = module.inclusive_chunks(lower, count, parts)
                    actual = [n for a, b in chunks for n in range(a, b+1)]
                    self.assertEqual(actual, list(range(lower, lower+count)))
                    if chunks:
                        sizes = [b-a+1 for a, b in chunks]
                        self.assertLessEqual(max(sizes)-min(sizes), 1)

    def test_reject_invalid_arguments(self):
        for args in ((0, 10, 2), (1, -1, 2), (1, 2, 0)):
            with self.assertRaises(ValueError):
                module.inclusive_chunks(*args)
        with self.assertRaises(TypeError):
            module.inclusive_chunks(1, 2.5, 2)


if __name__ == '__main__':
    unittest.main()
