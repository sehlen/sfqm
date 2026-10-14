# Finite Quadratic Modules and Simple Lattices

Algorithms and supplementary material for **Lattices with many Borcherds
products**, by Jan H. Bruinier, Stephan Ehlen and Eberhard Freitag.
The extended paper is [included here](bruinier_ehlen_freitag_extended.pdf).

## Installation

Use an isolated SageMath environment: this distribution includes the historical
`psage` dependencies and should not be installed over another PSAGE installation.

```sh
git clone https://github.com/sehlen/sfqm.git
cd sfqm
sage -pip install --no-build-isolation .
```

For development, use `sage -pip install --no-build-isolation -e .`.
The build uses Sage's installed Python, Cython and headers; installing with an
unrelated system Python is not supported. A C compiler is required. Linux builds
use OpenMP; other platforms build without that flag and have not been validated.

The modernization branch has been exercised on **Sage 10.7 / Python 3.12**.
CI additionally targets **Sage 10.10**; consult the PR checks before treating that
version as supported. The historical `master` implementation targeted Sage 7.2
and cannot be imported unchanged with Python 3.

The Cython files deliberately retain `language_level=2` to preserve their
integer-division semantics. The Python files have been ported to Python 3.

## Examples

```python
from sage.all import Integer
from sfqm.fqm.genus_symbol import GenusSymbol

s = GenusSymbol('3^+1')
assert s.dimension_cusp_forms(3) == 0
assert s.dimension_modular_forms(3) == 1

# The construction rules in Section 5 of the paper.
extensions = s.C(2)
assert all(t.order() == 4*s.order() for t in extensions)

# A small anisotropic search: inclusive levels 1 through 24, order <= 24.
from sfqm.simple.find_simple_c import _find_simple_anisotropic
simple, elapsed = _find_simple_anisotropic(1, 24, 24, weights=[Integer(4)])
```

The parallel wrapper's `num` is the number of levels to search, starting at
`lower`; its order cap is `lower + num - 1`. This historical API mixes level and
order bounds. For independent bounds, use the serial worker as above.

## Verification

```sh
OMP_NUM_THREADS=2 sage -python -m unittest discover -s tests -v
```

The suite includes direct enumeration of quadratic-form values and torsion,
the divisor sum (2.7), the class-number formula, Example 5.1, exact Weil-matrix
checks, low-weight corrections, and serial/parallel search comparisons.
`OMP_NUM_THREADS=2` avoids excessive native threading inside forked searches.

The old `test_global.py` files are historical computational notebooks in Python
form, not an automated regression suite. Some depend on Magma or missing session
globals; they are outside the supported test command above.

See the [review report](docs/review-2026-10.md) for confirmed defects, changes,
validation scope and remaining work. Passing these small-case checks is **not a
reproduction of the complete classification tables**.
