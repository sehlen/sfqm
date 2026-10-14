"""Build against the Python and Cython supplied by SageMath."""
import os
from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize
from sage.env import sage_include_directories

include_dirs = sage_include_directories()
# Keep legacy Cython integer division while the Python modules use Python 3.
# Do not enable fast-math for mathematical computations.
openmp = ["-fopenmp"] if os.name == "posix" and os.uname().sysname == "Linux" else []
extensions = [
    Extension("sfqm.simple.find_simple_c", ["sfqm/simple/find_simple_c.pyx"],
              include_dirs=include_dirs),
    Extension("psage.external.weil_invariants.weil_invariants",
              ["psage/psage/external/weil_invariants/weil_invariants.pyx"],
              include_dirs=include_dirs, extra_compile_args=openmp,
              extra_link_args=openmp),
]
setup(
    name="sfqm", version="0.3.0",
    description="Finite quadratic modules and simple lattices",
    author="Stephan Ehlen", url="https://github.com/sehlen/sfqm",
    packages=find_packages(include=["sfqm", "sfqm.*"]) + find_packages("psage"),
    package_dir={"psage": "psage/psage"},
    ext_modules=cythonize(extensions, compiler_directives={"language_level": 2}),
    python_requires=">=3.10",
)
