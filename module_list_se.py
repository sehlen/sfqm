r"""
List of extension modules
"""
#from distutils.extension import Extension
import os
import pkgconfig

# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = list(cblas_pc['libraries'])
cblas_library_dirs = list(cblas_pc['library_dirs'])
cblas_include_dirs = list(cblas_pc['include_dirs'])
# TODO: Remove Cygwin hack by installing a suitable cblas.pc
if os.path.exists('/usr/lib/libblas.dll.a'):
    cblas_libs = ['gslcblas']
# GNU Scientific Library
# Note we replace the built-in gslcblas with the above cblas
gsl_pc = pkgconfig.parse('gsl')
gsl_libs = list(set(gsl_pc['libraries']).difference(set(['gslcblas'])).union(set(cblas_libs)))
gsl_library_dirs = list(gsl_pc['library_dirs'])
gsl_include_dirs = list(gsl_pc['include_dirs'])

aliases = dict(
    GSL_LIBRARIES=gsl_libs,
    GSL_LIBDIR=gsl_library_dirs,
    GSL_INCDIR=gsl_include_dirs,
)

numpy_include_dirs = []
import setuptools
class Extension(setuptools.extension.Extension):
    def __init__(self, name, sources, include_dirs=[],
                  language="c", force=False, **kwds):
        #print "kwds=",kwds
        #print "module=",module
        setuptools.Extension.__init__(self, name, sources, language=language,
                                       include_dirs=include_dirs, **kwds)         
#     def __init__(self, module, sources, include_dirs,
#                  language="c", force=False, **kwds):
#         self.cython_cmds = []
#         for i in range(len(sources)):
#             f = sources[i]
#             if f.endswith('.pyx'):
#                 sources[i], cmds = cython(f) #, language, include_dirs, force)
#                 for c in cmds:
#                     self.cython_cmds.append(c)
#         setuptools.Extension.__init__(self, module, sources, language=language,
#                                       include_dirs=include_dirs, **kwds)


ext_modules = [
        Extension("psage.groups.dirichlet_conrey",
                                ['psage/groups/dirichlet_conrey.pyx'],
                                extra_compile_args = ['-w','-O2'])
]

## Stephan Ehlen's additional modules

sehlen_extensions = [
      Extension('psage.external.weil_invariants.weil_invariants',
              sources = ['psage/external/weil_invariants/weil_invariants.pyx'],
              libraries = ['m']
#              include_dirs = numpy_include_dirs),
     )
]

ext_modules.extend(sehlen_extensions)
