#################################################################################
#
# (c) Copyright 2010 William Stein, 2018 Stephan Ehlen
#
#  This file has been copied from PSAGE and modified to our needs.
#
#  PSAGE and all components in this code repository are free software: 
#  you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

import os, sys, tarfile
from sage.env import sage_include_directories,SAGE_INC,SAGE_LIB,SAGE_LOCAL
import subprocess 

if sys.maxint != 2**63 - 1:
    print "*"*70
    print "The PSAGE library only works on 64-bit computers.  Terminating build."
    print "*"*70
    sys.exit(1)

if '-ba' in sys.argv:
    print "Rebuilding all Cython extensions."
    sys.argv.remove('-ba')
    FORCE = True
else:
    FORCE = False

if '-np' in sys.argv:
    print 'Not including psage in build'
    sys.argv.remove('-np')
    INSTALL_PSAGE = False
else:
    INSTALL_PSAGE = True
    print 'Also installing psage dependencies...'
    #pt = tarfile.open('psage.tar.bz2', mode='r:bz2')
    #pt.extractall()

#from module_list import ext_modules,aliases

include_dirs = sage_include_directories(use_sources=True)
include_dirs = include_dirs + [os.path.join(SAGE_LIB,"cysignals")]
include_dirs = include_dirs + [os.path.join(SAGE_LIB,"sage/ext/")]

extra_compile_args = [ "-fno-strict-aliasing" ]
extra_link_args = [ ]

DEVEL = False
if DEVEL:
    extra_compile_args.append('-ggdb')

# Work around GCC-4.8.0 bug which miscompiles some sig_on() statements,
# as witnessed by a doctest in sage/libs/gap/element.pyx if the
# compiler flag -Og is used. See also
# * http://trac.sagemath.org/sage_trac/ticket/14460
# * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56982
if subprocess.call("""$CC --version | grep -i 'gcc.* 4[.]8' >/dev/null """, shell=True) == 0:
    extra_compile_args.append('-fno-tree-dominator-opts')

    
lib_headers = { "gmp":     [ os.path.join(SAGE_INC, 'gmp.h') ],   # cf. #8664, #9896
                "gmpxx":   [ os.path.join(SAGE_INC, 'gmpxx.h') ],
                "ntl":     [ os.path.join(SAGE_INC, 'NTL', 'config.h') ]
            }

print "include_dirs=",include_dirs

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

ext_modules = [
               Extension('sfqm.simple.find_simple_c',
                         sources = ['sfqm/simple/find_simple_c.pyx'],
                         libraries = ['m']
                         )
               ]

if INSTALL_PSAGE:
    ext_modules.extend(
                       [
                        Extension('psage.external.weil_invariants.weil_ivariants',
                                  sources = ['psage/external/weil_invariants/weil_invariants.pyx'],
                                  libraries = ['m'],
                                  extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
                                  extra_link_args=['-fopenmp']
                                  )
                        ]
                       )


for m in ext_modules:
    m.depends = m.depends + [__file__]

    # Add dependencies for the libraries
    for lib in lib_headers:
        if lib in m.libraries:
            m.depends += lib_headers[lib]

    m.extra_compile_args = m.extra_compile_args + extra_compile_args
    m.extra_link_args = m.extra_link_args + extra_link_args
    m.library_dirs = m.library_dirs + [os.path.join(SAGE_LOCAL, "lib")]
    m.include_dirs = m.include_dirs + include_dirs

print "include_dirs=",include_dirs

def run_cythonize():
    from Cython.Build import cythonize
    import Cython.Compiler.Options
    import Cython.Compiler.Main
    debug = False
    gdb_debug = True
    if os.environ.get('SAGE_DEBUG', None) == 'yes':
        print('Enabling Cython debugging support')
        debug = True
        Cython.Compiler.Main.default_options['gdb_debug'] = True
        Cython.Compiler.Main.default_options['output_dir'] = 'build'
        gdb_debug=True

    profile = False    
    if os.environ.get('SAGE_PROFILE', None) == 'yes':
        print('Enabling Cython profiling support')
        profile = True
   
    # Sage uses these directives (mostly for historical reasons).
    Cython.Compiler.Options.embed_pos_in_docstring = True
    Cython.Compiler.Options.get_directive_defaults()['autotestdict'] = False
    Cython.Compiler.Options.get_directive_defaults()['cdivision'] = True
    Cython.Compiler.Options.get_directive_defaults()['fast_getattr'] = True
    # The globals() builtin in Cython was fixed to return to the current scope,
    # but Sage relies on the broken behavior of returning to the nearest
    # enclosing Python scope (e.g. to perform variable injection).
    Cython.Compiler.Options.old_style_globals = True
    Cython.Compiler.Main.default_options['cache'] = False

    global ext_modules
    ext_modules = cythonize(
        ext_modules,
        gdb_debug=gdb_debug,
        nthreads=int(os.environ.get('SAGE_NUM_THREADS', 0)),
        #    build_dir=SAGE_CYTHONIZED,
        force=FORCE,
        include_path = include_dirs,
        aliases=aliases,
        compiler_directives={
            'embedsignature': True,
            'profile': profile,
        })
print("Updating Cython code....")
import time
t = time.time()
run_cythonize()
print("Finished Cythonizing, time: %.2f seconds." % (time.time() - t))
import distutils
#for m in ext_modules:
#    print m,isinstance(m,distutils.extension.Extension)
#from distutils.core import setup
from setuptools import setup

if '-np' in sys.argv:
    print 'Not including psage in build'
    sys.argv.remove('-np')
    INSTALL_PSAGE = False
else:
    INSTALL_PSAGE = True
    print 'Also installing psage dependencies...'
    pt = tarfile.open('psage.tar.bz2', mode='r:bz2')
    pt.extractall()



packages = [
             'sfqm',
             'sfqm.simple',
             'sfqm.fqm'
           ]

if INSTALL_PSAGE:
    packages.extend(
        [
          'psage.modules',
          'psage.modform.weilrep_tools'
        ]
    )

code = setup(
    name = 'sfqm',
    version = "0.2",
    description = "SFQM",
    author = 'Stephan Ehlen',
    author_email = 'stephan.j.ehlen@gmail.com',
    url = '',
    license = 'GPL v2+',
    packages = packages,
    platforms = ['any'],
    download_url = 'http://www.github.com/sehlen/sfqm',
    ext_modules = ext_modules,
)
