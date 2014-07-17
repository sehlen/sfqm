#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file has been copied from PSAGE and modified to our needs.
#
#  PSAGE is free software: you can redistribute it and/or modify
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


import os, sys


if sys.maxint != 2**63 - 1:
    print "*"*70
    print "The PSAGE library only works on 64-bit computers.  Terminating build."
    print "*"*70
    sys.exit(1)


import build_system

SAGE_ROOT = os.environ['SAGE_ROOT']
SAGE_LOCAL = os.environ['SAGE_LOCAL']

INCLUDES = ['%s/%s/'%(SAGE_ROOT,x) for x in
#             ('devel/sage/sage/ext', 'devel/sage', 'devel/sage/sage/gsl',
            ('src/sage/ext', 'src/sage', 'src/sage/gsl', 'src'
              )] \
         + ['%s/%s/'%(SAGE_LOCAL,x) for x in
             ('include/csage', 'include', 'include/python',
              'include/python2.7')]

if '-ba' in sys.argv:
    print "Rebuilding all Cython extensions."
    sys.argv.remove('-ba')
    FORCE = True
else:
    FORCE = False

def Extension(*args, **kwds):
    if not kwds.has_key('include_dirs'):
        kwds['include_dirs'] = INCLUDES
    else:
        kwds['include_dirs'] += INCLUDES
    if not kwds.has_key('force'):
        kwds['force'] = FORCE

    # Disable warnings when running GCC step -- cython has already parsed the code and
    # generated any warnings; the GCC ones are noise.
    if not kwds.has_key('extra_compile_args'):
        kwds['extra_compile_args'] = ['-w']
    else:
        kwds['extra_compile_args'].append('-w')

    E = build_system.Extension(*args, **kwds)
    E.libraries = ['csage'] + E.libraries
    return E


numpy_include_dirs = [os.path.join(SAGE_LOCAL,
                                   'lib/python/site-packages/numpy/core/include')]

ext_modules = [
      Extension('sfqm.simple.find_simple_c',
              sources = ['sfqm/simple/find_simple_c.pyx'],
              libraries = ['m']
     ),
     Extension('psage.modules.invariants',
              sources = ['psage.modules.invariants.pyx'],
              libraries = ['m']
     )
]

ext_modules.extend(ext_modules)

build_system.cythonize(ext_modules)

build_system.setup(
    name = 'sfqm',
    version = "2014.x.x",
    description = "SFQM",
    author = 'Stephan Ehlen',
    author_email = 'stephan.j.ehlen@gmail.com',
    url = '',
    license = 'GPL v2+',
    packages = ['sfqm',
                'sfqm.simple',
                'sfqm.fqm',
		'psage',
		'psage.modules',
		'psage.modules.finite_quadratic_module',
		'psage.modform',
		'psage.modform.weilrep_tools'
                ],
    platforms = ['any'],
    download_url = 'NA',
    ext_modules = ext_modules
)
