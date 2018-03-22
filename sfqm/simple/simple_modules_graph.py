"""
Tools to compute with finite quadratic modules and simple lattices.

AUTHOR: (c) Stephan Ehlen, 2014

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
#from Bsets import dict_to_genus_symbol_string, genus_symbol_string_to_dict, Bbf
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from sage.all import ZZ, Zmod, sys, parallel, is_prime, colors, cached_function, Integer, Partitions, Set, QQ, RR, is_prime_power, next_prime, prime_range, is_squarefree, uniq, MatrixSpace, kronecker, deepcopy, CC, exp, walltime, RealField, floor, pari, pi, ComplexField, sqrt, text, arrow, is_even, squarefree_part, polygon2d, line2d
from sage.parallel.decorate import *
from sage.misc.cachefunc import *
import itertools
from sage.all import DiGraph
from sage.misc.decorators import options

### to be removed ##
from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, Genus_Symbol_p_adic_ring, is_GlobalGenus
####################

import logging
import datetime

from sfqm.fqm.genus_symbol import GenusSymbol, anisotropic_symbols, prime_anisotropic_symbols

NCPUS0 = 4
NCPUS1 = 10

@parallel(ncpus=NCPUS1)
def check_simple(s, k, reduction = False, bound = 0, check_injectivity_criterion=False, p=0, q=0):
    simple = s.is_simple(k, reduction=reduction, bound=bound)
    if not simple and check_injectivity_criterion:
        return (not s.is_additive_lift_injective(p,q))
    else:
        return simple

def prime_pol(s, p, k):
    A = RR(s.order())
    A2 = RR(s.torsion(2))
    A3 = RR(s.torsion(3))
    p = RR(p)
    k = RR(k)
    
    return A * (p**2 - 1) * (k - 1) / 24 - A2 / 2 \
             - (2 + 2 * A3) / (3 * RR(3).sqrt()) - 3 * (p-1) * A / 2

def prime_pol_simple(p, k):
    p = RR(p)
    k = RR(k)
    return (p**2-1)*(k - 1)/24 - 0.5*p


class ColorFormatter(logging.Formatter):

    """
    This Formatter adds some colors.
    (Copied from the LMFDB.)
    """
    fmtString = '%(levelname)s:%(name)s: %(message)s'

    def __init__(self, *args, **kwargs):
        self._hl = kwargs.pop('hl', None)
        self._fmt_orig = kwargs.pop('fmt', None)
        logging.Formatter.__init__(self, self._fmt_orig, *args, **kwargs)

    def format(self, record):
        """modify the _mft string, call superclasses format method"""
        # reset fmt string
        self._fmt = self._fmt_orig or ColorFormatter.fmtString
        #fn = os.path.split(record.pathname)[-1]
        #record.pathname = "%s:%s" % (fn, record.lineno)

        # some colors for severity level
        if record.levelno >= logging.ERROR:
            self._fmt = '\033[91m' + self._fmt
        elif record.levelno >= logging.WARNING:
            self._fmt = '\033[93m' + self._fmt
        elif record.levelno <= logging.DEBUG:
            self._fmt = '\033[93m' + self._fmt
        elif record.levelno <= logging.INFO:
            self._fmt = '\033[92m' + self._fmt

        # bold, if module name matches
        if record.name == self._hl:
            self._fmt = "\033[1m" + self._fmt

        # reset, to unaffect the next line
        self._fmt += '\033[0m'

        return logging.Formatter.format(self, record)

formatter = ColorFormatter()

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)

fh = logging.FileHandler('simple_modules_graph-{0}.log'.format(os.getpid()))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)


def get_logger(name):
    logger = logging.getLogger('simple-{0}'.format(name))
    if len(logger.handlers) == 0:
        logger.setLevel(logging.DEBUG)
        logger.addHandler(ch)
        logger.addHandler(fh)
    return logger

logger = get_logger(os.getpid())


class SimpleModulesGraph(DiGraph):
    """
      SimpleModulesGraph is a class that can compute and plot
      a graph that contains all $k$-simple modules of signature $s$
      with bounded minimal number of generators.
    """

    def __init__(self, signature=0, weight=2, level_limit=34, rank_limit=4, primes=None, simple_color=None, nonsimple_color=None,
                 reduction=True, bound=0, check_injectivity_criterion=False, r=0, s=0):
        """
            Initialize a SimpleModulesGraph containing finite quadratic modules of signature ``signature``.
            They are checked for being ``weight``-simple if their minimal number of generators
            is at most `rank_limit`.

            INPUT:
            - ``signature``: the signature of the modules
            - ``weight``: check for cusp forms of weight ``weight``
            - ``level_limit``: only check for anisotropic modules with level smaller than ``level_limit``
            - ``rank_limit``: an upper bound for the minimal number of generators
            - ``bound``: upper bound for the dimension (for considered being simple), default=0

            OUTPUT:
            A SimpleModulesGraph object. No computations are done after initialization.
            Start the computation using the method ``compute()``.
        """
        ###########################
        # basic parameters
        ###########################
        self._level_limit = Integer(level_limit)
        self._rank_limit = Integer(rank_limit)
        self._signature = Integer(signature) % 8
        self._weight = QQ(weight)
        self._reduction = reduction
        self._bound = bound
        self._check_injectivity_criterion = check_injectivity_criterion
        self._r = r
        self._s = s
        logger.debug("Check: {0}".format(self._check_injectivity_criterion))
        #########################################################
        # Initialize the primes that need to be checked
        # According to Theorem 4.21 in [BEF],
        # in the worst case, we need to check primes p for which
        # prime_pol_simple(p, weight) <= bound.
        #########################################################
        if primes is None:
            p = 2
            while True:
                if prime_pol_simple(p, weight) > self._bound:
                    primes = list(prime_range(p))
                    break
                else:
                    p = next_prime(p)
        self._primes = primes
        self._simple_color = colors.darkred.rgb() if simple_color is None else simple_color
        self._nonsimple_color = colors.darkgreen.rgb() if nonsimple_color is None else nonsimple_color
        self._vertex_colors = dict()
        self._vertex_colors[self._simple_color] = list()
        self._vertex_colors[self._nonsimple_color] = list()
        self._heights = dict() # a height function for plotting
        self._simple = list() # will contain the list of k-simple modules
        super(SimpleModulesGraph, self).__init__()

    def compute_from_startpoints(self, points, p = None, cut_nonsimple_aniso = True, fast = 1, **kwds):
        self._reduction = kwds.get('reduction', self._reduction)
        if NCPUS0 == 1:
            for a in points:
                self._compute_simple_modules_graph_from_startpoint(a, p, cut_nonsimple_aniso, fast)
        else:
            computations = self._compute_simple_modules_graph_from_startpoint_parallel([(a, p, cut_nonsimple_aniso, fast) for a in points])
            for r in computations:
                if not isinstance(r[1], SimpleModulesGraph):
                    logger.error("Got something else... {0}".format(r))
                    # print r
                else:
                    G = r[1]
                    logger.debug("simple from G: {0}".format(len(G._simple)))
                    logger.debug(G._simple)
                    self.add_vertices(G.vertices())
                    self.add_edges(G.edges())
                    self._vertex_colors[
                        self._nonsimple_color] += G._vertex_colors[G._nonsimple_color]
                    self._vertex_colors[
                        self._simple_color] += G._vertex_colors[G._simple_color]
                    self._simple = uniq(self._simple + G._simple)
                    for h in G._heights.keys():
                        if not self._heights.has_key(h):
                            self._heights[h] = G._heights[h]
                        else:
                            self._heights[h] = uniq(self._heights[h] + G._heights[h])
                            logger.info("Found {0} {1}-simple module{2} so far.".format(
                                len(self._simple), self._weight, "s" if len(self._simple) != 1 else ""))
        logger.info("Found in total {0} {1}-simple module{3} with p-rank <= {2}".format(
            len(self._simple), self._weight, self._rank_limit, "s" if len(self._simple) != 1 else ""))
        return 0

    def compute(self, p=None, cut_nonsimple_aniso=True, fast=1):
        args = list()
        for N in range(1, self._level_limit):
            v2 = Integer(N).valuation(2)
            N2 = 2 ** v2
            if v2 in [0, 1, 2, 3] and is_squarefree(Integer(N) / N2):
                s = anisotropic_symbols(N, self._signature)
                if len(s) == 0:
                    continue
                args = args + s
        logger.debug('args = {0}'.format(args))
        logger.info('starting with {} anisotropic modules'.format(len(args)))
        self.compute_from_startpoints(args, p, cut_nonsimple_aniso, fast)
        return self._simple

    @parallel(ncpus=NCPUS0)
    def _compute_simple_modules_graph_from_startpoint_parallel(self, s, p=None, cut_nonsimple_aniso=True, fast=1):
        G = SimpleModulesGraph(
            self._signature, self._weight, self._level_limit, self._rank_limit, self._primes, reduction=self._reduction, bound=self._bound, check_injectivity_criterion=self._check_injectivity_criterion, r=self._r, s=self._s)
        G._compute_simple_modules_graph_from_startpoint(
            s, p, cut_nonsimple_aniso, fast)
        return G

    def _compute_simple_modules_graph_from_startpoint(self, s, p=None, cut_nonsimple_aniso=True, fast=1):
        # for forking this is necessary
        logger = get_logger(s)
        # print logger
        k = self._weight
        ###########################################################
        # Determine which primes need to be checked
        # According to the proof of Proposition XX in [BEF], we
        # only need to check primes not dividing the 6*level(s),
        # for which prime_pol(s,p,k) <= 0.
        # For those primes, we check if there is any
        # k-simple fqm in s.C(p) and if not, we do not have to
        # consider p anymore.
        ###########################################################
        if p == None:
            p = 2
            N = Integer(6) * s.level()
            slp = N.prime_factors()
            for q in prime_range(next_prime(N) + 1):
                if not q in slp:
                    logger.info(
                        "Smallest prime not dividing 6*level({0}) = {1} is p = {2}".format(s, Integer(6) * s.level(), q))
                    p = q
                    break
            while prime_pol(s, p, k) <= self._bound or p in slp:
                p = next_prime(p)
            p = uniq(prime_range(p) + slp)
        logger.info("Starting with s = {0} and primes = {1}".format(s, p))
        
        if isinstance(p, list):
            primes = p
        else:
            primes = [p]

        simple = s.is_simple(k, reduction = self._reduction, bound = self._bound)
        if self._check_injectivity_criterion and not simple:
            simple = not s.is_additive_lift_injective(self._r, self._s)
            logger.info("Checked injectivity for non-simple module {0} in signature ({2}, {3}): {1}".format(s,simple, self._r, self._s))
        
        if not simple:
            logger.info("{0} is not simple.".format(s))

        # print simple
        s = FQM_vertex(s)
        # print s
        self.add_vertex(s)
        # print "added vertex ", s
        
        ############################################################
        # Starting from the list of primes we generated,
        # we now compute which primes we actually need to consider.
        # That is, the primes such that there is a fqm in s.C(p)
        # which is k-simple.
        ############################################################
        np = list()
        if cut_nonsimple_aniso and simple:
            for i in range(len(primes)):
                p = primes[i]
                fs = False
                for t in s.genus_symbol().C(p, False):
                    if t.is_simple(k, bound = self._bound) \
                       or (self._check_injectivity_criterion and not t.is_additive_lift_injective(self._r, self._s)):
                        fs = True
                        logger.debug("Setting fs = True")
                        break
                if fs:
                    np.append(p)
        primes = np
        # print "here", primes
        logger.info("primes for graph for {0}: {1}".format(s, primes))
        # if len(primes) == 0:
        #    return
        heights = self._heights
        h = 0
        if not heights.has_key(h):
            heights[h] = [s]
        else:
            if heights[h].count(s) == 0:
                heights[h].append(s)
        vertex_colors = self._vertex_colors
        nonsimple_color = self._nonsimple_color
        simple_color = self._simple_color
        # Bs contains the modules of the current level (= height = h)
        Bs = [s]
        # set the correct color for the vertex s
        if simple:
            if vertex_colors[simple_color].count(s) == 0:
                vertex_colors[simple_color].append(s)
        elif vertex_colors[nonsimple_color].count(s) == 0:
            vertex_colors[nonsimple_color].append(s)
        
        ###################################################
        # MAIN LOOP
        # we loop until we haven't found any simple fqm's
        ###################################################
        while simple:
            h = h + 1
            if not heights.has_key(h):
                heights[h] = list()
            # the list Bss will contain the k-simple modules of the next height level
            # recall that Bs contains the modules of the current height level
            Bss = list()
            simple = False
            # checklist = modules that we need to check for with .is_simple(k)
            # we assemble this list because afterwards
            # we check them in parallel
            checklist = []
            for s1 in Bs:
                Bs2 = list()
                for p in primes:
                    # check if we really need to check p for s1
                    # otherwise none of the fqm's in s1.C(p) are simple
                    # and we will not consider them.
                    if prime_pol(s1.genus_symbol(), p, k) <= self._bound:
                        Bs2 = Bs2 + s1.genus_symbol().C(p, False)
                    else:
                        logger.info(
                            "Skipping p = {0} for s1 = {1}".format(p, s1))
                        # print "Skipping p = {0} for s1 = {1}".format(p, s1)
                # print "Bs2 = ", Bs2
                # now we actually check the symbols in Bs2
                for s2 in Bs2:
                    if s2.max_rank() > self._rank_limit:
                        # we skip s2 if its minimal number of generators
                        # is > than the given rank_limit.
                        continue
                    s2 = FQM_vertex(s2)
                    skip = False
                    for v in self._heights[h]:
                        # we skip symbols that correspond to isomorphic modules
                        if v.genus_symbol().defines_isomorphic_module(s2.genus_symbol()):
                            skip = True
                            logger.debug(
                                "skipping {0} b/c isomorphic to {1}".format(s2.genus_symbol(), v.genus_symbol()))
                            s2 = v
                            break
                    if skip:
                        continue
                    if not skip:
                        self.add_vertex(s2)
                        heights[h].append(s2)
                    self.update_edges(s2, h, fast=fast)
                    # before using the actual dimension formula
                    # we check if there is already a non-k-simple neighbor
                    # (an incoming edge from a non-k-simple module)
                    # which would imply that s2 is not k-simple.
                    has_nonsimple_neighbor = False
                    for e in self.incoming_edges(s2):
                        if vertex_colors[nonsimple_color].count(e[0]) > 0:
                            has_nonsimple_neighbor = True
                            logger.debug(
                                "Has nonsimple neighbor: {0}".format(s2.genus_symbol()))
                            break
                    if has_nonsimple_neighbor:
                        #not simple
                        if not vertex_colors[nonsimple_color].count(s2) > 0:
                            vertex_colors[nonsimple_color].append(s2)
                    else:
                        logger.debug("Check injectivity: {0} ({1},{2})".format(self._check_injectivity_criterion, self._r, self._s))
                        checklist.append((s2, k, self._reduction, self._bound, self._check_injectivity_criterion, self._r, self._s))
            logger.debug("checklist = {0}".format(checklist))
            # check the modules in checklist
            # for being k-simple
            # this is done in parallel
            # when a process returns
            # we add the vertex and give it its appropriate color
            if NCPUS1 == 1:
                checks = [([[s[0]]],check_simple(*s)) for s in checklist]
            else:
                checks = list(check_simple(checklist))
            logger.info("checks = {0}".format(checks))
            for check in checks:
                s2 = check[0][0][0]
                if check[1]:
                    simple = True
                    logger.info(
                        "Found simple module: {0}".format(s2.genus_symbol()))
                    Bss.append(s2)
                    if not vertex_colors[simple_color].count(s2) > 0:
                        vertex_colors[simple_color].append(s2)
                else:
                    if not vertex_colors[nonsimple_color].count(s2) > 0:
                        vertex_colors[nonsimple_color].append(s2)
            Bs = Bss
        simple = [v.genus_symbol() for v in vertex_colors[simple_color]]
        self._simple = uniq(simple)

    @options()
    def plot(self, textvertices=False, only_simple=False, fontsize=18, sort=False, fact = 0.5, thickness=4, edges_thickness=4, arrowshorten=8,
              linestyle_simple='solid', linestyle_nonsimple='dashed', arrowsize=2, arrows=False, **options):
        r"""
          Plots a SimpleModulesGraph.
        """
        vertex_colors = self._vertex_colors
        simple_color = self._simple_color
        nonsimple_color = self._nonsimple_color
        heights = self._heights
        if only_simple:
            for v in vertex_colors[nonsimple_color]:
                if self.has_vertex(v):
                    self.delete_vertex(v)
                # while vertex_colors[nonsimple_color].count(v)>0:
                #    vertex_colors[nonsimple_color].remove(v)
                for j in range(len(heights)):
                    while heights[j].count(v) > 0:
                        heights[j].remove(v)
            vertex_colors[nonsimple_color] = dict()

        pos = dict()
        labels = list()
        vertices = list()
        edges = list()
        min_d = 0.5
        widths = [float(sum([len(str(v)) / float(6) + min_d for v in heights[i]]))
                  for i in range(len(heights))]
        print widths
        max_w = max(widths)
        if len(widths) > 1:
            min_w = min([w for w in widths if w > 3])
        else:
            min_w = 0
        print min_w, max_w, widths
        max_vert = max([len(_) for _ in heights.values()])
        for i in range(len(heights)):
            if sort:
                heights[i] = sorted(heights[i])
            #print heights[i]
            real_w = widths[i]
            prev_w = min_w if i == 0 else widths[i - 1]
            next_w = min_w if i == len(heights) - 1 else widths[i + 1]
            cur_w = max(float(next_w) * 0.9, float(prev_w) * 0.9, real_w)
            d = max(2 * float(cur_w - real_w) / len(heights[i]), min_d)
            print real_w, cur_w
            print "d = ", d
            p = [-(cur_w), float(max_vert) * fact * i]
            w = 0
            for j in range(len(heights[i])):
                v = heights[i][j]
                p = [p[0] + w + 0.2, p[1]]
                w = float(len(str(v))) / float(6)
                p[0] = p[0] + w + d
                pos[heights[i][j]] = p
                c = simple_color if vertex_colors[
                    simple_color].count(v) > 0 else nonsimple_color
                if textvertices:
                    ct = colors.black.rgb()
                else:
                    ct = colors.white.rgb()
                labels.append(
                    text("$\mathbf{" + str(v)[1:len(str(v))-1] + "}$", (p[0] + 0.2, p[1]), rgbcolor=ct, zorder=8, fontsize=fontsize))
                print w
                if textvertices:
                    P = line2d([[p[0] - w, p[1] - 0.9], [p[0] - w, p[1] + 1.1]], rgbcolor=c, thickness=thickness,  linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple)
                    P += line2d([[p[0] - w, p[1] + 1.1], [p[0] + w + 0.2, p[1] + 1.1]], rgbcolor=c, thickness=thickness,  linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple)
                    P += line2d([[p[0] + w + 0.2, p[1] + 1.1], [p[0] + w + 0.2, p[1] - 0.9]], rgbcolor=c, thickness=thickness,  linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple)
                    P += line2d([[p[0] + w + 0.2, p[1] - 0.9], [p[0] - w, p[1] - 0.9]], rgbcolor=c, thickness=thickness,  linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple)
                else:
                    polygon2d([[p[0] - w, p[1] - 0.9], [p[0] - w, p[1] + 1.1], [
                              p[0] + w + 0.2, p[1] + 1.1], [p[0] + w + 0.2, p[1] - 0.9]], fill=(not textvertices), rgbcolor=c, thickness=thickness,  linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple)
                vertices.append(P)
        for e in self.edges():
            v = e[0]
            if not self.has_vertex(e[0]) or not self.has_vertex(e[1]):
                print "deleting edge ", e
                self.delete_edge(e[0], e[1])
            else:
                c = simple_color if vertex_colors[
                    simple_color].count(v) > 0 else nonsimple_color
                if arrows:
                    edges.append(arrow([pos[e[0]][0], pos[e[0]][
                             1] + 1.1], [pos[e[1]][0], pos[e[1]][1] - 0.9], rgbcolor=c, zorder=-1, arrowsize=arrowsize, arrowshorten=arrowshorten, width=edges_thickness, linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple))
                else:
                    edges.append(line2d([[pos[e[0]][0], pos[e[0]][1] + 1.1], [pos[e[1]][0], pos[e[1]][1] - 0.9]], rgbcolor=c, zorder=-1, thickness=edges_thickness, linestyle=linestyle_simple if c == simple_color else linestyle_nonsimple))
        print "calculation ended"
        gp = self.graphplot(dpi=300, pos=pos, vertex_size=2050, figsize=round(
            float(max_vert) * 1.5), vertex_colors=vertex_colors)
        gp._plot_components['vertex_labels'] = labels
        gp._plot_components['vertices'] = vertices
        gp._plot_components['edges'] = edges
        #self._pos = pos
        return gp.plot(**options)

    def update_edges(self, vertex=None, h=None, fast=0):
        r'''
           Update the edges.
           This function inserts missing edges between the vertices of self.

           INPUT:
           - ``vertex``: A FQM_vertex. Only edges which are incoming edges for vertex are inserted.
                         If vertex=None (Default) all edges are updated.
           - ``h``: Only insert edges coming from height h - 1.
                    Note that this means that ``vertex`` has height ``h``.
                    This parameter is given for efficiency so we do not need to look up the height of ``vertex``.
           - ``fast``: An Integer in [0,1,2]:
                       An optimization parameter that determines which edges we insert in the following way:
                       fast = 0: insert all edges (with the given constrains)
                       fast = 1: insert at most one edge to non-simple module then quit
                       fast = 2: as fast = 1 but don\'t check for isomorphy
        '''
        primes = self._primes
        if h != None:
            if not self._heights.has_key(h - 1):
                return False
            else:
                verts = self._heights[h - 1]
        else:
            verts = self.vertices()
        for v in verts:
            for p in primes:
                for s in v.genus_symbol().C(p, False):
                    if vertex != None:
                        if fast == 2:
                            check = FQM_vertex(s) == vertex
                        else:
                            check = s.defines_isomorphic_module(
                                vertex.genus_symbol())
                        if check and not self.has_edge(v, vertex):
                            self.add_edge(v, vertex)
                            logger.debug(
                                "Adding edge: {0} - {1}".format(v.genus_symbol(), s))
                            if fast > 0 and self._vertex_colors[self._nonsimple_color].count(v) > 0:
                                return True
                    else:
                        if self.has_vertex(FQM_vertex(s)) and not self.has_edge(v, FQM_vertex(s)):
                            self.add_edge(v, FQM_vertex(s))
                            logger.debug(
                                "Adding edge: {0} - {1}".format(v.genus_symbol(), s))
                            # print "Adding edge: ", v, s

    def simple_dict(self):
        r"""
          Return a dictionary containing the $k$-simple finite quadratic modules
          labeled by their level.
        """
        d = {}
        if self._simple != None:
            for s in Set(self._simple):
                d[s.level()] = list()
            for s in Set(self._simple):
                d[s.level()].append(s)
        return d

    #*********************************
    # Properties
    #*********************************
    @property
    def level_limit(self):
        return self._level_limit

    @level_limit.setter
    def level_limit(self, N):
        self._level_limit = N

    @property
    def signature(self):
        return self._signature

    @signature.setter
    def signature(self, s):
        self._signature = s

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, k):
        self._weight = k

    @property
    def primes(self):
        return self._primes

    @primes.setter
    def primes(self, primes, test=True):
        self._primes = primes
        if test:
            for p in primes:
                if not is_prime(p):
                    self._primes.remove(p)


class FQM_vertex(object):

    def __init__(self, genus_symbol):
        self._genus_symbol = genus_symbol

    def genus_symbol(self):
        return self._genus_symbol

    def FQM(self):
        return FiniteQuadraticModule(self._genus_symbol)

    @cached_method
    def is_simple(self, k, reduction=True, bound=0):
        return self._genus_symbol.is_simple(k, reduction=reduction, bound=bound)

    def is_additive_lift_injective(self, r, s):
        return self._genus_symbol.is_additive_lift_injective(r,s)

    def __eq__(self, o):
        if o.genus_symbol() == self.genus_symbol():
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.genus_symbol())

    def __repr__(self):
        return self.genus_symbol().latex()

    # def __str__(self):
    #    return str(self.genus_symbol())

def test_a5prime_formula(upper=100, lower=1):
    for N in range(lower, upper):
        for s in anisotropic_symbols(N):
            a = s.a5prime()
            b = s.a5prime_formula()
            if not abs(a - b) < 1e-10:
                print s, a, b


def gaussum(n, N, prec=53):
    CC = ComplexField(prec)
    return sum(CC(exp(2 * CC.pi() * CC(0, 1) * n * m ** 2 / N)) for m in range(N))


def test_gaussum(s):
    return [(n, s.gaussum(n) / gaussum(n, s.level()) ** s.p_rank(2)) for n in range(s.level())]


@parallel
def test_dimension_jacobi(N, k):
    from nils.jacobiforms.dimension_jac_forms import dimension_jac_cusp_forms
    N = Integer(N)
    dj = dimension_jac_cusp_forms(k + Integer(1) / Integer(2), N)
    M = FiniteQuadraticModule([2 * N], [1 / (4 * N)])
    s = GenusSymbol(M.jordan_decomposition().genus_symbol())
    dv = s.dimension_cusp_forms(k, reduction=True)
    if not dj == dv:
        print "Error: ", N, dj, dv
        return False
    return True

#######################################################
# HELPER functions
# Mainly for output
######################################################

def dict_to_latex_table(d, n=None, max_per_table=10, midrules=False):

    beginstr = r"\begin{table}" + "\n" + r"\centering"
    if not n == None:
        beginstr += r"\caption{$" + str(Integer(n + 2) / 2) + \
            r"$-simple discriminant forms of signature $" + \
            str(str((2 - n) % 8)) + r"$}"
    beginstr += r"\begin{tabularx}{lX} \toprule level  & $k$-simple modules\\\midrule" + \
        "\n"
    endstr = r"\bottomrule\end{tabularx}\end{table}" + "\n"
    s = beginstr
    dk = sorted(d.keys())
    for i in range(len(d.keys())):
        N = dk[i]
        if len(d[N]) == 0:
            continue
        s += "$" + str(N) + "$" + " & "
        l = sorted(
            d[N], lambda x, y: int(-1) if x.order() <= y.order() else int(1))
        for m in l:
            s += str(FQM_vertex(m)).strip() + ", "
        s = s[:-2]
        if midrules:
            s += r"\\\midrule" + "\n"
        else:
            s += r"\\" + "\n"
        if (i + 1) % max_per_table == 0:
            s += endstr + beginstr
    s += endstr
    return s


def aniso_dict_to_latex_table(d, max_per_table=10, include_weight=True):
    if include_weight:
        beginstr = r"\begin{tabularx}{llX} \toprule $k$ & signature & $k$-simple modules\\\midrule" + \
            "\n"
    else:
        beginstr = r"\begin{tabular}{llX} \toprule $k$ & signature & $k$-simple modules\\\midrule" + \
            "\n"
    endstr = r"\bottomrule\end{tabularx}" + "\n"
    s = beginstr
    ks = sorted(d.keys())
    for i in range(len(d.keys())):
        k = ks[i]
        sigs = uniq([_.signature() for _ in d[k]])
        for sig in sigs:
            if include_weight:
                s += "$" + str(k) + "$" + " & "
            s += " " + str(sig) + " & "
            ms = [m for m in d[k] if m.signature() == sig]
            ms = sorted(
                ms, lambda x, y: int(-1) if x.order() <= y.order() else int(1))
            for m in ms:
                s += str(FQM_vertex(m)).strip() + ", "
            s = s[:-2]
            s += r"\\\midrule" + "\n"
            if (i + 1) % max_per_table == 0:
                s += endstr + beginstr
    s += endstr
    return s

########################
# to be removed/replaced
@parallel
def is_quotient(M, sym, rank):
    symbol = GenusSymbol_global_ring(MatrixSpace(ZZ, rank, rank).one())
    symbol._local_symbols = [
        Genus_Symbol_p_adic_ring(p, syms) for p, syms in sym.iteritems()]
    s = get_symbol_string(symbol)
    print s
    N = FiniteQuadraticModule(s)
    t = N.order() / M.order()
    if not Integer(t).is_square():
        return False
    else:
        t = sqrt(t)
    for p in Integer(N.order()).prime_factors():
        if not N.signature(p) == M.signature(p):
            return False
    # if not N.signature() == M.signature():
    #    return False
    for G in N.subgroups():
        if G.order() == t and G.is_isotropic():
            Q = G.quotient()
            if Q.is_isomorphic(M):
                print Q
                return N
            else:
                del Q
    del N
    return False


def Bbf(M, m):  # brute force algorithm
    if isinstance(M, str):
        M = FiniteQuadraticModule(M)
    rank = len(M.gens()) + 2
    sign = M.signature()
    symbols = all_symbols(sign, rank, M.order() * m ** 2)
    B = list(is_quotient([(M, sym, rank) for sym in symbols]))
    B = filter(lambda x: x[1] != False, B)
    B = map(lambda x: x[1], B)
    modules = Set(B)
    modules_unique = list()
    for M in modules:
        seen = False
        for N in modules_unique:
            if N.is_isomorphic(M):
                seen = True
        if not seen:
            modules_unique.append(M)
    return [M.jordan_decomposition().genus_symbol() for M in modules_unique]


def all_symbols(sign, rank, D):
    symbols = list()
    # print D
    fac = Integer(D).factor()
    symbols = list()
    for p, v in fac:
        psymbols = list()
        parts = Partitions(v)
        exp_mult = list()
        for vs in parts:
            exponents = Set(list(vs))
            if len(vs) <= rank:
                mult = dict()
                for vv in exponents:
                    mult[vv] = list(vs).count(vv)
                exp_mult.append(mult)
        Dp = D // (p ** v)
        for vs in exp_mult:
            # print "partition:", vs
            ell = list()
            if p == 2:
                for vv, mult in vs.iteritems():
                    ll = list()
                    for t in [0, 1]:  # even(0) or odd(1) type
                        for det in [1, 3, 5, 7]:  # the possible determinants
                            if mult % 2 == 0 and t == 0:
                                ll.append([vv, mult, det, 0, 0])
                            if mult == 1:
                                odds = [det]
                            elif mult == 2:
                                if det in [1, 7]:
                                    odds = [0, 2, 6]
                                else:
                                    odds = [2, 4, 6]
                            else:
                                odds = [
                                    o for o in range(8) if o % 2 == mult % 2]
                            for oddity in odds:
                                if t == 1:
                                    ll.append([vv, mult, det, 1, oddity])
                    ell.append(ll)
            else:
                for vv, mult in vs.iteritems():
                    ll = [[vv, mult, 1], [vv, mult, -1]]
                    ell.append(ll)
            l = list(itertools.product(*ell))
            l = map(lambda x: {p: list(x)}, l)
            psymbols = psymbols + l
        if len(symbols) == 0:
            symbols = psymbols
        else:
            symbols_new = list()
            for sym in symbols:
                for psym in psymbols:
                    newsym = deepcopy(sym)
                    newsym.update(psym)
                    symbols_new.append(newsym)
            symbols = symbols_new
    return symbols


def get_symbol_string(sym):
    symstr = ''
    # print sym._local_symbols
    for lsym in sym._local_symbols:
        p = lsym._prime
        for s in lsym.symbol_tuple_list():
            if s[0] == 0:
                continue
            if len(symstr) != 0:
                symstr = symstr + '.'
            symstr = symstr + str(p ** s[0])
            if p == 2:
                sgn = '+' if (Integer(s[2]).kronecker(2) == 1) else '-'
                if s[3] == 1:
                    symstr = symstr + '_' + str(s[4])
                symstr = symstr + '^' + sgn + str(s[1])
                # else:
                #    symstr = symstr + '^' + str(s[1])
            else:
                sgn = '+' if (s[2] == 1) else '-'
                symstr = symstr + '^' + sgn + str(s[1])
    return symstr

def SimpleModulesGraph2n(n, aniso_level_limit, **kwds):
    return SimpleModulesGraph((2 - n) % 8, QQ(2 + n) / QQ(2), aniso_level_limit, 2 + n, **kwds)

def SimpleModulesGraphn2(n, aniso_level_limit, **kwds):
    return SimpleModulesGraph((n-2) % 8, QQ(2 + n) / QQ(2), aniso_level_limit, 2 + n, r=2, s=n, **kwds)
