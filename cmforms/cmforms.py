"""
Tools to check the dimension formulas for CM newforms.

AUTHOR: (c) Stephan Ehlen, 2013

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

from sage.rings.integer import Integer
from sage.modular.dirichlet import DirichletCharacter, DirichletGroup
from sage.modular.modform.constructor import Newforms
from sage.all import euler_phi, gcd, prime_range, squarefree_part, SageObject, sturm_bound, Gamma1, is_fundamental_discriminant, dimension_new_cusp_forms, load, save, copy, kronecker_character, is_even, is_odd, QuadraticField, ceil, floor, sqrt, RR, prod, log, is_squarefree, is_square, trivial_character

from nosqlite import client, server

import socket

from sage.parallel.decorate import *
from sage.misc.cachefunc import *

import sqlite3
import logging
import datetime

DEFAULT_PORT = 8212
dbclient = None
logger = None
dbserver = None
data_path = 'test'

class ColorFormatter(logging.Formatter):
    """
    This Formatter adds some colors for the moment.
    (Copied from the LMFDB.)
    """
    fmtString = '%(levelname)s:%(name)s@%(asctime)s: %(message)s'

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
            self._fmt = '\033[94m' + self._fmt
        elif record.levelno <= logging.INFO:
            self._fmt = '\033[92m' + self._fmt

        # bold, if module name matches
        if record.name == self._hl:
            self._fmt = "\033[1m" + self._fmt

        # reset, to unaffect the next line
        self._fmt += '\033[0m'

        return logging.Formatter.format(self, record)


def setup(path = 'test', port = DEFAULT_PORT):
    global logger, fh, ch, dbclient, dbserver, data_path

    if not os.path.exists(path):
        raise RuntimeError, "please create the data directory '%s'"%path
    data_path = path    
    
    logger = logging.getLogger('cm-{0}'.format(os.getpid()))
    logger.setLevel(logging.DEBUG)

    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter = ColorFormatter()

    if len(logger.handlers) == 0:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        fh = logging.FileHandler('cmforms-{0}.log'.format(os.getpid()))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    tries = 0
    maxtries = 10
    server_found = False
    while tries < maxtries:
        try:
            logger.info('Looking for running nosqlite server at port %d'%port)
            dbclient = client(port)
            dbclient.cmforms.collections()
            server_found = True
            break
        except:
            port += 1
            tries += 1
    if not server_found:
        port -= maxtries
        logger.info('No server found. Starting new.')
        dbserver = server(directory = os.path.join(path, 'db'), port = port)
        logger.info(dbserver)
        dbclient = client(dbserver.port)
        logger.info(dbclient)
        port = dbserver.port

    # sometimes I got timeouts when using the
    # nosqlite server/client
    socket.setdefaulttimeout(25)
    
    return 'server running on port %d'%port

class CmForms(SageObject):

    def __init__(self, Nchi, k, backend = None):
        self._k = k
        self._Nchi = Nchi
        if isinstance(Nchi, Integer):
            self._N = Nchi
            self._f = 1
            self._bound = sturm_bound(N, k)
            self._chi = trivial_character(N)
            self._logger = logging.getLogger("CmForms_{0}_{1}".format(N,k))
            if len(self._logger.handlers) == 0:
                self._logger.setLevel(logging.DEBUG)
                self._logger.addHandler(fh)
                self._logger.addHandler(ch)
        elif isinstance(Nchi, DirichletCharacter):
            self._chi = chi = Nchi
            self._N = N = chi.modulus()
            self._f = f = chi.conductor()
            self._logger = logging.getLogger("CmForms_{0}_{1}_{2}".format(N, character_index(chi) ,k))
            self._logger.setLevel(logging.DEBUG)
            if len(self._logger.handlers) == 0:
                self._logger.setLevel(logging.DEBUG)
                self._logger.addHandler(fh)
                self._logger.addHandler(ch)
            self._bound = sturm_bound(Gamma1(N), k) #take a lower bound?!!!
            self._logger.debug("in CmForms.init: N = {0}, f = {1}, bound = {2}".format(N,f,self._bound))
        else:
            raise ValueError("Nchi has to be an integer or a DirchichletCharacter")

        self._backend = backend
        if backend != None:
            self._logger.debug('Asking Backend for Newforms')
            if backend.has_newforms(chi,k):
                self._newforms = backend.load_newforms(chi,k)

    def newforms(self):
        self._logger.debug('Dimension Newforms: {0}'.format(dimension_new_cusp_forms(self._Nchi,self._k)))
        try:
            if self._newforms != None:
                return self._newforms
        except AttributeError:
            self._logger.debug("Started computing newforms for {0} and weight {1}".format(self._Nchi,self._k))
            self._newforms = Newforms(self._Nchi, self._k, names = 'a')
            self._logger.debug("Finished computing newforms for {0} and weight {1}".format(self._Nchi,self._k))
        return self._newforms

    def dimension(self, D = None, method = 'formula', ratio = Integer(1)/Integer(3)):
        N = self._N
        chi = self._chi
        bound = self._bound
        DG = DirichletGroup(N)
        self._logger.debug("in CmForms.dimension, chi = {0}, D = {1}".format(chi,D))
        #logger.debug(DG)
        if D==None:
            dr = [D for D in N.divisors() if is_fundamental_discriminant(-D)]
        else:
            dr = [D]
        if method == 'ratio' or method == 'exact':
            S = self.newforms()
            bound = max(bound, 100)
            #print 'Newforms:', S
            bound = len(S)*2*N # experimental hack
            test_primes = [p for p in prime_range(bound) if gcd(p,N) == 1]
            self._logger.debug('len(test_primes): {0}'.format(len(test_primes)))

        if method == 'ratio':
            logger.debug('Method: ratio, chi = {0}, D = {1}'.format(chi,D))
            ncm = 0
            for f in S:
                zerocount = 0
                for p in test_primes:
                    if f[p] == 0:
                        zerocount = zerocount + 1
                r = zerocount / Integer(len(test_primes))
                logger.debug('ratio: {0}'.format(r))
                if r > ratio:
                    ncm = ncm + f.hecke_eigenvalue_field().absolute_degree()/len(chi.galois_orbit())

        elif method == 'exact':
            self._logger.debug('Method: exact, chi = {0}, D = {1}'.format(chi,D))
            ncm = 0
            for d in dr:
                Sd = copy(S)
                psi = DG(kronecker_character(-d))
                for p in test_primes:
                    if psi(p) == -1:
                        for i,f in enumerate(Sd):
                            if f[p] != 0:
                                self._logger.debug('throwing away newform {0} because a_{1} != 0 and ({3}, {1}) = {4}'.format(i,p,f[p],-d,psi(p)))
                                Sd.pop(i)
                ncmd = sum([f.hecke_eigenvalue_field().absolute_degree()/len(chi.galois_orbit()) for f in Sd])
                self._logger.debug('CM forms with CM by {0}: {1}'.format(-d, ncmd))
                ncm = ncm + ncmd

        elif method == 'formula':
            self._logger.debug('Method: formula, chi = {0}, D = {1}'.format(chi,D))
            ncm = 0
            try:
                for d in dr:
                    ncm += self._dim_cmforms_formula_new_D(d)
            except NotImplementedError:
                raise
        return ncm

    def _dim_cmforms_formula_new_D(self, D):
        chi = self._chi
        k = self._k
        N = self._N
        f = self._f
        #print "D = {0}".format(D)
        #print 'got {0}'.format(chi)
        #check general compatibility of weight and resulting character
        if chi.is_odd() and is_even(k) or chi.is_even() and is_odd(k):
            return 0
        dimcmD = 0
        N = chi.modulus()
        DG = DirichletGroup(N)

        self._logger.debug("in _dim_cmforms_formula_new_D, D = {0}, chi = {1}".format(D, chi))

        if is_fundamental_discriminant(-D) and D.divides(N):
            K = QuadraticField(-D)
            psiD = kronecker_character(-D)
            fD = (chi*DG(psiD)).conductor()

            dimcmD = K.class_number()
            m = Integer(N//D)
            self._logger.debug("fD = {0}, m = {1}".format(fD, m))
            if not fD.divides(m):
                return 0
            # We also return the formula for D=3,4 for the moment
            # it should be an upper bound!
            #if D == 3 or D == 4:
                # TODO!
            #    raise NotImplementedError('D = 3 and D = 4 are not implemented at the moment.')
            for p,r in m.factor():
                self._logger.debug("p = {0}, r = {1}".format(p,r))
                dimcmp = 0
                if p.divides(D): # p ramified
                    self._logger.debug("p = {0} ramified for D = {1}".format(p,D))
                    s = ceil(Integer(r)/Integer(2))
                    t = floor(Integer(r)/Integer(2))
                    self._logger.debug("ordp(f_D) = {0}, s = {1}, t = {2}".format(fD.valuation(p), s, t))
                    if s < fD.valuation(p):
                        return 0
                    else:
                        if is_odd(r) and fD.valuation(p) == s:
                            dimcmp = p**t
                        elif is_even(r):
                            dimcmp = p**t-p**(t-1)
                        else:
                            dimcmp = 0
                elif psiD(p) == 1: # p split
                    self._logger.debug("p = {0} split for D = {1}".format(p,D))
                    s = fD.valuation(p)
                    if s > r:
                        return 0
                    elif r == s:
                        dimcmp = 2
                    elif s > r//2:
                        dimcmp = 2*p**(r-s-1)*(p-1)
                        if r-s >= 2:
                            dimcmp -= 2*p**(r-s-2)*(p-1)
                        elif r-s == 1:
                            dimcmp -= 2
                    elif is_odd(r):
                        dimcmp = 0
                    else:
                        dimcmp = euler_phi(p**(r//2)) - euler_phi(p**(r//2 -1))
                        if s == r//2:
                            dimcmp -= euler_phi(p**(r//2 -1))
                else: # p inert
                    self._logger.debug("p = {0} inert for D = {1}".format(p,D))
                    if is_odd(r):
                        return 0
                    else:
                        s = r//2
                        self._logger.debug("s = {0}".format(s))
                        if s < fD.valuation(p):
                            return 0
                        dimcmp += (p+1)*p**(s-1)
                        if s > fD.valuation(p):
                            if s >= 2:
                                dimcmp -= p**(s-2)*(p+1)
                            if s == 1:
                                dimcmp -= 1
                self._logger.debug("dimcmp = {0} for p = {1}".format(dimcmp,p))
                dimcmD = dimcmD*dimcmp
            self._logger.debug('CM forms with CM by {0}: {1}'.format(-D, dimcmD))
        return dimcmD



class CmFormsTester(SageObject):

    def __init__(self, k = 2, quadratic = False, save_tests = True, save_newforms= True, preserve = True, update_db = False, port = DEFAULT_PORT):
        self._backend = Backend(port = port, update_db = update_db)
        self.preserve = preserve
        self.save_newforms = save_newforms
        self.save_tests = save_tests
        self.quadratic = quadratic
        self.k = k
        logger.info('finished init of CmFormsTester')


    def characters(self, N):
        logger.debug('in CmFormsTester.characters')
        if self.quadratic:
            X = [chi for chi in characters(N) if chi.order() <= 2]
        else:
            X = characters(N)
        return X

    def run_test(self, r):
        l = [(chi, self.k) for N in r for chi in self.characters(N)]
        logger.debug("in run_test: l = {0}".format(l))
        self._params = l

        lr = list(self._comp_cmforms(l))
        #logger.debug(lr)
        errors = False
        for d in lr:
            logger.debug(d)
            if isinstance(d[1], str):
                logger.error(d)
                errors = True
                continue
            if d[1] != 'NO DATA' and len(d[1])>0:
                for D, td in d[1].iteritems():
                    dim_ex = td['dim']
                    dim_form = td['formula']
                    if dim_form != None:
                        if dim_ex > dim_form or (dim_ex != dim_form and D != 3 and D != 4):
                            logger.critical("Formula has an error for D = {0}; data: {1}".format(D,d))
                            raise RuntimeError((d, D))
        if not errors:
            return "Test successful"
        else:
            return "There were errors"

    @parallel(ncpus=15)
    def _comp_cmforms(self, chi, k):
        l = dict()
        N = chi.modulus()
        ci = character_index(chi)
        logger.debug("in _comp_cmforms: N = {0}, chi = {1}".format(N, ci))
        self._backend.add_running(N, ci, k)

        cm = CmForms(chi, k, self._backend)

        dr = [D for D in N.divisors() if is_fundamental_discriminant(-D)]
        #logger.debug(dr)

        for D in dr:
            fh = self._backend.formula_holds(chi, k, D)
            logger.debug("N = {0}, D = {1}, chi = {2}, formula_holds from backend: {3}".format(N, D, ci, fh))
            if not fh == True:
                fh = False
            if self.preserve and fh:
                logger.info('skipping N = %d, chi = %d, k= %d, D = %d'%(N, ci, k, D))
                continue
            else:
                logger.info('computing for N = %d, chi = %d, k= %d, D = %d'%(N, ci,k, D))

            dim_ex = cm.dimension(method = 'exact', D=D)
            logger.debug("We have for N = %d, chi = %d, D = %d that dim_ex = %d"%(N, ci, D, dim_ex))

            try:
                #if D != 3 and D != 4:
                logger.debug('running formula for N = %d, chi = %d, k= %d, D = %d'%(N, ci, k, D))
                dim_form = cm.dimension(method='formula', D=D)
                #else:
                #    dim_form = None
                #    logger.warn("No formula for D = {0}".format(D))
            except NotImplementedError:
                dim_form = None
                logger.exception("No formula for D = {0}".format(D))

            if not dim_form == None:
                logger.debug("We have for N = %d, D = %d that  dim_ex = %d, dim_form =%d"%(N, D, dim_ex, dim_form))

            if self.save_tests and dim_form != None:
                self._backend.save_test(chi, k, D, dim_ex, dim_form)

            l[D] = {'dim': dim_ex, 'formula': dim_form}

            if dim_form != None:
                if dim_ex > dim_form or (dim_ex != dim_form and D != 3 and D != 4):
                    logger.critical("Formula has an error for chi = {0}, k = {1}, D = {2} we have {3} CM forms but formula says {4}".format(chi, k, D, dim_ex, dim_form))

        logger.debug("Before saving newforms for N = %d, chi =%d"%(N, character_index(chi)))
        if self.save_newforms and not self._backend.has_newforms(chi,k):
            newforms = cm.newforms()
            if len(l) > 0:
                dimcm = sum([d['dim'] for d in l.values()])
            else:
                logger.warn("No suitable D")
                dimcm = 0
            self._backend.save_newforms(chi, k, newforms, dimension_new_cusp_forms(chi,k), dimcm)
        logger.debug("Finished: N = {0}, char = {1}, results(l) = {2}".format(N, character_index(chi), l))
        self._backend.remove_running(N, ci, k)
        return l

class Backend(object):

    def __init__(self, update_db = False, port = DEFAULT_PORT, auto_setup = True):
        global dbclient, data_path
        if dbclient == None and auto_setup:
            setup(port=port)
        self._client = dbclient
        self._path = data_path
        #, username='cm', password=r'sdjHHTESWQ~~kfb423895y2bBKHHAg.NbsTVDJKB/sdkjsdnmm/:::^#&_')
        self._db = dbclient.cmforms
        
        if update_db:
            self.update_db()

    def space_name(self, N, i, k):
        return os.path.join(self._path,'%05d-%04d-%02d'%(N,i,k))

    def cmforms_test_name(self, N, i, k, D, t):
        if t or t == '':
            t = ''
        else:
            t = '-error'
        return os.path.join(self.space_name(N, i, k), '%05d%s.sobj'%(D,t))

    def newforms_name(self, N, i, k):
        return os.path.join(self.space_name(N,i,k),'newforms.sobj')

    def save_test(self, chi, k, D, dim, form):
        N = chi.modulus()
        i = character_index(chi)
        if form != dim:
            t = 0
            if D == 3 or D == 4:
                if dim <= form:
                    t = 2 # status for upper bound holds
        else:
            t = 1

        f = self.space_name(N, i, k)
        if not os.path.exists(f):
            os.makedirs(f)

        f = self.cmforms_test_name(N, i, k, D, t)
        rd = {'char': chi, 'k': k, 'D': D, 'dim_cmforms': dim, 'formula': form}
        logger.debug('Saving test: {0}'.format(rd))

        save(rd, f)

        cmforms_data = self._db.cmforms_data
        try:
            r = cmforms_data.find_one(N = N, char = i, k = k, D = D)
            u = cmforms_data.update({'dim_cmforms': dim, 'formula_holds': t}, N = N, char = i, k = k, D = D)
        except ValueError:
            ins = cmforms_data.insert({'dim_cmforms': dim, 'formula_holds': t}, N = N, char = i, k = k, D = D)

        logger.debug('Finished saving test: {0}'.format(rd))

    def save_newforms(self, chi, k, newforms, dim_newforms, dim_cmforms):
        N = chi.modulus()
        i = character_index(chi)

        f = self.space_name(N, i, k)
        if not os.path.exists(f):
            os.makedirs(f)

        f = self.newforms_name(N, i, k)

        logger.debug("Saving newforms %s"%f)

        save(newforms,f)

        newforms_data = self._db.newforms_data
        try:
            r = newforms_data.find_one(N = N, char = i, k = k)
            u = newforms_data.update({'dim_newforms': dim_newforms, 'dim_cmforms': dim_cmforms, 'newforms_saved': 1}, N = N, char = i, k = k)
        except ValueError:
            ins = newforms_data.insert({'dim_newforms': dim_newforms, 'dim_cmforms': dim_cmforms, 'newforms_saved': 1}, N = N, char = i, k = k)
        logger.debug("Finished saving newforms %s"%f)

    def load_newforms(self, chi, k):
        N = chi.modulus()
        i = character_index(chi)

        f = self.space_name(N, i, k)
        if not os.path.exists(f):
            return None

        f = self.newforms_name(N, i, k)

        try:
            logger.debug("Loading newforms %s"%f)
            newforms = load(f)
        except:
            logger.exception("Loading newforms %s did not work"%f)
            raise RuntimeError('Newforms could not be loaded!')

        return newforms

    def formula_holds(self, chi, k, D):
        N = chi.modulus()
        i = character_index(chi)
        cmforms_data = self._db.cmforms_data
        try:
            r = cmforms_data.find_one(N = N, char = i, k = k, D = D)
            logger.debug(r)
            return (r['formula_holds'] == 1)
        except ValueError:
            logger.debug('in formula_holds: -1')
            return -1

    def has_newforms(self, chi, k):
        logger.debug("in has_newforms: chi = {0}, k = {1}".format(character_index(chi),k))
        N = chi.modulus()
        i = character_index(chi)
        newforms_data = self._db.newforms_data
        try:
            r = newforms_data.find_one(N = N, char = i, k = k)
            return r['newforms_saved']
        except ValueError:
            return False

    def find_data(self):
        """
        Return iterator ...

        Here N = level, k = weight, char = character,
        dim_newforms = dimension of newforms (cuspidal),
        D = discriminant,
        formula_holds = 1 if our formula is correct in this case (=0 otherwise - uh),
        newforms_saved = 1 if we saved the neforms for this space.
        """
        for Nik in os.listdir(self._path):
            z = Nik.split('-')
            if len(z) == 3:
                N, i, k = map(lambda x: int(x),z)
                #print N, i, k
                newforms_saved = os.path.exists(self.newforms_name(N, i, k))
                if newforms_saved:
                    newforms_saved = 1
                else:
                    newforms_saved = 0
                #print newforms_saved
                chi = character(N, i)
                dnew = int(dimension_new_cusp_forms(chi, k))

                d = [None]
                cmtot = 0

                Dl = [int(D) for D in Integer(N).divisors() if is_fundamental_discriminant(-D)]
                #print Dl

                for x in os.listdir(os.path.join(self._path, Nik)):
                    if x.endswith('newforms.sobj'):
                        continue
                    else:
                        #print x
                        Dt = x[:-5].split('-')
                        D = int(Dt[0])
                        t = ''
                        if len(Dt) == 2:
                            formula_holds = int(0)
                            t = Dt[1]
                        else:
                            formula_holds = int(1)
                        data = load(self.cmforms_test_name(N,i,k,D,t))
                        #print data
                        if data.has_key('dim_cmforms'):
                            dimcm = int(data['dim_cmforms'])
                        elif data.has_key('dim'):
                            dimcm = int(data['dim'])
                        else:
                            continue
                        d.append({'N': N, 'char': i, 'k': k, 'D':D, 'dim_cmforms': dimcm, 'formula_holds': formula_holds})
                        cmtot += dimcm
                        Dl.remove(D)
                #print 'nachher:', Dl
                if len(Dl)>0:
                    cmtot = int(-1)
                d[0] = {'N': N,'char': i, 'k': k, 'dim_newforms': dnew, 'dim_cmforms': cmtot, 'newforms_saved': newforms_saved}
                yield d

    def update_db(self):
        newforms_data = self._db.newforms_data
        cmforms_data = self._db.cmforms_data
        newforms_data.delete()
        cmforms_data.delete()

        # 2. iterate over known data
        for t in self.find_data():
            logger.debug(t)
            newforms_data.insert(t[0])
            t.remove(t[0])
            if len(t) > 0:
                cmforms_data.insert(t)

    def _migrate(self):
        for NiDkt in os.listdir(self._path):
            if not NiDkt.endswith('.sobj'):
                continue
            z = NiDkt[:-5].split('-')
            if len(z) == 4 or len(z) == 5:
                if len(z) == 4:
                    t = True
                    N,i,D,k = z
                    logger.debug(N, i, D, k)
                else:
                    N,i,D,k,t = z
                    logger.debug(N, i, D, k, t)
                    t = False
                sn = self.space_name(int(N), int(i), int(k))
                #snp = os.path.join(self._path, sn)
                #print snp
                if not os.path.exists(sn):
                    os.makedirs(sn)
                tn = self.cmforms_test_name(int(N),int(i),int(k),int(D),t)
                logger.debug(tn)
                os.rename(os.path.join(self._path,NiDkt),tn)
                
    def add_running(self, N, i, k):
        self._db.jobs.insert(N = N, char = i, k = k, timestamp=datetime.datetime.now())

    def remove_running(self, N, i, k):
        self._db.jobs.delete(N = N, char = i, k = k)

def character_index(chi):
    return characters(chi.modulus()).index(chi)

@cached_function
def characters(N):
    return  [X[0] for X in DirichletGroup(N).galois_orbits()]

@cached_function
def character(N, i):
    return characters(N)[i]

def newforms_estimate(N,k):
    return RR(Integer(k-1)/Integer(12)*euler_phi(N)*0.733957-(Integer(1)/Integer(2)+Integer(7)/Integer(12)*(k-1))*sqrt(N))

def cmforms_estimate(N,k):
    N = Integer(N)
    N0 = prod([p for p in N.prime_divisors() if p != 2])
    #N0 = 8*N0
    if Integer(4).divides(N):
        if Integer(8).divides(N):
            N0=8*N0
        else:
            N0=4*N0
    est = RR(6.87552932554442*sqrt(N)*N0**(Integer(1)/Integer(4))*log(N0)/RR.pi())
    return est 

def typeA_test(r, k, method='formula'):
    for N in range(1,r):
        if is_squarefree(N) or (N/squarefree_part(N) in [4,8]):
            DG = DirichletGroup(N)
            for chi in DG:
                if chi.order()<=2:
                    if is_even(k) and chi.is_odd():
                        continue
                    if is_odd(k) and chi.is_even():
                        continue
                    if not is_odd(gcd(N/chi.conductor(),chi.conductor())):
                        continue
                    logger.info(N)
                    d = dimension_new_cusp_forms(chi,k)
                    cm = CmForms(chi,k).dimension(method=method)
                    if d - cm <= 0:
                        logger.critical('0-Hit: N={0}, chi={1}, dim={2}, dimcm={3}'.format(N, chi, d, cm))



def type_trivial_test(r, k, method='formula'):
    '''
       Test for modules Witt-equivalent to the trivial
       module that are not tested by the typeA test.
    '''
    for p in prime_range(1,r):
        stop = False
        N = 1
        while not stop:
            N = N*p**2
            chi = trivial_character(N)
            d = dimension_new_cusp_forms(chi,k)
            cm = CmForms(chi,k).dimension(method=method)
            if d - cm <= 0:
                logger.critical('0-Hit: {0}, {1}, {2}'.format(N, d, cm))
            else:
                stop = True
