from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, Genus_Symbol_p_adic_ring, is_GlobalGenus
#from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from sage.matrix.matrix_space import MatrixSpace
from sage.all import ZZ, Zmod, sys, parallel, is_prime, colors, cached_function, Integer, Partitions, Set
import itertools

@parallel
def is_quotient(M,sym,rank):
    symbol=GenusSymbol_global_ring(MatrixSpace(ZZ,rank,rank).one())
    symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in sym.iteritems()]
    s = get_symbol_string(symbol)
    print s
    N = FiniteQuadraticModule(s)
    t = N.order()/M.order()
    if not Integer(t).is_square():
        return False
    else:
        t = sqrt(t)
    for p in N.order().prime_factors():
        if not N.signature(p) == M.signature(p):
            return False
    #if not N.signature() == M.signature():
    #    return False
    for G in N.subgroups():
        if G.order() == t and G.is_isotropic():
            Q = G.quotient()
            if Q.is_isomorphic(M):
                print Q
                return N
            else:
                del(Q)
    del(N)
    return False

@cached_function
def B(M,m, algorithm='theorem',rec=True,dict=False):
    if algorithm == 'bf':
        return Bbf(M,m)
    elif not algorithm == 'theorem':
        raise ValueError('algorithm has to be one of "theorem" and "bf" (for brute force)')
    Bs = list()
    if not is_prime(m):
        raise ValueError('m is supposed to be prime (general case: TODO)')
    p = m

    if isinstance(M,FiniteQuadraticModule_ambient):
        sym_M = genus_symbol_string_to_dict(M.jordan_decomposition().genus_symbol())
    elif isinstance(M, str):
        if len(M)>0:
            sym_M = genus_symbol_string_to_dict(M)
        else:
            sym_M = {}
    else:
        raise ValueError('M needs to be a "FiniteQuadraticModule" or a "str".')
    if not sym_M.has_key(p):
        sym_M[p]=list()
    #print sym_M
    #sym_M = get_genus_symbol_dict(M)

    if p == 2:
        sym_new = deepcopy(sym_M)
        sym_new[p].append([1,2,7,0,0])
        Bs.append(sym_new)
        print sym_new, dict_to_genus_symbol_string(sym_new)
        sym_new = deepcopy(sym_M)
        sym_new[p].append([1,2,7,1,0])
        Bs.append(sym_new)
        print sym_new, dict_to_genus_symbol_string(sym_new)
        for s in sym_M[p]:
            print s
            t=s[3]
            print t
            if rec:
                for ss in sym_M[p]:
                    if t==1 and ss[0]==s[0]+1:
                        sym_new=deepcopy(sym_M)
                        sym_new[p].remove(s)
                        if s[2] in [1,5]:
                            s2n = s[2]+2
                        else:
                            s2n = s[2]-2
                        sym_new[p].append([s[0],s[1],s2n,s[3],(s[4]+2)%8])
                        sym_new[p].remove(ss)
                        if ss[2] in [1,5]:
                            ss2n = ss[2]+2
                        else:
                            ss2n = ss[2]-2
                        sym_new[p].append([ss[0],ss[1],ss2n,ss[3],(ss[4]+2)%8])
                        Br = B(dict_to_genus_symbol_string(sym_new),p,rec=False,dict=True)
                        print "Br=", B(dict_to_genus_symbol_string(sym_new),p,rec=False,dict=True)
                        Bs = Bs+Br
            if s[1]>2 or t == 0 or t == None:
                # q^+2 type II -> (2q)^+2 type II
                if s[1] > 2 or s[4]==0:
                    print "rule q^+2 -> (2q)^+2"
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,7,0,0])
                    if s[1] > 2: sym_new[p].append([s[0],s[1]-2,-s[2],s[3],s[4]])
                    Bs.append(sym_new)
                    print dict_to_genus_symbol_string(sym_new)
                # q^-2 -> (2q)_4^-2
                if s[1] > 2 or s[4]==4:
                    print "rule q^-2 -> (2q)_4^-2"
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,3,1,4])
                    if s[1] > 2: sym_new[p].append([s[0],s[1]-2,(s[2]*3) % 8, s[3], (s[4]+4) % 8])
                    Bs.append(sym_new)
                    print sym_new, dict_to_genus_symbol_string(sym_new)
                # q^+2 -> (2q)_0^+2
                if s[1] > 2 or s[4]==0:
                    print "rule q^+2 -> (2q)_0^+2"
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,7,1,0])
                    if s[1] > 2: sym_new[p].append([s[0],s[1]-2,-s[2],s[3],s[4]])
                    Bs.append(sym_new)
                    print dict_to_genus_symbol_string(sym_new)
            if t == 1:
                print "t=1"
                #rule for odd type, mult 1 power up
                odds = [1,3,5,7]
                for o in odds:
                    if (s[1] == 1 and (s[4]==o or s[0]==1 and s[4] % 8 == (o+4) %8)) \
                           or s[1]>3 \
                           or (s[1]==3 and not ((s[4]-o) % 8 ) in [4,6] and kronecker(s[2]*o,2) % 8 == 1) \
                           or (s[1]==2 and (s[4]-o) % 8 in odds):
                        sym_new = deepcopy(sym_M)
                        sym_new[p].remove(s)
                        sym_new[p].append([s[0]+2,1,o,1,o])
                        if s[1] > 1:
                                sym_new[p].append([s[0],s[1]-1,s[2]*o,1,(s[4]-o) % 8])
                        print dict_to_genus_symbol_string(sym_new)
                        Bs.append(sym_new)
                print s[1]
                if s[1] >= 2:
                    print s
                    if s[4] == 0 or (s[0]==1 and s[4]==4):
                        print "1,4"
                        print sym_M
                        sym_new = deepcopy(sym_M)
                        sym_new[p].remove(s)
                        sym_new[p].append([s[0]+1,2,7,0,0])
                        if s[1]>2:
                            sym_new[p].append([s[0],s[1]-2,(-s[2]) % 8, 1, s[4]])
                        Bs.append(sym_new)
                        print sym_new
                        print dict_to_genus_symbol_string(sym_new)
                    if s[4] == 4:
                        sym_new = deepcopy(sym_M)
                        sym_new[p].remove(s)
                        sym_new[p].append([s[0]+1,2,3,0,4])
                        if s[1]>2:
                            sym_new[p].append([s[0],s[1]-2,(s[2]*3) % 8, 1, (s[4]-4)%8])
                        Bs.append(sym_new)
                        print sym_new
                if s[1]>2 and s[1] % 2 == 0:
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,s[2],0,s[4]])
                    sym_new[p].append([s[0],s[1]-2,0,0,0])
                    Bs.append(sym_new)
        
    else:
        ep = -1 if p % 4 == 3 else 1
        if not sym_M.has_key(p):
            sym_M[p]=list()
        #rule 1+2 for triv p-module
        sym_new = deepcopy(sym_M)
        sym_new[p].append([2,1,1])
        Bs.append(sym_new)
        sym_new = deepcopy(sym_M)
        sym_new[p].append([2,1,-1])
        Bs.append(sym_new)
        sym_new = deepcopy(sym_M)
        sym_new[p].append([1,2,ep])
        Bs.append(sym_new)
        for s in sym_M[p]:
            sym_new = deepcopy(sym_M)
            sym_new[p].remove(s)
            #rule 1 (power+2, mult 1)
            sym_new[p].append([s[0]+2,1,s[2]])
            if s[1] > 1:
                sym_new[p].append([s[0],s[1]-1,1])
                Bs.append(sym_new)
                sym_new = deepcopy(sym_M)
                sym_new[p].remove(s)
                sym_new[p].append([s[0]+2,1,-s[2]])
                sym_new[p].append([s[0],s[1]-1,-1])
                Bs.append(sym_new)
            else:
                Bs.append(sym_new)
            #rule 2 (power+1, mult 2)
            if s[1] >= 2:
                if s[1]==2 and s[2]==ep:
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,ep])
                    Bs.append(sym_new)
                elif s[1]>2:
                    sym_new = deepcopy(sym_M)
                    sym_new[p].remove(s)
                    sym_new[p].append([s[0]+1,2,ep])
                    sym_new[p].append([s[0],s[1]-2,ep*s[2]])
                    Bs.append(sym_new)
    if dict:
        return Bs
    Bss = list()
    for sym in Bs:
        sym = reduce_symbol(sym)
        sym = dict_to_genus_symbol_string(sym)
        Bss.append(sym)
    if p == 2:
        Bs = list()
        for s in Bss:
            try:
                FiniteQuadraticModule(s)
                Bs.append(s)
            except:
                print "ERROR: ", s
        Bss = Bs
        print "Bs = ", Bs
        #mods_uniq=list()
        #for M in mods:                               
        #    seen = False                                 
        #    for N in mods_uniq:
        #        if N.is_isomorphic(M):
        #            seen = True
        #    if not seen:
        #        mods_uniq.append(M)
        #Bss = [M.jordan_decomposition().genus_symbol() for M in mods_uniq]
    return Bss

@cached_function
def is_simple(symstr,k):
    if not symstr == '':
        M = FiniteQuadraticModule(symstr)
        V = VectorValuedModularForms(M)
        d = V.dimension_cusp_forms(k)
    else:
        d = dimension_cusp_forms(1,k)
    if d == 0:
        return True
    else:
        return False

def simple_modules_graph(N_limit, sig, k, primes, rank_limit, only_simple=False):
    graphs=list()
    for N in range(N_limit):
        if N==1 or N==2 or is_fundamental_discriminant(N) or is_fundamental_discriminant(-N):
            s = anisotropic_symbols(N,sig)
            print s
            if len(s)==0:
                continue
            for ss in s:
                print "ss=", ss
                g,data=simple_modules_graph_from_aniso(ss,primes,k,rank_limit)
                graphs.append([g,data])
    print graphs
    simple = list()
    for g in graphs:
        if isinstance(g[1][0],list):
            simple = simple + g[1][0]
    simple = reps(simple)
    print "Found in total %d %d-simple modules with p-rank <= %d"%(sum([len(_) for _ in simple.values()]),k, rank_limit)
    return graphs, simple

def reps(l):
    simple_dict = dict()
    mods_seen = list()
    for s in simple:
        if not s == '':
            M = FiniteQuadraticModule(s)
        else:
            M = FiniteQuadraticModule()
        if not simple_dict.has_key(M.level()):
            simple_dict[M.level()] = list()
        add = True
        for N in mods_seen:
            if M.level() == N.level() and M.is_isomorphic(N):
                add = False
        if add:
            simple_dict[M.level()].append(s)
            mods_seen.append(M)
    return simple_dict

def dict_to_latex_table(d, max_per_table=10):
    beginstr = r"\begin{tabular}{|l|r|} \hline level   & $2$-simple modules\\\hline" + "\n"
    endstr = r"\end{tabular}" + "\n"
    s = beginstr
    for i in range(len(d.keys())):
        N = d.keys()[i]
        s += "$" + str(N) + "$" + " & "
        for m in d[N]:
            s += str(FQM_vertex(m)).strip() + ", "
        s=s[:-2]
        s += r"\\\hline" + "\n"
        if (i+1) % max_per_table == 0:
            s += endstr + beginstr
    s += endstr
    return s


def prime_anisotropic_symbols(p):
    print p
    if not (p==1 or is_prime(p) or p == 4 or p == 8):
        raise ValueError
    else:
        if is_odd(p) and is_prime(p):
            if p % 4 == 3:
                return [str(p) + '^+1', str(p) + '^-1', str(p) + '^+2']
            else:
                return [str(p) + '^+1', str(p) + '^-1', str(p) + '^-2']
        else:
            if p==1:
                return ['']
            if p == 2:
                return ['2^-2']
            elif p==4:
                return ['2_1^+1', '2_7^+1', '2_2^+2', '2_0^+2', '2_3^-3', '2_5^-3']
            elif p==8:
                return ['4_1^+1', '4_3^-1', '4_5^-1', '4_7^+1', '2_1^+1.4_1^+1', '2_1^+1.4_3^-1', '2_1^+1.4_5^-1', '2_1^+1.4_7^+1']

def anisotropic_symbols(N,sig):
    syms = []
    if N == 1 and sig % 8 == 0: return prime_anisotropic_symbols(1)
    for p,n in Integer(N).factor():
        print p,n
        if len(syms)==0:
            syms=prime_anisotropic_symbols(p**n)
            print syms
        else:
            syms_old=deepcopy(syms)
            syms=list()
            for r in itertools.product(syms_old,prime_anisotropic_symbols(p**n)):
                syms.append(r[0]+'.'+r[1])
            print syms
        if len(syms)==0:
            return []
    print syms
    syms_new = list()
    for s in syms:
        print s
        if (FiniteQuadraticModule(s).signature()) % 8 == sig:
            syms_new.append(s)
    return syms_new

def max_rank(s):
    d = genus_symbol_string_to_dict(s)
    m = 0
    for p,sl in d.iteritems():
        r = 0
        for s in sl:
            r += s[1]
        m = max(r,m)
    return m
        

def simple_modules_graph_from_aniso(s,p,k, rank_limit, g=None, heights=None, vertex_colors = None, only_simple=False):
    #M = FiniteQuadraticModule(s)
    print "Starting with ", s, p,k
    if isinstance(p,list):
        primes = p
    else:
        primes = [p]
    print "primes in from aniso: ", primes
    simple = is_simple(s,k)
    if type(g) == type(None):
        g = DiGraph()
    s = FQM_vertex(s)
    g.add_vertex(s)
    if heights == None:
        heights = dict()
    h = 0
    if not heights.has_key(h):
        heights[h] = [s]
    else:
        if heights[h].count(s)==0:
            heights[h].append(s)
    simple_color = colors.red.rgb()
    nonsimple_color = colors.green.rgb()
    if vertex_colors == None or len(vertex_colors)==0:
        vertex_colors = dict()
        vertex_colors[simple_color]=list()
        vertex_colors[nonsimple_color]=list()
    Bs = [s]
    print Bs
    print type(s)
    if simple:
        if vertex_colors[simple_color].count(s)==0:
            vertex_colors[simple_color].append(s)
    elif vertex_colors[nonsimple_color].count(s)==0:
        vertex_colors[nonsimple_color].append(s)
    while simple:
        h = h+1
        if not heights.has_key(h):
            heights[h]=list()
        Bss = list()
        simple = False
        for s1 in Bs:
            Bs2 = list()
            for p in primes:
                Bs2 = Bs2 + B(s1.genus_symbol(),p)
            print "Bs2 = ", Bs2
            for s2 in Bs2:
                if not max_rank(s2) <= rank_limit:
                    print "skipping ", s2
                    continue
                s2 = FQM_vertex(s2)
                if not g.has_edge(s1,s2):
                    g.add_vertex(s2)
                    g.add_edge(s1,s2)
                    heights[h].append(s2)
                update_edges(g,primes)
                has_nonsimple_neighbor = False
                for e in g.incoming_edges(s2):
                    print s2, e
                    if vertex_colors[nonsimple_color].count(e[0])>0:
                        has_nonsimple_neighbor = True
                        print "Has nonsimple neighbor: ", s2
                        break
                if has_nonsimple_neighbor:
                    #not simple
                    if not vertex_colors[simple_color].count(s2)>0: vertex_colors[nonsimple_color].append(s2)
                else:
                    print 'checking for simple module: ', s2
                    if is_simple(s2.genus_symbol(),k):
                        simple = True
                        print "Found simple module: ", s2
                        Bss.append(s2)
                        if not vertex_colors[simple_color].count(s2)>0: vertex_colors[simple_color].append(s2)
                    else:
                        if not vertex_colors[nonsimple_color].count(s2)>0: vertex_colors[nonsimple_color].append(s2)
                    
        Bs = Bss
    print g
    print len(vertex_colors[nonsimple_color])
    if only_simple:
        for v in vertex_colors[nonsimple_color]:
            if g.has_vertex(v): g.delete_vertex(v)
            #while vertex_colors[nonsimple_color].count(v)>0:
            #    vertex_colors[nonsimple_color].remove(v)
            for j in range(len(heights)):
                while heights[j].count(v)>0:
                    heights[j].remove(v)
        vertex_colors[nonsimple_color]= dict()
    print g

    pos = dict()
    labels=list()
    vertices = list()
    edges = list()
    min_d = 0.5
    widths = [float(sum([len(str(v))/float(6) + min_d for v in heights[i]])) for i in range(len(heights))]
    print widths
    max_w = max(widths)
    if len(widths)>1:
        min_w = min([w for w in widths if w > 3])
    else:
        min_w = 0
    print min_w, max_w, widths
    max_vert = max([len(_) for _ in heights.values()])
    for i in range(len(heights)):
        real_w = widths[i]
        prev_w = min_w if i == 0 else widths[i-1]
        next_w = min_w if i == len(heights)-1 else widths[i+1]
        cur_w = max(float(next_w)*0.9, float(prev_w)*0.9, real_w)
        d = max(2*float(cur_w - real_w)/len(heights[i]),min_d)
        print real_w, cur_w
        print "d = ", d
        p = [-(cur_w), float(max_vert)/2*i]
        w = 0
        for j in range(len(heights[i])):
            v = heights[i][j]
            p = [p[0] + w, p[1]]
            w = float(len(str(v)))/float(6)
            p[0] = p[0] + w + d
            pos[heights[i][j]] = p
            c = simple_color if vertex_colors[simple_color].count(v)>0 else nonsimple_color
            labels.append(text(str(v), (p[0]+0.2,p[1]), rgbcolor=colors.white.rgb(), zorder=8, fontsize=14))
            print w
            P = polygon2d([[p[0]-w,p[1]-1],[p[0]-w,p[1]+1],[p[0]+w,p[1]+1], [p[0]+w,p[1]-1]], fill=True, rgbcolor=c)
            vertices.append(P)
    for e in g.edges():
        v = e[1]
        if not g.has_vertex(e[0]) or not g.has_vertex(e[1]):
            print "deleting edge ", e
            g.delete_edge(e[0],e[1])
            continue
        c = simple_color if vertex_colors[simple_color].count(v)>0 else nonsimple_color
        edges.append(arrow([pos[e[0]][0],pos[e[0]][1]+1],[pos[e[1]][0],pos[e[1]][1]-1], rgbcolor=c, zorder=-1, arrowshorten=8))
    gp = g.graphplot(dpi=300, pos=pos, vertex_size=2050, figsize=round(float(max_vert)*1.5), vertex_colors=vertex_colors)
    gp._plot_components['vertex_labels'] = labels
    gp._plot_components['vertices'] = vertices
    gp._plot_components['edges'] = edges
    simple = [v.genus_symbol() for v in vertex_colors[simple_color]]
    return gp,[simple,pos,g]

def update_edges(g,primes):
    for v in g.vertices():
        for p in primes:
            for s in B(v.genus_symbol(),p):
                if g.has_vertex(FQM_vertex(s)) and not g.has_edge(v,FQM_vertex(s)):
                    g.add_edge(v,FQM_vertex(s))
                    print "Adding edge: ", v, s

class FQM_vertex(object):

    def __init__ (self, genus_symbol):
        self._genus_symbol = genus_symbol

    def genus_symbol(self):
        return self._genus_symbol

    def __eq__(self, o):
        if o.genus_symbol() == self.genus_symbol():
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.genus_symbol())

    def __repr__(self):
        o = r'$'
        l = self.genus_symbol().split('.')
        for s in l:
            ss=s.split('^')
            if len(ss) == 2:
                s1,s2 = ss
            else:
                s1 = ss[0]
                s2 = '+1'
            o = o + s1 + r'^{' + s2 + r'}'
        o = o + r'$'
        return LatexExpr(o)
           

def reduce_symbol(sym):
    for p,sl in sym.iteritems():
        if p != 2:
            sl.sort()
            for s in sl:
                #print s
                ssl = copy(sl)
                for j in range(sl.index(s)+1,len(ssl)):
                    t = ssl[j]
                    if s[0] == t[0]:
                        s[1] = s[1]+t[1]
                        s[2] = s[2]*t[2]
                        sl.remove(t)
        else:
            sl.sort()
            for s in sl:
                #print s
                ssl = copy(sl)
                for j in range(sl.index(s)+1,len(ssl)):
                    t = ssl[j]
                    if s[0] == t[0]:
                        print s, t
                        if s[3] == None:
                            s[3] = 0
                        if t[3] == None:
                            t[3] = 0
                        s[1] = s[1]+t[1]
                        if s[3] == 0 and t[3] == 1:
                            s[2] = t[2]
                        elif s[3] == t[3]:
                            s[2] = (s[2]*t[2]) % 8
                        s[4] = (s[4]+t[4]) % 8
                        if (s[3] == 0  and t[3]==1) or (t[3] == 0  and s[3]==1):
                            s[3] = 1
                        print s
                        sl.remove(t)
            print sym, dict_to_genus_symbol_string(sym)
    return sym

def dict_to_genus_symbol_string(sym):
    symbol = GenusSymbol_global_ring(MatrixSpace(ZZ,1,1).one())
    symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in sym.iteritems()]
    s = get_symbol_string(symbol)
    return s

def genus_symbol_string_to_dict(s):
    sl = s.split('.')
    d = dict()
    for s in sl:
        L1 = s.split('^')
        if len(L1)>2:
            raise ValueError()
        elif len(L1) == 1:
            nn = 1
        else:
            nn = L1[1]
        n = Integer(nn)
        L1= L1[0].split("_")
        q = Integer(L1[0])
        if len(L1) > 2: # more than one _ in item
            raise ValueError
        elif len(L1) == 2:
            if Integer(L1[1]) in range(8):
                t = Integer(L1[1])
            else:
                raise ValueError, "Type given, which ist not in 0..7: %s"%(L1[1])
        else:
            t = None
        if not (n != 0 and q != 1 and q.is_prime_power()
                and ( None == t or (is_even(q) and t%2 == n%2))
                and ( not (None == t and is_even(q)) or 0 == n%2)
                ):
            raise ValueError,"{0} is not a valid signature!".format(s)
        p = q.prime_factors()[0]
        r = q.factor()[0][1]
        eps = sign(n)
        n = abs(n)
        if not d.has_key(p):
            d[p]=list()
        if p==2:
            if t == None:
                print "eps = ", eps
                if not is_even(n):
                    raise ValueError()
                d[p].append([r, n, 3*(-1)**(Integer(n-2)/2) % 8 if eps == -1 else (-1)**(Integer(n)),t,4 if eps == -1 else 0])
                print d
            else:
                if t.kronecker(2) == eps:
                    det = t % 8
                else:
                    if eps == -1:
                        det = 3
                    else:
                        det = 1
                d[p].append([r,n,det,1,t % 8])
        else:
            d[p].append([r,n,eps])
    return d
        

def Bbf(M,m): #brute force algorithm
    if isinstance(M,str):
        M = FiniteQuadraticModule(M)
    rank = len(M.gens()) + 2
    sign = M.signature()
    symbols = all_symbols(sign,rank,M.order()*m**2)
    B = list(is_quotient([(M,sym,rank) for sym in symbols]))
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


def all_symbols(sign,rank,D):
    symbols=list()
    #print D
    fac = Integer(D).factor()
    symbols=list()
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
        Dp=D//(p**v)
        for vs in exp_mult:
            #print "partition:", vs
            ell = list()
            if p==2:
                for vv,mult in vs.iteritems():
                    ll=list()
                    for t in [0,1]: # even(0) or odd(1) type
                        for det in [1,3,5,7]: # the possible determinants
                            if mult % 2 == 0 and t==0:
                                ll.append([vv,mult,det,0,0])
                            if mult==1:
                                odds=[det]
                            elif mult==2:
                                if det in [1,7]:
                                    odds=[0,2,6]
                                else:
                                    odds=[2,4,6]
                            else:
                                odds=[o for o in range(8) if o%2==mult%2]
                            for oddity in odds:
                                if t==1:
                                    ll.append([vv,mult,det,1,oddity])
                    ell.append(ll)
            else:
                for vv,mult in vs.iteritems():
                    ll=[[vv,mult,1],[vv,mult,-1]]
                    ell.append(ll)
            l = list(itertools.product(*ell))
            l = map(lambda x: {p: list(x)},l)
            psymbols=psymbols+l
        if len(symbols)==0:
            symbols=psymbols
        else:
            symbols_new=list()
            for sym in symbols:
                for psym in psymbols:
                    newsym=deepcopy(sym)
                    newsym.update(psym)
                    symbols_new.append(newsym)
            symbols=symbols_new
    return symbols

def get_symbol_string(sym):
    symstr = ''
    #print sym._local_symbols
    for lsym in sym._local_symbols:
        p=lsym._prime
        for s in lsym.symbol_tuple_list():
            if s[0]==0:
                continue
            if len(symstr)!=0:
                symstr=symstr + '.'
            symstr = symstr + str(p**s[0])
            if p == 2:
                sgn = '+' if (Integer(s[2]).kronecker(2) == 1) else '-' 
                if s[3]==1:
                    symstr = symstr + '_' + str(s[4])
                symstr = symstr + '^' + sgn + str(s[1])
                    #else:
                    #    symstr = symstr + '^' + str(s[1])
            else:
                sgn = '+' if (s[2] == 1)  else '-'
                symstr = symstr + '^' + sgn + str(s[1])
    return symstr
                
