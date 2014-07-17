from .fqm.genus_symbol import GenusSymbol
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from .simple.simple_modules_graph import SimpleModulesGraph, SimpleModulesGraph2n
from sage.all import load

def compute_simple_shimura_curves(path=''):
    l = load(path + '/lattices_aniso_1_notgamma0.sobj')
    simshim = []
    aniso_shim_syms=[]
    for Q in l:
        M = FiniteQuadraticModule(Q.matrix())
        s = GenusSymbol(M.jordan_decomposition().genus_symbol())
        aniso_shim_syms.append(s)
    G = SimpleModulesGraph2n(1,None)
    G.compute_from_startpoints(aniso_shim_syms)
    simshim.extend(G._simple)
    return simshim
