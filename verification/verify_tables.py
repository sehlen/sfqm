"""Verify the computational tables of the bundled August 2014 extended paper.

Run using Sage's Python, from outside the repository after installing its wheel.
See README.md for the distinction between entry checks and completeness.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
from sage.all import QQ, ZZ, QQbar, prime_range
from sage.env import SAGE_VERSION
from sfqm.fqm.genus_symbol import GenusSymbol, anisotropic_symbols
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from psage.external.weil_invariants.weil_invariants import cython_invariants_dim
from psage.modform.weilrep_tools.weight_one_half import weight_one_half_dim

ROOT = Path(__file__).resolve().parent


def module(entry):
    if 'N' in entry:
        n = ZZ(entry['N'])
        return GenusSymbol(FiniteQuadraticModule([2*n], [1/(4*n)]).jordan_decomposition().genus_symbol())
    return GenusSymbol(entry['symbol'])


def exact_dimension(g, k, correct=True):
    """General dimension formula in QQbar, without rounding to an integer.

    Independent of the library's floating-point dimension evaluator, but shares
    its genus-value/character routines and exact low-weight invariant machinery.
    correct=False returns the virtual dimension before nonnegative corrections.
    """
    k = QQ(k)
    sig, m = ZZ(g.signature()), ZZ(g.order())
    if (2*k+sig) % 2:
        return ZZ(0)
    eps = 1 if (2*k+sig) % 4 == 0 else -1
    vals, v2 = g.values(), g.two_torsion_values()
    d = QQ(m+eps*g.torsion(2))/2
    a3 = QQ(sum((1-a)*v for a, v in vals.items() if a)
            + eps*sum((1-a)*v for a, v in v2.items() if a))/2
    a4 = QQ(vals.get(0, 0)+eps*v2.get(0, 0))/2
    def ch(t):
        a, b = g.char_invariant(t)
        return QQbar(a)*QQbar(b)
    g1, g2, g3 = ch(1), ch(2), ch(-3)
    a1 = d/4-QQbar(m).sqrt()/4*QQbar.zeta(8)**ZZ(2*k+sig)*(g2+eps*g2.conjugate())/2
    a2 = d/3+QQbar(m).sqrt()/(3*QQbar(3).sqrt())*(QQbar.zeta(24)**ZZ(4*k+3*sig-10)*(g1+eps*g3)).real()
    result = ZZ(QQbar(d+d*k/12-a1-a2-a3-a4))
    if correct and k == 2:
        result += cython_invariants_dim(g.finite_quadratic_module(), False, proof=True)
    elif correct and k == QQ(3)/2:
        result += weight_one_half_dim(g.finite_quadratic_module(), False, proof=True)
    return ZZ(result)


def simple(g, k):
    # The correction is nonnegative; positive virtual dimension is a certified
    # exclusion, with no expensive invariant calculation needed.
    if exact_dimension(g, k, correct=False) > 0:
        return False
    return exact_dimension(g, k) == 0


def unique_add(items, g):
    if any(g.defines_isomorphic_module(h) for h in items):
        return False
    items.append(g)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--entries', choices=['4','5','7','8','9','cyclic'])
    ap.add_argument('--search', type=int, help='signature (2,n), using paper seeds')
    ap.add_argument('--seed-level', type=int, default=0, help='also enumerate anisotropic seeds below this level; zero uses only paper seeds')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    start = time.monotonic()
    out = {'sage_version': SAGE_VERSION, 'python': platform.python_version(),
           'arithmetic': 'QQbar and exact characteristic-zero invariant kernels',
           'complete': False, 'results': []}
    def save():
        out['elapsed_seconds'] = time.monotonic()-start
        p = Path(args.output)
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix('.tmp')
        tmp.write_text(json.dumps(out, indent=2))
        tmp.replace(p)
    def load(name):
        return json.loads((ROOT/('table'+name+'.json')).read_text())
    try:
        if args.entries:
            entries = load(args.entries)
            out['table'] = args.entries
            out['expected_rows'] = len(entries)
            for entry in entries:
                g = module(entry)
                print('CHECK', args.entries, str(g), entry['k'], flush=True)
                row = dict(entry, symbol=str(g), dimension=int(exact_dimension(g, QQ(entry['k']))),
                           actual_signature=int(g.signature()), actual_rank=int(g.max_rank()), actual_level=int(g.level()))
                row['pass'] = row['dimension'] == 0 and row['actual_signature'] == entry['signature']
                if 'level' in entry:
                    row['pass'] &= row['actual_level'] == entry['level']
                if 'rank_limit' in entry:
                    row['pass'] &= row['actual_rank'] <= entry['rank_limit']
                if 'n' in entry:
                    row['actual_global'] = bool(g.is_global(2, entry['n']))
                    row['pass'] &= row['actual_global'] == entry['global']
                if entry.get('anisotropic'):
                    row['isotropic_elements'] = int(g.values().get(0, 0))
                    row['pass'] &= row['isotropic_elements'] == 1
                out['results'].append(row)
                save()
                print('RESULT', json.dumps(row), flush=True)
        if args.search is not None:
            n = args.search
            k, sig = QQ(n+2)/2, (2-n) % 8
            # Table 6, no dynamic or numerical prime pruning.
            pk = 73 if k == QQ(3)/2 else 37 if k == 2 else 29 if k == QQ(5)/2 else 19 if k == 3 else 17 if k == QQ(7)/2 else 13 if k == 4 else 11 if k <= QQ(13)/2 else 7 if k <= 9 else 5
            primes = list(prime_range(pk))
            seed_entries = load('4' if n == 1 else '5')
            seeds = []
            for entry in seed_entries:
                if QQ(entry['k']) == k and entry['signature'] == sig:
                    g = module(entry)
                    if g.max_rank() <= n+2:
                        if not simple(g, k):
                            raise AssertionError(('paper seed is not simple', str(g)))
                        unique_add(seeds, g)
            if args.seed_level:
                scanned = []
                for level in range(1, args.seed_level):
                    for g in anisotropic_symbols(level, sig):
                        if g.max_rank() <= n+2 and simple(g, k):
                            unique_add(scanned, g)
                out['seed_scan_extra'] = [str(g) for g in scanned if not any(g.defines_isomorphic_module(h) for h in seeds)]
                out['seed_scan_missing'] = [str(g) for g in seeds if g.level() < args.seed_level and not any(g.defines_isomorphic_module(h) for h in scanned)]
                for g in scanned:
                    unique_add(seeds, g)
            out.update(n=n, weight=str(k), primes=[int(p) for p in primes],
                       seeds=[str(g) for g in seeds], seed_level_limit_exclusive=args.seed_level,
                       completeness='Closure conditional on completeness of the paper anisotropic seed lists; no independent full-bound seed exclusion.')
            found, seen, index = list(seeds), set(), 0
            save()
            while index < len(found):
                g = found[index]
                index += 1
                for p in primes:
                    for h in g.C(p):
                        key = str(h)
                        if key in seen or h.max_rank() > n+2:
                            continue
                        seen.add(key)
                        if simple(h, k):
                            unique_add(found, h)
                out.update(found=[str(h) for h in found], processed=index, candidates=len(seen))
                save()
                print('PROGRESS', n, index, len(found), len(seen), flush=True)
            entries = load('9')+load('cyclic') if n == 1 else load('7') if n == 2 else [e for e in load('8') if e['n'] == n]
            expected = [module(e) for e in entries]
            out['missing'] = [str(g) for g in expected if not any(g.defines_isomorphic_module(h) for h in found)]
            out['extra'] = [str(g) for g in found if not any(g.defines_isomorphic_module(h) for h in expected)]
            out['global_genus_count'] = sum(bool(g.is_global(2, n)) for g in found)
        out['complete'] = True  # completion of this job, NOT an unconditional classification proof
    except BaseException as e:
        out['error'] = repr(e)
        raise
    finally:
        save()
        print('SUMMARY', json.dumps(out), flush=True)
    if any(not r['pass'] for r in out['results']) or out.get('missing') or out.get('extra'):
        raise SystemExit(1)

if __name__ == '__main__':
    main()
