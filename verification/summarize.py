"""Summarize preserved verification JSON; never treat a partial file as success."""
import json
from pathlib import Path
root=Path(__file__).resolve().parent
rows=[]
for p in sorted((root/'results').glob('*.json')):
    if p.name == 'summary.json': continue
    d=json.loads(p.read_text())
    passed=d.get('complete',False) and not d.get('error') and all(r.get('pass',False) for r in d.get('results',[])) and not d.get('missing') and not d.get('extra')
    rows.append(dict(file=p.name,completed=d.get('complete',False),passed=passed,entries=len(d.get('results',[])),found=len(d.get('found',[])),seconds=round(d['elapsed_seconds'],3),error=d.get('error')))
# This is an accounting cross-check, conditional on Proposition 6.2, not a
# recomputation of the number of integral isometry classes in each genus.
counts={}
for table in ['7','8','9','cyclic']:
    for e in json.loads((root/('table'+table+'.json')).read_text()):
        if e['global']: counts[e['n']]=counts.get(e['n'],0)+1
counts[1]+=1
expected={1:256,2:67,3:15,4:5,5:3,6:4,7:2,8:3,9:3,10:2,18:1,26:1}
assert counts==expected and sum(counts.values())==362
print(json.dumps({'runs':rows,'table1_accounting':counts,'total_conditional_on_proposition_6_2':362},indent=2))
