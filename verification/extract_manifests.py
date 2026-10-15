"""Extract machine-readable entries; requires PyMuPDF and arXiv v1 dimest.tex.

Usage: python verification/extract_manifests.py /path/to/dimest.tex
Table 9 is from the bundled extended PDF, absent from the arXiv source.
"""
import hashlib,json,re,sys
from pathlib import Path
import fitz
ROOT=Path(__file__).resolve().parent
s=Path(sys.argv[1]).read_text()
def save(name,rows):
    (ROOT/('table'+name+'.json')).write_text(json.dumps(rows,indent=2)+'\n')
    print(name,len(rows))
def symbols(text):
    out=[]
    for token in re.findall(r'\$([^$]+)\$',text):
        if '^{' not in token: continue
        for t in token.split(','):
            t=t.strip()
            parts=re.findall(r'(\d+)(?:_(\d+))?\^\{([+-]\d+)\}',t)
            if parts: out.append({'symbol':'.'.join(a+('_'+b if b else '')+'^'+c for a,b,c in parts),'global':not t.startswith('[')})
    return out
part=s[s.index('\\caption{The table shows the $70$'):s.index('\\caption{The table shows all finite quadratic')]
rows=symbols(part)
for r in rows:r.update(k='2',signature=0,rank_limit=4,n=2)
assert len(rows)==70
save('7',rows)
part=s[s.index('\\caption{\\label{simpleanisotable}'):s.index('\\subsection{Differences of dimensions}')]
rows=[]
for row in part.split('\\\\'):
    cells=row.split('&')
    if len(cells)<3 or not symbols(cells[2]):continue
    weights=[]
    for x in re.findall(r'\$([^$]+)\$',cells[0]):
        m=re.fullmatch(r'\\frac\{(\d+)\}\{(\d+)\}',x)
        if m:weights.append('/'.join(m.groups()))
        elif x.isdigit():weights.append(x)
    for k in weights:
        for i in [1,3]:
            if len(cells)<=i+1:continue
            sig=re.sub(r'[^0-9]','',cells[i])
            if not sig:continue
            for r in symbols(cells[i+1]):r.update(k=k,signature=int(sig),anisotropic=True);rows.append(r)
assert len(rows)==79
save('5',rows)
part=s[s.index('\\caption{The 75'):s.index('\\caption{\\label{simpleanisotable}')]
rows=[]
for sig,start,end in [(1,'$2_3^{+3}3^{+1}','\\midrule\n   $5$'),(5,'$4_5^{-1}$','\\\\\\bottomrule')]:
    p=part[part.index(start):];p=p[:p.index(end)]
    for r in symbols(p):r.update(k='3/2',signature=sig,anisotropic=True);rows.append(r)
def squarefree(n):return all(n%(p*p) for p in range(2,int(n**.5)+1))
Ns=[n for n in range(1,37) if squarefree(n)]+[38,39,41,42,46,47,51,55,59,62,66,69,70,71,78,87,94,95,105,110,119]
rows += [dict(N=n,k='3/2',signature=1,anisotropic=True) for n in Ns]
assert len(rows)==75
save('4',rows)
Ns=list(range(1,37))+[38,39,40,41,42,44,45,46,47,48,49,50,51,52,54,55,56,59,60,62,63,64,66,68,69,70,71,72,75,76,78,80,81,84,87,90,94,95,96,98,100,104,105,108,110,119,120,126,132,140,144,150,168,180]
rows=[dict(N=n,k='3/2',signature=1,rank_limit=3,n=1) for n in Ns]
for r in rows:r['global']=True
assert len(rows)==90
save('cyclic',rows)
TABLE8={
3:'2_7^+1 2_7^+3 2_7^+1.4^+2 2_7^+5 2_7^+3.4^+2 2_7^+1.4^+4 4_7^+1 2^+2.4_7^+1 2^+4.4_7^+1 2_1^+1.3^+1 2_7^+1.3^-2 2_7^+1.3^+4 8_3^-1 8_7^+1 2^+2.8_3^-1 [2^+4.8_3^-1]',
4:'3^+1 3^-3 3^+5 2^+4.3^+1 [2^+6.3^+1] 2^+2.3^+1',
5:'4_5^-1 2^+2.4_5^-1 2^+4.4_5^-1 [2^+6.4_5^-1]',
6:'2^-2 2^-4 2^-6 [2^-8] 5^+1',
7:'4_3^-1 2_1^+1.3^-1',8:'3^-1 2_2^+2 7^+1',
9:'2_1^+1 4_1^+1 8_1^+1',10:'1^+1 2^+2',18:'1^+1',26:'1^+1'}
rows=[]
for n,entries in TABLE8.items():
    for e in entries.split():
        r=dict(symbol=e.strip('[]'),k=str((n+2)//2) if n%2==0 else str(n+2)+'/2',signature=(2-n)%8,rank_limit=n+2,n=n)
        r['global']=not e.startswith('[');rows.append(r)
save('8',rows)
D=fitz.open(ROOT.parent/'bruinier_ehlen_freitag_extended.pdf');rows=[];level=None
for page,ymin in [(35,375),(36,85)]:
    chars=[]
    for b in D[page].get_text('rawdict')['blocks']:
        for line in b.get('lines',[]):
            for span in line['spans']:
                for c in span['chars']:
                    x,y=c['origin'];t=c['c'].replace('−','-')
                    if ymin<y<730 and not t.isspace():chars.append((x,y,t,span['size']))
    ys=sorted(set(round(y,1) for x,y,t,size in chars if size>10 and (t.isdigit() or t in ',[]')))
    for y in ys:
        base=sorted(c for c in chars if c[3]>10 and abs(c[1]-y)<.15)
        lev=''.join(c[2] for c in base if c[0]<105)
        if lev.isdigit():level=int(lev)
        base=[c for c in base if c[0]>=105]
        if not base or any(c[2] not in '0123456789,[]' for c in base):continue
        small=[c for c in chars if c[3]<10 and -6<c[1]-y<4]
        tokens=[];i=0
        while i<len(base):
            x,yy,t,size=base[i]
            if t in ',[]':tokens.append(t);i+=1;continue
            num=t;end=x;i+=1
            while i<len(base) and base[i][2].isdigit() and base[i][0]-end<8:
                end=base[i][0];num+=base[i][2];i+=1
            stop=base[i][0] if i<len(base) else 550
            scripts=sorted(c for c in small if end+2<c[0]<stop)
            exp=''.join(c[2] for c in scripts if c[1]<y)
            sub=''.join(c[2] for c in scripts if c[1]>y)
            assert re.fullmatch('[+-][0-9]+',exp),(page,y,num,exp,sub)
            tokens.append(num+('_'+sub if sub else '')+'^'+exp)
        line=''.join(('.' if i and t[0].isdigit() and tokens[i-1][-1].isdigit() else '')+t for i,t in enumerate(tokens))
        for item in line.split(','):
            if item:
                r=dict(level=level,symbol=item.strip('[]'),printed_page=page+1,k='3/2',signature=1,rank_limit=3,n=1)
                r['global']=not item.startswith('[');rows.append(r)
assert len(rows)==169
save('9',rows)
(ROOT/'sources.json').write_text(json.dumps({'tex_url':'https://arxiv.org/src/1408.4148v1','tex_sha256':hashlib.sha256(s.encode()).hexdigest(),'pdf_sha256':hashlib.sha256((ROOT.parent/'bruinier_ehlen_freitag_extended.pdf').read_bytes()).hexdigest(),'note':'Tables 4,5,7 from v1 TeX; Table 8 manually transcribed; Table 9 extracted from bundled extended PDF.'},indent=2)+'\n')
