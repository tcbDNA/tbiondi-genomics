"""Step 1 of the number trace: flatten the manuscript (LaTeX) into [section, kind, text] entries. usage: extract.py MANUSCRIPT.tex"""
import re,sys,json,os
import paths
s=open(sys.argv[1],encoding='utf8').read()
def cl(t):
    for _ in range(4): t=re.sub(r'\\(textit|textbf|texttt|path|emph|seqsplit|mathrm)\{((?:[^{}]|\{[^{}]*\})*)\}',r'\2',t)
    t=re.sub(r'~?\\ref\{([^}]*)\}',lambda m:' [ref:'+m[1].replace('fig:fig','F').replace('tab:','T-').replace('sec:','S-')+']',t); t=re.sub(r'\\label\{[^}]*\}','',t)
    t=t.replace('\\ ',' ').replace('\\%','%').replace('\\&','&').replace('\\_','_').replace('---',' — ').replace('--','–').replace('\\S ','§').replace('\\S','§').replace('{,}',',')
    t=t.replace('\\ge','≥').replace('\\le','≤').replace('\\times','×').replace('\\sim','~').replace('\\cdot','·').replace('\\ne','≠').replace('\\to','→').replace('$','').replace('``','"').replace("''",'"')
    return re.sub(r'\s+',' ',t).strip()
doc=[]; sec='Header'
body=s[s.index('\\begin{document}'):]
# blocks: split by lines, track section, caption, table rows
intab=None
for line in body.split('\n'):
    if line.startswith('%') or not line.strip(): continue
    m=re.match(r'\\(sub)?section\*?\{((?:[^{}]|\{[^{}]*\})*)\}',line)
    if m: sec=cl(m[2]); doc.append((sec,'heading',cl(m[2]))); continue
    if '\\caption{' in line: doc.append((sec,'caption',cl(line[line.index('\\caption{')+9:]))); continue
    if re.search(r' & ',line) and line.rstrip().endswith('\\\\'):
        for k,c in enumerate(re.split(r'(?<!\\)&',line.rstrip()[:-2])): doc.append((sec,f'table cell {k+1}',cl(c)))
        continue
    if re.match(r'\\(begin|end|setlength|renewcommand|centering|toprule|midrule|bottomrule|vspace|clearpage|includegraphics|makebox|end\w+|noindent\\small)',line) and '\\parbox' not in line: 
        if not re.match(r'\\vspace.*',line): pass
        if re.match(r'\\(makebox|includegraphics|begin|end|setlength|renewcommand|centering|toprule|midrule|bottomrule|clearpage|vspace)',line): continue
    t=cl(line)
    for x in re.split(r'(?<=[.;:])\s+(?=[A-Z(\["≥0-9])',t): 
        if x: doc.append((sec,'text',x))
json.dump(doc,open(paths.WORK+'doc.json','w'))
import collections
c=collections.Counter(); n=collections.Counter()
for a,b,t in doc:
    k=len(re.findall(r'\d[\d,]*(?:\.\d+)?',t)); c[a]+=k; n[b.split()[0]]+=k
for k,v in c.items(): print(v,k)
print(n,sum(c.values()),len(doc))
