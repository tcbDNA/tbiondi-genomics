"""Find the number tokens of an entry and class them (number, label, accession, hash, version, date, cite, path, url, ref)."""
import re,json
import paths
doc=json.load(open(paths.WORK+'doc.json'))
# order matters: first pattern that covers a token decides its class
PAT=[('ref',r'\[ref:[^\]]*\]'),('url',r'\\url\{[^}]*\}|https?://\S+'),('path',r'(?:results|reports|src|env|data)/[\w./*\-]+|\b[\w]+\.(?:tsv|txt|py|md|c|cm|jsonl|lock)\b|\\path\{[^}]*\}'),
('hash',r'\b(?=[0-9a-f]*\d)(?=[0-9a-f]*[a-f])[0-9a-f]{7,64}\b'),('date',r'\b20\d\d-\d\d-\d\d(?: \d\d:\d\d)?\b'),
('accession',r'\b(?:NC|GCA|NP)_\d+(?:\.\d+)?|\b(?:SRR|ERR|DRR|PRJNA|PRJEB|SAMN)\d+\b|\b(?:CM|OZ|KU|PV|PX|PZ|QB|AP|AB|ON|OP|MN|MK|KY|KT|KJ|LR|PP|OQ|DQ|OR|FN|GQ|LDKF|JBQQXK|CDFYFP)\d{5,}(?:\.\d+)?|\bPV962\b|Pfam-A_enone'),
('version',r'\bv\d+(?:\.\d+)*(?:[-.]r?\d+)*\b|HMMER 3\.3\.2|HMMER v?3\.4|GPT-6\.2|\brefseq89m\b|RefSeq release 238|release 238|OrthoDB v13|Codetta archive|minimap2 v2\.31-r1302|\(v3\)|\bv3\b'),
('label',r'\b[ACGUTN]-?33C?\b|\bU36C\b|\b[TSQKI](?:216S?|203|178|191S?)\b|\bG34\b|\bCOX[123]\b|\bND(?:4L|[1-6])\b|\bATP[68]\b|\bcode-?s? ?\d+\b|\bcode 2\b|translation table \d+|table 5\b|\bT[0-6]\b|\b[AB]1\b|\bC[1-5]\b(?=,)|\bTier 1–3|\bTier \d(?: or \d)?|\bT[123]\b|GATE_\d A\d+|GATE_\d|Gate \d|round \d|Table S\d|Tables? S\d(?: and S\d)?|Note S\d|S\d and S\d|Figures? S\d|MTRF1|MRPL58|mtRF1a|\bRF1\b|RF1a|position 33|positions? (?:203|216|39|34–36)(?: and 216)?|position 216|\bRF2a\b|§ ?\d+[\w.]*|§11b|HI-2016|\b5\'|3\'|\^\{[a-z]\}|C1|fig\d|alpha|S-\w+|MITOS2|BLAST\+|TRNA'),
('cite',r'(?:et al\.|[A-Z][\w-]+) \(?(?:19|20)\d\d(?:, (?:19|20)\d\d)?\)?|^\(?(?:19|20)\d\d\)?[;.,):]?|\b(?:19|20)\d\d\)')]
def tokens(t):
    spans=[]
    for k,p in PAT:
        for m in re.finditer(p,t):
            if not any(a<m.end() and m.start()<b for a,b,_ in spans): spans.append((m.start(),m.end(),k))
    out=[]
    for m in re.finditer(r'\d[\d,]*(?:\.\d+)?',t):
        v=m[0].rstrip(',')
        k=next((k for a,b,k in spans if a<=m.start() and m.end()<=b+1 and m.start()<b),'number')
        pct='%' in t[m.end():m.end()+2] or bool(re.match(r'[–,\d.% ]*(?:–|, | or | and )?[\d.,]*%',t[m.end():m.end()+14]))
        out.append((v,k,m.start(),pct))
    return out
if __name__=='__main__':
    import collections
    c=collections.Counter()
    for a,b,t in doc:
        for x in tokens(t): c[x[1]]+=1
    print(c)
    seen=collections.Counter()
    for a,b,t in doc:
        if a in('References','Main Tables'): continue
        for v,k,s,p in tokens(t):
            if k=='cite' or (k=='number' and re.fullmatch(r'(19|20)\d\d',v)): seen[(k,t[max(0,s-18):s+8])]+=1
    print([x for x in seen if x[0]=='number'])
