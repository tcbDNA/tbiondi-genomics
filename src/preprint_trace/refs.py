"""Step 2: compare each reference with its bibliographic record. usage: refs.py MANUSCRIPT.tex
Europe PMC by DOI; Crossref for three DOIs and Zenodo for one that Europe PMC does not hold. Responses are kept in
results/work/preprint_trace/refs/ and only fetched when missing. Writes results/posthoc2/preprint_refs_check.tsv."""
import re,json,os,sys,time,csv,urllib.request,urllib.parse
import paths
D=paths.WORK+"refs/"; os.makedirs(D,exist_ok=True)
s=open(sys.argv[1],encoding="utf8").read(); b=s[s.index("\\section*{References}"):]; b=b[:b.index("\\section*{Main Tables}")] if "\\section*{Main Tables}" in b else b
items=re.findall(r"\\item\s+(.*?)(?=\\item|\\end\{)",b,re.S)
def get(u,f):
    if not os.path.exists(f): open(f,"wb").write(urllib.request.urlopen(urllib.request.Request(u,headers={"User-Agent":"python-urllib"}),timeout=40).read()); time.sleep(0.4)
    return json.load(open(f))
out=[]; bad=0
for k,it in enumerate(items):
    it=it.split("\\vspace")[0]; m=re.search(r"doi\.org/([^}\s]+)",it); doi=m[1].replace("\\_","_")
    j=get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&query="+urllib.parse.quote(f'DOI:"{doi}"'),D+f"{k:02d}.json")
    r=[x for x in (j.get("resultList",{}).get("result") or []) if (x.get("doi") or "").lower()==doi.lower()][:1]; src="Europe PMC"
    if r:
        r=r[0]; ji=r.get("journalInfo") or {}; year,vol,iss,page,title=r.get("pubYear"),ji.get("volume"),ji.get("issue"),r.get("pageInfo"),r.get("title")
    elif "zenodo" in doi:
        src="Zenodo"; x=get("https://zenodo.org/api/records/"+doi.split(".")[-1],D+f"c{k:02d}.json")["metadata"]; year,vol,iss,page,title=x["publication_date"][:4],None,None,None,x["title"]
    else:
        src="Crossref"; x=get("https://api.crossref.org/works/"+doi,D+f"c{k:02d}.json")["message"]; year,vol,iss,page,title=str(x["published"]["date-parts"][0][0]),x.get("volume"),x.get("issue"),x.get("page") or x.get("article-number"),x["title"][0]
    t=re.sub(r"\s+"," ",re.sub(r"\\url\{[^}]*\}|\\[a-z]+\{|\}","",it)).replace("--","-").replace("\\ "," ")
    miss=[n for n,v in (("year",year),("volume",vol),("pages",page)) if v and not re.search(r"(?<!\d)"+re.escape(str(v).lstrip("E") if n=="pages" and re.fullmatch(r"E\d{1,2}",str(v)) else str(v))+r"(?!\d)",t)]
    if iss and "("+str(iss)+")" in t: issm="yes"
    else: issm="not given in the manuscript" if iss else "none in the record"
    bad+=bool(miss); out.append([k+1,doi,src,year,vol or "",iss or "",page or "",title,"yes" if not miss else "no: "+", ".join(miss),issm])
with open(paths.ROOT+"results/posthoc2/preprint_refs_check.tsv","w",newline="",encoding="utf8") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow("item doi source year volume issue pages title manuscript_matches issue_in_manuscript".split()); w.writerows(out)
print(len(out),"references; not matching:",bad)
