"""Release check: README.md must hold every paragraph of the REPORT.md Abstract word for word (whitespace aside). Exit 1 if not."""
import re,sys,os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),".."); n=lambda t:re.sub(r"\s+"," ",t); s=open(os.path.join(R,"reports/REPORT.md")).read(); q=n(open(os.path.join(R,"README.md")).read())
ab=[x for x in s[s.index("## Abstract\n")+12:s.index("\n\n## 1. Introduction")].split("\n\n") if x.strip()]; bad=[x[:80] for x in ab if n(x) not in q]; print(len(ab),"abstract paragraphs;",len(bad),"not in README.md word for word",bad); sys.exit(1 if bad else 0)
