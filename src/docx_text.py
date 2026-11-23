"""Extract paragraph/table text from a .docx using the standard library only."""
import zipfile,re,sys,html
x=zipfile.ZipFile(sys.argv[1]).read("word/document.xml").decode()
for p in re.findall(r"<w:p[ >].*?</w:p>",x,flags=re.S):
    t=html.unescape(re.sub(r"<[^>]+>","",p)).strip()
    if t: print(t)
